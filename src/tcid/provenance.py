"""Hashing, code version and manifests for traceable simulations."""

from __future__ import annotations

import hashlib
import json
import platform
import subprocess
import sys
from datetime import datetime, timezone
from functools import lru_cache
from pathlib import Path

import numpy as np

from . import __version__


def _default(o):
    if isinstance(o, np.generic):
        return o.item()
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, Path):
        return str(o)
    raise TypeError(f"not serialisable: {type(o)}")


def canonical_json(obj) -> str:
    return json.dumps(obj, sort_keys=True, default=_default, separators=(",", ":"))


def config_hash(obj, n: int = 16) -> str:
    """Hash of a complete configuration (cache key)."""
    return hashlib.sha256(canonical_json(obj).encode()).hexdigest()[:n]


@lru_cache(maxsize=1)
def git_commit() -> str:
    root = Path(__file__).resolve().parents[2]
    try:
        c = subprocess.run(["git", "rev-parse", "--short", "HEAD"], cwd=root, capture_output=True, text=True, timeout=5)
        commit = c.stdout.strip() or "uncommitted"
        d = subprocess.run(["git", "status", "--porcelain", "src"], cwd=root, capture_output=True, text=True, timeout=5)
        return commit + ("+dirty" if d.stdout.strip() else "")
    except Exception:
        return "unknown"


def code_version() -> str:
    return f"tcid-{__version__}@{git_commit()}"


def environment() -> dict:
    import gmsh
    import pandas
    import scipy
    import shapely
    import skfem
    import sklearn

    return {
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "scikit-fem": skfem.__version__,
        "gmsh": gmsh.__version__,
        "shapely": shapely.__version__,
        "scikit-learn": sklearn.__version__,
        "pandas": pandas.__version__,
    }


def file_sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def write_manifest(directory: Path, config: dict, extra: dict | None = None) -> Path:
    """Manifest listing every result file with its checksum, the config hash and environment."""
    directory = Path(directory)
    files = []
    for p in sorted(directory.rglob("*")):
        if p.is_file() and p.name != "manifest.json":
            files.append({"path": str(p.relative_to(directory)), "bytes": p.stat().st_size, "sha256": file_sha256(p)})
    manifest = {
        "created": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "code_version": code_version(),
        "config_hash": config_hash(config),
        "config": config,
        "environment": environment(),
        "files": files,
        **(extra or {}),
    }
    out = directory / "manifest.json"
    out.write_text(json.dumps(manifest, indent=2, default=_default))
    return out
