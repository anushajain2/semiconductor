"""Compact reproducibility benchmark.

Run from a clean checkout of a tagged release: re-runs the solver verification gates and
recomputes a small fixed set of stored confirmation-batch solves **without the cache**, then
compares them with the committed results. Because geometry generation, meshing and the direct
solver are deterministic, the recomputed resistances should match the stored ones to rounding.
The benchmark records the exact code version it ran with; a dirty working tree is flagged.
"""

from __future__ import annotations

import json
import time
from dataclasses import replace
from pathlib import Path

import numpy as np
import pandas as pd

from .config import Design
from .provenance import code_version, config_hash, environment
from .simulate import make_jobs, run_jobs

TOLERANCE = 1e-9  # relative


def run_benchmark(study, out_dir: str | Path = "results/benchmark", n_seeds: int = 3) -> dict:
    from .verification import run_all

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    version = code_version()
    t0 = time.perf_counter()
    gates = run_all(out_dir)
    cands = json.loads((study.out / "final" / "candidates.json").read_text())
    designs = {k: Design(**cands[k]) for k in ("recommended", "baseline_random_orientation")}
    seeds = study.budget["confirm_seeds"][:n_seeds]
    ms_ref = replace(study.cfg.mesh, size_factor=float(study.budget["final_refinements"][0]))
    jobs = make_jobs(list(designs.values()), seeds, {**study.base, "mesh": ms_ref}, "benchmark")
    rows, fails = run_jobs(jobs, workers=study.workers, cache_dir=None, progress=False)  # no cache: recompute
    stored = pd.read_parquet(study.out / "final" / "final_rows.parquet")
    stored = stored[stored["tag"] == "confirm"].set_index(["design_id", "seed", "scenario"])["r_layer"]
    cmp = []
    for _, r in rows.iterrows():
        key = (r["design_id"], int(r["seed"]), r["scenario"])
        ref = stored.get(key)
        cmp.append({"design_id": r["design_id"], "seed": int(r["seed"]), "scenario": r["scenario"],
                    "r_layer_recomputed": r["r_layer"], "r_layer_stored": ref,
                    "rel_diff": abs(r["r_layer"] - ref) / ref if ref is not None else np.nan})
    cmp = pd.DataFrame(cmp)
    cmp.to_csv(out_dir / "benchmark_comparison.csv", index=False)
    max_diff = float(cmp["rel_diff"].max()) if len(cmp) else float("nan")
    res = {
        "code_version": version,
        "clean_tree": "dirty" not in version and "uncommitted" not in version,
        "verification_gates_pass": all(g["passed"] for g in gates),
        "designs": {k: config_hash(d.to_dict(), 12) for k, d in designs.items()},
        "seeds": seeds,
        "n_solves_recomputed": int(len(cmp)),
        "n_failures": int(len(fails)),
        "max_rel_diff_vs_stored": max_diff,
        "matches_stored_results": bool(len(cmp) and max_diff <= TOLERANCE and not len(fails)),
        "tolerance": TOLERANCE,
        "wall_time_s": time.perf_counter() - t0,
        "environment": environment(),
        "note": "stored results were produced from working trees recorded as dirty; this benchmark ties them "
                "to the tagged source by exact recomputation of a subset",
    }
    (out_dir / "benchmark.json").write_text(json.dumps(res, indent=2, default=str))
    return res
