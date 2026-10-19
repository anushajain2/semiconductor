"""Configuration objects for layers, geometry generation, meshing and designs.

Lengths are stored in metres. Every object serialises to a plain dict so that a
complete simulation configuration can be hashed for caching and traceability.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field, replace
from pathlib import Path

import yaml

from .units import to_si


@dataclass(frozen=True)
class Layer:
    """A finite adhesive layer: heat flows through ``thickness`` between isothermal faces."""

    thickness: float = 100e-6
    width: float = 200e-6
    margin: float = 3e-6  # resin-rich margin between particles and every boundary

    def to_dict(self) -> dict:
        return asdict(self)


@dataclass(frozen=True)
class GeometrySettings:
    min_gap: float = 1.0e-6  # minimum particle-particle surface separation
    d_min: float = 1.0e-6  # truncation of the equivalent-diameter distribution
    d_max: float = 60e-6
    max_long_axis_fraction: float = 0.9  # long axis <= fraction * (thickness - 2 margin)
    max_truncated_mass: float = 0.10  # reject designs whose size distribution is heavily truncated
    loading_tolerance: float = 0.005  # absolute tolerance on achieved area fraction
    max_attempts: int = 3000  # placement attempts per particle
    restarts: int = 3
    check_vertices: int = 64  # polygon used for (conservative) overlap checks

    def to_dict(self) -> dict:
        return asdict(self)


@dataclass(frozen=True)
class MeshSettings:
    n_across_gap: int = 3  # elements across the narrowest resolved gap
    n_perimeter: int = 48  # minimum vertices on a particle boundary
    h_bulk: float = 6e-6  # maximum element size
    h_min: float = 0.15e-6
    face_elements_in_margin: int = 2  # elements across the resin margin at the faces
    size_factor: float = 1.0  # global refinement factor (0.5 halves every size)
    min_quality: float = 0.05  # minimum gmsh SICN quality accepted

    def to_dict(self) -> dict:
        return asdict(self)

    def refined(self, factor: float) -> "MeshSettings":
        return replace(self, size_factor=self.size_factor * factor)


@dataclass(frozen=True)
class Design:
    """The five design descriptors (distribution-level, not per particle)."""

    loading: float  # filler area fraction of the whole layer (2D; NOT a 3D volume fraction)
    d50: float  # median equivalent-area diameter (number basis), m
    sigma_ln: float  # log-normal width of the diameter distribution
    aspect: float  # ellipse major/minor axis ratio (>= 1)
    order: float  # orientation order parameter S = <cos 2(theta - theta0)> in [0, 1)
    preferred_angle: float = 1.5707963267948966  # theta0, rad; pi/2 = along heat flow (y)

    def to_dict(self) -> dict:
        return asdict(self)

    def with_(self, **kw) -> "Design":
        return replace(self, **kw)

    def label(self) -> str:
        return (
            f"phi={self.loading:.2f} d50={self.d50 * 1e6:.1f}um s={self.sigma_ln:.2f} "
            f"AR={self.aspect:.2f} S={self.order:.2f}"
        )


def _lengths(d: dict, keys: tuple[str, ...]) -> dict:
    out = dict(d)
    for k in keys:
        if k in out:
            out[k] = to_si(out[k])
    return out


def layer_from_dict(d: dict) -> Layer:
    return Layer(**_lengths(d, ("thickness", "width", "margin")))


def geometry_from_dict(d: dict) -> GeometrySettings:
    return GeometrySettings(**_lengths(d, ("min_gap", "d_min", "d_max")))


def mesh_from_dict(d: dict) -> MeshSettings:
    return MeshSettings(**_lengths(d, ("h_bulk", "h_min")))


def design_from_dict(d: dict) -> Design:
    d = dict(d)
    d["d50"] = to_si(d["d50"])
    if "preferred_angle" in d:
        d["preferred_angle"] = to_si(d["preferred_angle"], "rad")
    return Design(**d)


@dataclass
class StudyConfig:
    name: str
    material_path: str
    scenarios_path: str
    layer: Layer
    geometry: GeometrySettings
    mesh: MeshSettings
    design_space: dict
    budget: dict
    raw: dict = field(default_factory=dict)
    root: Path = Path(".")

    def resolve(self, p: str) -> Path:
        q = Path(p)
        return q if q.is_absolute() else (self.root / q)


def load_study(path: str | Path) -> StudyConfig:
    path = Path(path)
    with open(path) as fh:
        raw = yaml.safe_load(fh)
    root = path.resolve().parent.parent if path.resolve().parent.name == "configs" else path.resolve().parent
    return StudyConfig(
        name=raw["name"],
        material_path=raw["material"],
        scenarios_path=raw["scenarios"],
        layer=layer_from_dict(raw.get("layer", {})),
        geometry=geometry_from_dict(raw.get("geometry", {})),
        mesh=mesh_from_dict(raw.get("mesh", {})),
        design_space=raw["design_space"],
        budget=raw.get("budget", {}),
        raw=raw,
        root=root,
    )
