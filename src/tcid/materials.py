"""Constituent properties, conductivity tensors and physical-input scenarios.

The matrix is isotropic. Each filler particle carries a 2x2 conductivity tensor
``diag(k_parallel, k_perp)`` expressed in the particle frame (major axis first)
and rotated into global coordinates with the particle orientation:

    K_global = Q(theta) diag(k_par, k_perp) Q(theta)^T

Rotating the particle geometry without rotating the tensor would be a physical
error, so the rotation is performed in exactly one place (:func:`rotated_tensor`).
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from pathlib import Path

import numpy as np
import yaml

from .units import to_si


def rotation(theta: float | np.ndarray) -> np.ndarray:
    """2D rotation matrix (or stack of them, shape (..., 2, 2))."""
    c, s = np.cos(theta), np.sin(theta)
    return np.stack([np.stack([c, -s], -1), np.stack([s, c], -1)], -2)


def rotated_tensor(k_par: float, k_perp: float, theta: float | np.ndarray) -> np.ndarray:
    """Global conductivity tensor of a particle whose major axis is at ``theta``."""
    q = rotation(theta)
    d = np.array([[k_par, 0.0], [0.0, k_perp]])
    return q @ d @ np.swapaxes(q, -1, -2)


@dataclass(frozen=True)
class MaterialCase:
    """A coherent material case: one matrix, one filler, density for fraction conversion."""

    name: str
    status: str
    k_matrix: float  # W/(m K)
    k_filler_par: float  # W/(m K), along particle major axis
    k_filler_perp: float  # W/(m K), across particle major axis
    rho_matrix: float  # kg/m^3
    rho_filler: float  # kg/m^3
    sources: dict = field(default_factory=dict)

    @property
    def filler_isotropic(self) -> bool:
        return self.k_filler_par == self.k_filler_perp

    def to_dict(self) -> dict:
        return asdict(self)


@dataclass(frozen=True)
class Scenario:
    """A physical-input stress scenario (not a probability-weighted sample)."""

    name: str
    r_interface: float  # m^2 K/W, particle-matrix area-specific resistance (0 = perfect contact)
    k_filler_scale: float = 1.0  # multiplies both filler conductivities of the material case
    k_matrix_scale: float = 1.0
    r_top: float = 0.0  # outer contact resistance, m^2 K/W
    r_bottom: float = 0.0
    role: str = "stress"  # "stress", "nominal" or "boundary" (outside supported range)

    def to_dict(self) -> dict:
        return asdict(self)


def material_from_dict(d: dict) -> MaterialCase:
    m, f = d["matrix"], d["filler"]
    k_par = to_si(f.get("k_parallel", f.get("k")))
    k_perp = to_si(f.get("k_perp", f.get("k")))
    return MaterialCase(
        name=d["name"],
        status=d.get("status", ""),
        k_matrix=to_si(m["k"]),
        k_filler_par=k_par,
        k_filler_perp=k_perp,
        rho_matrix=to_si(m.get("density", "1.2 g/cm3")),
        rho_filler=to_si(f.get("density", "3.95 g/cm3")),
        sources={"matrix": m.get("source", ""), "filler": f.get("source", "")},
    )


def load_material(path: str | Path) -> MaterialCase:
    with open(path) as fh:
        return material_from_dict(yaml.safe_load(fh))


def scenario_from_dict(d: dict, material: MaterialCase | None = None) -> Scenario:
    """Parse a scenario; ``k_filler``/``k_matrix`` (absolute) are converted to scale factors."""
    kfs = float(d.get("k_filler_scale", 1.0))
    kms = float(d.get("k_matrix_scale", 1.0))
    if "k_filler" in d:
        if material is None:
            raise ValueError("absolute k_filler in a scenario requires the material case")
        kfs = to_si(d["k_filler"]) / material.k_filler_par
    if "k_matrix" in d:
        if material is None:
            raise ValueError("absolute k_matrix in a scenario requires the material case")
        kms = to_si(d["k_matrix"]) / material.k_matrix
    return Scenario(
        name=d["name"],
        r_interface=to_si(d.get("r_interface", 0.0)),
        k_filler_scale=kfs,
        k_matrix_scale=kms,
        r_top=to_si(d.get("r_top", 0.0)),
        r_bottom=to_si(d.get("r_bottom", 0.0)),
        role=d.get("role", "stress"),
    )


def load_scenarios(path: str | Path, material: MaterialCase | None = None) -> list[Scenario]:
    with open(path) as fh:
        data = yaml.safe_load(fh)
    return [scenario_from_dict(s, material) for s in data["scenarios"]]


def params_for(material: MaterialCase, scenario: Scenario):
    """Solver material parameters for a scenario (outer contacts handled by the caller)."""
    from .solver import MaterialParams

    return MaterialParams(
        k_matrix=material.k_matrix * scenario.k_matrix_scale,
        k_filler_par=material.k_filler_par * scenario.k_filler_scale,
        k_filler_perp=material.k_filler_perp * scenario.k_filler_scale,
        r_interface=scenario.r_interface,
    )
