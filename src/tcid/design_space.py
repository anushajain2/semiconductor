"""Bounded design space over the distribution descriptors (loading held at fixed levels)."""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from scipy.stats import qmc

from .config import Design
from .units import to_si

# name, log-scaled?
VARIABLES = (("d50", True), ("sigma_ln", False), ("aspect", False), ("order", False))


@dataclass
class DesignSpace:
    bounds: dict  # name -> (lo, hi) in SI
    loading: float
    preferred_angle: float

    @classmethod
    def from_config(cls, ds: dict, loading: float | None = None) -> "DesignSpace":
        b = {}
        for name, _ in VARIABLES:
            lo, hi = ds["bounds"][name]
            b[name] = (to_si(lo), to_si(hi))
        return cls(
            bounds=b,
            loading=float(loading if loading is not None else ds["loading_levels"][0]),
            preferred_angle=to_si(ds.get("preferred_angle", "90 deg"), "rad"),
        )

    @property
    def dim(self) -> int:
        return len(VARIABLES)

    def names(self) -> list[str]:
        return [n for n, _ in VARIABLES]

    def to_unit(self, d: Design) -> np.ndarray:
        x = []
        for name, logv in VARIABLES:
            lo, hi = self.bounds[name]
            v = getattr(d, name)
            x.append((math.log(v) - math.log(lo)) / (math.log(hi) - math.log(lo)) if logv else (v - lo) / (hi - lo))
        return np.array(x)

    def unit_from_columns(self, df) -> np.ndarray:
        """Unit-cube coordinates for a table with columns d50, sigma_ln, aspect, order."""
        cols = []
        for name, logv in VARIABLES:
            lo, hi = self.bounds[name]
            v = np.asarray(df[name], float)
            cols.append((np.log(v) - math.log(lo)) / (math.log(hi) - math.log(lo)) if logv else (v - lo) / (hi - lo))
        return np.column_stack(cols)

    def from_unit(self, x) -> Design:
        x = np.clip(np.asarray(x, float), 0.0, 1.0)
        kw = {}
        for (name, logv), xi in zip(VARIABLES, x):
            lo, hi = self.bounds[name]
            kw[name] = math.exp(math.log(lo) + xi * (math.log(hi) - math.log(lo))) if logv else lo + xi * (hi - lo)
        return Design(loading=self.loading, preferred_angle=self.preferred_angle, **kw)

    def contains(self, d: Design, tol: float = 1e-9) -> bool:
        x = self.to_unit(d)
        return bool(np.all(x >= -tol) and np.all(x <= 1 + tol) and abs(d.loading - self.loading) < 1e-9)

    def sobol(self, n: int, seed: int) -> np.ndarray:
        s = qmc.Sobol(d=self.dim, scramble=True, seed=seed)
        m = int(math.ceil(math.log2(max(n, 2))))
        return s.random_base2(m)[:n]

    def uniform(self, n: int, seed: int) -> np.ndarray:
        return np.random.default_rng(seed).uniform(size=(n, self.dim))

    def grid(self, levels: int = 3) -> np.ndarray:
        g = np.linspace(0.0, 1.0, levels)
        mesh = np.meshgrid(*[g] * self.dim, indexing="ij")
        return np.column_stack([m.ravel() for m in mesh])

    def baseline_of(self, d: Design) -> Design:
        """Conventional reference: same loading, sizes and shape, random orientation."""
        return d.with_(order=0.0)
