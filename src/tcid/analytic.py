"""Closed-form references used by the verification gates."""

from __future__ import annotations

import numpy as np


def series_resistance(thicknesses, conductivities, interface_resistances=()) -> float:
    """Area-specific resistance of a 1D stack: sum t_i/k_i + sum R_j."""
    return float(np.sum(np.asarray(thicknesses) / np.asarray(conductivities)) + np.sum(interface_resistances))


def parallel_resistance(widths, conductivities, thickness: float) -> float:
    """Strips parallel to the heat flow: conductances add (interfaces carry no normal flux)."""
    w = np.asarray(widths, float)
    k = np.asarray(conductivities, float)
    return float(thickness * w.sum() / np.sum(k * w))


def wiener_bounds(phi: float, k_m: float, k_f: float) -> tuple[float, float]:
    """Harmonic (lower) and arithmetic (upper) bounds on effective conductivity.

    Rigorous for any two-phase microstructure between isothermal faces with insulated
    sidewalls and perfect contact; not applicable with interface resistance.
    """
    lo = 1.0 / (phi / k_f + (1.0 - phi) / k_m)
    hi = phi * k_f + (1.0 - phi) * k_m
    return lo, hi


def hashin_shtrikman_2d(phi: float, k_m: float, k_f: float) -> tuple[float, float]:
    """2D Hashin-Shtrikman bounds (statistically isotropic, infinite medium, perfect contact).

    Informational only for a finite layer with resin-rich margins and anisotropic
    arrangements; such layers may legitimately fall outside them.
    """
    k1, k2 = min(k_m, k_f), max(k_m, k_f)
    phi2 = phi if k_f >= k_m else 1.0 - phi
    phi1 = 1.0 - phi2
    lo = k1 + phi2 / (1.0 / (k2 - k1) + phi1 / (2.0 * k1)) if k2 > k1 else k1
    hi = k2 + phi1 / (1.0 / (k1 - k2) + phi2 / (2.0 * k2)) if k2 > k1 else k2
    return lo, hi
