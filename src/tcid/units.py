"""Unit handling.

All internal quantities are SI: metres, W/(m K), m^2 K/W, kelvin, radians.
Configuration files may give quantities either as plain numbers (taken as SI)
or as strings such as ``"100 um"`` / ``"3e-6 m2K/W"`` / ``"0.05 cm2K/W"``.
"""

from __future__ import annotations

import math
import re

# factor converting the given unit to SI
_UNITS: dict[str, float] = {
    # length
    "m": 1.0,
    "mm": 1e-3,
    "um": 1e-6,
    "µm": 1e-6,
    "micron": 1e-6,
    "nm": 1e-9,
    # conductivity
    "W/mK": 1.0,
    "W/(m K)": 1.0,
    "W/m/K": 1.0,
    # area-specific resistance
    "m2K/W": 1.0,
    "m^2K/W": 1.0,
    "m2 K/W": 1.0,
    "K m2/W": 1.0,
    "mm2K/W": 1e-6,
    "K mm2/W": 1e-6,
    "cm2K/W": 1e-4,
    "K cm2/W": 1e-4,
    # angle
    "rad": 1.0,
    "deg": math.pi / 180.0,
    # dimensionless
    "": 1.0,
    "-": 1.0,
    "1": 1.0,
    # density
    "kg/m3": 1.0,
    "g/cm3": 1000.0,
    # temperature difference / absolute
    "K": 1.0,
}

_QTY = re.compile(r"^\s*([-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?)\s*(.*?)\s*$")


def to_si(value, unit: str | None = None) -> float:
    """Convert ``value`` to SI.

    ``value`` may be a number (with optional ``unit``), a string ``"<number> <unit>"``
    or a mapping ``{"value": x, "unit": u}``.
    """
    if isinstance(value, dict):
        return to_si(value["value"], value.get("unit", ""))
    if isinstance(value, (int, float)):
        u = "" if unit is None else unit
        if u not in _UNITS:
            raise ValueError(f"unknown unit {u!r}")
        return float(value) * _UNITS[u]
    if isinstance(value, str):
        m = _QTY.match(value)
        if not m:
            raise ValueError(f"cannot parse quantity {value!r}")
        num, u = float(m.group(1)), m.group(2)
        if unit is not None and u == "":
            u = unit
        if u not in _UNITS:
            raise ValueError(f"unknown unit {u!r} in {value!r}")
        return num * _UNITS[u]
    raise TypeError(f"unsupported quantity type {type(value)!r}")


def from_si(value: float, unit: str) -> float:
    """Express an SI ``value`` in ``unit``."""
    return value / _UNITS[unit]


def weight_to_volume_fraction(w_filler: float, rho_filler: float, rho_matrix: float) -> float:
    """Convert filler weight fraction to volume fraction (void-free two-phase mixture)."""
    if not 0.0 <= w_filler <= 1.0:
        raise ValueError("weight fraction must be in [0, 1]")
    vf = w_filler / rho_filler
    vm = (1.0 - w_filler) / rho_matrix
    return vf / (vf + vm)


def volume_to_weight_fraction(phi_filler: float, rho_filler: float, rho_matrix: float) -> float:
    """Convert filler volume fraction to weight fraction (void-free two-phase mixture)."""
    if not 0.0 <= phi_filler <= 1.0:
        raise ValueError("volume fraction must be in [0, 1]")
    mf = phi_filler * rho_filler
    mm = (1.0 - phi_filler) * rho_matrix
    return mf / (mf + mm)
