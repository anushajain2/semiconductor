"""Constrained 2D microstructure generator.

Particles are non-overlapping ellipses (circles when aspect = 1) placed inside a
finite layer with a resin-rich margin on every boundary. Requests that cannot be
packed are rejected with an explicit reason; the achieved loading is recorded and
never silently changed.

Distribution-level design descriptors (see :class:`tcid.config.Design`):

* equivalent-area diameter ~ log-normal(median ``d50``, width ``sigma_ln``),
  truncated to ``[d_min, d_max_eff]`` where ``d_max_eff`` keeps the long axis
  within the layer;
* every particle has the same aspect ratio;
* orientations follow an axial von Mises distribution around the preferred
  direction with order parameter ``S = <cos 2(theta - theta0)>`` (S = 0 is random,
  S -> 1 is perfectly aligned).

A 2D area fraction is *not* an experimentally established 3D volume fraction.
"""

from __future__ import annotations

import hashlib
import json
import math
from dataclasses import dataclass, field

import numpy as np
import shapely
from scipy import optimize, special, stats

from .config import Design, GeometrySettings, Layer

GENERATOR_VERSION = "geomgen-1.0"


# ----------------------------------------------------------------------------- orientation


def kappa_from_order(order: float) -> float:
    """Concentration of the doubled-angle von Mises with mean resultant ``order``.

    For theta' = 2 (theta - theta0) ~ vonMises(0, kappa), S = I1(kappa) / I0(kappa).
    """
    if order <= 0.0:
        return 0.0
    if order >= 0.999:
        raise ValueError("orientation order must be < 0.999")
    f = lambda k: special.i1e(k) / special.i0e(k) - order  # noqa: E731
    return optimize.brentq(f, 1e-9, 2000.0)


def sample_orientations(rng: np.random.Generator, n: int, order: float, theta0: float) -> np.ndarray:
    """Major-axis angles in [0, pi)."""
    if order <= 0.0:
        th = rng.uniform(0.0, math.pi, n)
    else:
        k = kappa_from_order(order)
        th = theta0 + 0.5 * rng.vonmises(0.0, k, n)
    return np.mod(th, math.pi)


def realized_order(theta: np.ndarray, theta0: float) -> float:
    if theta.size == 0:
        return float("nan")
    return float(np.mean(np.cos(2.0 * (theta - theta0))))


# ----------------------------------------------------------------------------- shapes


def semi_axes(d_eq: np.ndarray, aspect: float) -> tuple[np.ndarray, np.ndarray]:
    """Semi-axes (a >= b) of an ellipse with equivalent-area diameter ``d_eq``."""
    r = 0.5 * np.asarray(d_eq)
    s = math.sqrt(aspect)
    return r * s, r / s


def half_extents(a, b, theta) -> tuple[np.ndarray, np.ndarray]:
    c, s = np.cos(theta), np.sin(theta)
    ex = np.sqrt((a * c) ** 2 + (b * s) ** 2)
    ey = np.sqrt((a * s) ** 2 + (b * c) ** 2)
    return ex, ey


def ellipse_perimeter(a, b):
    """Ramanujan's second approximation."""
    h = ((a - b) / (a + b)) ** 2
    return math.pi * (a + b) * (1 + 3 * h / (10 + np.sqrt(4 - 3 * h)))


def ellipse_vertices(cx, cy, a, b, theta, n: int, scale: float = 1.0) -> np.ndarray:
    """``n`` vertices on the ellipse boundary, uniformly spaced in the parametric angle.

    Parametric spacing concentrates vertices at the high-curvature tips.
    ``scale > 1`` produces a circumscribing polygon (used for conservative overlap checks).
    """
    t = np.linspace(0.0, 2.0 * math.pi, n, endpoint=False)
    x = scale * a * np.cos(t)
    y = scale * b * np.sin(t)
    c, s = math.cos(theta), math.sin(theta)
    return np.column_stack([cx + c * x - s * y, cy + s * x + c * y])


def circumscribed_scale(n: int) -> float:
    """Scale of an n-gon (vertices uniform in parametric angle) that contains the ellipse.

    The affine image of a circle-circumscribing regular n-gon contains the ellipse.
    """
    return 1.0 / math.cos(math.pi / n)


# ----------------------------------------------------------------------------- result


@dataclass
class Geometry:
    design: dict
    layer: dict
    settings: dict
    seed: int
    generator_version: str
    status: str  # "ok" or "rejected"
    reason: str = ""
    cx: np.ndarray = field(default_factory=lambda: np.zeros(0))
    cy: np.ndarray = field(default_factory=lambda: np.zeros(0))
    a: np.ndarray = field(default_factory=lambda: np.zeros(0))
    b: np.ndarray = field(default_factory=lambda: np.zeros(0))
    theta: np.ndarray = field(default_factory=lambda: np.zeros(0))
    stats: dict = field(default_factory=dict)

    @property
    def ok(self) -> bool:
        return self.status == "ok"

    @property
    def n(self) -> int:
        return int(self.cx.size)

    def key(self) -> str:
        payload = {
            "design": self.design,
            "layer": self.layer,
            "settings": self.settings,
            "seed": self.seed,
            "gen": self.generator_version,
        }
        return hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest()[:16]

    def arrays(self) -> dict:
        return {"cx": self.cx, "cy": self.cy, "a": self.a, "b": self.b, "theta": self.theta}

    def rotated90(self) -> "Geometry":
        """The same particles rotated by +90 degrees about the origin then shifted into x >= 0.

        Used by the orientation-implementation check: layer thickness becomes the x extent.
        """
        t = self.layer["thickness"]
        g = Geometry(**{**self.__dict__})
        g.cx, g.cy = t - self.cy, self.cx.copy()
        g.theta = np.mod(self.theta + math.pi / 2, math.pi)
        return g


def _design_seed(design: Design, seed: int) -> np.random.SeedSequence:
    h = hashlib.sha256(json.dumps(design.to_dict(), sort_keys=True).encode()).hexdigest()
    return np.random.SeedSequence([int(seed), int(h[:8], 16)])


def max_equivalent_diameter(design: Design, layer: Layer, gs: GeometrySettings) -> float:
    avail = gs.max_long_axis_fraction * (min(layer.thickness, layer.width) - 2 * layer.margin)
    # long axis = d * sqrt(aspect)
    return min(gs.d_max, avail / math.sqrt(design.aspect))


def truncated_mass(design: Design, layer: Layer, gs: GeometrySettings) -> float:
    dist = stats.lognorm(s=max(design.sigma_ln, 1e-12), scale=design.d50)
    dmax = max_equivalent_diameter(design, layer, gs)
    return float(dist.cdf(gs.d_min) + dist.sf(dmax))


def _sample_diameters(rng, design: Design, dmin: float, dmax: float, n: int) -> np.ndarray:
    lo = (math.log(dmin) - math.log(design.d50)) / max(design.sigma_ln, 1e-12)
    hi = (math.log(dmax) - math.log(design.d50)) / max(design.sigma_ln, 1e-12)
    if design.sigma_ln <= 1e-9:
        return np.full(n, np.clip(design.d50, dmin, dmax))
    z = stats.truncnorm.rvs(lo, hi, size=n, random_state=rng)
    return np.exp(math.log(design.d50) + design.sigma_ln * z)


def _sample_population(rng, design: Design, layer: Layer, gs: GeometrySettings, tries: int = 50):
    """Draw a particle population whose total area matches the target within tolerance.

    Populations are redrawn (rejection sampling on the total area) when the target cannot
    be matched; a narrow distribution of large particles in a small domain may make the
    requested loading unattainable, in which case the request is rejected.
    """
    reason = ""
    for _ in range(tries):
        ds, reason = _sample_population_once(rng, design, layer, gs)
        if ds is not None or "attainable" not in reason:
            return ds, reason
    return None, reason


def _sample_population_once(rng, design: Design, layer: Layer, gs: GeometrySettings):
    """Draw diameters until the target area is matched within tolerance."""
    target = design.loading * layer.width * layer.thickness
    # sample against a tighter tolerance: the meshed (inscribed) polygons lose up to ~0.3%
    # of the ellipse area, and the polygon loading is checked against the full tolerance
    tol = 0.6 * gs.loading_tolerance * layer.width * layer.thickness
    dmax = max_equivalent_diameter(design, layer, gs)
    if dmax <= gs.d_min:
        return None, "particles cannot fit in the layer at this aspect ratio"
    ds: list[float] = []
    area = 0.0
    while area < target - tol:
        d = float(_sample_diameters(rng, design, gs.d_min, dmax, 1)[0])
        ds.append(d)
        area += math.pi * d * d / 4
        if len(ds) > 200000:
            return None, "too many particles requested"
    # redraw the last particle until total area is within tolerance
    for _ in range(500):
        if abs(area - target) <= tol:
            return np.array(ds), ""
        if not ds:
            break
        area -= math.pi * ds[-1] ** 2 / 4
        if abs(area - target) <= tol:
            ds.pop()
            return np.array(ds), ""
        d = float(_sample_diameters(rng, design, gs.d_min, dmax, 1)[0])
        ds[-1] = d
        area += math.pi * d * d / 4
    return None, "requested loading not attainable with discrete particles from this size distribution"


def _place(rng, d_eq, design: Design, layer: Layer, gs: GeometrySettings):
    a_all, b_all = semi_axes(d_eq, design.aspect)
    th_all = sample_orientations(rng, d_eq.size, design.order, design.preferred_angle)
    order = np.argsort(-d_eq)  # largest first
    a_all, b_all, th_all = a_all[order], b_all[order], th_all[order]
    n = d_eq.size
    cx = np.empty(n)
    cy = np.empty(n)
    scale = circumscribed_scale(gs.check_vertices)
    polys: list = []
    g = gs.min_gap
    m = layer.margin
    for i in range(n):
        a, b, th = a_all[i], b_all[i], th_all[i]
        ex, ey = half_extents(a, b, th)
        xlo, xhi = m + ex, layer.width - m - ex
        ylo, yhi = m + ey, layer.thickness - m - ey
        if xlo > xhi or ylo > yhi:
            return None, f"particle {i} (long axis {2 * a * 1e6:.2f} um) does not fit inside the margins"
        placed = False
        for _ in range(gs.max_attempts):
            x = rng.uniform(xlo, xhi)
            y = rng.uniform(ylo, yhi)
            if i > 0:
                dist = np.hypot(cx[:i] - x, cy[:i] - y)
                if np.any(dist < b + b_all[:i] + g):
                    continue  # certain overlap (inscribed circles)
                amb = np.nonzero(dist < scale * (a + a_all[:i]) + g)[0]
                if amb.size:
                    cand = shapely.Polygon(ellipse_vertices(x, y, a, b, th, gs.check_vertices, scale))
                    dd = shapely.distance(cand, [polys[j] for j in amb])
                    if np.any(dd < g):
                        continue
            cx[i], cy[i] = x, y
            polys.append(shapely.Polygon(ellipse_vertices(x, y, a, b, th, gs.check_vertices, scale)))
            placed = True
            break
        if not placed:
            return None, f"packing infeasible: particle {i + 1}/{n} could not be placed"
    return (cx, cy, a_all, b_all, th_all), ""


def pairwise_gaps(
    cx, cy, a, b, theta, cutoff: float, n_vertices: int = 128, conservative: bool = False
) -> np.ndarray:
    """Per-particle distance to the nearest neighbouring particle surface (capped at cutoff).

    With ``conservative=True`` circumscribing polygons are used, giving a lower bound on
    the true ellipse-ellipse distance; otherwise inscribed polygons (an upper bound that
    converges to the true distance as ``n_vertices`` grows).
    """
    n = cx.size
    gaps = np.full(n, cutoff)
    if n < 2:
        return gaps
    sc = circumscribed_scale(n_vertices) if conservative else 1.0
    polys = [
        shapely.Polygon(ellipse_vertices(cx[i], cy[i], a[i], b[i], theta[i], n_vertices, sc)) for i in range(n)
    ]
    for i in range(n):
        dist = np.hypot(cx - cx[i], cy - cy[i])
        near = np.nonzero((dist < a + a[i] + cutoff) & (np.arange(n) != i))[0]
        if near.size:
            d = shapely.distance(polys[i], [polys[j] for j in near])
            gaps[i] = min(cutoff, float(np.min(d)))
    return gaps


def generate(design: Design, layer: Layer, gs: GeometrySettings, seed: int) -> Geometry:
    """Generate one feasible realization or a rejected record with the reason."""
    base = Geometry(
        design=design.to_dict(),
        layer=layer.to_dict(),
        settings=gs.to_dict(),
        seed=int(seed),
        generator_version=GENERATOR_VERSION,
        status="rejected",
    )
    if not (0.0 < design.loading < 1.0):
        base.reason = "loading must be in (0, 1)"
        return base
    if design.aspect < 1.0:
        base.reason = "aspect ratio must be >= 1"
        return base
    if not (0.0 <= design.order < 0.999):
        base.reason = "orientation order must be in [0, 0.999)"
        return base
    tm = truncated_mass(design, layer, gs)
    base.stats["truncated_mass"] = tm
    if tm > gs.max_truncated_mass:
        base.reason = (
            f"size distribution truncated by {tm:.1%} (> {gs.max_truncated_mass:.0%}) to fit the layer; "
            "requested distribution is incompatible with the layer dimensions"
        )
        return base

    ss = _design_seed(design, seed)
    last_reason = ""
    for attempt, child in enumerate(ss.spawn(gs.restarts)):
        rng = np.random.default_rng(child)
        d_eq, reason = _sample_population(rng, design, layer, gs)
        if d_eq is None:
            base.reason = reason
            return base
        placed, reason = _place(rng, d_eq, design, layer, gs)
        if placed is None:
            last_reason = reason
            continue
        cx, cy, a, b, th = placed
        g = Geometry(**{**base.__dict__})
        g.status = "ok"
        g.reason = ""
        g.cx, g.cy, g.a, g.b, g.theta = cx, cy, a, b, th
        g.stats = {**base.stats, **describe(g), "placement_restarts": attempt}
        return g
    base.reason = last_reason + f" after {gs.restarts} restarts"
    return base


def describe(g: Geometry) -> dict:
    """Particle statistics recorded with every run."""
    lay = g.layer
    d_eq = 2.0 * np.sqrt(g.a * g.b)
    area = math.pi * g.a * g.b
    w = np.sort(d_eq)
    cum = np.cumsum(np.pi * w**2 / 4)
    d50_area = float(w[np.searchsorted(cum, 0.5 * cum[-1])]) if w.size else float("nan")
    return {
        "n_particles": int(g.n),
        "loading_ellipse": float(area.sum() / (lay["width"] * lay["thickness"])),
        "d10": float(np.percentile(d_eq, 10)) if d_eq.size else float("nan"),
        "d50_realized": float(np.median(d_eq)) if d_eq.size else float("nan"),
        "d90": float(np.percentile(d_eq, 90)) if d_eq.size else float("nan"),
        "d50_area_weighted": d50_area,
        "sigma_ln_realized": float(np.std(np.log(d_eq))) if d_eq.size > 1 else 0.0,
        "order_realized": realized_order(g.theta, g.design["preferred_angle"]),
        "long_axis_max": float(2 * g.a.max()) if g.n else 0.0,
    }


def validate(g: Geometry, n_vertices: int = 256) -> dict:
    """Independent validity check: no overlap below min gap and all particles inside margins."""
    lay, gs = g.layer, g.settings
    ex, ey = half_extents(g.a, g.b, g.theta)
    tol = 1e-12
    inside = (
        (g.cx - ex >= lay["margin"] - tol)
        & (g.cx + ex <= lay["width"] - lay["margin"] + tol)
        & (g.cy - ey >= lay["margin"] - tol)
        & (g.cy + ey <= lay["thickness"] - lay["margin"] + tol)
    )
    # circumscribing polygons give a lower bound on the true surface separation
    gaps = pairwise_gaps(
        g.cx, g.cy, g.a, g.b, g.theta, cutoff=10 * gs["min_gap"], n_vertices=n_vertices, conservative=True
    )
    min_gap = float(gaps.min()) if g.n > 1 else float("inf")
    # polygon error of a 256-gon is < 1e-4 of the semi-axis; allow 0.5% of the gap
    gap_ok = bool(min_gap >= gs["min_gap"] * (1 - 5e-3))
    return {
        "all_inside_margins": bool(inside.all()),
        "min_gap_lower_bound": min_gap,
        "gap_ok": gap_ok,
        "valid": bool(inside.all() and gap_ok),
    }
