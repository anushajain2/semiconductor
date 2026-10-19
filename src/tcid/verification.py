"""Solver acceptance gates (guide section 5).

Each check returns a record ``{name, criterion, value, tolerance, passed, detail}``.
The tolerances are the guide's proposed engineering acceptance criteria; they are
tightened automatically where a design conclusion is of similar magnitude
(see :func:`tcid.studies.final_verification`).
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np

from . import analytic
from .config import Design, GeometrySettings, Layer, MeshSettings
from .design_space import DesignSpace
from .geometry import generate, validate
from .mesh import FILLER, MATRIX, mesh_particles, mesh_strips
from .solver import LayerModel, MaterialParams

KM, KF = 0.2, 30.0


def _rec(name, criterion, value, tol, passed, **detail):
    return {
        "name": name,
        "criterion": criterion,
        "value": float(value),
        "tolerance": float(tol),
        "passed": bool(passed),
        "detail": detail,
    }


def check_homogeneous(tol=0.005) -> dict:
    worst = 0.0
    cases = []
    for t, k, h in [(100e-6, 0.2, 10e-6), (50e-6, 1.3, 4e-6), (150e-6, 30.0, 20e-6)]:
        m = mesh_strips([0.0, t], [MATRIX], 2 * t, h)
        s = LayerModel(m).solve(MaterialParams(k, KF, KF, 0.0))
        err = abs(s.r_layer - t / k) / (t / k)
        worst = max(worst, err)
        cases.append({"t": t, "k": k, "R": s.r_layer, "R_exact": t / k, "rel_err": err})
    return _rec("homogeneous_layer", "R = t/k within 0.5%", worst, tol, worst <= tol, cases=cases)


def check_layered(tol=0.01) -> dict:
    worst = 0.0
    cases = []
    # series stack matrix/filler/matrix/filler/matrix with finite interface resistance
    lv = [0, 20e-6, 35e-6, 60e-6, 70e-6, 100e-6]
    kinds = [MATRIX, FILLER, MATRIX, FILLER, MATRIX]
    m = mesh_strips(lv, kinds, 200e-6, 4e-6)
    lm = LayerModel(m)
    th = np.diff(lv)
    for r in [0.0, 1e-7, 1e-6, 1e-5]:
        k = [KM if kk == MATRIX else KF for kk in kinds]
        exact = analytic.series_resistance(th, k, [r] * 4)
        s = lm.solve(MaterialParams(KM, KF, KF, r))
        err = abs(s.r_layer - exact) / exact
        worst = max(worst, err)
        cases.append({"layout": "series", "r_interface": r, "R": s.r_layer, "R_exact": exact, "rel_err": err})
    # parallel strips: interfaces parallel to the flow carry no normal flux
    xl = [0, 50e-6, 80e-6, 150e-6, 200e-6]
    kinds = [MATRIX, FILLER, MATRIX, FILLER]
    m = mesh_strips(xl, kinds, 100e-6, 5e-6, orientation="parallel")
    lm = LayerModel(m)
    for r in [0.0, 1e-6]:
        k = [KM if kk == MATRIX else KF for kk in kinds]
        exact = analytic.parallel_resistance(np.diff(xl), k, 100e-6)
        s = lm.solve(MaterialParams(KM, KF, KF, r))
        err = abs(s.r_layer - exact) / exact
        worst = max(worst, err)
        cases.append({"layout": "parallel", "r_interface": r, "R": s.r_layer, "R_exact": exact, "rel_err": err})
    # outer contacts as Robin conditions on a series stack
    lv = [0, 30e-6, 60e-6, 100e-6]
    kinds = [MATRIX, FILLER, MATRIX]
    m = mesh_strips(lv, kinds, 200e-6, 5e-6)
    s = LayerModel(m).solve(MaterialParams(KM, KF, KF, 1e-6, r_hot=1e-5, r_cold=3e-5))
    exact = analytic.series_resistance(np.diff(lv), [KM, KF, KM], [1e-6, 1e-6]) + 4e-5
    err = abs(s.r_assembly - exact) / exact
    worst = max(worst, err)
    cases.append({"layout": "series+robin contacts", "R": s.r_assembly, "R_exact": exact, "rel_err": err})
    return _rec("layered_interfaces", "summed layer + interface resistances within 1%", worst, tol, worst <= tol, cases=cases)


def representative_geometry(seed=3, design: Design | None = None, layer: Layer | None = None):
    design = design or Design(loading=0.30, d50=10e-6, sigma_ln=0.3, aspect=2.5, order=0.6)
    layer = layer or Layer()
    g = generate(design, layer, GeometrySettings(), seed)
    if not g.ok:
        raise RuntimeError(g.reason)
    return g


def check_energy_balance(tol=0.001, size_factor=2.0) -> dict:
    g = representative_geometry()
    m = mesh_particles(g, MeshSettings(size_factor=size_factor))
    lm = LayerModel(m)
    worst, cases = 0.0, []
    for r in [0.0, 1e-8, 1e-6, 1e-5]:
        for kp, kq in [(KF, KF), (300.0, 5.0)]:
            s = lm.solve(MaterialParams(KM, kp, kq, r))
            worst = max(worst, s.heat_balance_error)
            cases.append(
                {
                    "r_interface": r,
                    "k_par": kp,
                    "k_perp": kq,
                    "heat_balance_error": s.heat_balance_error,
                    "flux_consistency_error": s.flux_consistency_error,
                    "residual": s.residual,
                }
            )
    return _rec("energy_balance", "|Q_in - Q_out| / Q < 0.1%", worst, tol, worst <= tol, cases=cases)


def check_mesh_refinement(tol=0.02, factors=(1.0, 0.7, 0.5)) -> dict:
    g = representative_geometry()
    res = {}
    for f in factors:
        lm = LayerModel(mesh_particles(g, MeshSettings(size_factor=f)))
        res[f] = {r: lm.solve(MaterialParams(KM, KF, KF, r)).r_layer for r in (1e-8, 1e-6)}
    changes = []
    for a, b in zip(factors[:-1], factors[1:]):
        for r in (1e-8, 1e-6):
            changes.append(abs(res[b][r] - res[a][r]) / res[b][r])
    worst = max(changes)
    return _rec(
        "mesh_refinement",
        "resistance change < 2% on the next refinement",
        worst,
        tol,
        worst <= tol,
        resistances={str(k): v for k, v in res.items()},
    )


def check_orientation(tol_exact=1e-8, tol_remesh=0.01) -> dict:
    """Rotating geometry, tensor and boundary problem together preserves the result.

    (a) rigid rotation of the mesh by an arbitrary angle with an anisotropic filler (exact);
    (b) geometry rotated by 90 degrees and re-meshed with heat flow along x (mesh tolerance);
    (c) sensitivity: rotating geometry but *not* the tensor must change the result.
    """
    g = representative_geometry()
    ms = MeshSettings(size_factor=2.0)
    mesh = mesh_particles(g, ms)
    mp_ = MaterialParams(KM, 300.0, 5.0, 1e-7)
    r0 = LayerModel(mesh).solve(mp_).r_layer
    errs = []
    for ang in (math.radians(37.0), math.radians(90.0), math.radians(-123.0)):
        r1 = LayerModel(mesh.rotated(ang)).solve(mp_).r_layer
        errs.append(abs(r1 - r0) / r0)
    exact_err = max(errs)
    g90 = g.rotated90()
    m90 = mesh_particles(g90, ms, flow_axis="x")
    r90 = LayerModel(m90).solve(mp_).r_layer
    remesh_err = abs(r90 - r0) / r0
    # deliberate error: rotate the geometry but keep the unrotated material frames
    bad = mesh.rotated(math.radians(90.0))
    bad.region_theta = mesh.region_theta.copy()
    r_bad = LayerModel(bad).solve(mp_).r_layer
    bad_change = abs(r_bad - r0) / r0
    passed = exact_err <= tol_exact and remesh_err <= tol_remesh and bad_change > 2 * tol_remesh
    return _rec(
        "orientation_implementation",
        "rigid rotation exact (<1e-8); 90-degree re-meshed rotation within 1%; un-rotated tensor detected",
        max(exact_err, remesh_err),
        tol_remesh,
        passed,
        rigid_rotation_errors=errs,
        remeshed_90deg_error=remesh_err,
        unrotated_tensor_change=bad_change,
    )


def check_geometry_validity(n_designs=16, seed=11) -> dict:
    space = DesignSpace(
        bounds={"d50": (4e-6, 24e-6), "sigma_ln": (0.05, 0.5), "aspect": (1.0, 4.0), "order": (0.0, 0.95)},
        loading=0.30,
        preferred_angle=math.pi / 2,
    )
    layer, gs = Layer(), GeometrySettings()
    stats = {"accepted": 0, "rejected": 0, "invalid_accepted": 0, "reasons": [], "loading_error_max": 0.0}
    min_gap = math.inf
    for x in space.sobol(n_designs, seed):
        d = space.from_unit(x)
        g = generate(d, layer, gs, 0)
        if not g.ok:
            stats["rejected"] += 1
            stats["reasons"].append(g.reason)
            continue
        v = validate(g)
        stats["accepted"] += 1
        if not v["valid"]:
            stats["invalid_accepted"] += 1
        min_gap = min(min_gap, v["min_gap_lower_bound"])
        stats["loading_error_max"] = max(stats["loading_error_max"], abs(g.stats["loading_ellipse"] - d.loading))
    passed = stats["invalid_accepted"] == 0 and stats["loading_error_max"] <= gs.loading_tolerance + 1e-12
    return _rec(
        "geometry_validity",
        "no overlap or margin violation in accepted geometries; achieved loading within tolerance and recorded",
        stats["invalid_accepted"],
        0,
        passed,
        min_gap_lower_bound=min_gap,
        **stats,
    )


def check_limits_and_bounds() -> dict:
    """Pure-phase limits, monotonic response to interface resistance, mixture bounds."""
    g = representative_geometry()
    mesh = mesh_particles(g, MeshSettings(size_factor=2.0))
    lm = LayerModel(mesh)
    t = g.layer["thickness"]
    phi = mesh.filler_fraction()
    out = {}
    # pure phase: filler identical to matrix and perfect contact -> t / k_m
    r_pure = lm.solve(MaterialParams(KM, KM, KM, 0.0)).r_layer
    out["pure_phase_error"] = abs(r_pure - t / KM) / (t / KM)
    # monotonic increase with a positive interface resistance
    rs = [0.0, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-2]
    vals = [lm.solve(MaterialParams(KM, KF, KF, r)).r_layer for r in rs]
    out["monotonic_in_r_interface"] = bool(np.all(np.diff(vals) > 0))
    out["r_vs_r_interface"] = dict(zip(map(str, rs), vals))
    # large interface resistance -> particles act like insulating holes
    r_hole = lm.solve(MaterialParams(KM, 1e-6, 1e-6, 0.0)).r_layer
    out["huge_resistance_vs_holes_rel"] = abs(vals[-1] - r_hole) / r_hole
    # Wiener bounds (rigorous, perfect contact)
    k_eff = t / vals[0]
    lo, hi = analytic.wiener_bounds(phi, KM, KF)
    out["wiener"] = {"k_eff": k_eff, "lower": lo, "upper": hi, "inside": bool(lo <= k_eff <= hi)}
    hlo, hhi = analytic.hashin_shtrikman_2d(phi, KM, KF)
    out["hashin_shtrikman_2d_informational"] = {"lower": hlo, "upper": hhi, "inside": bool(hlo <= k_eff <= hhi)}
    passed = (
        out["pure_phase_error"] < 1e-8
        and out["monotonic_in_r_interface"]
        and out["huge_resistance_vs_holes_rel"] < 0.01
        and out["wiener"]["inside"]
    )
    return _rec(
        "limits_and_bounds",
        "pure-phase exact; monotone in R_interface; R->inf matches holes within 1%; inside Wiener bounds",
        out["pure_phase_error"],
        1e-8,
        passed,
        **out,
    )


CHECKS = {
    "homogeneous_layer": check_homogeneous,
    "layered_interfaces": check_layered,
    "energy_balance": check_energy_balance,
    "mesh_refinement": check_mesh_refinement,
    "orientation_implementation": check_orientation,
    "geometry_validity": check_geometry_validity,
    "limits_and_bounds": check_limits_and_bounds,
}


def run_all(out_dir: str | Path | None = None, only: list[str] | None = None) -> list[dict]:
    results = []
    for name, fn in CHECKS.items():
        if only and name not in only:
            continue
        print(f"  verifying {name} ...", flush=True)
        results.append(fn())
        print(f"    {'PASS' if results[-1]['passed'] else 'FAIL'}  value={results[-1]['value']:.3e}", flush=True)
    if out_dir:
        out = Path(out_dir)
        out.mkdir(parents=True, exist_ok=True)
        (out / "verification.json").write_text(json.dumps(results, indent=2, default=float))
        lines = ["| Check | Proposed acceptance criterion | Worst value | Pass |", "|---|---|---|---|"]
        for r in results:
            lines.append(f"| {r['name']} | {r['criterion']} | {r['value']:.2e} | {'yes' if r['passed'] else '**NO**'} |")
        (out / "verification.md").write_text("\n".join(lines) + "\n")
    return results
