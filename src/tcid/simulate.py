"""Traceable simulation runner: design + seed -> geometry -> mesh -> solves over scenarios.

Every solve produces one row recording its full configuration hash, geometry seed,
achieved loading, particle statistics, constituent properties, interface
assumptions, boundary conditions, mesh size, solver residual, heat-balance error,
integrated flux, resistance, elapsed time, code version and status. Geometry and
mesh rejections and solver failures are stored with explicit reasons; they are
excluded from optimisation but their frequency is reported.

Results are cached by the *complete* configuration (not geometry alone).
"""

from __future__ import annotations

import json
import multiprocessing as mp
import os
import time
import traceback
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

from .config import Design, GeometrySettings, Layer, MeshSettings
from .geometry import GENERATOR_VERSION, classify_failure, generate, validate
from .materials import MaterialCase, Scenario, params_for
from .mesh import MESHER_VERSION, MeshError, mesh_particles
from .provenance import code_version, config_hash
from .solver import SOLVER_VERSION, LayerModel

BOUNDARY_CONDITIONS = "isothermal hot (y=0, 1 K) and cold (y=t, 0 K) faces; insulated sidewalls"


@dataclass
class Job:
    design: Design
    seed: int
    layer: Layer
    geometry: GeometrySettings
    mesh: MeshSettings
    material: MaterialCase
    scenarios: list[Scenario]
    tag: str = ""
    contact_model: str = "series"  # "series" (isothermal faces + series contacts) or "robin"
    save_fields: bool = False
    meta: dict = field(default_factory=dict)

    def geometry_config(self) -> dict:
        return {
            "design": self.design.to_dict(),
            "seed": int(self.seed),
            "layer": self.layer.to_dict(),
            "geometry": self.geometry.to_dict(),
            "generator": GENERATOR_VERSION,
        }

    def job_key(self) -> str:
        return config_hash({**self.geometry_config(), "mesh": self.mesh.to_dict(), "mesher": MESHER_VERSION})

    def run_key(self, sc: Scenario) -> str:
        return config_hash(
            {
                **self.geometry_config(),
                "mesh": self.mesh.to_dict(),
                "mesher": MESHER_VERSION,
                "material": self.material.to_dict(),
                "scenario": sc.to_dict(),
                "contact_model": self.contact_model,
                "solver": SOLVER_VERSION,
            }
        )


def _design_cols(d: Design) -> dict:
    return {
        "loading_target": d.loading,
        "d50": d.d50,
        "sigma_ln": d.sigma_ln,
        "aspect": d.aspect,
        "order": d.order,
        "preferred_angle": d.preferred_angle,
    }


def _base_row(job: Job) -> dict:
    return {
        "tag": job.tag,
        "design_id": config_hash(job.design.to_dict(), 12),
        "seed": int(job.seed),
        **_design_cols(job.design),
        "thickness": job.layer.thickness,
        "width": job.layer.width,
        "margin": job.layer.margin,
        "min_gap_setting": job.geometry.min_gap,
        "size_factor": job.mesh.size_factor,
        "h_bulk": job.mesh.h_bulk,
        "n_across_gap": job.mesh.n_across_gap,
        "material": job.material.name,
        "k_matrix_nominal": job.material.k_matrix,
        "generator_version": GENERATOR_VERSION,
        "mesher_version": MESHER_VERSION,
        "solver_version": SOLVER_VERSION,
        "code_version": code_version(),
        "geometry_key": config_hash(job.geometry_config()),
        "job_key": job.job_key(),
        **{f"meta_{k}": v for k, v in job.meta.items()},
    }


def _failure(job: Job, stage: str, reason: str, extra: dict | None = None) -> dict:
    return {**_base_row(job), "status": "failed", "stage": stage, "reason": reason,
            "failure_kind": classify_failure(stage, reason), **(extra or {})}


def run_job(job: Job, cache_dir: str | None = None, fields_dir: str | None = None) -> tuple[list[dict], list[dict]]:
    """Run one geometry realisation over all scenarios. Returns (rows, failures)."""
    cache = Path(cache_dir) if cache_dir else None
    rows: list[dict] = []
    todo: list[Scenario] = []
    if cache is not None:
        fcache = cache / "failures" / f"{job.job_key()}.json"
        if fcache.exists():
            f = json.loads(fcache.read_text())
            f["tag"] = job.tag
            f.setdefault("failure_kind", classify_failure(f.get("stage", ""), f.get("reason", "")))
            return [], [f]
        for sc in job.scenarios:
            p = cache / "runs" / f"{job.run_key(sc)}.json"
            if p.exists() and not (job.save_fields and fields_dir and not (Path(fields_dir) / f"{job.run_key(sc)}.npz").exists()):
                r = json.loads(p.read_text())
                r["tag"] = job.tag
                r.update({f"meta_{k}": v for k, v in job.meta.items()})
                r["cached"] = True
                rows.append(r)
            else:
                todo.append(sc)
    else:
        todo = list(job.scenarios)
    if not todo:
        return rows, []

    t_job = time.perf_counter()
    failures: list[dict] = []

    def fail(stage, reason, extra=None):
        # CPU spent on rejected work counts toward a campaign's cost
        f = _failure(job, stage, reason, {**(extra or {}), "elapsed_s": time.perf_counter() - t_job})
        if cache is not None:
            (cache / "failures").mkdir(parents=True, exist_ok=True)
            (cache / "failures" / f"{job.job_key()}.json").write_text(json.dumps(f, default=float))
        return rows, [f]

    t0 = time.perf_counter()
    g = generate(job.design, job.layer, job.geometry, job.seed)
    geo_time = time.perf_counter() - t0
    if not g.ok:
        return fail("geometry", g.reason, {"truncated_mass": g.stats.get("truncated_mass", np.nan)})
    v = validate(g)
    if not v["valid"]:
        return fail("geometry_validation", f"independent validity check failed: {v}")
    try:
        m = mesh_particles(g, job.mesh)
    except MeshError as exc:
        return fail("mesh", str(exc))
    if not m.info["mesh_ok"]:
        return fail("mesh_quality", m.info["mesh_reason"], {k: m.info[k] for k in ("bridge_elements", "margin_bridges", "min_quality")})
    if abs(m.info["loading_polygon"] - job.design.loading) > job.geometry.loading_tolerance:
        return fail(
            "loading",
            f"simulated (polygon) loading {m.info['loading_polygon']:.4f} differs from target "
            f"{job.design.loading:.4f} by more than {job.geometry.loading_tolerance}",
        )
    lm = LayerModel(m)

    geo_cols = {
        **{f"geo_{k}": val for k, val in g.stats.items()},
        "loading_achieved": m.info["loading_polygon"],
        "min_gap_meshed": m.info["min_gap_meshed"],
        "min_gap_lower_bound": v["min_gap_lower_bound"],
        "n_elements": m.info["n_elements"],
        "n_nodes": m.info["n_nodes"],
        "min_quality": m.info["min_quality"],
        "bridge_elements": m.info["bridge_elements"],
        "h_particle_min": m.info["h_particle_min"],
        "geometry_time_s": geo_time,
        "mesh_time_s": m.info["mesh_time_s"],
    }
    for sc in todo:
        mp_ = params_for(job.material, sc)
        robin = job.contact_model == "robin" and (sc.r_top > 0 or sc.r_bottom > 0)
        if robin:
            mp_.r_hot, mp_.r_cold = sc.r_bottom, sc.r_top
        try:
            s = lm.solve(mp_)
        except Exception as exc:  # record and continue with other scenarios
            failures.append(_failure(job, "solve", f"{type(exc).__name__}: {exc}", {"scenario": sc.name}))
            continue
        r_total_series = sc.r_top + s.r_layer + sc.r_bottom
        row = {
            **_base_row(job),
            **geo_cols,
            "scenario": sc.name,
            "scenario_role": sc.role,
            "k_matrix": mp_.k_matrix,
            "k_filler_par": mp_.k_filler_par,
            "k_filler_perp": mp_.k_filler_perp,
            "r_interface": sc.r_interface,
            "interface_model": "duplicated DOFs + interface conductance" if sc.r_interface > 0 else "tied DOFs (perfect contact)",
            "r_top": sc.r_top,
            "r_bottom": sc.r_bottom,
            "contact_model": "robin" if robin else "series",
            "boundary_conditions": BOUNDARY_CONDITIONS,
            "n_dof": s.n_dof,
            "solver_residual": s.residual,
            "heat_balance_error": s.heat_balance_error,
            "flux_consistency_error": s.flux_consistency_error,
            "q_hot": s.q_hot,
            "q_cold": s.q_cold,
            "r_layer": s.r_layer,
            "k_eff": job.layer.thickness / s.r_layer,
            "r_total": s.r_assembly if robin else r_total_series,
            "r_total_series": r_total_series,
            "solve_time_s": s.solve_time_s,
            "status": "ok",
            "cached": False,
        }
        row["run_key"] = job.run_key(sc)
        rows.append(row)
        if job.save_fields and fields_dir:
            Path(fields_dir).mkdir(parents=True, exist_ok=True)
            d = lm.disc(mp_.r_interface > 0)
            np.savez_compressed(
                Path(fields_dir) / f"{row['run_key']}.npz",
                T=s.T,
                points=d.basis.mesh.p.T,
                triangles=d.dof_of_elem,
                region=m.region,
                region_kind=m.region_kind,
                cx=g.cx, cy=g.cy, a=g.a, b=g.b, theta=g.theta,
                r_layer=s.r_layer,
                width=job.layer.width, thickness=job.layer.thickness,
            )
    wall = time.perf_counter() - t_job
    for r in rows:
        if not r.get("cached"):
            r["job_time_s"] = wall
            if cache is not None:
                (cache / "runs").mkdir(parents=True, exist_ok=True)
                (cache / "runs" / f"{r['run_key']}.json").write_text(json.dumps(r, default=float))
    return rows, failures


def _run_job_safe(args):
    job, cache_dir, fields_dir = args
    try:
        return run_job(job, cache_dir, fields_dir)
    except Exception as exc:  # never lose the whole batch
        return [], [_failure(job, "internal", f"{type(exc).__name__}: {exc}\n{traceback.format_exc(limit=3)}")]


def run_jobs(
    jobs: list[Job],
    workers: int = 1,
    cache_dir: str | os.PathLike | None = "cache",
    fields_dir: str | os.PathLike | None = None,
    progress: bool = True,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    # one thread per worker process: parallelism comes from the process pool
    for var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ.setdefault(var, "1")
    args = [(j, str(cache_dir) if cache_dir else None, str(fields_dir) if fields_dir else None) for j in jobs]
    rows, fails = [], []
    t0 = time.perf_counter()
    if workers <= 1 or len(jobs) <= 1:
        results = map(_run_job_safe, args)
    else:
        ctx = mp.get_context("spawn")
        ex = ProcessPoolExecutor(max_workers=workers, mp_context=ctx, max_tasks_per_child=8)
        results = ex.map(_run_job_safe, args)
    for i, (r, f) in enumerate(results):
        rows.extend(r)
        fails.extend(f)
        if progress and (i + 1) % max(1, len(jobs) // 10) == 0:
            print(f"  [{i + 1}/{len(jobs)} jobs, {time.perf_counter() - t0:.0f}s]", flush=True)
    if workers > 1 and len(jobs) > 1:
        ex.shutdown()
    return pd.DataFrame(rows), pd.DataFrame(fails)


def make_jobs(designs, seeds, base: dict, tag: str, meta_fn=None) -> list[Job]:
    jobs = []
    for i, d in enumerate(designs):
        for s in seeds:
            meta = meta_fn(i, d) if meta_fn else {}
            jobs.append(Job(design=d, seed=s, tag=tag, meta=meta, **base))
    return jobs
