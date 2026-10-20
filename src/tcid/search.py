"""Closed design loop: robust Bayesian optimisation and matched-budget baselines.

Objective (to minimise) for a design x over the finite stress-scenario set S:

    J_robust(x) = max_{s in S} E_seeds[ R_layer(x, s) ]       (worst scenario mean)
    J_nominal(x) = E_seeds[ R_layer(x, nominal) ]

The scenarios are stress cases, not probabilities: J_robust is a worst-case statement,
never a reliability estimate. Outer contact resistance is added after the search (Q2):
an identical series contact cannot change the ordering of designs.
"""

from __future__ import annotations

import time

import numpy as np
import pandas as pd

from .design_space import DesignSpace
from .surrogate import GPSurrogate, aggregate, features, physical_features


def design_objectives(rows: pd.DataFrame, n_stress: int | None = None) -> pd.DataFrame:
    """Per-design robust/nominal objectives from solver rows (designs with failures excluded)."""
    agg = aggregate(rows)
    if agg.empty:
        return pd.DataFrame()
    stress = agg[agg["scenario_role"] == "stress"]
    nominal = agg[agg["scenario_role"] == "nominal"]
    g = stress.groupby("design_id")
    out = pd.DataFrame(
        {
            "robust": g["r_mean"].max(),
            "worst_scenario": g.apply(lambda d: d.loc[d["r_mean"].idxmax(), "scenario"], include_groups=False),
            "mean_stress": g["r_mean"].mean(),
            "best_stress": g["r_mean"].min(),
            "n_stress": g["scenario"].nunique(),
            "seed_cv_max": g.apply(lambda d: float((d["r_std"] / d["r_mean"]).max()), include_groups=False),
        }
    )
    out["nominal"] = nominal.set_index("design_id")["r_mean"]
    first = agg.groupby("design_id").first()
    for c in ("d50", "sigma_ln", "aspect", "order", "loading_target", "tag", "n_seeds"):
        out[c] = first[c]
    if n_stress is not None:
        out = out[out["n_stress"] == n_stress]
    return out.reset_index()


def job_cost(rows: pd.DataFrame) -> pd.DataFrame:
    """Simulator cost per design: solver calls and CPU seconds (geometry + mesh + solves)."""
    ok = rows[rows["status"] == "ok"]
    per_job = ok.groupby("job_key").agg(
        design_id=("design_id", "first"), t=("job_time_s", "first"), calls=("scenario", "count")
    )
    return per_job.groupby("design_id").agg(calls=("calls", "sum"), cpu_s=("t", "sum"))


class RobustEI:
    """Monte-Carlo expected improvement of the worst-scenario objective."""

    def __init__(self, model: GPSurrogate, space: DesignSpace, stress_phys: np.ndarray, n_mc: int = 256, seed: int = 0):
        self.model = model
        self.space = space
        self.phys = stress_phys  # (S, 2) physical features of the stress scenarios
        self.z = np.random.default_rng(seed).standard_normal((n_mc, stress_phys.shape[0]))

    def predict(self, U: np.ndarray):
        n, s = U.shape[0], self.phys.shape[0]
        X = np.hstack([np.repeat(U, s, axis=0), np.tile(self.phys, (n, 1))])
        mu, sd = self.model.predict(X)
        return mu.reshape(n, s), sd.reshape(n, s)

    def __call__(self, U: np.ndarray, best: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        mu, sd = self.predict(U)
        # samples of worst-scenario resistance (scenarios treated as independent given x)
        samples = np.exp(mu[:, None, :] + sd[:, None, :] * self.z[None, :, :]).max(axis=2)
        ei = np.maximum(best - samples, 0.0).mean(axis=1)
        robust_mean = np.exp(mu).max(axis=1)
        return ei, robust_mean, sd.mean(axis=1)


def propose_batch(
    acq: RobustEI,
    space: DesignSpace,
    evaluated_u: np.ndarray,
    failed_u: np.ndarray,
    best: float,
    top_u: np.ndarray,
    batch: int,
    explore: int,
    seed: int,
    min_dist: float = 0.12,
):
    """Constrained acquisition search over a space-filling pool plus local perturbations."""
    rng = np.random.default_rng(seed)
    pool = [space.sobol(4096, seed)]
    for u in top_u:
        pool.append(np.clip(u + 0.06 * rng.standard_normal((256, space.dim)), 0, 1))
    pool = np.vstack(pool)

    def far(u, ref, d):
        return ref.size == 0 or np.min(np.linalg.norm(ref - u, axis=1)) >= d

    keep = np.array([far(u, evaluated_u, 0.02) and far(u, failed_u, 0.08) for u in pool])
    pool = pool[keep]
    ei, rmean, sd = acq(pool, best)
    chosen: list[np.ndarray] = []
    kinds: list[str] = []
    for i in np.argsort(-ei):
        if len(chosen) >= batch - explore:
            break
        if far(pool[i], np.array(chosen).reshape(-1, space.dim), min_dist):
            chosen.append(pool[i])
            kinds.append("ei")
    for i in np.argsort(-sd):
        if len(chosen) >= batch:
            break
        if far(pool[i], np.array(chosen).reshape(-1, space.dim), min_dist):
            chosen.append(pool[i])
            kinds.append("explore")
    return np.array(chosen), kinds, {"pool": int(pool.shape[0]), "max_ei": float(ei.max())}


def bayes_opt(
    run, rows: pd.DataFrame, space: DesignSpace, scenarios, material, iterations: int, batch: int, explore: int, seed: int,
    tag_prefix: str = "bo",
):
    """``run(designs, tag, kinds) -> rows`` evaluates designs with the physical solver."""
    from .provenance import config_hash

    stress = [s for s in scenarios if s.role == "stress"]
    phys = physical_features(
        np.array([s.r_interface for s in stress]),
        np.array([material.k_filler_par * s.k_filler_scale for s in stress]),
    )
    history = []
    failed_u = np.zeros((0, space.dim))
    for it in range(iterations):
        t0 = time.perf_counter()
        agg = aggregate(rows)
        obj = design_objectives(rows, n_stress=len(stress))
        X = features(agg, space)
        y = agg["log_r_mean"].to_numpy()
        nv = (agg["log_r_var"] / agg["n_seeds"].clip(lower=1)).to_numpy()
        model = GPSurrogate(n_restarts=2, seed=seed + it).fit(X, y, nv)
        t_train = time.perf_counter() - t0
        t1 = time.perf_counter()
        best = float(obj["robust"].min())
        top = obj.nsmallest(5, "robust")
        top_u = space.unit_from_columns(top)
        evaluated_u = space.unit_from_columns(agg.drop_duplicates("design_id"))
        acq = RobustEI(model, space, phys, seed=seed + it)
        U, kinds, info = propose_batch(acq, space, evaluated_u, failed_u, best, top_u, batch, explore, seed + 100 * it)
        t_acq = time.perf_counter() - t1
        designs = [space.from_unit(u) for u in U]
        new = run(designs, f"{tag_prefix}-{it + 1}", kinds)
        fails = new.attrs.get("failed_design_ids", set())
        for u, d in zip(U, designs):
            if config_hash(d.to_dict(), 12) in fails:
                failed_u = np.vstack([failed_u, u])
        rows = pd.concat([rows, new], ignore_index=True)
        obj_new = design_objectives(rows, n_stress=len(stress))
        history.append(
            {
                "iteration": it + 1,
                "train_time_s": t_train,
                "acquisition_time_s": t_acq,
                "kernel": str(model.model.kernel_),
                "best_robust": float(obj_new["robust"].min()),
                "design_ids": [config_hash(d.to_dict(), 12) for d in designs],
                "kinds": kinds,
                "n_designs": int(obj_new.shape[0]),
                "n_failed_designs": int(failed_u.shape[0]),
                **info,
            }
        )
        print(f"  BO iteration {it + 1}: best robust R = {history[-1]['best_robust']:.4e} m2K/W", flush=True)
    return rows, history


def best_so_far(obj: pd.DataFrame, cost: pd.DataFrame, order: list[str], col: str = "robust", overhead_s: float = 0.0):
    """Best-so-far trajectory against cumulative solver calls and CPU time."""
    o = obj.set_index("design_id")
    c = cost.reindex(order).fillna(0)
    calls = np.cumsum(c["calls"].to_numpy())
    cpu = np.cumsum(c["cpu_s"].to_numpy())
    vals = o[col].reindex(order).to_numpy()
    best = np.fmin.accumulate(np.where(np.isnan(vals), np.inf, vals))
    return pd.DataFrame({"design_id": order, "calls": calls, "cpu_s": cpu + overhead_s, "best": best})
