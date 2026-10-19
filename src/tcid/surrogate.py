"""Surrogates over design descriptors and uncertain physical inputs.

Target: log of the positive layer resistance. Random microstructure realizations are
treated as variability: repeated seeds are aggregated and their mean is modelled with
an observation-noise estimate (per-point GP noise = seed variance / n_seeds). The model
never tries to predict the identity of an unrecorded particle arrangement.

Two kinds of uncertainty are kept separate:
* surrogate uncertainty (GP predictive sd) - what the approximation does not know about
  its simulator; more data reduces it;
* physical-input uncertainty - the scenario inputs (interface resistance, filler
  conductivity) that are unknown for the real material; more data does not reduce it.
"""

from __future__ import annotations

import math
import time
import warnings

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.ensemble import ExtraTreesRegressor
from sklearn.exceptions import ConvergenceWarning
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import ConstantKernel, Matern, WhiteKernel
from sklearn.model_selection import GroupKFold

from .design_space import DesignSpace

# fixed (not data-fitted) scaling of the physical inputs -> no preprocessing leakage
LOG_R_BOUNDS = (math.log10(1e-9), math.log10(1e-5))
KF_BOUNDS = (10.0, 50.0)


def physical_features(r_interface: np.ndarray, k_filler: np.ndarray) -> np.ndarray:
    lr = np.log10(np.maximum(np.asarray(r_interface, float), 1e-12))
    a = (lr - LOG_R_BOUNDS[0]) / (LOG_R_BOUNDS[1] - LOG_R_BOUNDS[0])
    b = (np.asarray(k_filler, float) - KF_BOUNDS[0]) / (KF_BOUNDS[1] - KF_BOUNDS[0])
    return np.column_stack([a, b])


def aggregate(rows: pd.DataFrame) -> pd.DataFrame:
    """Aggregate seeds: one row per (design, scenario) with mean/var of log resistance."""
    ok = rows[rows["status"] == "ok"].copy()
    ok["log_r"] = np.log(ok["r_layer"])
    keys = ["design_id", "scenario"]
    first = ok.groupby(keys).first()
    agg = ok.groupby(keys).agg(
        r_mean=("r_layer", "mean"),
        r_std=("r_layer", "std"),
        log_r_mean=("log_r", "mean"),
        log_r_var=("log_r", "var"),
        n_seeds=("seed", "nunique"),
        loading_achieved=("loading_achieved", "mean"),
        order_realized=("geo_order_realized", "mean"),
        sim_time_s=("solve_time_s", "sum"),
    )
    cols = ["d50", "sigma_ln", "aspect", "order", "loading_target", "r_interface", "k_filler_par", "scenario_role", "tag"]
    out = agg.join(first[cols]).reset_index()
    out["log_r_var"] = out["log_r_var"].fillna(0.0)
    out["r_std"] = out["r_std"].fillna(0.0)
    return out


def features(agg: pd.DataFrame, space: DesignSpace) -> np.ndarray:
    xd = space.unit_from_columns(agg)
    xp = physical_features(agg["r_interface"].to_numpy(), agg["k_filler_par"].to_numpy())
    return np.hstack([xd, xp])


class GPSurrogate:
    name = "gaussian_process"

    def __init__(self, n_restarts: int = 2, seed: int = 0):
        self.n_restarts = n_restarts
        self.seed = seed
        self.model = None
        self.fit_time_s = 0.0

    def fit(self, X: np.ndarray, y: np.ndarray, noise_var: np.ndarray | None = None) -> "GPSurrogate":
        t0 = time.perf_counter()
        d = X.shape[1]
        kernel = ConstantKernel(1.0, (1e-3, 1e3)) * Matern(
            length_scale=np.full(d, 0.5), length_scale_bounds=(1e-2, 1e2), nu=2.5
        ) + WhiteKernel(1e-5, (1e-9, 1e-1))
        alpha = 1e-8 if noise_var is None else np.maximum(noise_var, 1e-8)
        self.model = GaussianProcessRegressor(
            kernel=kernel,
            alpha=alpha,
            normalize_y=True,
            n_restarts_optimizer=self.n_restarts,
            random_state=self.seed,
        )
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", ConvergenceWarning)
            self.model.fit(X, y)
        self.fit_time_s = time.perf_counter() - t0
        return self

    def predict(self, X: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        mu, sd = self.model.predict(X, return_std=True)
        return mu, sd

    def describe(self) -> dict:
        return {"model": self.name, "kernel": str(self.model.kernel_), "fit_time_s": self.fit_time_s}


class TreeSurrogate:
    """Small extra-trees ensemble: practical alternative when the GP becomes expensive."""

    name = "extra_trees"

    def __init__(self, n_estimators: int = 300, seed: int = 0):
        self.model = ExtraTreesRegressor(n_estimators=n_estimators, min_samples_leaf=2, random_state=seed, n_jobs=1)
        self.fit_time_s = 0.0

    def fit(self, X, y, noise_var=None):
        t0 = time.perf_counter()
        self.model.fit(X, y)
        self.fit_time_s = time.perf_counter() - t0
        return self

    def predict(self, X):
        preds = np.stack([t.predict(X) for t in self.model.estimators_])
        return preds.mean(0), preds.std(0)

    def describe(self) -> dict:
        return {"model": self.name, "fit_time_s": self.fit_time_s}


def cross_validate(agg: pd.DataFrame, space: DesignSpace, make_model, n_splits: int = 5) -> dict:
    """Held-out error with folds split by design (all scenarios/seeds of a design together)."""
    X = features(agg, space)
    y = agg["log_r_mean"].to_numpy()
    nv = (agg["log_r_var"] / agg["n_seeds"].clip(lower=1)).to_numpy()
    groups = agg["design_id"].to_numpy()
    pred = np.full_like(y, np.nan)
    psd = np.full_like(y, np.nan)
    for tr, te in GroupKFold(n_splits=min(n_splits, len(np.unique(groups)))).split(X, y, groups):
        m = make_model().fit(X[tr], y[tr], nv[tr])
        pred[te], psd[te] = m.predict(X[te])
    rel = np.abs(np.exp(pred) - np.exp(y)) / np.exp(y)
    low = y <= np.quantile(y, 0.25)
    z = (y - pred) / np.sqrt(psd**2 + nv)
    # ordering of designs within each scenario
    rhos = []
    for _, g in agg.assign(pred=pred).groupby("scenario"):
        if len(g) > 3:
            rhos.append(spearmanr(g["log_r_mean"], g["pred"]).statistic)
    return {
        "n_points": int(len(y)),
        "n_designs": int(len(np.unique(groups))),
        "median_rel_error": float(np.median(rel)),
        "p90_rel_error": float(np.quantile(rel, 0.9)),
        "max_rel_error": float(np.max(rel)),
        "median_rel_error_low_resistance_quartile": float(np.median(rel[low])),
        "coverage_95": float(np.mean(np.abs(z) < 1.96)),
        "spearman_within_scenario_mean": float(np.mean(rhos)) if rhos else float("nan"),
        "predictions": pred.tolist(),
        "pred_sd": psd.tolist(),
    }
