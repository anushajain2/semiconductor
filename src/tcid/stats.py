"""Small-sample statistics used to compare finalists.

Comparisons between designs use independent geometry realizations on each side, so
they are two-sample (Welch) comparisons. Families of comparisons are Holm-corrected.
"Not distinguishable" means only that no difference was resolved; it is never read as
equivalence.
"""

from __future__ import annotations

import math

import numpy as np
from scipy import stats


def welch(a, b, alpha: float = 0.05) -> dict:
    """Difference of means a - b with Welch standard error, df, CI and two-sided p-value."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    na, nb = a.size, b.size
    va, vb = a.var(ddof=1) / na, b.var(ddof=1) / nb
    d = a.mean() - b.mean()
    se = math.sqrt(va + vb)
    if se == 0:
        return {"diff": d, "se": 0.0, "df": float("inf"), "ci_low": d, "ci_high": d, "p": 0.0 if d else 1.0,
                "n_a": na, "n_b": nb}
    df = (va + vb) ** 2 / (va**2 / (na - 1) + vb**2 / (nb - 1))
    tcrit = stats.t.ppf(1 - alpha / 2, df)
    p = 2 * stats.t.sf(abs(d) / se, df)
    return {"diff": d, "se": se, "df": df, "ci_low": d - tcrit * se, "ci_high": d + tcrit * se, "p": p,
            "n_a": na, "n_b": nb}


def holm(pvalues: dict) -> dict:
    """Holm-Bonferroni adjusted p-values for a family of tests."""
    keys = sorted(pvalues, key=pvalues.get)
    m = len(keys)
    out, running = {}, 0.0
    for i, k in enumerate(keys):
        running = max(running, min(1.0, (m - i) * pvalues[k]))
        out[k] = running
    return out


def bootstrap(samples: dict, fn, n: int = 4000, seed: int = 0, alpha: float = 0.05) -> dict:
    """Percentile bootstrap of ``fn(**means)``; each sample array is resampled independently.

    ``samples`` maps a name to the per-realization values of one design/condition: shape (n,)
    or (n, k) for whole curves (rows are realizations). ``fn`` receives the resampled means
    (over rows) by name and returns a float (NaN allowed).
    """
    rng = np.random.default_rng(seed)
    arrs = {k: np.asarray(v, float) for k, v in samples.items()}
    point = fn(**{k: v.mean(axis=0) for k, v in arrs.items()})
    draws = np.empty(n)
    for i in range(n):
        draws[i] = fn(**{k: v[rng.integers(0, len(v), len(v))].mean(axis=0) for k, v in arrs.items()})
    finite = draws[np.isfinite(draws)]
    lo, hi = (np.quantile(finite, [alpha / 2, 1 - alpha / 2]) if finite.size else (np.nan, np.nan))
    return {"estimate": float(point), "ci_low": float(lo), "ci_high": float(hi),
            "fraction_undefined": float(1 - finite.size / n), "n_boot": n}


def first_crossing(x: np.ndarray, d: np.ndarray) -> float:
    """First sign change of d(x) on a log-spaced grid, by log-linear interpolation (NaN if none)."""
    s = np.sign(d)
    for k in range(len(s) - 1):
        if s[k] != 0 and s[k + 1] != 0 and s[k] != s[k + 1]:
            x0, x1 = math.log10(x[k]), math.log10(x[k + 1])
            return 10 ** (x0 + (x1 - x0) * d[k] / (d[k] - d[k + 1]))
    return float("nan")
