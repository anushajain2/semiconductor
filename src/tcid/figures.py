"""Figures for the technical report (guide section 8).

Colour use: categorical slots in a fixed order (identity follows the candidate, never
its rank), a single-hue sequential ramp for temperature, recessive grids, text in ink
colours, direct labels where there are few series.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.collections import PolyCollection  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap, ListedColormap  # noqa: E402
from matplotlib.patches import Ellipse, Patch  # noqa: E402

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK2 = "#52514e"
GRID = "#e4e3df"
SLOTS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
BLUE_RAMP = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]
TEMP_CMAP = LinearSegmentedColormap.from_list("temp_blue", BLUE_RAMP[::-1])

# fixed identity -> colour assignment for every figure
IDENTITY = {
    "baseline_random_orientation": SLOTS[0],
    "recommended": SLOTS[1],
    "small_particles_aligned": SLOTS[2],
    "large_particles_aligned": SLOTS[3],
    "medium_particles_aligned": SLOTS[4],
    "bayesian_optimisation": SLOTS[1],
    "random_search": SLOTS[0],
    "sweep": SLOTS[2],
}
NICE = {
    "baseline_random_orientation": "Random-orientation baseline",
    "recommended": "Recommended (verified best)",
    "small_particles_aligned": "Aligned, small particles",
    "large_particles_aligned": "Aligned, large particles",
    "medium_particles_aligned": "Aligned, medium particles",
    "bayesian_optimisation": "Bayesian optimisation",
    "random_search": "Random search",
    "sweep": "3-level sweep",
}


def _style():
    plt.rcParams.update(
        {
            "figure.facecolor": SURFACE,
            "axes.facecolor": SURFACE,
            "axes.edgecolor": GRID,
            "axes.labelcolor": INK2,
            "axes.titlecolor": INK,
            "axes.titlesize": 11,
            "axes.labelsize": 9.5,
            "xtick.color": INK2,
            "ytick.color": INK2,
            "xtick.labelsize": 8.5,
            "ytick.labelsize": 8.5,
            "axes.grid": True,
            "grid.color": GRID,
            "grid.linewidth": 0.6,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "legend.frameon": False,
            "legend.fontsize": 8.5,
            "lines.linewidth": 2.0,
            "font.size": 9.5,
            "savefig.dpi": 160,
            "savefig.bbox": "tight",
        }
    )


def _nice(label: str) -> str:
    return NICE.get(label, label.replace("_", " "))


def plot_geometry(g, path, ax=None):
    _style()
    own = ax is None
    if own:
        fig, ax = plt.subplots(figsize=(7, 3.8))
    lay = g.layer
    ax.add_patch(plt.Rectangle((0, 0), lay["width"] * 1e6, lay["thickness"] * 1e6, fc="#f3f2ee", ec=INK2, lw=0.8))
    for i in range(g.n):
        ax.add_patch(Ellipse((g.cx[i] * 1e6, g.cy[i] * 1e6), 2 * g.a[i] * 1e6, 2 * g.b[i] * 1e6,
                             angle=math.degrees(g.theta[i]), fc=SLOTS[0], ec=SURFACE, lw=0.6))
    ax.set_xlim(0, lay["width"] * 1e6)
    ax.set_ylim(0, lay["thickness"] * 1e6)
    ax.set_aspect("equal")
    ax.grid(False)
    ax.set_xlabel("x (µm)")
    ax.set_ylabel("y, heat-flow direction (µm)")
    if own:
        fig.savefig(path)
        plt.close(fig)


def _field_panel(ax, z, title, vmin=0.0, vmax=1.0):
    p, t = z["points"] * 1e6, z["triangles"]
    tc = ax.tripcolor(p[:, 0], p[:, 1], t, z["T"], cmap=TEMP_CMAP, vmin=vmin, vmax=vmax, shading="gouraud", rasterized=True)
    for i in range(len(z["cx"])):
        ax.add_patch(Ellipse((z["cx"][i] * 1e6, z["cy"][i] * 1e6), 2 * z["a"][i] * 1e6, 2 * z["b"][i] * 1e6,
                             angle=math.degrees(z["theta"][i]), fc="none", ec="white", lw=0.5))
    ax.set_xlim(0, float(z["width"]) * 1e6)
    ax.set_ylim(0, float(z["thickness"]) * 1e6)
    ax.set_aspect("equal")
    ax.grid(False)
    ax.set_title(title, loc="left")
    ax.set_xlabel("x (µm)")
    return tc


def fig_fields(study, out: Path):
    rows = pd.read_parquet(study.out / "final" / "final_rows.parquet")
    f = rows[rows["tag"] == "final-fields"]
    cands = json.loads((study.out / "final" / "candidates.json").read_text())
    from .provenance import config_hash

    want = [("baseline_random_orientation", "Random-orientation baseline"), ("recommended", "Recommended design")]
    fig, axes = plt.subplots(1, 2, figsize=(12, 3.9), sharey=True, layout="constrained")
    tc = None
    for ax, (lab, nice) in zip(axes, want):
        did = config_hash(cands[lab], 12)
        r = f[f["design_id"] == did].iloc[0]
        z = np.load(study.out / "fields" / f"{r['run_key']}.npz")
        tc = _field_panel(ax, z, f"{nice}\nR_layer = {r['r_layer'] * 1e6:.0f} mm²K/W, k_eff = {r['k_eff']:.2f} W/mK")
    axes[0].set_ylabel("y (µm); hot face y = 0")
    cb = fig.colorbar(tc, ax=axes, shrink=0.85, pad=0.02)
    cb.set_label("Temperature (K above cold face)")
    cb.outline.set_visible(False)
    fig.suptitle("Temperature fields at the same scale (nominal scenario, fresh seed, refined mesh)", x=0.02, ha="left", color=INK)
    fig.savefig(out / "fig1_fields.png")
    plt.close(fig)


def fig_design_map(study, out: Path):
    m = pd.read_parquet(study.out / "questions" / "q1_map.parquet")
    g = pd.read_parquet(study.out / "questions" / "q1_summary.parquet")
    labels = [k for k in IDENTITY if k in set(m["winner"]) | set(m["runner_up"]) | set(g["label"])]
    labels = [l for l in labels if l in g["label"].unique()]
    cmap = ListedColormap([IDENTITY[l] for l in labels])
    rvals = sorted(m["r_interface"].unique())
    fr = sorted(m["order_fraction"].unique())
    kfs = sorted(m["k_filler"].unique())
    fig = plt.figure(figsize=(12, 6.6))
    gs = fig.add_gridspec(2, len(kfs), height_ratios=[1, 1.15], hspace=0.45, wspace=0.12)
    for j, kf in enumerate(kfs):
        ax = fig.add_subplot(gs[0, j])
        sub = m[m["k_filler"] == kf]
        Z = np.full((len(fr), len(rvals)), np.nan)
        for _, r in sub.iterrows():
            Z[fr.index(r["order_fraction"]), rvals.index(r["r_interface"])] = labels.index(r["winner"])
        ax.imshow(Z, cmap=cmap, vmin=-0.5, vmax=len(labels) - 0.5, aspect="auto", origin="lower")
        for _, r in sub.iterrows():
            if not r["distinguishable"]:
                y, x = fr.index(r["order_fraction"]), rvals.index(r["r_interface"])
                ax.add_patch(plt.Rectangle((x - 0.5, y - 0.5), 1, 1, fill=False, hatch="///", ec=SURFACE, lw=0))
        ax.axvspan(-0.5, rvals.index(1e-8) - 0.5, color=SURFACE, alpha=0.55, lw=0)
        ax.axvspan(rvals.index(1e-6) + 0.5, len(rvals) - 0.5, color=SURFACE, alpha=0.55, lw=0)
        ax.set_xticks(range(len(rvals)), [f"{r:.0e}".replace("e-0", "e-") for r in rvals], rotation=60)
        ax.set_yticks(range(len(fr)), [f"{f:.2f}" for f in fr] if j == 0 else [])
        ax.grid(False)
        ax.set_title(f"k_filler = {kf:g} W/mK", loc="left")
        ax.set_xlabel("Interface resistance (m²K/W)")
        if j == 0:
            ax.set_ylabel("Achieved / target orientation order")
    handles = [Patch(fc=IDENTITY[l], label=_nice(l)) for l in labels]
    handles.append(Patch(fc="#999999", hatch="///", ec=SURFACE, label="winner not resolved from runner-up (Welch p ≥ 0.05)"))
    handles.append(Patch(fc="#dddddd", label="faded: outside supported range"))
    fig.legend(handles=handles, loc="upper center", ncol=4, bbox_to_anchor=(0.5, 1.02))
    ax = fig.add_subplot(gs[1, :])
    sub = g[(g["k_filler_par"] == 30.0) & (g["order_fraction"] == 1.0)]
    ends = []
    for lab in labels:
        s = sub[sub["label"] == lab].sort_values("r_interface")
        if s.empty:
            continue
        ax.errorbar(s["r_interface"], s["r_mean"] * 1e6, yerr=2 * s["se"] * 1e6, color=IDENTITY[lab], marker="o", ms=4,
                    capsize=2, label=_nice(lab))
        ends.append([s["r_mean"].iloc[-1] * 1e6, lab, s["r_interface"].iloc[-1]])
    # direct labels at the line ends, nudged apart so they never overlap
    ends.sort()
    gap = 0.045 * (sub["r_mean"].max() - sub["r_mean"].min()) * 1e6
    for i in range(1, len(ends)):
        ends[i][0] = max(ends[i][0], ends[i - 1][0] + gap)
    for y, lab, x in ends:
        ax.annotate(_nice(lab), (x, y), xytext=(6, 0), textcoords="offset points", color=INK2, fontsize=8, va="center")
    ax.axvspan(1e-8, 1e-6, color=GRID, alpha=0.35, lw=0)
    ax.text(1.1e-8, ax.get_ylim()[1], "supported interface range", va="top", color=INK2, fontsize=8)
    ax.set_xscale("log")
    ax.set_xlabel("Particle-matrix interface resistance (m²K/W)")
    ax.set_ylabel("Layer resistance (mm²K/W)")
    ax.set_title("Layer resistance vs interface resistance (k_filler = 30 W/mK, full orientation order; bars = ±2 SE over seeds)", loc="left")
    ax.set_xlim(right=ax.get_xlim()[1] * 4)
    fig.savefig(out / "fig2_design_map.png")
    plt.close(fig)


def fig_contacts(study, out: Path):
    q = pd.read_parquet(study.out / "questions" / "q2_benefit.parquet")
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(12, 4.2), gridspec_kw={"width_ratios": [1.25, 1], "wspace": 0.25})
    nominal = study.nominal.name
    th = sorted(q["thickness"].unique())
    for (pair, ls), k0 in zip((("recommended", "-"), ("medium_particles_aligned", "--")), (0, 0)):
        for k, t in enumerate(th):
            s = q[(q["thickness"] == t) & (q["pair"] == pair)]
            if s.empty:
                continue
            col = SLOTS[k % 3]
            nom = s[s["scenario"] == nominal].sort_values("r_contact_total")
            x = nom["r_contact_total"].replace(0, 3e-7)
            a1.fill_between(x, nom["benefit_ci_low"] * 100, nom["benefit_ci_high"] * 100, color=col, alpha=0.12, lw=0)
            a1.plot(x, nom["benefit_rel_total"] * 100, color=col, ls=ls, marker="o", ms=3.5,
                    label=f"{'recommended' if pair == 'recommended' else 'medium pair'}, t = {t * 1e6:.0f} µm")
    a1.axhline(10, color=INK2, lw=1, ls="--")
    a1.text(1.4e-3, 10.3, "10% benefit", color=INK2, fontsize=8, ha="right")
    a1.axvspan(3e-7, 3e-5 * 2, color=GRID, alpha=0.3, lw=0)
    a1.set_xscale("log")
    a1.set_xlim(2.5e-7, 1.5e-3)
    a1.set_xlabel("Total outer contact resistance R_top + R_bottom (m²K/W; leftmost point = 0)")
    a1.set_ylabel("Assembly benefit vs baseline (%)")
    a1.set_title("Assembly benefit vs outer contacts, nominal scenario\n(band: 95% bootstrap interval over 12 seeds per design)",
                 loc="left", fontsize=10)
    a1.legend(loc="upper right", fontsize=7.5)
    s = q[np.isclose(q["thickness"], study.cfg.layer.thickness) & (q["scenario"] == nominal) & (q["pair"] == "recommended")]
    cs = [0, 1e-5, 3e-5, 1e-4]
    xs = np.arange(len(cs))
    w = 0.36
    for off, col_layer, lab in ((-w / 2 - 0.01, SLOTS[0], "baseline"), (w / 2 + 0.01, SLOTS[1], "candidate")):
        lay = [float(s[s["r_contact_total"] == c][f"r_layer_{lab}"].iloc[0]) * 1e6 for c in cs]
        a2.bar(xs + off, lay, w, color=col_layer, label=f"layer ({'baseline' if lab == 'baseline' else 'recommended'})")
        a2.bar(xs + off, [c * 1e6 for c in cs], w, bottom=lay, color=SLOTS[3], edgecolor=SURFACE, lw=2,
               label="outer contacts" if lab == "candidate" else None)
    a2.set_xticks(xs, [f"{c * 1e6:g}" for c in cs])
    a2.set_xlabel("Total contact resistance (mm²K/W), t = 100 µm, nominal scenario")
    a2.set_ylabel("Area-specific resistance (mm²K/W)")
    a2.set_title("Layer vs contact contributions\n(t = 100 µm)", loc="left", fontsize=10)
    a2.legend(loc="upper left")
    fig.savefig(out / "fig3_contacts.png")
    plt.close(fig)


def fig_realizations(study, out: Path):
    rows = pd.read_parquet(study.out / "final" / "final_rows.parquet")
    fin = json.loads((study.out / "final" / "final.json").read_text())
    labels = fin["labels"]
    f = rows[rows["tag"] == "final"].copy()
    f["label"] = f["design_id"].map(labels)
    worst = study.stress
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(12, 4.6), gridspec_kw={"width_ratios": [2.2, 1]})
    order = sorted(f["label"].unique(), key=lambda l: fin["verified"][l]["robust"])
    order = [l for l in order if l != "conventional_circles"] + ["conventional_circles"]
    scen = [study.nominal.name] + [s.name for s in worst]
    rng = np.random.default_rng(0)
    width = 0.8 / len(scen)
    for j, sc in enumerate(scen):
        col = SLOTS[j]
        for i, lab in enumerate(order):
            v = f[(f["label"] == lab) & (f["scenario"] == sc)]["r_layer"].to_numpy() * 1e6
            x = i - 0.4 + width * (j + 0.5)
            a1.scatter(x + rng.uniform(-width / 4, width / 4, v.size), v, s=12, color=col, alpha=0.8, lw=0,
                       label=sc if i == 0 else None)
            a1.hlines(v.mean(), x - width / 2.5, x + width / 2.5, color=INK, lw=1.2)
    a1.set_xticks(range(len(order)), [l.replace("+", " + ").replace("__S0", " (S=0)") for l in order], rotation=40,
                  ha="right", fontsize=7.5)
    a1.set_ylabel("Layer resistance (mm²K/W)")
    a1.set_title("Fresh realizations (6 new seeds per design; bar = mean), refined mesh", loc="left")
    a1.legend(ncol=3, loc="upper left", fontsize=7.5)
    ref = pd.read_parquet(study.out / "final" / "mesh_refinement.parquet")
    sfs = {"loop": 1.0, "refined": float(study.budget["final_refinements"][0]), "fine": float(study.budget["final_refinements"][1])}
    for role, sub in ref.groupby("role"):
        key = "recommended" if role == "recommended" else "baseline_random_orientation"
        for sc, ss in sub.groupby("scenario"):
            y = [ss[c].iloc[0] / ss["fine"].iloc[0] for c in sfs]
            a2.plot(list(sfs.values()), y, color=IDENTITY[key], marker="o", ms=4, alpha=0.8,
                    label=_nice(key) if sc == study.nominal.name else None)
    a2.axhspan(0.98, 1.02, color=GRID, alpha=0.4, lw=0)
    a2.invert_xaxis()
    a2.set_xlabel("Mesh size factor (smaller = finer)")
    a2.set_ylabel("R / R(finest)")
    a2.set_title("Mesh refinement (2 seeds, all scenarios)", loc="left")
    a2.legend(loc="best")
    fig.savefig(out / "fig4_realizations_mesh.png")
    plt.close(fig)


def fig_search(study, out: Path):
    c = pd.read_parquet(study.out / "search" / "best_so_far.parquet")
    fin = json.loads((study.out / "final" / "final.json").read_text())
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(12, 4.2), sharey=True)
    for meth, s in c.groupby("method"):
        col = IDENTITY[meth]
        dy = {"bayesian_optimisation": -7, "sweep": 7, "random_search": 0}[meth]
        for ax, x in ((a1, "calls"), (a2, "cpu_s")):
            xv = s[x] / (3600 if x == "cpu_s" else 1)
            ax.step(xv, s["best"] * 1e6, where="post", color=col, label=_nice(meth))
            ax.annotate(_nice(meth), (xv.iloc[-1], s["best"].iloc[-1] * 1e6), xytext=(40, dy), textcoords="offset points",
                        color=INK2, fontsize=8, va="center")
    for meth, lab in fin["method_picks"].items():
        v = fin["confirmed"][lab]
        s = c[c["method"] == meth]
        a1.errorbar([s["calls"].iloc[-1]], [v["robust"] * 1e6], yerr=[2 * v["robust_se"] * 1e6], fmt="D", ms=6,
                    color=IDENTITY[meth], mec=SURFACE, ecolor=IDENTITY[meth], capsize=2, zorder=5)
    a1.set_xlabel("Solver calls (one call = one geometry × scenario solve)")
    a1.set_ylabel("Best robust objective (mm²K/W)")
    a1.set_title("Best design found vs simulator calls (one campaign per method)\n"
                 "◆ = method's best-of-4, confirmed on 12 independent seeds (±2 SE)", loc="left", fontsize=10)
    a2.set_xlabel("CPU hours (geometry + meshing + solves + surrogate training/acquisition)")
    a2.set_title("…and vs elapsed CPU time\n", loc="left", fontsize=10)
    for ax in (a1, a2):
        ax.set_xlim(right=ax.get_xlim()[1] * 1.3)
    fig.savefig(out / "fig5_search_efficiency.png")
    plt.close(fig)


def fig_surrogate(study, out: Path):
    p = study.out / "surrogate" / "cv_predictions_loop.parquet"
    if not p.exists():
        p = study.out / "surrogate" / "cv_predictions_initial.parquet"
    d = pd.read_parquet(p)
    fig, ax = plt.subplots(figsize=(4.6, 4.3))
    y, yp = np.exp(d["log_r_mean"]) * 1e6, np.exp(d["pred"]) * 1e6
    ax.errorbar(y, yp, yerr=[yp - np.exp(d["pred"] - 2 * d["pred_sd"]) * 1e6, np.exp(d["pred"] + 2 * d["pred_sd"]) * 1e6 - yp],
                fmt="o", ms=3.5, color=SLOTS[0], ecolor=GRID, elinewidth=1)
    lo, hi = min(y.min(), yp.min()), max(y.max(), yp.max())
    ax.plot([lo, hi], [lo, hi], color=INK2, lw=1)
    ax.set_xlabel("Solver (seed mean), mm²K/W")
    ax.set_ylabel("GP held-out prediction, mm²K/W")
    ax.set_title("Design-grouped 5-fold CV (±2 sd)", loc="left")
    fig.savefig(out / "fig6_surrogate_parity.png")
    plt.close(fig)


def fig_q3(study, out: Path):
    a = pd.read_parquet(study.out / "questions" / "q3_advantage.parquet")
    a = a[a["scenario"] == study.nominal.name].copy()
    order = list(dict.fromkeys(a["variant"]))
    fig, ax = plt.subplots(figsize=(8, 3.8))
    y = a.set_index("variant").loc[order]
    ax.barh(range(len(order)), y["advantage_rel"] * 100, xerr=2 * y["se"] / (y["advantage"] / y["advantage_rel"]) * 100,
            color=[SLOTS[1] if v == "reference" else SLOTS[0] for v in order], height=0.6, capsize=2, ecolor=INK2)
    ax.set_yticks(range(len(order)), order)
    ax.invert_yaxis()
    ax.set_xlabel("Resistance advantage over baseline (%), nominal scenario, ±2 SE")
    ax.set_title("Does the optimum rely on ideal geometry? (fresh seeds per variant)", loc="left")
    fig.savefig(out / "fig7_perturbations.png")
    plt.close(fig)


def fig_mesh_detail(study, out: Path):
    """Near-contact and margin resolution for visual inspection."""
    from .config import Design
    from .geometry import generate
    from .mesh import FILLER, mesh_particles

    c = json.loads((study.out / "final" / "candidates.json").read_text())
    d = Design(**c["recommended"])
    g = generate(d, study.cfg.layer, study.cfg.geometry, study.budget["final_seeds"][0])
    m = mesh_particles(g, study.cfg.mesh)
    from .geometry import pairwise_gaps

    gaps = pairwise_gaps(g.cx, g.cy, g.a, g.b, g.theta, cutoff=20e-6)
    i = int(np.argmin(gaps))
    cx, cy = g.cx[i] * 1e6, g.cy[i] * 1e6
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.2))
    for ax, (x0, y0, half, title) in zip(
        axes,
        [(cx, cy, max(3 * g.a[i] * 1e6 / 2, 6), "Narrowest particle gap"),
         (cx, g.layer["margin"] * 1e6, 8, "Resin margin at the hot face")],
    ):
        p = m.points * 1e6
        tri = m.triangles
        cen = p[tri].mean(axis=1)
        sel = (np.abs(cen[:, 0] - x0) < half * 1.4) & (np.abs(cen[:, 1] - y0) < half)
        kind = m.region_kind[m.region[sel]]
        pc = PolyCollection(p[tri[sel]], facecolors=np.where(kind == FILLER, "#b7d3f6", "#f3f2ee"), edgecolors=INK2, linewidths=0.25)
        ax.add_collection(pc)
        ax.set_xlim(x0 - half * 1.4, x0 + half * 1.4)
        ax.set_ylim(max(0, y0 - half), y0 + half)
        ax.set_aspect("equal")
        ax.grid(False)
        ax.set_title(title, loc="left")
        ax.set_xlabel("x (µm)")
    axes[0].set_ylabel("y (µm)")
    fig.suptitle(f"Mesh detail, recommended design, seed {study.budget['final_seeds'][0]}: min gap {gaps.min() * 1e6:.2f} µm, "
                 f"{m.n_elements} elements, 0 bridge elements = {m.info['bridge_elements'] == 0}", x=0.07, ha="left")
    fig.savefig(out / "fig8_mesh_detail.png")
    plt.close(fig)


def fig_crossings(study, out: Path):
    cur = pd.read_parquet(study.out / "questions" / "q1_crossing_curves.parquet")
    cr = json.loads((study.out / "questions" / "q1_crossings.json").read_text())
    pairs = [k for k, p in cr["pairs"].items() if "estimate" in p]
    fig, axes = plt.subplots(1, len(pairs), figsize=(13, 3.9), sharey=False)
    for ax, (pair, col) in zip(axes, zip(pairs, SLOTS)):
        a, b = pair.split("|")
        s = cur[cur["pair"] == pair].sort_values("r_interface")
        ax.fill_between(s["r_interface"], s["ci_low"] * 1e6, s["ci_high"] * 1e6, color=col, alpha=0.18, lw=0)
        ax.plot(s["r_interface"], s["diff"] * 1e6, color=col, marker="o", ms=3)
        ax.axhline(0, color=INK2, lw=1)
        ax.axvspan(s["r_interface"].min(), 1e-6, color=GRID, alpha=0.35, lw=0)
        p = cr["pairs"][pair]
        if np.isfinite(p["ci_low"]):
            ax.axvspan(p["ci_low"], p["ci_high"], color=INK2, alpha=0.12, lw=0)
            ax.axvline(p["estimate"], color=INK2, lw=1, ls="--")
        ax.set_xscale("log")
        ax.set_xlabel("Interface resistance (m²K/W)")
        ax.set_title(f"{_nice(a)} − {_nice(b)}\ncrossing ≈ {p['estimate']:.1e} [{p['ci_low']:.1e}, {p['ci_high']:.1e}]",
                     loc="left", fontsize=9)
        ax.text(0.02, 0.96, p["status"], transform=ax.transAxes, va="top", fontsize=7.5, color=INK2)
    axes[0].set_ylabel("Difference in layer resistance (mm²K/W)")
    fig.suptitle("Particle-size crossings: mean difference ±95% bootstrap band (12 seeds per design); grey column = "
                 "supported range; dark band = 95% interval of the crossing", x=0.01, ha="left", fontsize=10, color=INK)
    fig.tight_layout()
    fig.savefig(out / "fig9_crossings.png")
    plt.close(fig)


def make_all(study):
    _style()
    out = Path("docs/figures")
    out.mkdir(parents=True, exist_ok=True)
    for fn in (fig_fields, fig_design_map, fig_contacts, fig_realizations, fig_search, fig_surrogate, fig_q3, fig_mesh_detail,
               fig_crossings):
        try:
            fn(study, out)
            print(f"  figure {fn.__name__}: ok")
        except FileNotFoundError as exc:
            print(f"  figure {fn.__name__}: skipped ({exc})")
