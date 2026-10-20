"""Minimal user workflow (guide section 9).

The user chooses the material case, layer thickness, loading constraint and scenario set.
Inside the verified domain the tool returns a few verified candidate distributions,
downloadable configuration files, a baseline comparison and a statement of the modelled
domain. Requests outside the verified domain return an explicit unsupported-case
response, or (``run_if_unsupported``) trigger new solver runs of the verified candidates
at the requested conditions - clearly labelled as not re-optimised.
"""

from __future__ import annotations

import json
import math
from dataclasses import replace
from pathlib import Path

import pandas as pd
import yaml

from .config import Design
from .units import to_si

DOMAIN_TEXT = (
    "Modelled domain: 2D cross-section of a cured particle-filled epoxy layer between isothermal faces with "
    "insulated sidewalls; steady linear conduction; isotropic alumina-like development filler (NOT boron nitride); "
    "area fraction (not a 3D volume fraction); particle-matrix interface resistance and filler conductivity treated "
    "as stress scenarios, not probabilities; outer contacts added in series. Not modelled: 3D connectivity, BN "
    "anisotropy, viscosity/processing, adhesion, cure, dielectric breakdown, stress, ageing. Not a qualified adhesive."
)


def _verified_domain(study) -> dict:
    return {
        "material": study.material.name,
        "thickness_m": [study.cfg.layer.thickness],
        "loading_levels": [float(x) for x in study.cfg.design_space["loading_levels"]],
        "scenarios": sorted(s.name for s in study.scenarios),
    }


def recommend(out, material, thickness, loading, scenarios, out_dir, run_if_unsupported=False,
              config="configs/study_phase1.yaml", cache="cache", workers=None) -> dict:
    from .materials import load_scenarios
    from .studies import Study

    study = Study(config, out, cache, workers)
    t = to_si(thickness)
    dom = _verified_domain(study)
    try:
        req_scen = sorted(s.name for s in load_scenarios(scenarios, study.material))
    except Exception as exc:  # noqa: BLE001
        req_scen = [f"unreadable scenario file: {exc}"]
    problems = []
    if material != dom["material"]:
        problems.append(f"material case '{material}' is not supported (verified: '{dom['material']}')")
    if not any(abs(t - x) < 1e-9 for x in dom["thickness_m"]):
        problems.append(f"thickness {t * 1e6:.0f} um outside the verified set {[round(x * 1e6, 3) for x in dom['thickness_m']]} um")
    if not any(abs(loading - x) < 1e-9 for x in dom["loading_levels"]):
        problems.append(f"loading {loading} not among the verified loading levels {dom['loading_levels']}")
    if not set(req_scen) <= set(dom["scenarios"]):
        problems.append(f"scenario set {req_scen} is not a subset of the verified scenarios {dom['scenarios']}")
    final_path = Path(out) / "final" / "final.json"
    if not final_path.exists():
        problems.append("no verified study results found; run `tcid all` first")
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    if problems:
        text = "UNSUPPORTED REQUEST\n  - " + "\n  - ".join(problems) + "\n" + DOMAIN_TEXT
        res = {"supported": False, "problems": problems, "text": text}
        if run_if_unsupported and final_path.exists() and material == dom["material"]:
            res.update(_run_outside_domain(study, t, loading, out_dir))
            res["text"] += "\n\n" + res["run_text"]
        (out_dir / "response.json").write_text(json.dumps(res, indent=2, default=str))
        return res

    fin = json.loads(final_path.read_text())
    cands = json.loads((Path(out) / "final" / "candidates.json").read_text())
    summ = pd.read_parquet(Path(out) / "final" / "final_summary.parquet")
    q1 = Path(out) / "questions" / "q1_map.parquet"
    q1m = pd.read_parquet(q1) if q1.exists() else None
    lo = fin["label_of"]
    conf = fin["confirmed"]
    base_lab = lo[f"{fin['recommended']}__S0"]
    base = {**fin["verified"][base_lab], "confirmed": conf.get(base_lab)}
    sep = fin["separability_from_recommended"]
    # recommended design first, then every optimised candidate re-measured on the confirmation batch
    picks = [{"label": lo[fin["recommended"]], "role": "recommended", "resolved": None,
              **fin["verified"][lo[fin["recommended"]]], "confirmed": conf.get(lo[fin["recommended"]])}]
    for key, v in sorted(sep.items(), key=lambda kv: kv[1]["diff"]):
        lab = lo[key]
        if lab in conf and lab not in [p["label"] for p in picks]:
            picks.append({"label": lab, "role": "alternative", "resolved": v["resolved"], **fin["verified"][lab],
                          "confirmed": conf[lab]})
    lines = ["Candidate distributions. Robust R = mean resistance in the worst stress scenario (mm2K/W).",
             f"Verification: {len(fin['selection_batch_seeds'])} seeds (used for selection); confirmation: "
             f"{len(fin['confirmation_seeds'])} independent seeds (quote these).", ""]
    lines.append(f"{'candidate':34s} {'d50 um':>7s} {'sigma':>6s} {'AR':>5s} {'S':>5s} {'verif.':>8s} {'confirmed':>15s} {'vs ref. base':>12s}")
    rows_csv = []
    bconf = base["confirmed"]["robust"] if base.get("confirmed") else base["robust"]
    for p in picks + [{"label": "baseline (random orientation)", "role": "baseline", **base}]:
        d = p["design"]
        c = p.get("confirmed")
        r_c = c["robust"] if c else float("nan")
        rel = 1 - r_c / bconf if c else float("nan")
        tag = {"recommended": " *", "alternative": "" if p.get("resolved") else " ~"}.get(p["role"], "")
        cstr = f"{r_c * 1e6:6.1f} ± {c['robust_se'] * 1e6:4.1f}" if c else "–"
        lines.append(f"{(p['label'] + tag)[:34]:34s} {d['d50'] * 1e6:7.2f} {d['sigma_ln']:6.2f} {d['aspect']:5.2f} {d['order']:5.2f} "
                     f"{p['robust'] * 1e6:8.1f} {cstr:>15s} {rel * 100:+10.1f}%")
        sc = summ[summ["design_id"] == p["design_id"]].set_index("scenario")["r_mean"]
        rows_csv.append({"candidate": p["label"], **d, "r_robust_verification_m2K_W": p["robust"],
                         "r_robust_confirmed_m2K_W": r_c, "r_robust_confirmed_se": c["robust_se"] if c else float("nan"),
                         "r_nominal_verification_m2K_W": p["nominal"], "confirmed_benefit_vs_baseline": rel,
                         **{f"r_{k}_verification": v for k, v in sc.items()}})
        if p["label"] != "baseline (random orientation)":
            safe = p["label"].split("+")[0]
            cfg = {
                "material": str(study.cfg.resolve(study.cfg.material_path)),
                "scenarios": str(study.cfg.resolve(study.cfg.scenarios_path)),
                "design": {**d, "d50": f"{d['d50'] * 1e6:.4f} um"},
                "seed": int(fin["seeds"][0]),
                "layer": {k: f"{v * 1e6:g} um" for k, v in study.cfg.layer.to_dict().items()},
                "mesh": {**study.cfg.mesh.to_dict(), "size_factor": fin["mesh_size_factor"]},
            }
            (out_dir / f"{safe}.yaml").write_text(yaml.safe_dump(cfg, sort_keys=False))
    pd.DataFrame(rows_csv).to_csv(out_dir / "comparison.csv", index=False)
    lines.append("")
    lines.append("vs ref. base = confirmed robust R relative to the recommended design's own random-orientation baseline.")
    lines.append("* recommended;  ~ difference from the recommended design not resolved (Holm-corrected Welch test on the")
    lines.append("  verification batch). Not resolved is not the same as equivalent; with no resolved difference, choose on")
    lines.append("  other grounds (cost, viscosity, supply).")
    lines.append("Sensitivity: worst-case values above are the maximum over the stress scenarios "
                 f"{', '.join(s.name for s in study.stress)} (stress cases, not probabilities).")
    if q1m is not None:
        sup = q1m[(q1m["r_interface"] >= 1e-8) & (q1m["r_interface"] <= 1e-6)]
        frac = 1 - sup["distinguishable"].mean()
        lines.append(f"Design map (Q1, 3 seeds per cell): within the supported interface range the most frequent leader is "
                     f"{sup['winner'].mode().iloc[0]}; {frac:.0%} of map cells cannot distinguish the top two.")
    lines.append(f"Configuration files written to {out_dir}/ (solve one with `tcid solve <file>`).")
    lines.append("")
    lines.append(DOMAIN_TEXT)
    text = "\n".join(lines)
    (out_dir / "domain.txt").write_text(DOMAIN_TEXT + "\n")
    res = {"supported": True, "text": text, "candidates": picks, "baseline": base}
    (out_dir / "response.json").write_text(json.dumps(res, indent=2, default=str))
    return res


def _run_outside_domain(study, t: float, loading: float, out_dir: Path) -> dict:
    """Direct solver runs of verified candidates at new conditions (no re-optimisation)."""
    cands = json.loads((study.out / "final" / "candidates.json").read_text())
    designs = {k: Design(**v).with_(loading=loading) for k, v in cands.items()
               if k in ("recommended", "baseline_random_orientation")}
    layer = replace(study.cfg.layer, thickness=t, width=2 * t)
    rows = study.run(list(designs.values()), [5000, 5001, 5002], f"request-t{t * 1e6:.0f}-phi{loading:.2f}", layer=layer)
    if rows.empty:
        return {"run_text": "New solver runs: every candidate was infeasible at the requested conditions."}
    from .provenance import config_hash

    lab = {config_hash(d.to_dict(), 12): k for k, d in designs.items()}
    rows["label"] = rows["design_id"].map(lab)
    g = rows.groupby(["label", "scenario"])["r_layer"].mean().unstack()
    g.to_csv(out_dir / "new_runs.csv")
    lines = [f"New solver runs at t = {t * 1e6:.0f} um, loading = {loading:.2f} (3 seeds; candidates NOT re-optimised):"]
    for k, r in g.iterrows():
        lines.append(f"  {k:30s} worst stress R = {r[[s.name for s in study.stress]].max() * 1e6:.1f} mm2K/W")
    missing = set(designs) - set(g.index)
    if missing:
        lines.append(f"  infeasible at these conditions: {sorted(missing)}")
    return {"run_text": "\n".join(lines), "new_runs": g.reset_index().to_dict(orient="records")}
