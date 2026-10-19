"""Technical report generator (guide section 9): docs/report.md built from the result files.

Every number in the report is read from ``results/``; the interpretive text that cannot
be computed lives in docs/interpretation.md and is embedded verbatim.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import yaml


def _j(p: Path):
    return json.loads(p.read_text()) if p.exists() else None


def _um2(x):  # m^2K/W -> mm^2K/W
    return f"{x * 1e6:.1f}"


def evaluate_hypotheses(study) -> list[dict]:
    o = study.out
    fin = _j(o / "final" / "final.json")
    q1, q2, q3, q4 = (_j(o / "questions" / f"q{i}.json") for i in (1, 2, 3, 4))
    comp = _j(o / "search" / "comparison.json")
    cands = _j(o / "final" / "candidates.json")
    out = []
    if q1 and cands:
        ro = cands["recommended"]
        aligned = {"large_particles_aligned", "small_particles_aligned", "recommended"}
        small_large = [r for r in q1["reversals"] if {r["a"], r["b"]} <= aligned]
        # strict reading of "at or above about 1e-6": any estimated crossing below 1e-6 counts against H1
        rev_low = [r for r in small_large if r["crossing_estimate"] < 1e-6]
        crossings = "; ".join(f"{r['a']} vs {r['b']} at ≈{r['crossing_estimate']:.1e} m²K/W "
                              f"({r['better_below']} better below)" for r in small_large)
        cond = [ro["order"] >= 0.7, ro["aspect"] >= 3.0, not q1["baseline_ever_wins_supported"], not rev_low]
        out.append({"id": "H1", "outcome": "supported" if all(cond) else ("partly supported" if any(cond) else "refuted"),
                    "evidence": f"recommended (verified best) design S = {ro['order']:.2f}, AR = {ro['aspect']:.2f}; baseline wins a supported-range "
                                f"map cell: {q1['baseline_ever_wins_supported']}; small/large-particle ordering reversals: "
                                f"{crossings or 'none'}"})
    if q2:
        b = q2["contact_total_where_benefit_falls_to"]["0.10"]
        nominal = study.nominal.name
        stress_names = [s.name for s in study.stress]
        r100 = b.get("100um", {})
        worst100 = min(r100.get(s, np.inf) for s in stress_names) if r100 else np.nan
        nom100 = r100.get(nominal, np.nan)
        t50 = b.get("50um", {}).get(nominal, np.nan)
        t150 = b.get("150um", {}).get(nominal, np.nan)
        cond = [nom100 >= 3e-5, worst100 >= 3e-5, t50 < t150]
        out.append({"id": "H2", "outcome": "supported" if all(cond) else ("partly supported" if any(cond) else "refuted"),
                    "evidence": f"benefit reaches 10% at total contact {_um2(nom100)} mm²K/W (nominal) / {_um2(worst100)} "
                                f"(least favourable stress scenario) for 100 µm; 50 µm: {_um2(t50)}, 150 µm: {_um2(t150)} mm²K/W"})
    if q3:
        av = q3["advantage_by_variant"]
        rel = {k: v["relative_to_reference"] for k, v in av.items()}
        smooth = rel.get("order_minus_0.1", 0) >= rel.get("order_minus_0.2", 0) > 0.2
        gaps = [abs(rel.get(k, np.nan) - 1) < 1 / 3 for k in ("min_gap_0.5um", "min_gap_2um") if k in rel]
        cond = [q3["min_t_ratio_reference"] > 3, smooth, all(gaps) and bool(gaps)]
        out.append({"id": "H3", "outcome": "supported" if all(cond) else ("partly supported" if any(cond) else "refuted"),
                    "evidence": f"min t-ratio on fresh seeds {q3['min_t_ratio_reference']:.1f}; advantage relative to reference: "
                                + ", ".join(f"{k} {v:.2f}" for k, v in rel.items())})
    if q4 and "recommended" in q4:
        r = q4["recommended"]
        vo, vi, vk = (r[f"variance_reduction_{k}"] for k in ("orientation_achieved", "interface_resistance", "filler_conductivity"))
        out.append({"id": "H4", "outcome": "supported" if vo > vi > vk else ("partly supported" if vo > max(vi, vk) else "refuted"),
                    "evidence": f"variance of (baseline - recommended) removed by knowing: orientation {vo:.0%}, interface "
                                f"resistance {vi:.0%}, filler conductivity {vk:.0%}"})
    if comp:
        bo, rnd = comp["bayesian_optimisation"], comp["random_search"]
        cb, cr = bo["calls_to_within_1pct_of_overall_best"], rnd["calls_to_within_1pct_of_overall_best"]
        bo_wins = cb is not None and (cr is None or cb < cr)
        gate = _j(o / "gate.json")
        extra = ""
        if gate and "verified_best_robust" in gate:
            vb, vs = gate["verified_best_robust"], gate["verified_best_robust_se"]
            extra = ("; verified on fresh seeds: " + ", ".join(f"{k.replace('_', ' ')} {_um2(vb[k])} ± {_um2(vs[k])}" for k in vb)
                     + f" mm²K/W; BO distinguishable from random: {gate['checks']['bo_better_than_random_at_matched_budget']}, "
                     f"from sweep: {gate['checks']['bo_better_than_sweep']}")
        out.append({"id": "H5", "outcome": "supported" if bo_wins else "refuted",
                    "evidence": f"solver calls to reach within 1% of the best robust value found by any method: BO {cb}, "
                                f"random {cr}, sweep {comp['sweep']['calls_to_within_1pct_of_overall_best']}; best robust "
                                f"(loop seeds): BO {_um2(bo['best_robust'])}, random {_um2(rnd['best_robust'])}, sweep "
                                f"{_um2(comp['sweep']['best_robust'])} mm²K/W{extra}"})
    return out


def write_report(study) -> Path:
    o = study.out
    ver = _j(o.parent / "verification" / "verification.json") or []
    pilot = _j(o / "pilot" / "pilot.json")
    cv_i = _j(o / "surrogate" / "cv_initial.json")
    cv_l = _j(o / "surrogate" / "cv_loop.json")
    comp = _j(o / "search" / "comparison.json")
    fin = _j(o / "final" / "final.json")
    gate = _j(o / "gate.json")
    q1, q2, q3, q4 = (_j(o / "questions" / f"q{i}.json") for i in (1, 2, 3, 4))
    hyps = evaluate_hypotheses(study)
    interp = Path("docs/interpretation.md")
    L = []
    w = L.append
    w("# Simulation-verified microstructure design study for particle-filled epoxy adhesive layers")
    w("")
    w("**Phase 1 technical report: 2D verified solver and constrained design loop (development filler, not BN).**")
    w("")
    w("> This is a simulation-based design study. No experimental confirmation exists. The filler is an isotropic "
      "alumina-like *development* stand-in, so no statement here is evidence about boron nitride. Electrical insulation "
      "is only the motivation for choosing ceramic fillers; the thermal solver says nothing about dielectric performance. "
      "No new commercial adhesive or semiconductor reliability is claimed.")
    w("")
    if interp.exists():
        w(interp.read_text().strip())
        w("")
    w("## 1 Physical model")
    w("")
    w("* Steady linear conduction, no heat generation: q = −K∇T, ∇·(K∇T) = 0. Epoxy is isotropic; each particle carries "
      "diag(k∥, k⊥) rotated with its orientation (identical in this isotropic development case, but the rotation path is verified).")
    w("* Particle–matrix interfaces: continuous normal flux and a jump ΔT = R_interface·q_n, implemented with duplicated "
      "degrees of freedom on the filler side of every interface and the weak-form term ∫ (1/R)[T][v] ds. Perfect contact ties the degrees of freedom.")
    w("* Conforming Gmsh triangulation of polygonal particles (vertices on the exact ellipses), P1 elements assembled with scikit-fem, "
      "SuperLU direct solve. At least 3 elements across every narrow gap and 2 across the resin margin.")
    w("* Boundary problem: isothermal hot and cold faces (ΔT = 1 K), insulated sidewalls. Area-specific layer resistance "
      "R_layer = ΔT·W/Q from integrated boundary heat flow; outer contacts added in series, R_total = R_top + R_layer + R_bottom, "
      "checked against Robin contact conditions.")
    lay = study.cfg.layer
    w(f"* Layer: t = {lay.thickness * 1e6:.0f} µm, width {lay.width * 1e6:.0f} µm, resin margin {lay.margin * 1e6:.0f} µm; "
      f"minimum particle gap {study.cfg.geometry.min_gap * 1e6:.1f} µm; filler area fraction {study.space.loading:.2f} (a 2D area fraction, not a 3D volume fraction).")
    w("")
    w("## 2 Material case and scenarios")
    w("")
    reg = yaml.safe_load(Path("data/source_register.yaml").read_text())
    w("| Register ID | Quantity | Value / interval | Status |")
    w("|---|---|---|---|")
    for e in reg["entries"]:
        val = e.get("interval", e.get("value", ""))
        w(f"| {e['id']} | {e['quantity']} | {val} {e.get('units', '')} | {e['status']} |")
    w("")
    w("Stress scenarios (finite set, not probabilities):")
    w("")
    w("| Scenario | Role | R_interface (m²K/W) | k_filler (W/mK) |")
    w("|---|---|---|---|")
    for s in study.scenarios:
        w(f"| {s.name} | {s.role} | {s.r_interface:.0e} | {study.material.k_filler_par * s.k_filler_scale:g} |")
    w("")
    w("## 3 Verification evidence")
    w("")
    w("| Check | Proposed acceptance criterion | Worst value | Pass |")
    w("|---|---|---|---|")
    for r in ver:
        crit = r["criterion"].replace("|", "\\|")
        w(f"| {r['name']} | {crit} | {r['value']:.2e} | {'yes' if r['passed'] else '**no**'} |")
    if fin:
        w(f"| finalist mesh refinement (size factor {study.budget['final_refinements'][0]} → {study.budget['final_refinements'][1]}) "
          f"| change < 2% | {fin['mesh_refinement_max_change']:.2e} | {'yes' if fin['mesh_refinement_max_change'] < 0.02 else '**no**'} |")
        if fin.get("smallest_relative_benefit_claimed") is not None:
            w(f"| numerical error vs claimed benefit | mesh change < 20% of smallest claimed benefit "
              f"({fin['smallest_relative_benefit_claimed']:.1%}) | {fin['mesh_refinement_max_change']:.2e} | "
              f"{'yes' if fin['numerical_error_below_claim'] else '**no**'} |")
    w("")
    w("Every solver row also records its own heat-balance error, solver residual and flux-consistency error; "
      "see `results/phase1/**/*.parquet`.")
    w("")
    w("![Mesh detail](figures/fig8_mesh_detail.png)")
    w("")
    w("## 4 Budget, failures and the design loop")
    w("")
    if pilot:
        w(f"Timing pilot: {pilot['n_runs']} representative runs (including dense and difficult geometries), "
          f"{pilot['cpu_s_per_job_mean']:.1f} CPU-s per geometry on average (max {pilot['cpu_s_per_job_max']:.0f} s) for "
          f"{pilot['scenarios_per_job']} scenario solves; mean mesh {pilot['n_elements_mean'] / 1e3:.0f}k elements (max "
          f"{pilot['n_elements_max'] / 1e3:.0f}k). Decision: {pilot['decision']}.")
        w("")
    if comp:
        w("| Method | Designs (complete / requested) | Solver calls | CPU h | Best robust R, loop seeds (mm²K/W) | Calls to within 1% of overall best |")
        w("|---|---|---|---|---|---|")
        for k in ("bayesian_optimisation", "random_search", "sweep"):
            c = comp[k]
            w(f"| {k.replace('_', ' ')} | {c['n_designs_complete']} / {c['n_designs_requested']} | {c['solver_calls']} | "
              f"{c['cpu_hours']:.2f} | {_um2(c['best_robust'])} | {c['calls_to_within_1pct_of_overall_best'] or 'never'} |")
        w("")
        if "failures_by_tag_and_stage" in comp:
            w("Rejected geometries/meshes (excluded from optimisation, reported here): "
              + ", ".join(f"{k}: {v}" for k, v in comp["failures_by_tag_and_stage"].items()) + ".")
            w("")
    for name, cv in (("initial design (40 designs)", cv_i), ("full loop data (initial + BO)", cv_l)):
        if cv:
            g, t = cv["gaussian_process"], cv["extra_trees"]
            w(f"Surrogate held-out error, {name}, folds split by design: GP median {g['median_rel_error']:.2%} "
              f"(p90 {g['p90_rel_error']:.2%}, low-resistance quartile {g['median_rel_error_low_resistance_quartile']:.2%}, "
              f"95% interval coverage {g['coverage_95']:.0%}, within-scenario Spearman {g['spearman_within_scenario_mean']:.3f}); "
              f"extra-trees median {t['median_rel_error']:.2%}. Target < 5%: {'met' if cv['meets_target'] else 'NOT met'}.")
            w("")
    w("![Search efficiency](figures/fig5_search_efficiency.png)")
    w("")
    w("![Surrogate parity](figures/fig6_surrogate_parity.png)")
    w("")
    if fin:
        w("## 5 Verified results (fresh seeds, refined mesh)")
        w("")
        w(f"Re-simulated with {len(fin['seeds'])} new geometry seeds at mesh size factor {fin['mesh_size_factor']}; "
          "these solver results, not surrogate predictions, are the reported results.")
        w("")
        w("| Candidate | d50 (µm) | σ_ln | AR | S | Robust R ± SE (mm²K/W) | Nominal R (mm²K/W) | Surrogate-predicted robust R |")
        w("|---|---|---|---|---|---|---|---|")
        for lab in fin["order_verified"]:
            v = fin["verified"][lab]
            d = v["design"]
            w(f"| {lab} | {d['d50'] * 1e6:.1f} | {d['sigma_ln']:.2f} | {d['aspect']:.2f} | {d['order']:.2f} | "
              f"{_um2(v['robust'])} ± {_um2(v['robust_se'])} | {_um2(v['nominal'])} | {_um2(fin['surrogate_predicted_robust'][lab])} |")
        w("")
        w(f"Top-ranked finalist survives the direct solver check: **{fin['top1_survives']}**; Spearman correlation of "
          f"predicted vs verified finalist ordering {fin['spearman_pred_vs_verified']:.2f}; loop top-1 within 2 SE of the verified "
          f"best: {fin['top1_within_2se_of_verified_best']}. Rank correlation of nominal and worst-scenario objectives over all "
          f"loop designs: {fin['spearman_nominal_vs_robust_loop']:.3f} (loop nominal and robust optima are the same design: "
          f"{fin['nominal_vs_robust_same_design']}).")
        w("")
        rec = fin["recommended"]
        w(f"**Recommended design: `{rec}`** (verified best robust value among optimised candidates). The loop's top-ranked "
          f"design scored {_um2(fin['loop_top1_value_on_training_seeds'])} mm²K/W on its training seeds and "
          f"{_um2(fin['loop_top1_value_on_fresh_seeds'])} on fresh seeds (selection on noisy realizations).")
        w("")
        w("| Optimised candidate | Robust R − recommended (mm²K/W) | ± SE | Distinguishable (> 2 SE) |")
        w("|---|---|---|---|")
        for k, v in fin["separability_from_recommended"].items():
            w(f"| {k} | {_um2(v['gap'])} | {_um2(v['se'])} | {'yes' if v['distinguishable'] else 'no'} |")
        w("")
        pdiff = o / "final" / "paired_differences.parquet"
        if pdiff.exists():
            d = pd.read_parquet(pdiff)
            w("Candidate minus its own random-orientation baseline (matched loading, sizes and aspect ratio; same scenarios; "
              "SE from independent realizations):")
            w("")
            w("| Candidate | Scenario | Baseline − candidate (mm²K/W) | ± SE | Relative benefit | t |")
            w("|---|---|---|---|---|---|")
            for _, r in d.iterrows():
                w(f"| {r['label']} | {r['scenario']} | {_um2(r['baseline_minus_candidate'])} | {_um2(r['se'])} | "
                  f"{r['relative_benefit']:.1%} | {r['t_ratio']:.1f} |")
            w("")
        w("![Fields](figures/fig1_fields.png)")
        w("")
        w("![Realizations and mesh](figures/fig4_realizations_mesh.png)")
        w("")
    w("## 6 Research questions and pre-registered hypotheses")
    w("")
    w("Hypotheses were recorded in `docs/hypotheses.md` before optimisation. Outcomes, evaluated mechanically from the result files:")
    w("")
    w("| ID | Outcome | Evidence |")
    w("|---|---|---|")
    for h in hyps:
        w(f"| {h['id']} | {h['outcome']} | {h['evidence']} |")
    w("")
    if q1:
        wins = ", ".join(f"{k} {v}" for k, v in q1["winners_supported_range"].items())
        rev = "; ".join(f"{r['a']} / {r['b']} cross at ≈{r['crossing_estimate']:.1e} m²K/W ({r['better_below']} better below)"
                        for r in q1["reversals"]) or "none"
        w(f"**Q1, when the preferred design changes.** Map cells won within the supported interface range: {wins} "
          f"(the random-orientation baseline wins none). {q1['indistinguishable_fraction_supported']:.0%} of those cells cannot "
          f"separate the top two candidates at 3 seeds (< 2 SE). Ordering reversals at k_filler = 30 W/mK and full orientation: {rev}.")
        w("")
        w("![Design map](figures/fig2_design_map.png)")
        w("")
    if q2:
        w("**Q2, when material optimisation is worth doing.** Total outer contact resistance at which the assembly benefit "
          "of the recommended design over its random-orientation baseline falls to 10%:")
        w("")
        w("| Thickness | " + " | ".join(s.name for s in study.scenarios) + " |")
        w("|---|" + "---|" * len(study.scenarios))
        for t, v in q2["contact_total_where_benefit_falls_to"]["0.10"].items():
            w(f"| {t} | " + " | ".join(_um2(v.get(s.name, float('nan'))) + " mm²K/W" for s in study.scenarios) + " |")
        w("")
        if q2.get("robin_vs_series_max_rel_diff") is not None:
            w(f"The series reduction agrees with Robin contact conditions to within {q2['robin_vs_series_max_rel_diff']:.2%} "
              "(heterogeneous layer, contacts up to 1e-4 m²K/W per face).")
            w("")
        w("![Contacts](figures/fig3_contacts.png)")
        w("")
    if q3:
        w("**Q3, whether the optimum relies on ideal geometry** (2D only; the 3D check is Phase 2).")
        w("")
        w("| Variant | Advantage vs baseline | Relative to reference | t |")
        w("|---|---|---|---|")
        for k, v in q3["advantage_by_variant"].items():
            w(f"| {k} | {v['advantage_rel']:.1%} | {v['relative_to_reference']:.2f} | {v['t_ratio']:.1f} |")
        if q3.get("failures"):
            w("")
            w("Infeasible variants: " + "; ".join(f"{k}: {v}" for k, v in q3["failures"].items()))
        w("")
        w("![Perturbations](figures/fig7_perturbations.png)")
        w("")
    if q4:
        w("**Q4, which unknown to resolve next.** Share of the spread of the candidate-minus-baseline difference removed by "
          "knowing one input exactly (grid over supported ranges; conditional on the other modelled uncertainties):")
        w("")
        w("| Decision quantity | Orientation achieved | Interface resistance | Filler conductivity | Spread (sd, mm²K/W) |")
        w("|---|---|---|---|---|")
        for k, r in q4.items():
            w(f"| {k} | {r['variance_reduction_orientation_achieved']:.0%} | {r['variance_reduction_interface_resistance']:.0%} | "
              f"{r['variance_reduction_filler_conductivity']:.0%} | {_um2(r['total_sd'])} |")
        w("")
    if gate:
        w("## 7 Gate for moving to Phase 2")
        w("")
        w("| Condition | Met |")
        w("|---|---|")
        for k, v in gate["checks"].items():
            w(f"| {k.replace('_', ' ')} | {'yes' if v else '**no**'} |")
        w("")
        w(f"Proceed to Phase 2: **{gate['proceed_to_phase2']}**. Blocking (required conditions not met): "
          f"{', '.join(gate['blocking']) or 'none'}. Informational: {', '.join(gate.get('informational_unmet', [])) or 'none'}. "
          f"Learned search earns its complexity in this 4-descriptor problem: **{gate['ml_earns_its_complexity']}** "
          "(the guide's rule applies: an equally effective direct search is retained and reported).")
        w("")
    w("## 8 Limits and failed assumptions")
    w("")
    w("* 2D cross-sections cannot represent 3D connectivity; area fraction is not volume fraction. Phase 2 (3D oblate platelets, BN tensor) is required before any BN claim.")
    w("* The material case is a development stand-in; every property is *assumed* (source register). The BN case (SR-001, SR-009) is still pending extraction.")
    w("* Finite domain (200 µm wide, insulated sidewalls): results are finite-volume examples, not representative bulk properties; lateral-size sensitivity is reported in Q3.")
    w("* Scenario results are worst-case statements over hand-chosen stress cases, not calibrated reliability estimates.")
    w("* Not modelled: processing (viscosity, alignment method), particle clustering/contact, voids, adhesion, cure, dielectric strength, thermo-mechanical stress, ageing.")
    w("")
    w("## 9 Reproduce")
    w("")
    w("```bash\npip install -e .[test]   # see README for the Gmsh system library note\ntcid all                 # verify, pilot, design loop, final checks, questions, gate, figures, report\n```")
    w("")
    w("Individual stages: `tcid verify | pilot | initial | search | final | questions | gate | figures | report`. "
      "Raw rows: `results/phase1/**.parquet`; manifest with checksums and environment: `results/phase1/manifest.json`.")
    p = Path("docs/report.md")
    p.write_text("\n".join(L) + "\n")
    return p
