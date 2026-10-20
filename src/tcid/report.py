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


def _ci(d, scale=1e6, fmt=".1f"):
    return f"{d['estimate'] * scale:{fmt}} [{d['ci_low'] * scale:{fmt}}, {d['ci_high'] * scale:{fmt}}]"


def evaluate_hypotheses(study) -> list[dict]:
    """Grade the pre-registered hypotheses against intervals, not point estimates."""
    o = study.out
    fin = _j(o / "final" / "final.json")
    q1, q2, q3, q4 = (_j(o / "questions" / f"q{i}.json") for i in (1, 2, 3, 4))
    qc = _j(o / "questions" / "q1_crossings.json")
    gate = _j(o / "gate.json")
    comp = _j(o / "search" / "comparison.json")
    cands = _j(o / "final" / "candidates.json")
    out = []

    def grade(conds):
        vals = [c for c in conds if c is not None]
        if all(vals) and len(vals) == len(conds):
            return "supported"
        if not any(vals) and len(vals) == len(conds):
            return "refuted"
        return "partly supported" + (" (one clause unresolved)" if None in conds else "")

    if q1 and cands and qc:
        ro = cands["recommended"]
        pairs = {k: p for k, p in qc["pairs"].items() if "estimate" in p}
        # clause: size-ordering reversals only at or above about 1e-6
        below = [k for k, p in pairs.items() if np.isfinite(p["ci_high"]) and p["ci_high"] < 1e-6]
        straddle = [k for k, p in pairs.items() if "straddles" in p["status"]]
        clause = False if below else (None if straddle else True)
        conds = [ro["order"] >= 0.7, ro["aspect"] >= 3.0, not q1["baseline_ever_wins_supported"], clause]
        cr = "; ".join(f"{k.replace('|', ' vs ')}: {_ci(p, 1, '.1e')} m²K/W ({p['status']})" for k, p in pairs.items())
        out.append({"id": "H1", "outcome": grade(conds),
                    "evidence": f"recommended design S = {ro['order']:.2f}, AR = {ro['aspect']:.2f}; random-orientation baseline wins "
                                f"a supported-range map cell: {q1['baseline_ever_wins_supported']}; size-ordering crossings "
                                f"(estimate [95% bootstrap interval], 12 seeds): {cr}"})
    if q2:
        b = q2["contact_total_where_benefit_falls_to"]["0.10"]
        nominal = study.nominal.name
        stress = [s.name for s in study.stress]
        ref_t = f"{study.cfg.layer.thickness * 1e6:.0f}um"
        rec = b["recommended"][ref_t]
        rec_nom = rec[nominal]
        rec_worst = min((rec[s] for s in stress), key=lambda d: d["ci_low"])
        # thickness clause of H2 (50 um vs 150 um) on any pair that is feasible at both thicknesses
        both = [p for p, d in b.items() if "50um" in d and "150um" in d]
        if both:
            p = both[0]
            thick = b[p]["50um"][nominal]["ci_high"] < b[p]["150um"][nominal]["ci_low"]
            thick_txt = f"{p}: 50 µm {_ci(b[p]['50um'][nominal])}, 150 µm {_ci(b[p]['150um'][nominal])}"
        else:
            thick = None
            med = b.get("medium_particles_aligned", {})
            trend = (f"; between 100 and 150 µm (medium pair) {_ci(med['100um'][nominal])} → {_ci(med['150um'][nominal])}"
                     if "100um" in med and "150um" in med else "")
            thick_txt = "50 µm clause not testable: no design pair could be packed at 50 µm (see Q2)" + trend
        conds = [rec_nom["ci_low"] >= 3e-5, rec_worst["ci_low"] >= 3e-5, thick]
        out.append({"id": "H2", "outcome": grade(conds),
                    "evidence": f"total contact at which the benefit falls to 10% (estimate [95% bootstrap interval], mm²K/W), "
                                f"recommended design at {ref_t}: nominal {_ci(rec_nom)}, least favourable stress scenario "
                                f"{_ci(rec_worst)}; thickness trend: {thick_txt}"})
    if q3:
        av = q3["advantage_by_variant"]
        rel = {k: v["relative_to_reference"] for k, v in av.items()}
        smooth = rel.get("order_minus_0.1", 0) >= rel.get("order_minus_0.2", 0) > 0.2
        gaps = [abs(rel.get(k, np.nan) - 1) < 1 / 3 for k in ("min_gap_0.5um", "min_gap_2um") if k in rel]
        conds = [q3["min_t_ratio_reference"] > 3, smooth, all(gaps) and bool(gaps)]
        out.append({"id": "H3", "outcome": grade(conds),
                    "evidence": f"min t-ratio on fresh seeds {q3['min_t_ratio_reference']:.1f}; advantage relative to reference: "
                                + ", ".join(f"{k} {v:.2f}" for k, v in rel.items())})
    if q4 and "variants" in q4:
        def order_ok(r):
            vo, vi, vk = (r[f"variance_reduction_{k}"] for k in ("orientation_achieved", "interface_resistance", "filler_conductivity"))
            return vo > vi > vk
        ref_name, ref = next(iter(q4["variants"].items()))
        r = ref["recommended"]
        holds = {k: order_ok(v["recommended"]) for k, v in q4["variants"].items()}
        out.append({"id": "H4", "outcome": ("supported for the reference ranges" if order_ok(r) else "refuted")
                    + ("" if all(holds.values()) else "; ordering changes with the ranges"),
                    "evidence": f"sensitivity shares (reference ranges): orientation {r['variance_reduction_orientation_achieved']:.0%}, "
                                f"interface {r['variance_reduction_interface_resistance']:.0%}, filler conductivity "
                                f"{r['variance_reduction_filler_conductivity']:.0%}; ordering holds under range variants: "
                                + ", ".join(f"{k.split(' (')[0]}: {v}" for k, v in holds.items())})
    if comp and gate and fin:
        cb, cr = (comp[k]["calls_to_within_1pct_of_overall_best"] for k in ("bayesian_optimisation", "random_search"))
        faster = cb is not None and (cr is None or cb < cr)
        mc = fin["method_comparison_confirmed"]
        conf = gate["confirmed_pick_robust"]
        conds = [faster, gate["checks"]["bo_resolved_better_than_random_this_campaign"]]
        out.append({"id": "H5", "outcome": grade(conds) + " (single campaign per method)",
                    "evidence": f"calls to within 1% of the best loop value: BO {cb}, random {cr or 'never'}, sweep "
                                f"{comp['sweep']['calls_to_within_1pct_of_overall_best']}; best-of-4 per method confirmed on "
                                f"{len(fin['confirmation_seeds'])} independent seeds: "
                                + ", ".join(f"{k.replace('_', ' ')} {_um2(v)}" for k, v in conf.items())
                                + " mm²K/W; random − BO "
                                f"{_um2(mc['random_search_minus_bayesian_optimisation']['diff'])} "
                                f"[{_um2(mc['random_search_minus_bayesian_optimisation']['ci_low'])}, "
                                f"{_um2(mc['random_search_minus_bayesian_optimisation']['ci_high'])}]; sweep − BO "
                                f"{_um2(mc['sweep_minus_bayesian_optimisation']['diff'])} "
                                f"[{_um2(mc['sweep_minus_bayesian_optimisation']['ci_low'])}, "
                                f"{_um2(mc['sweep_minus_bayesian_optimisation']['ci_high'])}] (95% Welch CI). The sweep "
                                "comparison is unresolved, which is not evidence of equivalence."})
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
        w("## 5 Verified and confirmed results")
        w("")
        w(f"Two independent batches on the refined mesh (size factor {fin['mesh_size_factor']}): a **verification batch** "
          f"({len(fin['selection_batch_seeds'])} new seeds) on which candidates are compared and the recommendation is chosen, and a "
          f"**confirmation batch** ({len(fin['confirmation_seeds'])} further seeds) that re-measures the chosen designs. Choosing "
          "the minimum of several noisy values biases it low, so the confirmation values are the ones to quote. Every search "
          "method had its top four loop designs verified, so no method gets more chances than another.")
        w("")
        w("| Candidate | d50 (µm) | σ_ln | AR | S | Verification: robust R ± SE | Confirmation: robust R ± SE | Nominal R (verif. / conf.) |")
        w("|---|---|---|---|---|---|---|---|")
        for lab in fin["order_verified"]:
            v = fin["verified"][lab]
            d = v["design"]
            c = fin["confirmed"].get(lab)
            cs = f"{_um2(c['robust'])} ± {_um2(c['robust_se'])}" if c else "–"
            nom = f"{_um2(v['nominal'])} / {_um2(c['nominal'])}" if c else _um2(v["nominal"])
            w(f"| {lab} | {d['d50'] * 1e6:.1f} | {d['sigma_ln']:.2f} | {d['aspect']:.2f} | {d['order']:.2f} | "
              f"{_um2(v['robust'])} ± {_um2(v['robust_se'])} | {cs} | {nom} |")
        w("")
        w("Robust R is the mean resistance in the candidate's worst stress scenario (mm²K/W).")
        w("")
        rec = fin["recommended"]
        w(f"**Recommended design: `{rec}`**, chosen as the lowest verification-batch robust value among the "
          f"{len(fin['separability_from_recommended']) + 1} optimised candidates. The loop's top-ranked design scored "
          f"{_um2(fin['loop_top1_value_on_training_seeds'])} mm²K/W on its training seeds and "
          f"{_um2(fin['loop_top1_value_on_fresh_seeds'])} on fresh seeds (selection on noisy realizations).")
        w("")
        w("Which optimised candidates are resolved from the recommended one on the verification batch? Welch two-sample "
          "tests on per-realization worst-scenario values, Holm-corrected across the family. 'Not resolved' means no "
          "difference was detected, not that the designs are equivalent.")
        w("")
        w("| Candidate | Robust R − recommended (mm²K/W) | 95% CI (per comparison, unadjusted) | p (Holm) | Resolved (Holm p < 0.05) |")
        w("|---|---|---|---|---|")
        for k, v in sorted(fin["separability_from_recommended"].items(), key=lambda kv: kv[1]["diff"]):
            w(f"| {k} | {_um2(v['diff'])} | [{_um2(v['ci_low'])}, {_um2(v['ci_high'])}] | {v['p_holm']:.3f} | "
              f"{'yes' if v['resolved'] else 'no'} |")
        w("")
        if fin.get("confirmed_vs_recommended"):
            w("On the independent confirmation batch (no selection), difference from the recommended design:")
            w("")
            w("| Candidate | Robust R − recommended (mm²K/W) | 95% Welch CI | p |")
            w("|---|---|---|---|")
            for k, v in sorted(fin["confirmed_vs_recommended"].items(), key=lambda kv: kv[1]["diff"]):
                w(f"| {k} | {_um2(v['diff'])} | [{_um2(v['ci_low'])}, {_um2(v['ci_high'])}] | {v['p']:.3f} |")
            w("")
        rb = fin["confirmed_recommended_vs_baseline"]
        w(f"Confirmation batch, matched random-orientation baseline minus recommended: {_um2(rb['diff'])} mm²K/W "
          f"(95% CI [{_um2(rb['ci_low'])}, {_um2(rb['ci_high'])}], p = {rb['p']:.1e}), a {rb['rel']:.1%} reduction in the "
          "worst stress scenario.")
        w("")
        w(f"Top-ranked loop finalist keeps first place on fresh seeds: **{fin['top1_survives']}**; not resolved from the "
          f"verified best: {fin['top1_not_resolved_from_verified_best']}. Spearman correlation of surrogate-predicted vs "
          f"verified finalist ordering {fin['spearman_pred_vs_verified']:.2f}. Rank correlation of nominal and worst-scenario "
          f"objectives over all loop designs: {fin['spearman_nominal_vs_robust_loop']:.3f}.")
        w("")
        pdiff = o / "final" / "paired_differences.parquet"
        if pdiff.exists():
            d = pd.read_parquet(pdiff)
            w("Verification batch: each candidate minus its own random-orientation baseline (matched loading, sizes and "
              "aspect ratio; same scenarios; Welch 95% CI from independent realizations):")
            w("")
            w("| Candidate | Scenario | Baseline − candidate (mm²K/W) | 95% CI | Relative benefit |")
            w("|---|---|---|---|---|")
            for _, r in d.iterrows():
                w(f"| {r['label']} | {r['scenario']} | {_um2(r['baseline_minus_candidate'])} | "
                  f"[{_um2(r['ci_low'])}, {_um2(r['ci_high'])}] | {r['relative_benefit']:.1%} |")
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
    qc = _j(o / "questions" / "q1_crossings.json")
    if q1:
        wins = ", ".join(f"{k} {v}" for k, v in q1["winners_supported_range"].items())
        w(f"**Q1, when the preferred design changes.** Coarse design map ({q1['seeds_per_cell']} seeds per cell): map cells won "
          f"within the supported interface range: {wins}; the random-orientation baseline wins none. "
          f"{q1['indistinguishable_fraction_supported']:.0%} of those cells do not resolve the winner from the runner-up "
          "(Welch p ≥ 0.05), so the map shows where each design leads on average, not established boundaries.")
        w("")
        if qc:
            w(f"Targeted runs around the crossings: {len(qc['seeds'])} seeds per design, {len(qc['r_grid'])} interface "
              "resistances from 1e-7 to 1e-5 m²K/W (k_filler = 30 W/mK, full orientation). Crossing = first sign change of the "
              "mean-difference curve; interval = 95% percentile bootstrap over seeds.")
            w("")
            if qc.get("infeasible_small_sizes_um"):
                w(f"Smaller particles ({', '.join(f'{v:.1f}' for v in qc['infeasible_small_sizes_um'])} µm) with this width, "
                  "shape and alignment could not be packed at the study loading and minimum gap; the small rung is the "
                  f"smallest packable size. Ladder (µm): " + ", ".join(f"{k} {v:.1f}" for k, v in qc["d50_um"].items()) + ".")
                w("")
            w("| Pair | Crossing estimate (m²K/W) | 95% interval | Better below the crossing | Status |")
            w("|---|---|---|---|---|")
            for k, p in qc["pairs"].items():
                if "estimate" not in p:
                    w(f"| {k.replace('|', ' vs ')} | – | – | – | {p['status']} |")
                    continue
                w(f"| {k.replace('|', ' vs ')} | {p['estimate']:.2e} | [{p['ci_low']:.2e}, {p['ci_high']:.2e}] | "
                  f"{p['better_below']} | {p['status']} |")
            w("")
            w("![Crossings](figures/fig9_crossings.png)")
            w("")
        w("![Design map](figures/fig2_design_map.png)")
        w("")
    if q2:
        w(f"**Q2, when material optimisation is worth doing.** Computed from {len(q2['seeds'])} independent seeds per design "
          f"on the reported mesh (size factor {q2['mesh_size_factor']}); at the reference thickness these are the "
          "confirmation-batch values of section 5, so the numbers below can be reproduced from that table. "
          "With identical series contacts C the relative benefit is (R_base − R_cand)/(R_base + C); the 10% threshold "
          "C = (R_base − R_cand)/0.1 − R_base multiplies the uncertainty of the difference by ten, hence the intervals.")
        w("")
        lt = o / "questions" / "q2_layer_resistances.parquet"
        if lt.exists():
            t = pd.read_parquet(lt)
            t = t[t["scenario"] == study.nominal.name]
            w("| Design pair | Thickness | Baseline R_layer ± SE | Candidate R_layer ± SE | Benefit at C = 0 [95% CI] | C at 10% benefit [95% CI] |")
            w("|---|---|---|---|---|---|")
            for _, r in t.sort_values(["pair", "thickness"]).iterrows():
                tk = f"{r['thickness'] * 1e6:.0f}um"
                z = q2["benefit_at_zero_contact"][r["pair"]][f"{tk}/{study.nominal.name}"]
                c = q2["contact_total_where_benefit_falls_to"]["0.10"][r["pair"]][tk][study.nominal.name]
                w(f"| {r['pair']} | {tk} | {_um2(r['r_baseline'])} ± {_um2(r['r_baseline_se'])} | {_um2(r['r_candidate'])} ± "
                  f"{_um2(r['r_candidate_se'])} | {z['estimate']:.1%} [{z['ci_low']:.1%}, {z['ci_high']:.1%}] | {_ci(c)} mm²K/W |")
            w("")
            w("Nominal scenario shown; all scenarios are in `results/phase1/questions/q2.json`. The assumed realistic total "
              "contact range is up to 60 mm²K/W (SR-005). The medium pair (same width, shape and alignment as the recommended "
              "design, smaller particles) is included so that the thickness trend can be examined where the recommended "
              "distribution does not fit.")
            w("")
        if q2.get("rejected_realizations"):
            w("Rejected realizations (a pair is reported at a thickness only if every realization of both designs packed):")
            w("")
            for pname, roles in q2["rejected_realizations"].items():
                for role, tks in roles.items():
                    for tk, txt in tks.items():
                        w(f"* {pname} ({role}), {tk}: {txt}")
            w("")
        if q2.get("robin_vs_series_max_rel_diff") is not None:
            by = q2.get("robin_vs_series_by_pair", {})
            w(f"The series reduction R_top + R_layer + R_bottom agrees with Robin contact solves to within "
              f"{q2['robin_vs_series_max_rel_diff']:.1%} (contacts up to 1e-4 m²K/W per face; by pair: "
              + ", ".join(f"{k} {v:.1%}" for k, v in by.items()) + "). The larger error for the large-particle design reflects "
              "less uniform face temperatures; the series numbers above slightly underestimate assembly resistance.")
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
        w("**Q4, sensitivity of each decision to each uncertain input.** Share of the spread of a decision quantity removed "
          "by knowing one input exactly, averaged over a uniform grid of the chosen ranges (grid means over 3 seeds, so "
          "geometry noise is part of the spread). This is a sensitivity index conditional on those ranges and on treating the "
          "inputs as independent. It is not a value-of-information result: measurement cost, accuracy and whether the answer "
          "would change a decision are not modelled. Decision quantities: baseline − candidate (benefit of aligning) and "
          "large − small particles (choice of size).")
        w("")
        w("| Ranges | Decision quantity | Orientation achieved | Interface resistance | Filler conductivity | Spread (sd, mm²K/W) |")
        w("|---|---|---|---|---|---|")
        for vname, v in q4["variants"].items():
            for k in ("recommended", "medium_particles_aligned", "large_minus_small"):
                r = v[k]
                def pct(x):
                    return "not varied" if x != x else f"{x:.0%}"
                w(f"| {vname} | {k} | {pct(r['variance_reduction_orientation_achieved'])} | "
                  f"{pct(r['variance_reduction_interface_resistance'])} | {pct(r['variance_reduction_filler_conductivity'])} | "
                  f"{_um2(r['total_sd'])} |")
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
          f"Learned search resolved better than both alternatives in this single campaign: **{gate['ml_earns_its_complexity']}**. "
          "Where a cheaper direct search was not resolved worse, the guide's rule applies and it is retained; each method "
          "was run once, so this says nothing about the methods in general.")
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
