# Simulation-verified microstructure design study for particle-filled epoxy adhesive layers

**Phase 1 technical report: 2D verified solver and constrained design loop (development filler, not BN).**

> This is a simulation-based design study. No experimental confirmation exists. The filler is an isotropic alumina-like *development* stand-in, so no statement here is evidence about boron nitride. Electrical insulation is only the motivation for choosing ceramic fillers; the thermal solver says nothing about dielectric performance. No new commercial adhesive or semiconductor reliability is claimed.

## Summary of findings (Phase 1, 2D development case)

All numbers below come from the result files and are repeated in the tables later in this report.
They describe a 2D, 100 µm, area-fraction-0.30 layer with an *assumed* alumina-like filler.
They are not BN results.

1. **Alignment and elongation give a robust, well-resolved benefit.**
   * The recommended distribution has median diameter 6.7 µm, log-normal width 0.42, aspect ratio 4 and orientation order 0.87 about the heat-flow direction.
   * Verified on 6 fresh realizations with a refined mesh, its worst-scenario layer resistance is 180.8 ± 1.2 mm²K/W.
   * Its own random-orientation counterpart (same loading, sizes and shape) gives 232.8 mm²K/W. That is a 22–23% reduction in every scenario, with t ≈ 15.
   * Conventional random circles at the same loading reach 279.9 mm²K/W.
   * The benefit survives fresh seeds, a ±20% change in particle size, a doubled lateral domain, a halved or doubled minimum gap, and a halved or doubled resin margin. Its relative size stays between 0.93 and 1.22 of the reference.
   * A wider size distribution (width +0.1) keeps 89% of the benefit.
   * It is most sensitive to the orientation actually achieved: an orientation order 0.2 lower keeps 74% of the benefit.
2. **The evidence cannot choose the particle size.**
   * Every leading candidate sits at the upper bound of the supported aspect ratio (4) with order 0.84–0.95. Their median diameters span 6.7–21 µm.
   * Their verified worst-case resistances span only 180.8–186.9 mm²K/W. Three of the five optimised candidates are within 2 standard errors of the recommended one.
   * The tool therefore reports that size cannot be decided from the current evidence, rather than claiming a size optimum.
   * The optimum on the aspect-ratio bound means the search wants more elongated particles than the Phase 1 envelope admits. Whether that is achievable is a geometry and processing question for Phase 2.
3. **The optimiser's own winner did not survive re-simulation, which is selection bias, not a solver error.**
   * The loop's top-ranked design used large, nearly monodisperse particles (21 µm equivalent diameter). Their long axes are about 42 µm, almost half the 100 µm layer.
   * It scored 175.7 mm²K/W on its three training seeds and 186.9 on six fresh seeds.
   * Those designs vary by up to 6% between realizations, so ranking on 3 seeds favours lucky draws.
   * The verified best came from the same loop (its nominal optimum). Its ordering survives, and the loop's winner is statistically indistinguishable from it.
   * Phase 2 shortlists should use at least 6 seeds, or a variance-penalised objective.
4. **Interface uncertainty within the supported range does not change the decision. Interface knowledge would decide the size.** (Q1, Q4)
   * Random orientation never wins any cell of the design map, even at 1e-5 m²K/W.
   * Among aligned designs the size ordering changes with interface resistance:
     * smaller particles (4 µm) are better below about 6e-7 m²K/W;
     * the recommended 6.7 µm design is better between about 6e-7 and 1.5e-6;
     * larger particles (20 µm) are better above about 1.5e-6, because they carry less interface area per unit filler.
   * Nominal and worst-scenario objectives rank the loop designs almost identically (Spearman 0.992). Accounting for the stress scenarios therefore did not change the recommendation.
   * Q4 splits the uncertainty by decision:
     * For "align or not", knowing the achieved orientation removes 98% of the spread in the benefit, and interface resistance only 2%. Characterise achieved orientation first.
     * For "which size", interface resistance explains 92% of the spread. Measure it next if particle size must be chosen.
5. **Outer contacts do not make the microstructure irrelevant in the assumed range.** (Q2)
   * With identical series contacts the ordering cannot change, so only the size of the benefit matters.
   * For a 100 µm layer the assembly benefit is 25% with ideal contacts. It falls to 10% only when the two contacts together reach about 330 mm²K/W, roughly twice the layer resistance itself and far above the assumed contact range (≤ 60 mm²K/W in total).
   * The threshold scales with thickness: about 99 mm²K/W at 50 µm and about 463 mm²K/W at 150 µm.
   * Thin layers also gain less at zero contact (19% at 50 µm), because particles are no longer small relative to the thickness.
   * Robin contact solves agree with the series reduction to within 0.6%.
6. **The learned search did not earn its complexity here.** (H5)
   * Bayesian optimisation found better designs than matched-budget random search: verified 180.8 ± 1.2 vs 196.2 ± 3.2 mm²K/W.
   * A plain 3-level factorial sweep was statistically equivalent (183.1 ± 2.6). It also got within 1% of the best loop value in 585 solver calls versus 1020 for BO.
   * The GP surrogate met the 5% held-out target (4.1% median), but its intervals were too narrow (85% coverage of nominal 95%). A small extra-trees ensemble was more accurate (1.8%) and better calibrated.
   * Following the guide, the simple sweep is retained as the Phase 1 search. BO may pay off in Phase 2, where 3D solves are far more expensive and more descriptors are free.

### Phase 2 gate

The gate is **not passed**, for two reasons:

* the BN material case and BN tensor (source register SR-001 and SR-009) are still pending extraction;
* the loop's top-ranked finalist did not keep first place on fresh seeds, although it is within 2 standard errors of the verified best.

The solver gates, the held-out error target and the numerical-error-versus-claim check all pass.

### Notes on the hypothesis evaluation

* H1's clause "reversal only at or above about 1e-6 m²K/W" is graded strictly. The recommended-versus-small-particle crossing is estimated at 6e-7 m²K/W, inside the supported range, so H1 is only *partly supported*. Its other clauses (strong alignment, high aspect ratio, no reversal against the random baseline) hold.
* The scenarios are hand-chosen stress cases, so every "worst case" above is a stress statement, not a reliability estimate.

## 1 Physical model

* Steady linear conduction, no heat generation: q = −K∇T, ∇·(K∇T) = 0. Epoxy is isotropic; each particle carries diag(k∥, k⊥) rotated with its orientation (identical in this isotropic development case, but the rotation path is verified).
* Particle–matrix interfaces: continuous normal flux and a jump ΔT = R_interface·q_n, implemented with duplicated degrees of freedom on the filler side of every interface and the weak-form term ∫ (1/R)[T][v] ds. Perfect contact ties the degrees of freedom.
* Conforming Gmsh triangulation of polygonal particles (vertices on the exact ellipses), P1 elements assembled with scikit-fem, SuperLU direct solve. At least 3 elements across every narrow gap and 2 across the resin margin.
* Boundary problem: isothermal hot and cold faces (ΔT = 1 K), insulated sidewalls. Area-specific layer resistance R_layer = ΔT·W/Q from integrated boundary heat flow; outer contacts added in series, R_total = R_top + R_layer + R_bottom, checked against Robin contact conditions.
* Layer: t = 100 µm, width 200 µm, resin margin 3 µm; minimum particle gap 1.0 µm; filler area fraction 0.30 (a 2D area fraction, not a 3D volume fraction).

## 2 Material case and scenarios

| Register ID | Quantity | Value / interval | Status |
|---|---|---|---|
| SR-001 | BN-epoxy formulation, cure, BN grades and through-plane conductivity data |   | pending |
| SR-002 | epoxy matrix thermal conductivity | [0.17, 0.23] W/(m K) | assumed |
| SR-003 | alumina-like isotropic filler conductivity (development stand-in) | [20, 35] W/(m K) | assumed |
| SR-004 | particle-matrix area-specific interface resistance (effective) | [1e-08, 1e-06] m^2 K/W | assumed |
| SR-005 | outer adhesive-to-substrate contact resistance (each face) | [0, 3e-05] m^2 K/W | assumed |
| SR-006 | densities for weight/volume conversion | {'epoxy': 1.2, 'alumina': 3.95, 'hBN': '2.1-2.28'} g/cm^3 | assumed |
| SR-007 | adhesive layer thickness | ['50e-6', '150e-6'] m | assumed |
| SR-008 | ceramic composite microstructure and conductivity dataset |   | pending |
| SR-009 | hBN platelet conductivity tensor (in-plane / through-plane) |   | pending |

Stress scenarios (finite set, not probabilities):

| Scenario | Role | R_interface (m²K/W) | k_filler (W/mK) |
|---|---|---|---|
| nominal | nominal | 1e-07 | 30 |
| good_iface_low_kf | stress | 1e-08 | 20 |
| good_iface_high_kf | stress | 1e-08 | 35 |
| poor_iface_low_kf | stress | 1e-06 | 20 |
| poor_iface_high_kf | stress | 1e-06 | 35 |

## 3 Verification evidence

| Check | Proposed acceptance criterion | Worst value | Pass |
|---|---|---|---|
| homogeneous_layer | R = t/k within 0.5% | 1.76e-15 | yes |
| layered_interfaces | summed layer + interface resistances within 1% | 1.01e-13 | yes |
| energy_balance | \|Q_in - Q_out\| / Q < 0.1% | 8.64e-12 | yes |
| mesh_refinement | resistance change < 2% on the next refinement | 2.50e-04 | yes |
| orientation_implementation | rigid rotation exact (<1e-8); 90-degree re-meshed rotation within 1%; un-rotated tensor detected | 1.26e-05 | yes |
| geometry_validity | no overlap or margin violation in accepted geometries; achieved loading within tolerance and recorded | 0.00e+00 | yes |
| limits_and_bounds | pure-phase exact; monotone in R_interface; R->inf matches holes within 1%; inside Wiener bounds | 1.34e-14 | yes |
| finalist mesh refinement (size factor 0.7 → 0.5) | change < 2% | 4.84e-04 | yes |
| numerical error vs claimed benefit | mesh change < 20% of smallest claimed benefit (22.0%) | 4.84e-04 | yes |

Every solver row also records its own heat-balance error, solver residual and flux-consistency error; see `results/phase1/**/*.parquet`.

![Mesh detail](figures/fig8_mesh_detail.png)

## 4 Budget, failures and the design loop

Timing pilot: 20 representative runs (including dense and difficult geometries), 11.3 CPU-s per geometry on average (max 56 s) for 5 scenario solves; mean mesh 141k elements (max 615k). Decision: proceed with planned grid.

| Method | Designs (complete / requested) | Solver calls | CPU h | Best robust R, loop seeds (mm²K/W) | Calls to within 1% of overall best |
|---|---|---|---|---|---|
| bayesian optimisation | 77 / 80 | 1155 | 0.40 | 175.7 | 1020 |
| random search | 80 / 80 | 1200 | 0.54 | 193.0 | never |
| sweep | 72 / 81 | 1080 | 0.71 | 176.5 | 585 |

Rejected geometries/meshes (excluded from optimisation, reported here): bo/geometry: 9, sweep/geometry: 25.

Surrogate held-out error, initial design (40 designs), folds split by design: GP median 3.93% (p90 6.37%, low-resistance quartile 4.33%, 95% interval coverage 62%, within-scenario Spearman 0.901); extra-trees median 2.04%. Target < 5%: met.

Surrogate held-out error, full loop data (initial + BO), folds split by design: GP median 4.09% (p90 12.78%, low-resistance quartile 3.41%, 95% interval coverage 85%, within-scenario Spearman 0.892); extra-trees median 1.83%. Target < 5%: met.

![Search efficiency](figures/fig5_search_efficiency.png)

![Surrogate parity](figures/fig6_surrogate_parity.png)

## 5 Verified results (fresh seeds, refined mesh)

Re-simulated with 6 new geometry seeds at mesh size factor 0.7; these solver results, not surrogate predictions, are the reported results.

| Candidate | d50 (µm) | σ_ln | AR | S | Robust R ± SE (mm²K/W) | Nominal R (mm²K/W) | Surrogate-predicted robust R |
|---|---|---|---|---|---|---|---|
| nominal_opt | 6.7 | 0.42 | 4.00 | 0.87 | 180.8 ± 1.2 | 167.5 | 180.1 |
| sweep_best | 9.8 | 0.28 | 4.00 | 0.95 | 183.1 ± 2.6 | 172.1 | 176.8 |
| shortlist_2 | 9.6 | 0.40 | 4.00 | 0.95 | 183.6 ± 1.0 | 172.7 | 178.3 |
| shortlist_3 | 16.6 | 0.17 | 4.00 | 0.95 | 185.5 ± 1.8 | 177.1 | 179.0 |
| robust_opt | 21.1 | 0.05 | 4.00 | 0.84 | 186.9 ± 3.3 | 179.4 | 176.0 |
| random_search_best | 17.3 | 0.22 | 3.94 | 0.56 | 196.2 ± 3.2 | 187.5 | 219.2 |
| shortlist_3__S0 | 16.6 | 0.17 | 4.00 | 0.00 | 230.8 ± 5.1 | 220.0 | 236.7 |
| shortlist_2__S0 | 9.6 | 0.40 | 4.00 | 0.00 | 231.4 ± 3.1 | 217.8 | 234.1 |
| sweep_best__S0 | 9.8 | 0.28 | 4.00 | 0.00 | 232.1 ± 3.3 | 217.0 | 237.9 |
| nominal_opt__S0 | 6.7 | 0.42 | 4.00 | 0.00 | 232.8 ± 3.0 | 215.7 | 241.2 |
| robust_opt__S0 | 21.1 | 0.05 | 4.00 | 0.00 | 233.8 ± 7.4 | 224.5 | 235.0 |
| random_search_best__S0 | 17.3 | 0.22 | 3.94 | 0.00 | 237.6 ± 6.4 | 226.9 | 235.6 |
| conventional_circles | 10.0 | 0.30 | 1.00 | 0.00 | 279.9 ± 1.5 | 267.0 | 293.3 |

Top-ranked finalist survives the direct solver check: **False**; Spearman correlation of predicted vs verified finalist ordering 0.81; loop top-1 within 2 SE of the verified best: True. Rank correlation of nominal and worst-scenario objectives over all loop designs: 0.992 (loop nominal and robust optima are the same design: False).

**Recommended design: `nominal_opt`** (verified best robust value among optimised candidates). The loop's top-ranked design scored 175.7 mm²K/W on its training seeds and 186.9 on fresh seeds (selection on noisy realizations).

| Optimised candidate | Robust R − recommended (mm²K/W) | ± SE | Distinguishable (> 2 SE) |
|---|---|---|---|
| robust_opt | 6.1 | 3.5 | no |
| shortlist_2 | 2.8 | 1.5 | no |
| shortlist_3 | 4.7 | 2.2 | yes |
| nominal_opt | 0.0 | 1.7 | no |
| random_search_best | 15.4 | 3.4 | yes |
| sweep_best | 2.3 | 2.9 | no |

Candidate minus its own random-orientation baseline (matched loading, sizes and aspect ratio; same scenarios; SE from independent realizations):

| Candidate | Scenario | Baseline − candidate (mm²K/W) | ± SE | Relative benefit | t |
|---|---|---|---|---|---|
| robust_opt | good_iface_high_kf | 45.0 | 8.1 | 20.2% | 5.6 |
| robust_opt | good_iface_low_kf | 44.7 | 7.9 | 19.7% | 5.6 |
| robust_opt | nominal | 45.1 | 8.1 | 20.1% | 5.6 |
| robust_opt | poor_iface_high_kf | 47.3 | 8.3 | 20.5% | 5.7 |
| robust_opt | poor_iface_low_kf | 46.9 | 8.1 | 20.1% | 5.8 |
| shortlist_2 | good_iface_high_kf | 44.9 | 3.5 | 20.8% | 13.0 |
| shortlist_2 | good_iface_low_kf | 44.7 | 3.4 | 20.4% | 13.2 |
| shortlist_2 | nominal | 45.1 | 3.4 | 20.7% | 13.2 |
| shortlist_2 | poor_iface_high_kf | 48.0 | 3.3 | 21.0% | 14.5 |
| shortlist_2 | poor_iface_low_kf | 47.9 | 3.2 | 20.7% | 14.7 |
| shortlist_3 | good_iface_high_kf | 42.7 | 5.5 | 19.5% | 7.8 |
| shortlist_3 | good_iface_low_kf | 42.4 | 5.4 | 19.1% | 7.9 |
| shortlist_3 | nominal | 42.9 | 5.4 | 19.5% | 7.9 |
| shortlist_3 | poor_iface_high_kf | 45.6 | 5.5 | 20.0% | 8.3 |
| shortlist_3 | poor_iface_low_kf | 45.3 | 5.4 | 19.6% | 8.4 |
| nominal_opt | good_iface_high_kf | 47.9 | 3.1 | 22.5% | 15.3 |
| nominal_opt | good_iface_low_kf | 47.7 | 3.1 | 22.0% | 15.4 |
| nominal_opt | nominal | 48.3 | 3.1 | 22.4% | 15.4 |
| nominal_opt | poor_iface_high_kf | 52.3 | 3.3 | 22.8% | 16.1 |
| nominal_opt | poor_iface_low_kf | 52.0 | 3.2 | 22.4% | 16.2 |
| random_search_best | good_iface_high_kf | 39.2 | 7.4 | 17.4% | 5.3 |
| random_search_best | good_iface_low_kf | 38.9 | 7.3 | 17.0% | 5.4 |
| random_search_best | nominal | 39.3 | 7.4 | 17.3% | 5.3 |
| random_search_best | poor_iface_high_kf | 41.6 | 7.3 | 17.7% | 5.7 |
| random_search_best | poor_iface_low_kf | 41.3 | 7.2 | 17.4% | 5.8 |
| sweep_best | good_iface_high_kf | 44.5 | 4.3 | 20.7% | 10.3 |
| sweep_best | good_iface_low_kf | 44.4 | 4.2 | 20.4% | 10.6 |
| sweep_best | nominal | 44.9 | 4.3 | 20.7% | 10.5 |
| sweep_best | poor_iface_high_kf | 49.1 | 4.3 | 21.4% | 11.3 |
| sweep_best | poor_iface_low_kf | 49.0 | 4.2 | 21.1% | 11.5 |

![Fields](figures/fig1_fields.png)

![Realizations and mesh](figures/fig4_realizations_mesh.png)

## 6 Research questions and pre-registered hypotheses

Hypotheses were recorded in `docs/hypotheses.md` before optimisation. Outcomes, evaluated mechanically from the result files:

| ID | Outcome | Evidence |
|---|---|---|
| H1 | partly supported | recommended (verified best) design S = 0.87, AR = 4.00; baseline wins a supported-range map cell: False; small/large-particle ordering reversals: large_particles_aligned vs recommended at ≈1.5e-06 m²K/W (recommended better below); large_particles_aligned vs small_particles_aligned at ≈1.1e-06 m²K/W (small_particles_aligned better below); recommended vs small_particles_aligned at ≈6.0e-07 m²K/W (small_particles_aligned better below) |
| H2 | supported | benefit reaches 10% at total contact 330.9 mm²K/W (nominal) / 322.0 (least favourable stress scenario) for 100 µm; 50 µm: 99.0, 150 µm: 463.2 mm²K/W |
| H3 | supported | min t-ratio on fresh seeds 12.4; advantage relative to reference: d50_x0.8 1.09, d50_x1.2 1.22, margin_1.5um 1.07, margin_6um 0.96, min_gap_0.5um 0.95, min_gap_2um 1.03, order_minus_0.1 0.89, order_minus_0.2 0.74, reference 1.00, sigma_plus_0.1 0.89, width_x2 0.93 |
| H4 | supported | variance of (baseline - recommended) removed by knowing: orientation 98%, interface resistance 2%, filler conductivity 0% |
| H5 | supported | solver calls to reach within 1% of the best robust value found by any method: BO 1020, random None, sweep 585; best robust (loop seeds): BO 175.7, random 193.0, sweep 176.5 mm²K/W; verified on fresh seeds: bayesian optimisation 180.8 ± 1.2, random search 196.2 ± 3.2, sweep 183.1 ± 2.6 mm²K/W; BO distinguishable from random: True, from sweep: False |

**Q1, when the preferred design changes.** Map cells won within the supported interface range: small_particles_aligned 24, recommended 21 (the random-orientation baseline wins none). 100% of those cells cannot separate the top two candidates at 3 seeds (< 2 SE). Ordering reversals at k_filler = 30 W/mK and full orientation: large_particles_aligned / recommended cross at ≈1.5e-06 m²K/W (recommended better below); large_particles_aligned / small_particles_aligned cross at ≈1.1e-06 m²K/W (small_particles_aligned better below); recommended / small_particles_aligned cross at ≈6.0e-07 m²K/W (small_particles_aligned better below).

![Design map](figures/fig2_design_map.png)

**Q2, when material optimisation is worth doing.** Total outer contact resistance at which the assembly benefit of the recommended design over its random-orientation baseline falls to 10%:

| Thickness | nominal | good_iface_low_kf | good_iface_high_kf | poor_iface_low_kf | poor_iface_high_kf |
|---|---|---|---|---|---|
| 50um | 99.0 mm²K/W | 94.8 mm²K/W | 98.2 mm²K/W | 112.0 mm²K/W | 115.3 mm²K/W |
| 100um | 330.9 mm²K/W | 322.0 mm²K/W | 328.9 mm²K/W | 358.9 mm²K/W | 365.9 mm²K/W |
| 150um | 463.2 mm²K/W | 451.4 mm²K/W | 461.2 mm²K/W | 495.2 mm²K/W | 505.1 mm²K/W |

The series reduction agrees with Robin contact conditions to within 0.58% (heterogeneous layer, contacts up to 1e-4 m²K/W per face).

![Contacts](figures/fig3_contacts.png)

**Q3, whether the optimum relies on ideal geometry** (2D only; the 3D check is Phase 2).

| Variant | Advantage vs baseline | Relative to reference | t |
|---|---|---|---|
| d50_x0.8 | 23.6% | 1.09 | 25.1 |
| d50_x1.2 | 25.6% | 1.22 | 9.5 |
| margin_1.5um | 23.7% | 1.07 | 13.8 |
| margin_6um | 20.3% | 0.96 | 11.9 |
| min_gap_0.5um | 20.9% | 0.95 | 10.1 |
| min_gap_2um | 22.5% | 1.03 | 12.7 |
| order_minus_0.1 | 19.3% | 0.89 | 11.1 |
| order_minus_0.2 | 16.1% | 0.74 | 5.9 |
| reference | 21.8% | 1.00 | 13.1 |
| sigma_plus_0.1 | 19.9% | 0.89 | 13.6 |
| width_x2 | 20.3% | 0.93 | 12.2 |

![Perturbations](figures/fig7_perturbations.png)

**Q4, which unknown to resolve next.** Share of the spread of the candidate-minus-baseline difference removed by knowing one input exactly (grid over supported ranges; conditional on the other modelled uncertainties):

| Decision quantity | Orientation achieved | Interface resistance | Filler conductivity | Spread (sd, mm²K/W) |
|---|---|---|---|---|
| recommended | 98% | 2% | 0% | 8.5 |
| small_particles_aligned | 97% | 2% | 0% | 11.8 |
| large_particles_aligned | 89% | 11% | 0% | 11.8 |
| large_minus_small | 8% | 92% | 0% | 6.0 |

## 7 Gate for moving to Phase 2

| Condition | Met |
|---|---|
| solver checks pass | yes |
| source supported bn case defined | **no** |
| heldout median error below 5pct | yes |
| finalist top1 survives direct check | **no** |
| loop top1 within 2se of verified best | yes |
| bo better than random at matched budget | yes |
| bo better than sweep | **no** |
| numerical error below claimed benefit | yes |

Proceed to Phase 2: **False**. Blocking (required conditions not met): source_supported_bn_case_defined, finalist_top1_survives_direct_check. Informational: bo_better_than_sweep. Learned search earns its complexity in this 4-descriptor problem: **False** (the guide's rule applies: an equally effective direct search is retained and reported).

## 8 Limits and failed assumptions

* 2D cross-sections cannot represent 3D connectivity; area fraction is not volume fraction. Phase 2 (3D oblate platelets, BN tensor) is required before any BN claim.
* The material case is a development stand-in; every property is *assumed* (source register). The BN case (SR-001, SR-009) is still pending extraction.
* Finite domain (200 µm wide, insulated sidewalls): results are finite-volume examples, not representative bulk properties; lateral-size sensitivity is reported in Q3.
* Scenario results are worst-case statements over hand-chosen stress cases, not calibrated reliability estimates.
* Not modelled: processing (viscosity, alignment method), particle clustering/contact, voids, adhesion, cure, dielectric strength, thermo-mechanical stress, ageing.

## 9 Reproduce

```bash
pip install -e .[test]   # see README for the Gmsh system library note
tcid all                 # verify, pilot, design loop, final checks, questions, gate, figures, report
```

Individual stages: `tcid verify | pilot | initial | search | final | questions | gate | figures | report`. Raw rows: `results/phase1/**.parquet`; manifest with checksums and environment: `results/phase1/manifest.json`.
