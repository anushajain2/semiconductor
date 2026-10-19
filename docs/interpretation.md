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
