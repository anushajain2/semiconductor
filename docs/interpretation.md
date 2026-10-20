## Summary of findings (Phase 1, 2D development case)

These findings describe a 2D layer, 100 µm thick, with filler area fraction 0.30 and an *assumed* alumina-like filler. They are not boron nitride results.

**How the evidence was produced.** Every number below comes from the result files and is repeated in the tables later in this report.

* Comparisons between designs use Welch tests on independent random particle arrangements ("seeds"), with Holm correction where a family of comparisons is made.
* Locations and thresholds are given with 95% bootstrap intervals.
* "Not resolved" means no difference was detected. It does not mean the designs are equivalent.
* The candidates were chosen on a 6-seed verification batch. Their quoted values come from a separate 12-seed confirmation batch, because the minimum of several noisy values is biased low.

**Revision note.** This revision follows an external review. Several earlier statements were weaker than they read: the size crossovers, the contact-resistance threshold, the claimed "statistical equivalence" of BO and the sweep, the uncertainty table, and the measurement-priority wording. Section 8 lists what changed.

1. **Aligning elongated particles is a large, well-resolved benefit. It is the one firm design conclusion.**
   * Every optimised candidate is elongated (aspect ratio 3.3–4) and at least partly aligned (orientation order 0.56–0.95). The eight BO and sweep candidates all have aspect ratio 4 and order 0.84–0.95.
   * Each one has 17–27% lower worst-case layer resistance than its own random-orientation counterpart with the same loading, sizes and shape (verification batch).
   * On the independent confirmation batch, the recommended design beats its matched baseline by 45 mm²K/W (95% CI 36–54), a 20% reduction.
   * The same comparison on the batch it was selected from gave 26%. The gap is a measured example of best-of-N selection bias.
   * Random circles at the same loading reach 280 mm²K/W. The confirmed BO and sweep candidates lie at 178–185.
2. **The evidence cannot choose among the aligned aspect-ratio-4 designs.**
   * The mechanical recommendation (lowest verification value) is median size 24 µm, width 0.05, aspect ratio 4, orientation order 0.95. It confirms at 178.5 ± 3.0 mm²K/W.
   * Seven of the eleven other optimised candidates are not resolved from it. These include the previous recommendation (6.7 µm), which confirms at 182.9 ± 0.8.
   * On the confirmation batch, every aligned candidate's 95% CI for its difference from the recommended design includes zero. Only the random-search pick is resolved worse (+12.9 mm²K/W, 95% CI 6.2–19.6).
   * The recommended design sits at a corner of the design space (largest size, narrowest width, highest aspect ratio and alignment). It also has the largest realization-to-realization spread. It cannot be used in a 50 µm layer (its particles do not fit), and its advantage falls to 0.44× when particle size is scaled up by 20%.
   * Among the unresolved designs, a ~10 µm design such as `sweep_best` (179.3 ± 0.9) is less variable and fits thinner layers. These are legitimate grounds for choosing between designs the simulation cannot separate.
3. **Particle size versus interface resistance: the direction is robust, the location is not.** (Q1; 12 seeds; dense grid)
   * The size ladder is 5, 11 and 20 µm, all with the recommended width, shape and alignment. Smaller sizes (4 and 4.5 µm) could not be packed at this loading and minimum gap.
   * At low interface resistance, smaller particles are better. As interface resistance rises, their penalty grows steadily. This trend is consistent across seeds and well resolved.
   * The estimated crossings, with 95% intervals:
     * small vs medium: 7.4e-7 m²K/W [4.9e-7, 1.0e-6]
     * small vs large: 7.9e-7 [4.9e-7, 1.1e-6]
     * medium vs large: 8.9e-7 [1.7e-7, 2.2e-6]; 9% of bootstrap resamples show no crossing at all.
   * All three intervals straddle the 1e-6 upper limit of the supported interface range. So the crossings sit near the top of that range but are **not established design boundaries**.
   * Random orientation never wins a cell of the design map. On the 3-seed map, 93% of supported-range cells do not resolve the winner from the runner-up.
4. **Outer contacts do not make the microstructure irrelevant in the assumed range, though the threshold is uncertain.** (Q2; same 12 seeds and mesh as the confirmation batch)
   * At 100 µm, the recommended design's assembly benefit is 20.3% [16.8, 23.8] with ideal contacts.
   * The benefit falls to 10% at a total contact resistance of 220 mm²K/W [144, 299] (nominal scenario), or 214 [138, 292] in the least favourable stress scenario. The threshold multiplies noise roughly tenfold, hence the wide interval.
   * Even its lower bound is more than twice the assumed realistic range (≤ 60 mm²K/W in total).
   * The earlier figure of 331 came from a separate 4-seed batch. It is not reproducible from the main table, and it has been replaced.
   * The threshold rises with thickness. For the medium pair it goes from 292 [238, 348] at 100 µm to 478 [438, 519] at 150 µm.
   * The pre-registered 50 µm comparison could not be made. The recommended distribution does not fit, and the medium pair failed to pack in 3 of 12 (aligned) and 7 of 12 (random) realizations.
   * The series formula R_top + R_layer + R_bottom underestimates Robin contact solves by up to 3.0% for the large-particle design (1.2% for the medium pair).
5. **Robustness of the recommended design to geometric perturbations.** (Q3; 12 seeds)
   * The reference advantage is 17.2% (t = 7.6).
   * It is insensitive to orientation: lowering the order by 0.1 or 0.2 gives 0.99× and 1.05× of the reference. The pre-registered "smooth decrease" was not observed, but neither was a collapse.
   * It is insensitive to minimum gap (0.82–1.26×), lateral domain width (0.95×) and a thinner margin (1.10×).
   * It is sensitive to larger particles (size ×1.2 gives 0.44×) and to a wider size distribution (0.69×).
   * A 6 µm resin margin could not be packed for this design.
6. **Which input each decision is sensitive to, conditionally.** (Q4)
   * Within the reference ranges, the spread of "baseline − candidate" (the value of aligning) is almost entirely explained by the achieved orientation: 98–100%. The interface and filler-conductivity shares are both ≈ 0%, so their relative order is not meaningful.
   * Widening the interface range to 1e-5 m²K/W moves the interface share to 18% (recommended design) and 39% (medium design).
   * For the size decision, the large − small difference varies only with interface resistance. But orientation was *not varied* there, because the small design could not be packed at reduced orientation.
   * These are sensitivity shares that depend on the chosen ranges, grid weighting and independence assumptions. They indicate which uncertainty matters to which decision. They do not say which measurement is worth buying, because measurement cost, accuracy and decision impact are not modelled.
7. **Search methods: one campaign each.**
   * Each method's best-of-four design was verified, then confirmed on 12 independent seeds: BO 182.9 ± 0.8, random search 191.4 ± 1.0, sweep 178.5 ± 3.0 mm²K/W.
   * BO is resolved better than random search in this campaign: random − BO = 8.4 mm²K/W [5.8, 11.1], Holm p < 0.001.
   * The sweep comparison is unresolved: sweep − BO = −4.5 [−11.1, +2.2]. The sweep also reached within 1% of the best loop value in fewer solver calls (585 vs 1020).
   * Because each method ran once, none of this generalises to the methods themselves, and "unresolved" is not equivalence. Following the guide, the cheaper sweep is retained for Phase 1.
   * The GP surrogate met the 5% held-out target (4.1% median), but its 95% intervals covered only 85% of points. A small extra-trees ensemble was more accurate (1.8%).

### Phase 2 gate

The gate is **not passed**, for two reasons:

* the BN material case and BN tensor (SR-001, SR-009) are still pending extraction;
* the loop's top-ranked finalist did not keep first place on fresh seeds, although it is not resolved from the verified best.

### Notes on the hypothesis evaluation

The hypotheses are graded mechanically against intervals.

* **H1** has one unresolved clause, because every crossing interval straddles 1e-6.
* **H2**'s thickness clause is untestable at 50 µm. The trend between 100 and 150 µm is resolved.
* **H3** is only partly supported because the predicted decrease with orientation did not appear, even though no collapse occurred.
* **H4** holds for the reference ranges. The interface-versus-filler part is trivially ordered, since both shares are ≈ 0.
* **H5** holds for this single campaign only.

The scenarios are hand-chosen stress cases, so every "worst case" above is a stress statement, not a reliability estimate.

### 8 What changed after the external review

| Review point | Change |
|---|---|
| Crossovers stated more strongly than the evidence | A dense 17-point grid with 12 seeds per design and bootstrap intervals on each crossing location. Crossings are reported as estimates with intervals; points above the supported range are labelled as stress cases. |
| Contact threshold not reproducible from the main table | Q2 recomputed from the confirmation seeds and mesh, with layer resistances shown, bootstrap intervals on the benefit and the threshold, and a second design pair for the thickness trend. |
| "Statistically equivalent" claims | Removed. Methods are compared on an independent confirmation batch, and each verified the same number (4) of candidates. Statements are limited to this single campaign. |
| Self-comparison row and the "2 SE" rule | The self-comparison is removed. Welch tests with Holm correction replace the 2 SE screen; confirmation-batch values are quoted to remove selection bias. |
| Measurement-priority wording | Reframed as conditional sensitivity shares, with range variants that show how the shares move; no measurement recommendation is made. |
| Verification ahead of validation | Unchanged and agreed: all properties remain assumed, and Phase 2 (BN case, 3D platelets) is required before any material claim. |
