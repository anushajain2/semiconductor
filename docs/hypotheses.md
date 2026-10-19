# Pre-registered hypotheses (Phase 1)

Written on 2026-10-19, **after** the solver verification gates had passed and **before** any
design-set, optimisation or research-question simulation was run (guide section 8). They
exist so that a disappointing result cannot quietly turn into a different success claim. The
outcome of each one is reported in `docs/report.md` against the text below, which has not
been edited since.

Scope: 2D development case (alumina-like isotropic filler in epoxy, area fraction 0.30,
100 µm layer, isothermal faces, insulated sidewalls). These hypotheses say nothing about BN
or 3D connectivity; Phase 2 tests whether conclusions survive there.

| ID | Question | Hypothesis (falsifiable statement) | What would refute it |
|---|---|---|---|
| H1 | Q1: when does the preferred design change? | Within the supported interface range (1e-8 to 1e-6 m²K/W) the best robust design is strongly aligned (S ≥ 0.7) with aspect ratio ≥ 3, and its ranking against the random-orientation baseline does not reverse. The ordering of two aligned designs that differ mainly in particle size can reverse only at or above about 1e-6 m²K/W, because the interface penalty scales with the interface area per filler area (∝ 1/d). | A reversal against the baseline inside the supported range, or a best robust design with S < 0.5 |
| H2 | Q2: when is material optimisation worth doing? | For a 100 µm layer the optimised design keeps at least 10% lower total resistance than the baseline for total outer contact resistance up to 3e-5 m²K/W. The contact level at which the benefit falls to 10% is lower for a 50 µm layer than for a 150 µm layer. | Benefit below 10% inside that contact range at 100 µm, or the thickness trend reversed |
| H3 | Q3: does the optimum rely on ideal geometry? | The advantage over the baseline survives fresh realisations (paired difference larger than 3× its standard error) and decreases smoothly, not abruptly, when the achieved orientation order drops by 0.2. Halving or doubling the minimum gap changes the advantage by less than a third. | Advantage indistinguishable from zero on fresh seeds, a collapse with small perturbations, or strong gap dependence |
| H4 | Q4: which unknown should be resolved next? | Narrowing the achieved orientation order reduces the spread of the candidate-minus-baseline difference more than narrowing the interface resistance within the supported range, which in turn helps more than narrowing filler conductivity. | A different ordering of variance reduction |
| H5 | Does the learned search earn its complexity? | GP-based Bayesian optimisation reaches a robust objective within 1% of its final best with fewer solver calls than matched-budget random search. A full-factorial sweep at three levels may be equally effective; if it is, that will be reported and retained. | Random search matching or beating BO at equal solver calls |
