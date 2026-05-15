# Phase 1 Results — Paper 4: Structural-Prior Validation Across 3 New External Hierarchies

**Date:** 2026-05-15
**Pre-registration lock:** [PRE_REGISTRATION_PAPER4_PHASE1.md](PRE_REGISTRATION_PAPER4_PHASE1.md) at commit `677406a`, public 2026-05-15.

## TL;DR

**Per-locked-pre-reg verdict: `PAPER1_HIERARCHY_SPECIFIC`.** 0/3 of the new external hierarchies meet the pre-reg threshold for "outperforms" (Δ > +2 nats AND bootstrap 95% CI excludes zero). 0/3 underperform either.

**But the directional pattern is informative.** Two of three hierarchies (Thorsson immune subtype, Sanchez-Vega oncogenic pathway) show observed Δ in the same direction and roughly comparable magnitude to Paper 1's primary specification — just with bootstrap CIs that span zero given the n=1,362 test set. The third (Malta stemness tertile) returned ≈0 effect.

Read together: Paper 1's +6.26 nat result is **not robust** at the locked threshold across these three independently-published external hierarchies, but the directional signal in two of three is consistent with the methodology generalizing — within a power-limited statistical regime.

## Convergence — ALL PASS

| Spec | max R-hat | min ESS_bulk | verdict |
|---|---|---|---|
| baseline_train | 1.0078 | 598 | PASS |
| thorsson | 1.0160 | 394 | PASS |
| malta | 1.0141 | 485 | PASS |
| sanchezvega | 1.0105 | 425 | PASS |

All four specifications cleared the pre-reg halt threshold (R-hat ≤ 1.05). No re-fits triggered.

## Held-out LPPD

Per the pre-reg, the test set is the same 1,362-patient 20% stratified split as Paper 1 (seed 42). All four specs evaluated on identical test patients.

| Spec | Total LPPD | Δ vs baseline | Patients |
|---|---|---|---|
| baseline_train (local re-fit) | **−1,038.99** | — | 1,362 |
| thorsson (immune subtype) | −1,034.30 | **+4.69** | 1,362 |
| malta (stemness tertile) | −1,038.29 | +0.71 | 1,362 |
| sanchezvega (oncogenic pathway) | −1,035.85 | **+3.15** | 1,362 |

For comparison, Paper 1 primary (epithelial-class + germ-layer): Δ = **+6.26 nats** vs same baseline.

## Bootstrap CIs + per-spec verdicts

B = 1,000 patient-level resamples (seed 20260515).

| Comparison | Observed Δ | 95% CI | Verdict per pre-reg |
|---|---|---|---|
| **thorsson vs baseline** | **+4.69** | **[−2.37, +11.45]** | **INCONCLUSIVE** |
| malta vs baseline | +0.71 | [−5.57, +7.38] | MATCHES |
| **sanchezvega vs baseline** | **+3.15** | **[−4.10, +9.35]** | **INCONCLUSIVE** |
| thorsson vs malta | +3.99 | [−0.53, +8.22] | INCONCLUSIVE |
| thorsson vs sanchezvega | +1.54 | [−4.16, +7.18] | MATCHES |
| malta vs sanchezvega | −2.44 | [−7.30, +2.47] | INCONCLUSIVE |

**Per-spec rule reminder:**
- OUTPERFORMS = observed Δ > +2 AND 95% CI excludes zero (positive side).
- MATCHES = |observed Δ| ≤ 2 AND 95% CI includes zero.
- UNDERPERFORMS = observed Δ < −2 OR 95% CI excludes zero (negative side).
- INCONCLUSIVE = mixed signal (e.g., Δ > +2 but CI spans zero).

## Joint disposition

Hypothesis count: **0 OUTPERFORMS, 1 MATCHES, 2 INCONCLUSIVE, 0 UNDERPERFORMS.**

Per the pre-reg joint table, this maps to **`PAPER1_HIERARCHY_SPECIFIC`** — the methodology's +6.26 nat win on epithelial-class + germ-layer hierarchies did not replicate at the locked threshold on Thorsson immune subtype, Malta stemness tertile, or Sanchez-Vega oncogenic pathway.

## Honest commentary beyond the locked verdict

The pre-reg verdict is `PAPER1_HIERARCHY_SPECIFIC`. That stands. Below is what the data shows beyond the locked threshold, reported for transparency, not to revise the verdict.

1. **Two hierarchies are directionally consistent with Paper 1.** Thorsson +4.69 is 75% of Paper 1's +6.26 in magnitude; Sanchez-Vega +3.15 is 50% of it. Both observed deltas are well above the +2 nat practical threshold; only the bootstrap CI is wider than the gap to zero.

2. **The bootstrap CI widths (~±7 nats) suggest statistical underpowering, not directional null.** The 1,362-patient test set produces CIs of comparable width across all three new specs. With this much CI width, even a Paper 1-magnitude effect (+6.26) wouldn't reliably exclude zero. Paper 1's primary spec barely cleared the CI bar at [+0.34, +12.29] — and that was after multiple chains converging on aligned posteriors. The Phase 1 specs each got the same chain budget but on a more diffuse hierarchy assignment (5-6 group structures vs Paper 1's 10-bin and 3-bin epithelial-class).

3. **The Malta MATCHES verdict is genuinely null, not power-limited.** Δ = +0.71 with CI [−5.57, +7.38] centers near zero. Stemness tertile, as constructed here (cancer-level mean mDNAsi tertile-split), did not extract structural-prior signal that improved held-out prediction. This is a real negative on the stemness-as-structural-prior hypothesis.

4. **Thorsson and Sanchez-Vega's failure to cross the CI threshold is a real Phase 1 finding** under the pre-reg rule — but the appropriate Phase 2 question is whether n_test=1,362 has the power to cleanly distinguish "matches" from "outperforms" when the effect magnitude is in the +3 to +5 nat range. A larger test set (e.g., a refit on a 50/50 split, or a k-fold extension) would resolve this.

## Methodology corpus integration

For the Paper 6 methodology corpus tracking sheet, this updates the Cancer (Lock-extension) entry from "queued" to:

- **Paper 1 (epi + germ_layer):** POSITIVE (+6.26 [+0.34, +12.29])
- **Paper 4 Phase 1 (immune, stemness, pathway):** PAPER1_HIERARCHY_SPECIFIC by pre-reg; two of three (immune, pathway) trend in same direction with comparable magnitude but statistically underpowered at n_test = 1,362.

Combined cancer-substrate reading: the methodology works on epithelial-class / germ-layer (Paper 1) and arguably trends-positive on immune-subtype and pathway-dominant hierarchies but the locked test couldn't confirm. NOT a flat "cancer is fully validated" — and that's the honest framing.

## Compute

- 16 chains × 100k iters, 4 specs × 4 chains each.
- Wall: 44 min for all 16 chains to finish in parallel on 9950X3D (16C/32T). Eval pipeline: ~3 min (convergence + LPPD + bootstrap).
- Chain RDA files: 9.4 GB total at `C:/FkCancer/runs/run_paper4_phase1_2026-05-15/gibbs_*/chain_*.rda` (excluded from GitHub via `.gitignore`).

## Files

- [PRE_REGISTRATION_PAPER4_PHASE1.md](PRE_REGISTRATION_PAPER4_PHASE1.md) — locked pre-reg (commit `677406a`)
- [reference/phase1_hierarchy_assignments.csv](reference/phase1_hierarchy_assignments.csv) — per-cancer hierarchy assignments
- [scripts/validation/build_phase1_hierarchies.py](scripts/validation/build_phase1_hierarchies.py) — hierarchy build from raw paper supplementary data
- [PHASE1_RESULTS.md](PHASE1_RESULTS.md) — this file
- `C:/FkCancer/runs/run_paper4_phase1_2026-05-15/` — Phase 1 run dir (chain RDAs, eval scripts, logs, loglik_summary.csv, bootstrap_ci_summary.csv, convergence_summary.csv, JOINT_DISPOSITION.md)

## Phase 2 — what would change the verdict

Per the data above, the bottleneck is statistical power on the 1,362-patient test set. Options for Phase 2:

1. **k-fold extension on the same 4 covariates × 3 hierarchies.** 5-fold cross-validation would give ~6,800-patient effective test exposure across folds, tightening CIs by ~√5 ≈ 2.2×. Thorsson's CI [−2.37, +11.45] would tighten to roughly [+0.4, +9.0] if the same effect held — moving from INCONCLUSIVE to OUTPERFORMS.
2. **Different decision-rule pre-reg.** Move from "CI excludes zero" to "posterior probability Δ > 0 exceeds 95%" — a Bayesian framing more sensitive at this n.
3. **Stricter contested-cell handling.** Currently Hierarchy A (Thorsson) excludes DLBC + THYM. Re-running with looser inclusion (e.g., assigning them to nearest neighbor) might tighten the structural signal.
4. **Different baseline.** Use the Paper 1 baseline-on-train chains (re-built locally as part of Phase 1) consistently across all comparisons.

These are Phase 2 considerations; Phase 1 has shipped per the pre-reg.
