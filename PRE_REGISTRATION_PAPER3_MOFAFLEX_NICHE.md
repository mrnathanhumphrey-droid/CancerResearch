# PRE_REGISTRATION_NICHE_JOINT_REFIT.md

**Status:** LOCKED 2026-05-11 on "send it" trigger. Final "exhaust the trail" pre-reg for the MOFA-FLEX niche operator family. Tests the joint VI refit with class-conditional prior MEAN shift, using a paired-comparison framework that cancels VI re-init drift between specs.

**Scope:** Methodology probe — joint VI refit of MOFA-FLEX with the existing `structural_prior_patch.py` machinery. The patch injects a Horseshoe-prior with shifted mean: `W[k, j] ~ Normal(λ · group_mean_c(k)[j], local_scale_horseshoe)`. Each refit is initialized from the reconverged baseline `W_train` via `weights_init_dict`. The PAIRED test framework compares each spec's refit to a same-seed, same-init λ=0 SANITY refit, cancelling the VI re-init drift that shelved this code path in 2026-05-10.

**Decision favors NULL by default. Expected outcome: FALSIFIED.** Per the pre-existing Option B shelving (λ=0 sanity drifted +31 nats/cell vs reconverged baseline), and per options 1+3 findings (λ=0 optimal in frozen-Z framework), the most likely result is that PRIMARY's prior-mean-shifted refit produces W that fits training data better (because it's pulled toward an average) but generalizes worse to held-out cells. This pre-reg is explicit exhaustion: we run for completeness so future investigators don't ask "did they try joint refit?"

---

## §0. Retro context

- **Option B (shelved 2026-05-10):** λ=0 sanity drifted +31 nats/cell vs on-disk reconverged baseline → operator interpretation ambiguous in absolute terms.
- **Option 1 (niche-conditional eval, 2026-05-11):** all 5 (seed × spec) cells found λ*=(0,...) — no positive class-uniform shrinkage helps. Identity is optimal in the frozen-Z eval grid.
- **Option 3 (niche-precision closed-form, 2026-05-11):** closed-form Gaussian posterior W produces ~17 nats drift vs VI baseline regardless of prior class structure — the closed-form/VI framework gap dominates any class signal.
- **This probe (joint refit, exhaustion):** does joint VI refit with class-conditional prior MEAN shift, paired against same-seed same-init λ=0 sanity refit, produce positive Δ_paired_LL? Predicted no, but corpus needs the explicit test.

---

## §1. Primary hypothesis

- **H_NULL:** Joint VI refit with class-conditional prior mean shift (λ=1) does NOT produce a higher paired-LL on held-out cells than the same-seed same-init λ=0 sanity refit.
- **H_ALT:** PRIMARY paired Δ_LL > 0 with 95% bootstrap CI strictly above 0, AND PRIMARY > NULL by z ≥ 2.

**Disposition ladder:**

| Disposition | Criteria |
|---|---|
| `VALIDATED_JOINT_REFIT` | PRIMARY paired Δ > +0.5 nats/cell with CI > 0 AND beats NULL by z ≥ 2 AND cross-seed reproducibility |
| `MATCHES_SANITY` | \|PRIMARY paired Δ\| ≤ 0.5 nats/cell AND CI includes 0 |
| `FALSIFIED_JOINT_REFIT` | PRIMARY paired Δ < −0.5 nats/cell with CI < 0 |
| `INCONCLUSIVE` | Cross-seed variance dominates signal OR sanity refit fails to converge |

**Predicted disposition:** `FALSIFIED_JOINT_REFIT`. The pre-reg makes this explicit; we run for completeness.

---

## §2. Data + baseline state

- **Reconverged baseline:** ON DISK at `runs/seed_<s>_reconverged/{model.h5, weights_rna.npy, factors_*.npy}` for s ∈ {2511021635, 2511021636, 2511021637, 2511021638}.
- **Held-out split:** Phase 5 split (seed 20260510, per-niche/celltype stratified 20% holdout — used by Option B's existing `reconstruction_ll_per_cell_subsetted`). NOTE: this is a DIFFERENT split than Option C1's 80/20 (seed 42); not directly comparable to options 1/3 numbers. The paired-Δ vs SANITY refit makes the comparison meaningful within this framework.
- **Niche labels per factor:** computed from baseline_Z_xen + xenium niche labels (Option B's `assign_niche_labels`).
- **No new fits to non-MOFA-FLEX components.**

---

## §3. Operator specification (LOCKED — reuses Option B machinery)

Existing `structural_prior_patch.py` injects:

```
W[k, j] ~ Normal(λ · prior_mean[k, j], local_scale_horseshoe[k, j])
```

For dense (uninformed) factors: `prior_mean[k, :] = 0` → no shift.
For informed factor k with class c(k): `prior_mean[k, :] = mean_{k' ∈ c(k)} W_train[k', :]` (per-class group mean from baseline W).

VI is initialized at `weights_init_dict = {"rna": baseline_W}` (init from reconverged baseline). All other hyperparameters from `run_option_b_refit.py` — `MAX_EPOCHS=1000`, `PATIENCE=50`, same Horseshoe `WEIGHT_PRIOR`, same `LIKELIHOOD`, same `BATCH_SIZE`, same `LR`.

The λ scalar controls operator strength. λ=0 is "no prior shift" (sanity). λ=1 is "full mean shift" (matches Option C1's full-replacement spec semantics).

---

## §4. Specs (LOCKED — 4 specs × 4 seeds = 16 refits)

| Spec ID | λ | hierarchy | shuffle | Description |
|---|---|---|---|---|
| **SANITY** | 0 | niche_3cat | 0 | No prior shift; same VI init + seed as PRIMARY. Reference for paired Δ. |
| **PRIMARY** | 1 | niche_3cat | 0 | Joint refit toward niche_3cat group means. |
| **SECONDARY** | 1 | niche_5cat | 0 | Joint refit toward niche_5cat group means. |
| **NULL** | 1 | niche_3cat | 42 | Joint refit toward SHUFFLED-label group means. |

Seeds: same 4 as Option C1 + Option B: {2511021635, 2511021636, 2511021637, 2511021638}. VI seed = MOFA-FLEX `training_opts.seed = args.seed` (per Option B's convention). So for given seed s, all 4 specs use the same VI init randomization, but their prior_mean differs.

---

## §5. Paired-test framework (LOCKED)

For each seed s, compute per-cell held-out LL for each of the 4 spec refits. Then per-cell paired Δ:

```
Δ_PRIMARY_paired[cell] = LL_PRIMARY_refit[cell] − LL_SANITY_refit[cell]
Δ_SECONDARY_paired[cell] = LL_SECONDARY_refit[cell] − LL_SANITY_refit[cell]
Δ_NULL_paired[cell] = LL_NULL_refit[cell] − LL_SANITY_refit[cell]
```

The SANITY refit's VI re-init drift cancels in the difference (both runs face same VI noise from same init + same seed).

Aggregate across seeds: mean paired Δ per spec + 95% CI via cell-level bootstrap (B = 1000).

---

## §6. Decision rule (LOCKED)

PRIMARY passes only if:

1. **PRIMARY paired Δ > +0.5 nats/cell** with 95% bootstrap CI lower bound > 0.
2. **PRIMARY paired Δ > NULL paired Δ** with z ≥ 2 (per-seed standard error).
3. **Cross-seed reproducibility:** standard deviation of PRIMARY paired Δ across 4 seeds ≤ 5 nats/cell (drift-allowing).

If (1) and (2) pass → `VALIDATED_JOINT_REFIT`.
If |PRIMARY Δ| ≤ 0.5 → `MATCHES_SANITY`.
If PRIMARY Δ < −0.5 with CI < 0 → `FALSIFIED_JOINT_REFIT`.

---

## §7. Adversarial cascade (LOCKED, all 3 required for VALIDATED)

A. **SANITY-vs-on-disk drift bound:** |LL_SANITY_refit − LL_reconverged_baseline| should be < 50 nats/cell (loose bound). If SANITY refit produces >50 nats deviation, the VI refit framework is unusable and disposition = `INVALID_REFIT_DIVERGENT`.

B. **Cross-seed reproducibility for SANITY:** σ(LL_SANITY) across 4 seeds ≤ 10 nats/cell. Documents VI re-init noise floor.

C. **NULL adversarial:** NULL paired Δ should not beat PRIMARY paired Δ (z ≥ 2 separation).

---

## §8. Pre-committed deliverables

| File | Contents |
|---|---|
| `PRE_REGISTRATION_NICHE_JOINT_REFIT.md` | this document |
| `run_niche_joint_refit_orchestrate.sh` | sequential firing of 16 refits |
| `run_niche_joint_refit_finalize.py` | paired-Δ aggregator + cascade + verdict |
| `runs/niche_joint_refit/seed<s>_<spec>/` | per-refit artifacts (W, Z, LL, summary) |
| `runs/niche_joint_refit/disposition_table.csv` | per-spec paired Δ + CI |
| `runs/niche_joint_refit/check_diagnostics.csv` | cascade pass/fail |
| `runs/niche_joint_refit/NICHE_JOINT_REFIT_RESULTS.md` | top-level disposition + narrative |

---

## §9. Compute

- 16 refits × ~45 min each = ~12 hr GPU wall (cuda:0, sequential).
- No CPU/disk contention with hydrology Stan runs.

---

## §10. Lock-in protocol

This pre-registration is **immutable** once committed to disk. Pre-reg violations during compute → `PRE_REG_DEVIATION_LOG.md`. **Compute fires on locked file; no further confirmation required.**

---

## §11. Explicit exhaustion framing

This is the fourth and final operator probe of the MOFA-FLEX niche corpus:
1. **Option C1 uniform-1 full replacement (FALSIFIED):** Δ ≈ −32 nats/cell in frozen-Z eval.
2. **Option 1 class-conditional λ grid (in-flight result: BEATS_UNIFORM_FAILS_BASELINE):** λ*=0 identity-optimal across all specs and NULL.
3. **Option 3 closed-form Gaussian posterior (INVALID_IMPLEMENTATION):** framework gap with VI baseline produces ~17 nats noise that drowns out class signal.
4. **This pre-reg — joint VI refit (predicted FALSIFIED):** completeness exhaustion.

Once this lands, the corpus statement is: "MOFA-FLEX's per-factor LL-optimal W substrate does not admit class-conditional structural priors via linear-mean-shift operators under either frozen-W eval, frozen-Z closed-form, or joint VI refit. The η²=0.31 per-class signal in damage-under-replacement is real but does not translate to a positive operator in this operator family."

Future investigators wanting to test class-conditional priors on factor-analysis substrates should consider: (a) different operator families (multiplicative scale, tensor-product, etc.), (b) different comparators (uninformed W baseline à la Lock 2022), or (c) different likelihood frameworks (NB or Beta-Binomial rather than Gaussian).
