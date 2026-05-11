#!/usr/bin/env bash
# Post-chain analysis pipeline for the 18 extension specs.
#
# Pre-condition: all 18 spec dirs gibbs_<spec>/ contain chain_1..4.rda.
# Verify via:
#   ls runs/run_extension_2026-05-10/gibbs_*/chain_4.rda | wc -l   # should = 18
#
# Steps fire sequentially because step k+1 reads csvs written by step k.

set -uo pipefail
RD="C:/FkCancer/runs/run_extension_2026-05-10"
RSCRIPT="C:/Program Files/R/R-4.6.0/bin/Rscript.exe"

cd "$RD"

# Pre-flight: verify all 18 spec dirs are complete
need=18
got=$(ls -1 gibbs_*/chain_4.rda 2>/dev/null | wc -l)
if [ "$got" -lt "$need" ]; then
  echo "ABORT: expected $need spec dirs with chain_4.rda, got $got"
  echo "Run 'ls -1 gibbs_*/chain_4.rda' to check"
  exit 1
fi
echo "Pre-flight OK: $got / $need spec dirs have all chains."

START_TIME=$(date)
echo "=========================================================="
echo "Extension analysis pipeline"
echo "Start: $START_TIME"
echo "=========================================================="

echo ""; echo "[1/5] Convergence diagnostics ..."
"$RSCRIPT" 03_convergence_18.R 2>&1 | tee convergence_extension.log

echo ""; echo "[2/5] Held-out LPPD ..."
"$RSCRIPT" 04_compute_loglik_18.R 2>&1 | tee loglik_extension.log

echo ""; echo "[3/5] Bootstrap CI + Bonferroni ..."
"$RSCRIPT" 05_bootstrap_bonferroni_18.R 2>&1 | tee bootstrap_extension.log

echo ""; echo "[4/5] Per-cancer LPPD breakdown ..."
"$RSCRIPT" 06_per_cancer_18.R 2>&1 | tee per_cancer_extension.log

echo ""; echo "[5/5] Synthesize results + apply decision rule ..."
"$RSCRIPT" 07_synthesize_extension_results.R 2>&1 | tee synthesize_extension.log

END_TIME=$(date)
echo ""
echo "=========================================================="
echo "Pipeline complete."
echo "Start: $START_TIME"
echo "End:   $END_TIME"
echo "=========================================================="
echo ""
echo "Outputs:"
echo "  convergence_summary_extension.csv"
echo "  loglik_summary_extension.csv"
echo "  bootstrap_extension_summary.csv"
echo "  per_cancer_extension_long.csv / per_cancer_extension_wide.csv"
echo "  operator_composition_summary.csv"
echo "  extension_decision_log.txt"
echo "  EXTENSION_VALIDATION_RESULTS.md"
