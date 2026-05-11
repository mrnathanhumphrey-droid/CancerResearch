#!/usr/bin/env bash
# Master runner: fires all 18 extension Gibbs specifications sequentially,
# with 4 chains in parallel within each spec.
#
# Total wall time estimate: 18 specs * ~37 min = ~11 hours.

set -uo pipefail
RD="C:/FkCancer/runs/run_extension_2026-05-10"
RSCRIPT="C:/Program Files/R/R-4.6.0/bin/Rscript.exe"

SPECS=(
  spec_01_cov8_2_rule1_pairA
  spec_02_cov8_2_rule1_pairB
  spec_03_cov8_2_rule2_pairA
  spec_04_cov8_2_rule2_pairB
  spec_05_cov8_2_rule3_pairA
  spec_06_cov8_2_rule3_pairB
  spec_07_cov10_1_rule1_pairA
  spec_08_cov10_1_rule1_pairB
  spec_09_cov10_1_rule2_pairA
  spec_10_cov10_1_rule2_pairB
  spec_11_cov10_1_rule3_pairA
  spec_12_cov10_1_rule3_pairB
  spec_13_cov1_1_rule1_pairA
  spec_14_cov1_1_rule1_pairB
  spec_15_cov1_1_rule2_pairA
  spec_16_cov1_1_rule2_pairB
  spec_17_cov1_1_rule3_pairA
  spec_18_cov1_1_rule3_pairB
)

cd "$RD"
START_TIME=$(date)
echo "=========================================================="
echo "Extension Gibbs runner — 18 specs sequential"
echo "Start: $START_TIME"
echo "=========================================================="

for spec in "${SPECS[@]}"; do
  echo ""
  echo "=========================================================="
  echo "[$(date)] Starting $spec"
  echo "=========================================================="
  for chain in 1 2 3 4; do
    "$RSCRIPT" run_extension_chain.R "$spec" "$chain" \
      > "${spec}_chain${chain}.log" 2>&1 &
  done
  wait
  echo "[$(date)] Finished $spec"
  ls -la "gibbs_${spec}/" 2>&1 || echo "MISSING: gibbs_${spec}/ directory not found"
done

END_TIME=$(date)
echo ""
echo "=========================================================="
echo "All 18 specs complete."
echo "Start: $START_TIME"
echo "End:   $END_TIME"
echo "=========================================================="
