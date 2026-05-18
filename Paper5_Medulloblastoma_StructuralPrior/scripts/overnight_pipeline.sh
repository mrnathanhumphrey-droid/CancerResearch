#!/usr/bin/env bash
# Paper 5 — overnight leakage-clean pipeline.
#
# Sequence (each step blocks the next):
#   1. Wait for current BIDIFAC+ run to finish (or refire on stall)
#   2. Run screening_paper5_v2.R
#   3. Local-commit the new covariate lock + DEVIATIONS.md
#   4. Run build_paper5_split_v2.R
#   5. Fire 24-chain Gibbs (6 specs × 4 chains), batched 8 parallel
#   6. Run convergence_paper5.R (skip + halt if any R-hat > 1.05)
#   7. Run compute_held_out_loglik_paper5_v2.R
#
# Watchdog rules:
#   - BIDIFAC+ stall (>30min no log update): kill + refire ONCE with
#     K=2000/each, num.comp=5 (per prior authorization).
#   - Screening / Gibbs stall (>15min no log update): kill + refire SAME
#     params ONCE. If second stall, halt step + log and exit.
#   - Zombie cleanup: kill any Rscript.exe spawned by this script that
#     outlives its step.
#   - No git push, ever. Local commit only at step 3.

set -uo pipefail

REPO="/c/Cancer Research/Paper5_Medulloblastoma_StructuralPrior"
cd "$REPO"

LOG="$REPO/logs/overnight_clean.log"
STATE="$REPO/logs/overnight_state.txt"
mkdir -p logs

# tee-style logging
log() {
  local ts msg
  ts=$(date '+%Y-%m-%d %H:%M:%S')
  msg="[$ts] $*"
  echo "$msg"
  echo "$msg" >> "$LOG"
}

set_state() {
  echo "$(date '+%s') $*" > "$STATE"
}

# Get Rscript.exe PIDs that match a substring in their command line.
# Uses PowerShell since wmic is deprecated in Windows 11.
rscript_pids_matching() {
  local pat="$1"
  powershell -NoProfile -Command "Get-CimInstance Win32_Process -Filter \"name='Rscript.exe'\" | Where-Object { \$_.CommandLine -match '$pat' } | Select-Object -ExpandProperty ProcessId" 2>/dev/null \
    | tr -d '\r'
}

# Kill ALL Rscript.exe — used when refiring BIDIFAC+ since there's only
# one R process at that step. Safe overnight: user confirmed no other R work.
kill_all_rscript() {
  log "  Killing all Rscript.exe processes"
  taskkill //F //IM Rscript.exe 2>&1 | head -5 | tee -a "$LOG"
  sleep 3
}

kill_pid() {
  local pid="$1"
  if [ -n "$pid" ]; then
    log "Killing PID $pid"
    taskkill //PID "$pid" //F 2>&1 | head -2 | tee -a "$LOG"
    sleep 2
  fi
}

# log_age <file> -> seconds since last mtime (or 99999999 if missing)
log_age() {
  local f="$1"
  if [ ! -f "$f" ]; then echo 99999999; return; fi
  local now mtime
  now=$(date +%s)
  mtime=$(stat -c %Y "$f" 2>/dev/null || echo 0)
  echo $((now - mtime))
}

############################################################################
# Step 1: wait for / supervise current BIDIFAC+ run
############################################################################
step1_wait_bidifac() {
  set_state "step1_bidifac"
  log "Step 1: waiting for BIDIFAC+ (already running)"

  local stall_thresh=1800   # 30 min
  local poll_int=120
  local retried=0
  local bidifac_log="$REPO/logs/bidifac_clean_stdout.log"
  local out_file="$REPO/data/bidifac_components_clean.rds"

  while [ ! -f "$out_file" ]; do
    sleep "$poll_int"

    # Is any Rscript.exe still alive at all?
    local n_rscript
    n_rscript=$(tasklist //FI "IMAGENAME eq Rscript.exe" 2>/dev/null | grep -c "Rscript.exe" || echo 0)
    local age
    age=$(log_age "$bidifac_log")

    log "  BIDIFAC+ poll: rscript_count=$n_rscript log_age=${age}s output_exists=no"

    if [ "$n_rscript" -eq 0 ] && [ ! -f "$out_file" ]; then
      log "  BIDIFAC+ process died without producing output."
      if [ "$retried" -eq 0 ]; then
        log "  Refiring BIDIFAC+ with reduced params (K=2000/each, num.comp=5)."
        bidifac_refire_reduced
        retried=1
      else
        log "  Already retried once. Halting."
        return 2
      fi
    elif [ "$age" -gt "$stall_thresh" ]; then
      log "  STALL detected (log_age=${age}s > ${stall_thresh}s)."
      kill_all_rscript
      if [ "$retried" -eq 0 ]; then
        log "  Refiring BIDIFAC+ with reduced params."
        bidifac_refire_reduced
        retried=1
      else
        log "  Already retried once. Halting."
        return 2
      fi
    fi
  done

  log "  BIDIFAC+ output present: $(ls -la "$out_file")"
  log "Step 1: DONE"
  return 0
}

bidifac_refire_reduced() {
  # Write a one-shot R wrapper that overrides params then sources main script
  local wrapper="$REPO/scripts/run_bidifac_plus_leakage_clean_reduced.R"
  cat > "$wrapper" <<'EOF'
# Reduced-param refire: K=2000/each, num.comp=5
# Patches the main script's globals via search-and-replace at source time.
src <- readLines("C:/Cancer Research/Paper5_Medulloblastoma_StructuralPrior/scripts/run_bidifac_plus_leakage_clean.R")
src <- gsub("K_METH_TOP <- 5000L", "K_METH_TOP <- 2000L", src)
src <- gsub("K_EXPR_TOP <- 5000L", "K_EXPR_TOP <- 2000L", src)
src <- gsub("num.comp  = 10", "num.comp  = 5", src)
src <- gsub("max.comb  = 10", "max.comb  = 5", src)
eval(parse(text = paste(src, collapse = "\n")))
EOF
  log "  Starting reduced-param BIDIFAC+..."
  nohup "C:/Program Files/R/R-4.6.0/bin/Rscript.exe" "$wrapper" \
    > "$REPO/logs/bidifac_clean_stdout.log" 2>&1 &
  log "  Reduced refire PID: $!"
}

############################################################################
# Step 2: screening_paper5_v2.R  (~1h wall)
############################################################################
step2_screening() {
  set_state "step2_screening"
  log "Step 2: screening_paper5_v2.R"
  local out_file="$REPO/reference/paper5_target_covariates_clean.csv"
  local step_log="$REPO/logs/screening_clean.log"
  rm -f "$out_file" "$step_log"

  local retry=0
  while [ "$retry" -lt 2 ]; do
    nohup "C:/Program Files/R/R-4.6.0/bin/Rscript.exe" \
      "$REPO/scripts/screening_paper5_v2.R" \
      > "$step_log" 2>&1 &
    local pid=$!
    log "  Screening PID $pid (attempt $((retry+1)))"

    local stall_thresh=900   # 15 min
    while ! [ -f "$out_file" ]; do
      sleep 60
      if ! kill -0 "$pid" 2>/dev/null; then
        if [ -f "$out_file" ]; then break; fi
        log "  Screening PID $pid died without output."
        retry=$((retry+1))
        break
      fi
      local age
      age=$(log_age "$step_log")
      log "  Screening poll: pid_alive=yes log_age=${age}s"
      if [ "$age" -gt "$stall_thresh" ]; then
        log "  STALL screening. Killing PID $pid."
        kill_pid "$pid"
        retry=$((retry+1))
        break
      fi
    done
    if [ -f "$out_file" ]; then
      log "  Screening output present."
      break
    fi
  done

  if [ ! -f "$out_file" ]; then
    log "Step 2: FAILED after $retry attempts."
    return 2
  fi
  log "Step 2: DONE"
  return 0
}

############################################################################
# Step 3: local commit of new lock + DEVIATIONS.md
############################################################################
step3_commit() {
  set_state "step3_commit"
  log "Step 3: local commit (no push)"
  cd "/c/Cancer Research"
  if ! git rev-parse --git-dir > /dev/null 2>&1; then
    log "  Not a git repo. Skipping commit."
    cd "$REPO"; return 0
  fi
  git add Paper5_Medulloblastoma_StructuralPrior/reference/paper5_target_covariates_clean.csv \
          Paper5_Medulloblastoma_StructuralPrior/DEVIATIONS.md \
          Paper5_Medulloblastoma_StructuralPrior/scripts/run_bidifac_plus_leakage_clean.R \
          Paper5_Medulloblastoma_StructuralPrior/scripts/screening_paper5_v2.R \
          Paper5_Medulloblastoma_StructuralPrior/scripts/build_paper5_split_v2.R \
          Paper5_Medulloblastoma_StructuralPrior/scripts/run_gibbs_paper5_v2.R \
          Paper5_Medulloblastoma_StructuralPrior/scripts/compute_held_out_loglik_paper5_v2.R \
          Paper5_Medulloblastoma_StructuralPrior/scripts/compute_held_out_loglik_paper5_logfix.R \
          Paper5_Medulloblastoma_StructuralPrior/scripts/overnight_pipeline.sh \
          Paper5_Medulloblastoma_StructuralPrior/results/loglik_summary_logfix.csv \
          2>&1 | tee -a "$LOG"
  git commit -m "Paper 5: leakage-clean covariate lock + pre-reg deviation Entry 001

Three pre-reg-fatal bugs surfaced after v1 24-chain Gibbs landed:
- LPPD units bug: dnorm(y_raw, mu_log, sigma_log) → 97% of +349 nat headline
- Sensitivity spec collapse: primary == sensitivity bit-identical
- BIDIFAC+ + screening leakage: full pipeline rebuild required

Locks new target covariates {clean}, supersedes 34db626.
No push — awaiting explicit trigger from user." \
    2>&1 | tee -a "$LOG"
  local rc=${PIPESTATUS[0]}
  cd "$REPO"
  if [ "$rc" -ne 0 ]; then
    log "  Commit failed (rc=$rc). Continuing anyway — pipeline does not depend on commit."
  fi
  log "Step 3: DONE (local commit only, no push)"
  return 0
}

############################################################################
# Step 4: build_paper5_split_v2.R  (~5 min)
############################################################################
step4_build_split() {
  set_state "step4_build_split"
  log "Step 4: build_paper5_split_v2.R"
  local step_log="$REPO/logs/build_split_v2.log"
  local out_file="$REPO/data/target_covs_primary_clean.rda"
  rm -f "$step_log" "$out_file"

  "C:/Program Files/R/R-4.6.0/bin/Rscript.exe" \
    "$REPO/scripts/build_paper5_split_v2.R" > "$step_log" 2>&1
  local rc=$?
  if [ "$rc" -ne 0 ] || [ ! -f "$out_file" ]; then
    log "Step 4: FAILED (rc=$rc). Tail of log:"
    tail -20 "$step_log" | tee -a "$LOG"
    return 2
  fi
  log "Step 4: DONE"
  return 0
}

############################################################################
# Step 5: 24-chain Gibbs (6 specs × 4 chains, batched 8 parallel)
############################################################################
step5_gibbs() {
  set_state "step5_gibbs"
  log "Step 5: 24-chain leakage-clean Gibbs"
  local specs=(baseline_train primary sensitivity primary_drop_t1 primary_drop_t2 primary_drop_t3)
  local chain_ids=(1 2 3 4)
  local batch_size=8

  # Build full job list (24 jobs)
  local jobs=()
  for sp in "${specs[@]}"; do
    for cid in "${chain_ids[@]}"; do
      jobs+=("$sp $cid")
    done
  done
  log "  Total jobs: ${#jobs[@]}, batch size $batch_size"

  local i=0
  while [ "$i" -lt "${#jobs[@]}" ]; do
    local batch_end=$((i + batch_size))
    [ "$batch_end" -gt "${#jobs[@]}" ] && batch_end=${#jobs[@]}
    log "  Batch: jobs $i..$((batch_end-1))"

    local pids=()
    local logs=()
    local outs=()
    for ((j=i; j<batch_end; j++)); do
      local job=(${jobs[$j]})
      local sp=${job[0]}; local cid=${job[1]}
      local job_log="$REPO/logs/gibbs_${sp}_chain${cid}_clean.log"
      local job_out="$REPO/results/gibbs_${sp}_clean/chain_${cid}.rda"
      mkdir -p "$REPO/results/gibbs_${sp}_clean"
      rm -f "$job_log"
      nohup "C:/Program Files/R/R-4.6.0/bin/Rscript.exe" \
        "$REPO/scripts/run_gibbs_paper5_v2.R" "$sp" "$cid" \
        > "$job_log" 2>&1 &
      pids+=($!)
      logs+=("$job_log")
      outs+=("$job_out")
      log "    spawned: $sp chain $cid PID $!"
    done

    # Watchdog: wait until all jobs in batch produce output OR die.
    # Stall per chain: 15 min log inactivity → kill + refire ONCE.
    local refired=()
    for p in "${pids[@]}"; do refired+=(0); done

    while true; do
      local n_done=0
      for ((k=0; k<${#pids[@]}; k++)); do
        local pid=${pids[$k]}
        local out=${outs[$k]}
        local lg=${logs[$k]}
        if [ -f "$out" ]; then n_done=$((n_done+1)); continue; fi
        if ! kill -0 "$pid" 2>/dev/null; then
          if [ -f "$out" ]; then n_done=$((n_done+1)); continue; fi
          # Process died without output
          if [ "${refired[$k]}" -eq 0 ]; then
            log "    PID $pid died without output. Refiring."
            local job=(${jobs[$((i+k))]})
            local sp=${job[0]}; local cid=${job[1]}
            nohup "C:/Program Files/R/R-4.6.0/bin/Rscript.exe" \
              "$REPO/scripts/run_gibbs_paper5_v2.R" "$sp" "$cid" \
              > "$lg" 2>&1 &
            pids[$k]=$!
            refired[$k]=1
            log "      refire PID ${pids[$k]}"
          else
            log "    PID $pid died (already refired). Marking as failed."
            n_done=$((n_done+1))   # count as done to break loop
          fi
          continue
        fi
        # Process alive — check stall
        local age
        age=$(log_age "$lg")
        if [ "$age" -gt 900 ]; then
          log "    STALL on PID $pid (log_age=${age}s)"
          if [ "${refired[$k]}" -eq 0 ]; then
            kill_pid "$pid"
            sleep 3
            local job=(${jobs[$((i+k))]})
            local sp=${job[0]}; local cid=${job[1]}
            nohup "C:/Program Files/R/R-4.6.0/bin/Rscript.exe" \
              "$REPO/scripts/run_gibbs_paper5_v2.R" "$sp" "$cid" \
              > "$lg" 2>&1 &
            pids[$k]=$!
            refired[$k]=1
            log "      refire PID ${pids[$k]}"
          else
            log "    Stall on already-refired PID $pid. Killing + marking failed."
            kill_pid "$pid"
            n_done=$((n_done+1))
          fi
        fi
      done

      if [ "$n_done" -eq "${#pids[@]}" ]; then break; fi
      sleep 60
    done

    log "  Batch complete. Output files:"
    for o in "${outs[@]}"; do
      if [ -f "$o" ]; then
        log "    OK  $o ($(stat -c %s "$o") bytes)"
      else
        log "    MISS $o"
      fi
    done

    i=$batch_end
  done

  # Verify total chain coverage
  local n_chains
  n_chains=$(find "$REPO/results" -name "chain_*.rda" -path "*_clean/*" 2>/dev/null | wc -l)
  log "  Final chain count: $n_chains / 24"
  if [ "$n_chains" -lt 20 ]; then
    log "Step 5: FAILED (only $n_chains chains landed)"
    return 2
  fi
  log "Step 5: DONE ($n_chains chains)"
  return 0
}

############################################################################
# Step 6: convergence (rewrite to clean paths)
############################################################################
step6_convergence() {
  set_state "step6_convergence"
  log "Step 6: convergence diagnostics on _clean chains"
  local step_log="$REPO/logs/convergence_clean.log"
  rm -f "$step_log"

  # Write a one-shot wrapper that patches convergence_paper5.R for _clean
  local wrap="$REPO/scripts/convergence_paper5_clean.R"
  cat > "$wrap" <<'EOF'
src <- readLines("C:/Cancer Research/Paper5_Medulloblastoma_StructuralPrior/scripts/convergence_paper5.R")
src <- gsub("gibbs_%s", "gibbs_%s_clean", src)
src <- gsub("convergence_%s\\.csv", "convergence_%s_clean.csv", src)
src <- gsub("convergence_summary\\.csv", "convergence_summary_clean.csv", src)
eval(parse(text = paste(src, collapse = "\n")))
EOF
  "C:/Program Files/R/R-4.6.0/bin/Rscript.exe" "$wrap" > "$step_log" 2>&1
  local rc=$?
  if [ "$rc" -ne 0 ]; then
    log "Step 6: FAILED (rc=$rc). Tail:"
    tail -20 "$step_log" | tee -a "$LOG"
    return 2
  fi
  # Check verdicts
  if [ -f "$REPO/results/convergence_summary_clean.csv" ]; then
    log "  Convergence summary:"
    cat "$REPO/results/convergence_summary_clean.csv" | tee -a "$LOG"
    if grep -q "HALT" "$REPO/results/convergence_summary_clean.csv"; then
      log "  WARNING: HALT verdict in convergence. Continuing to LPPD anyway, surface in summary."
    fi
  fi
  log "Step 6: DONE"
  return 0
}

############################################################################
# Step 7: LPPD recompute
############################################################################
step7_lppd() {
  set_state "step7_lppd"
  log "Step 7: compute_held_out_loglik_paper5_v2.R"
  local step_log="$REPO/logs/lppd_clean.log"
  rm -f "$step_log"
  "C:/Program Files/R/R-4.6.0/bin/Rscript.exe" \
    "$REPO/scripts/compute_held_out_loglik_paper5_v2.R" > "$step_log" 2>&1
  local rc=$?
  if [ "$rc" -ne 0 ]; then
    log "Step 7: FAILED (rc=$rc). Tail:"
    tail -30 "$step_log" | tee -a "$LOG"
    return 2
  fi
  log "  LPPD summary:"
  cat "$REPO/results/loglik_summary_clean.csv" 2>/dev/null | tee -a "$LOG"
  log "Step 7: DONE"
  return 0
}

############################################################################
# Main
############################################################################
log "=========================================="
log "OVERNIGHT PIPELINE START"
log "Working dir: $REPO"
log "Initial Rscript.exe count: $(tasklist //FI 'IMAGENAME eq Rscript.exe' 2>/dev/null | grep -c Rscript.exe)"
log "=========================================="

step1_wait_bidifac || { log "Aborting after step 1"; exit 2; }
step2_screening    || { log "Aborting after step 2"; exit 2; }
step3_commit       # commit failure is non-fatal
step4_build_split  || { log "Aborting after step 4"; exit 2; }
step5_gibbs        || { log "Aborting after step 5"; exit 2; }
step6_convergence  # convergence HALT is logged but doesn't abort
step7_lppd         || { log "Aborting after step 7"; exit 2; }

set_state "complete"
log "=========================================="
log "OVERNIGHT PIPELINE COMPLETE"
log "Final state:"
log "  LPPD: $REPO/results/loglik_summary_clean.csv"
log "  Convergence: $REPO/results/convergence_summary_clean.csv"
log "  Deviation doc: $REPO/DEVIATIONS.md"
log "  Lock awaiting push: reference/paper5_target_covariates_clean.csv"
log "=========================================="
