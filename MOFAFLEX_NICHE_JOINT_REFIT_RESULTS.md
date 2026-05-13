# NICHE_JOINT_REFIT_RESULTS.md

**Verdict: `FALSIFIED_JOINT_REFIT`**

## Disposition (paired Δ vs SANITY refit)

     spec  n_seeds  mean_paired_delta  sd_paired_delta  min_paired_delta  max_paired_delta
  PRIMARY        3          -0.844846         1.220186         -1.782879          0.534623
SECONDARY        3          -0.306129         0.929464         -0.913831          0.763834
     NULL        3          -0.631880         1.933498         -2.842664          0.743215

## Sanity refit drift vs on-disk reconverged baseline

      seed  delta_xen_vs_reconverged  delta_chrom_vs_reconverged
2511021636                 -0.124091                  -11.696605
2511021637                 -0.343109                    1.353778
2511021638                 -0.113100                  -21.967176

## Diagnostics

- **primary_paired_delta**: -0.8448461492856344
- **primary_sd**: 1.220185547193043
- **null_paired_delta**: -0.6318801591793696
- **null_sd**: 1.9334977945293501
- **null_z_separation**: -0.09314787099339486
- **sanity_drift_max_abs**: 21.96717643737793
- **cascade_A_sanity_drift_under_50**: True
- **passes_1_paired_delta_above_0p5**: False
- **passes_2_null_z_above_2**: False
- **passes_3_repro_sd_under_5**: True

## Per-seed paired Δ

      seed      spec  paired_delta_mean  paired_delta_xen_mean  paired_delta_chrom_mean  n_paired
2511021636   PRIMARY          -1.782879              -0.189198               -10.507591     39333
2511021636 SECONDARY          -0.768389              -0.167944                -4.055567     39333
2511021636      NULL           0.203809              -0.074619                 1.728076     39333
2511021637   PRIMARY          -1.286282              -0.063845                -7.978598     39333
2511021637 SECONDARY          -0.913831              -0.025449                -5.777338     39333
2511021637      NULL          -2.842664               0.059051               -18.728302     39333
2511021638   PRIMARY           0.534623              -0.167786                 4.380005     39333
2511021638 SECONDARY           0.763834              -0.039001                 5.159009     39333
2511021638      NULL           0.743215              -0.210286                 5.963221     39333