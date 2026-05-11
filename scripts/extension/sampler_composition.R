# Operator-composition Gibbs sampler — Paper 2 extension.
#
# Identical to sampler_structural_prior.R EXCEPT for the target-covariate
# update step, which now applies one of three composition rules over a pair
# of hierarchies (h1, h2):
#
#   RULE 1 (sum):       beta_c,k = (beta_g(P_h1(c)) + beta_g(P_h2(c))) / 2
#                       Each hierarchy's per-group beta is sampled separately
#                       from pooled within-h-group residuals; the per-cancer
#                       prior is the equal-weight average. PIP = 1 iff cancer
#                       is in non-contested groups under BOTH hierarchies.
#
#   RULE 2 (sequential): beta_c,k = beta_(P_h1(c), P_h2(c))_joint
#                       Per-(h1, h2) joint group beta sampled from pooled
#                       within-joint-group residuals. PIP = 1 iff joint group
#                       has >= 2 cancers (with covariate present).
#
#   RULE 3 (intersect-with-fallback):
#                       Use joint beta if joint group has >= 2 cancers; else
#                       fall back to h1 group beta. PIP = 1 iff cancer is in
#                       non-contested h1-group.
#
# `target_covs` is a list of:
#   list(cov_id, rule (1|2|3), pair = list(h1_label, gpc_h1, h2_label, gpc_h2))
# where gpc_h1 / gpc_h2 are length-n_cancer character vectors of group labels
# (NA for contested cancers).

library(MASS)
library(truncnorm)
library(EnvStats)

HierarchicalLogNormalSpikeSlab_Composition = function(
    Covariates, Survival, Censored,
    starting_values, iters, covariates_in_model,
    covariates_by_cancer, priors, pi_generation,
    target_covs = list(),
    progress = "normal") {

  logSum <- function(l) max(l) + log(sum(exp(l - max(l))))

  Insert0sIntoVec = function(vec, avail, add_NA) {
    if (add_NA) { new_vec = rep(NA, length(avail)); new_vec[avail] = vec }
    else        { new_vec = rep(0,  length(avail)); new_vec[avail] = vec }
    new_vec
  }

  NumIncludedAndAvail = function(covariates_by_cancer, covariates_in_model, gamma_i) {
    out = sapply(covariates_in_model, function(predictor) {
      sum(unlist(sapply(1:n_cancer, function(type) {
        cov_for_type = covariates_by_cancer[[type]]
        if (predictor %in% cov_for_type) gamma_i[[type]][cov_for_type == predictor] == 1
      })))
    })
    names(out) = covariates_in_model
    out
  }

  Beta_Posterior_Helper = function(x_i, y_i, gamma_ij, sigma2_i, lambda2_i, beta_tilde_i) {
    x_i = as.matrix(x_i); y_i = as.matrix(y_i)
    alpha_i = (1/spike_priorvar) * lambda2_i
    alpha_i[gamma_ij == 0] = 1
    Variance = spike_priorvar * alpha_i
    if (p == 1) VarCovarMatrix = 1/Variance else VarCovarMatrix = diag(1/Variance)
    beta_tilde_after_adding_0s = beta_tilde_i * gamma_ij
    B = solve((1/sigma2_i)*t(x_i)%*%x_i + VarCovarMatrix, tol = 1e-35)
    b = (1/sigma2_i) * t(x_i) %*% y_i + VarCovarMatrix %*% beta_tilde_after_adding_0s
    list(mean = B %*% b, var = B)
  }

  Lambda2_Posterior_Helper = function(current_betas, current_beta_tilde, gamma_i, n_cancer) {
    tot = rep(0, p)
    for (type in 1:n_cancer) {
      ps_i = covariates_by_cancer[[type]]
      avail = covariates_in_model %in% ps_i
      gamma_ij = gamma_i[[type]]
      summand_i = (current_betas[[type]] - current_beta_tilde[avail])^2 * gamma_ij
      tot = tot + Insert0sIntoVec(summand_i, avail, add_NA = FALSE)
    }
    tot
  }

  Sigma2_Posterior_Rate_Helper = function(X, Y, current_betas, n_cancer) {
    tot = 0
    for (i in 1:n_cancer) {
      x_i = as.matrix(X[[i]]); y_i = as.matrix(Y[[i]])
      tot = tot + sum((y_i - x_i %*% current_betas[[i]])^2)
    }
    tot
  }

  Beta_Tilde_Posterior_Helper = function(current_betas_mean, current_lambda2,
                                          num_inc_avail, tau2) {
    post_mean = (num_inc_avail*tau2*current_betas_mean) / (current_lambda2 + num_inc_avail*tau2)
    post_var  = (current_lambda2*tau2)                  / (current_lambda2 + num_inc_avail*tau2)
    list(mean = post_mean, var = post_var)
  }

  X = Covariates; Y = Survival
  p = length(covariates_in_model)
  n_vec = c(unlist(lapply(X, nrow)))
  n_cancer = length(n_vec)
  I_p = sapply(covariates_in_model,
               function(i) sum(sapply(covariates_by_cancer, function(cancer) i %in% cancer)))

  tau2_intercept   = priors$betatilde_priorvar_intercept
  tau2_coefficient = priors$betatilde_priorvar_coefficient
  lambda2_priorshape_intercept   = priors$lambda2_priorshape_intercept
  lambda2_priorrate_intercept    = priors$lambda2_priorrate_intercept
  lambda2_priorshape_coefficient = priors$lambda2_priorshape_coefficient
  lambda2_priorrate_coefficient  = priors$lambda2_priorrate_coefficient
  sigma2_priorshape = priors$sigma2_priorshape
  sigma2_priorrate  = priors$sigma2_priorrate
  spike_priorvar    = priors$spike_priorvar

  sigma2_start     = starting_values$sigma2_start
  beta_tilde_start = starting_values$beta_tilde_start
  lambda2_start    = starting_values$lambda2_start
  gamma_start      = starting_values$gamma_start
  pi_start         = starting_values$pi_start

  total_obs = sum(n_vec)
  num_censored  = sapply(lapply(Y, is.na), sum)
  which_censored = sapply(Y, is.na)

  betas = lapply(1:n_cancer, function(i) list())
  beta_tilde = matrix(ncol = p, nrow = iters + 1); beta_tilde[1,] = beta_tilde_start
  lambda2    = matrix(ncol = p, nrow = iters + 1); lambda2[1,]    = lambda2_start
  sigma2 = c(); sigma2[1] = sigma2_start
  gamma  = gamma_start
  pi     = matrix(ncol = p, nrow = iters + 1); pi[1,] = pi_start

  for (k in 1:n_cancer) {
    Y[[k]][which_censored[[k]]] = Censored[[k]][which_censored[[k]]]
    Y[[k]] = log(Y[[k]])
  }

  # ---- Pre-compute target-cov index lookups + group-membership tables ----
  target_meta = lapply(target_covs, function(t) {
    k_id = which(covariates_in_model == t$cov_id)
    stopifnot(length(k_id) == 1)
    k_idx = sapply(1:n_cancer, function(c) {
      m = which(covariates_by_cancer[[c]] == t$cov_id)
      if (length(m) == 0) NA_integer_ else as.integer(m)
    })
    # Joint group label: h1__h2, NA if either side is contested
    gpc_h1 = t$pair$gpc_h1
    gpc_h2 = t$pair$gpc_h2
    joint_label = ifelse(is.na(gpc_h1) | is.na(gpc_h2), NA_character_,
                          paste(gpc_h1, gpc_h2, sep = "__"))
    list(
      cov_id = t$cov_id,
      rule = t$rule,
      h1_label = t$pair$h1_label,
      h2_label = t$pair$h2_label,
      k_id_in_model = k_id,
      k_idx_in_cancer = k_idx,
      gpc_h1 = gpc_h1,
      gpc_h2 = gpc_h2,
      joint = joint_label,
      h1_groups = unique(gpc_h1[!is.na(gpc_h1)]),
      h2_groups = unique(gpc_h2[!is.na(gpc_h2)]),
      joint_groups = unique(joint_label[!is.na(joint_label)])
    )
  })

  if (length(target_covs) > 0) {
    cat(sprintf("[COMPOSITION SAMPLER ACTIVE] %d target covariates:\n", length(target_covs)))
    for (tm in target_meta) {
      n_h1 = sum(!is.na(tm$gpc_h1)); n_h2 = sum(!is.na(tm$gpc_h2))
      n_jt = sum(!is.na(tm$joint))
      cat(sprintf("  cov %s, RULE %d, %s x %s: h1 n_groups=%d (n=%d), h2 n_groups=%d (n=%d), joint n_groups=%d (n=%d)\n",
                  tm$cov_id, tm$rule, tm$h1_label, tm$h2_label,
                  length(tm$h1_groups), n_h1,
                  length(tm$h2_groups), n_h2,
                  length(tm$joint_groups), n_jt))
    }
  }

  # Sample group-shared beta given (cancers in group, pooled residuals + x for cov k)
  sample_group_beta = function(cancers_in_g, kc_per_cancer, prior_mean, prior_var, sig2) {
    if (length(cancers_in_g) == 0) return(NA_real_)
    x_pool = numeric(0); r_pool = numeric(0)
    for (c in cancers_in_g) {
      kc = kc_per_cancer[c]
      if (is.na(kc)) next
      X_c = as.matrix(X[[c]]); y_c = as.numeric(Y[[c]])
      beta_c = betas[[c]][[i]]
      full_resid = y_c - as.numeric(X_c %*% beta_c)
      r_for_k = full_resid + X_c[, kc] * beta_c[kc]
      x_pool = c(x_pool, X_c[, kc])
      r_pool = c(r_pool, r_for_k)
    }
    if (length(x_pool) == 0) return(NA_real_)
    post_prec = (1/prior_var) + (sum(x_pool^2)/sig2)
    post_var  = 1/post_prec
    post_mean = post_var * ((prior_mean/prior_var) + sum(x_pool*r_pool)/sig2)
    rnorm(1, post_mean, sqrt(post_var))
  }

  for (i in 1:iters) {
    if (progress == "normal") svMisc::progress(i/(iters/100))
    else if (progress == "model_comparison") {
      if (i %% 10000 == 0) cat(paste("Iter", i), "\n")
    }

    gamma_i = sapply(gamma, `[`, i)
    number_included_and_available = NumIncludedAndAvail(covariates_by_cancer,
                                                        covariates_in_model, gamma_i)

    # ---- Step 1: standard per-cancer beta update ----
    for (j in 1:n_cancer) {
      ps_j = covariates_by_cancer[[j]]
      avail = covariates_in_model %in% ps_j
      lambda2_i = lambda2[i,][avail]
      beta_tilde_i = beta_tilde[i,][avail]
      gamma_ij = gamma[[j]][[i]]
      beta_post_params = Beta_Posterior_Helper(X[[j]], Y[[j]], gamma_ij,
                                                sigma2[i], lambda2_i, beta_tilde_i)
      beta_j_gen = mvrnorm(1, mu = beta_post_params$mean, Sigma = beta_post_params$var)
      betas[[j]][[i]] = beta_j_gen
    }

    # ---- Step 1.5 [COMPOSITION]: override target covariates ----
    if (length(target_covs) > 0) {
      sig2 = sigma2[i]
      for (tm in target_meta) {
        k_id = tm$k_id_in_model
        prior_mean = beta_tilde[i, k_id]
        prior_var  = lambda2[i, k_id]

        # Sample group-shared betas as needed by the rule
        beta_h1 = setNames(numeric(length(tm$h1_groups)), tm$h1_groups)
        beta_h2 = setNames(numeric(length(tm$h2_groups)), tm$h2_groups)
        beta_joint = setNames(numeric(length(tm$joint_groups)), tm$joint_groups)
        joint_size = setNames(integer(length(tm$joint_groups)), tm$joint_groups)

        if (tm$rule %in% c(1L, 3L)) {
          for (g in tm$h1_groups) {
            cancers_in_g = which(tm$gpc_h1 == g & !is.na(tm$gpc_h1))
            cancers_in_g = cancers_in_g[!is.na(tm$k_idx_in_cancer[cancers_in_g])]
            beta_h1[g] = sample_group_beta(cancers_in_g, tm$k_idx_in_cancer,
                                            prior_mean, prior_var, sig2)
          }
        }
        if (tm$rule == 1L) {
          for (g in tm$h2_groups) {
            cancers_in_g = which(tm$gpc_h2 == g & !is.na(tm$gpc_h2))
            cancers_in_g = cancers_in_g[!is.na(tm$k_idx_in_cancer[cancers_in_g])]
            beta_h2[g] = sample_group_beta(cancers_in_g, tm$k_idx_in_cancer,
                                            prior_mean, prior_var, sig2)
          }
        }
        if (tm$rule %in% c(2L, 3L)) {
          for (g in tm$joint_groups) {
            cancers_in_g = which(tm$joint == g & !is.na(tm$joint))
            cancers_in_g = cancers_in_g[!is.na(tm$k_idx_in_cancer[cancers_in_g])]
            joint_size[g] = length(cancers_in_g)
            if (length(cancers_in_g) >= 2) {
              beta_joint[g] = sample_group_beta(cancers_in_g, tm$k_idx_in_cancer,
                                                 prior_mean, prior_var, sig2)
            } else {
              beta_joint[g] = NA_real_
            }
          }
        }

        # Apply per-cancer overrides per the rule
        for (c in 1:n_cancer) {
          kc = tm$k_idx_in_cancer[c]
          if (is.na(kc)) next   # cancer doesn't have this cov
          new_b = NA_real_
          if (tm$rule == 1L) {
            g1 = tm$gpc_h1[c]; g2 = tm$gpc_h2[c]
            if (!is.na(g1) && !is.na(g2)) {
              new_b = (beta_h1[g1] + beta_h2[g2]) / 2
            } else {
              new_b = 0
            }
          } else if (tm$rule == 2L) {
            gj = tm$joint[c]
            if (!is.na(gj) && joint_size[gj] >= 2) {
              new_b = beta_joint[gj]
            } else {
              new_b = 0
            }
          } else if (tm$rule == 3L) {
            gj = tm$joint[c]; g1 = tm$gpc_h1[c]
            if (!is.na(gj) && joint_size[gj] >= 2) {
              new_b = beta_joint[gj]
            } else if (!is.na(g1)) {
              new_b = beta_h1[g1]
            } else {
              new_b = 0
            }
          }
          if (is.na(new_b)) new_b = 0   # belt-and-suspenders
          betas[[c]][[i]][kc] = new_b
        }
      }
    }

    # ---- Step 2: lambda^2 update ----
    current_betas = sapply(betas, `[`, i)
    current_beta_tilde = beta_tilde[i,]
    W = Lambda2_Posterior_Helper(current_betas, current_beta_tilde, gamma_i, n_cancer)
    post_lambda2_shape = (number_included_and_available/2) +
      c(lambda2_priorshape_intercept, rep(lambda2_priorshape_coefficient, p-1))
    post_lambda2_rate = c(lambda2_priorrate_intercept,
                          rep(lambda2_priorrate_coefficient, p-1)) + 0.5*W
    lambda2[i+1,] = 1/rgamma(p, shape = post_lambda2_shape, rate = post_lambda2_rate)

    # ---- Step 3: beta_tilde update ----
    current_betas_augmented = lapply(1:n_cancer, function(k) {
      betas_i = current_betas[[k]]
      gamma_ij = gamma_i[[k]]
      betas_i[gamma_ij == 0] = NA
      ps_k = covariates_by_cancer[[k]]
      avail = c(covariates_in_model %in% ps_k)
      Insert0sIntoVec(betas_i, avail, add_NA = TRUE)
    })
    current_betas_mean = colMeans(do.call(rbind, current_betas_augmented), na.rm = TRUE)
    current_betas_mean[is.nan(current_betas_mean)] = 0
    current_lambda2 = lambda2[i+1,]
    post_beta_tilde = Beta_Tilde_Posterior_Helper(
      current_betas_mean, current_lambda2, number_included_and_available,
      tau2 = c(tau2_intercept, rep(tau2_coefficient, p-1)))
    beta_tilde[i+1,] = rnorm(p, post_beta_tilde$mean, sqrt(post_beta_tilde$var))

    # ---- Step 4: sigma^2 update ----
    post_sigma2_shape = (total_obs/2) + sigma2_priorshape
    post_sigma2_rate  = 0.5*Sigma2_Posterior_Rate_Helper(X, Y, current_betas, n_cancer) +
                        sigma2_priorrate
    sigma2[i+1] = 1/rgamma(1, post_sigma2_shape, rate = post_sigma2_rate)

    # ---- Step 5: pi update ----
    if (pi_generation == "shared_across_cancers") {
      n_inc = NumIncludedAndAvail(covariates_by_cancer, covariates_in_model, gamma_i)
      pi[i+1,] = rbeta(p, 1 + n_inc, 1 + I_p - n_inc)
    } else if (pi_generation == "fixed_at_0.5") {
      pi[i+1,] = rep(0.5, p)
    } else if (pi_generation == "fixed_at_1.0") {
      pi[i+1,] = rep(1.0, p)
    } else if (pi_generation == "shared_across_betas") {
      n_inc_all = sum(NumIncludedAndAvail(covariates_by_cancer, covariates_in_model, gamma_i)[-1])
      sum_I_p = sum(I_p[-1])
      pi[i+1,] = rbeta(1, 1 + n_inc_all, 1 + sum_I_p - n_inc_all)
    } else if (pi_generation == "fixed_at_0.0") {
      pi[i+1,] = c(1, rep(0, p-1))
    }

    # ---- Step 6: gamma update (standard) ----
    for (j in 1:n_cancer) {
      ps_j = covariates_by_cancer[[j]]
      avail = covariates_in_model %in% ps_j
      pi_j = pi[i+1,][avail]
      beta_tilde_i = beta_tilde[i+1,][avail]
      lambda2_i = lambda2[i+1,][avail]

      like_slab  = dnorm(current_betas[[j]], mean = beta_tilde_i, sd = sqrt(lambda2_i),    log = TRUE)
      like_spike = dnorm(current_betas[[j]], mean = 0,            sd = sqrt(spike_priorvar), log = TRUE)
      probs_j = sapply(1:length(ps_j), function(k) {
        x = log(pi_j[k]) + like_slab[k]
        y = log(1 - pi_j[k]) + like_spike[k]
        exp(x - logSum(c(x,y)))
      })
      new_gammas = rbinom(length(ps_j), 1, probs_j)
      new_gammas[1] = 1
      gamma[[j]][[i+1]] = new_gammas
    }

    # ---- Step 6.5 [COMPOSITION]: force gamma for target covariates ----
    if (length(target_covs) > 0) {
      for (tm in target_meta) {
        for (c in 1:n_cancer) {
          kc = tm$k_idx_in_cancer[c]
          if (is.na(kc)) next
          new_g = 0L
          if (tm$rule == 1L) {
            if (!is.na(tm$gpc_h1[c]) && !is.na(tm$gpc_h2[c])) new_g = 1L
          } else if (tm$rule == 2L) {
            gj = tm$joint[c]
            if (!is.na(gj)) {
              joint_size = sum(tm$joint == gj & !is.na(tm$joint))
              if (joint_size >= 2) new_g = 1L
            }
          } else if (tm$rule == 3L) {
            if (!is.na(tm$gpc_h1[c])) new_g = 1L
          }
          gamma[[c]][[i+1]][kc] = new_g
        }
      }
    }

    # ---- Step 7: censored Y imputation ----
    for (k in 1:n_cancer) {
      n_gens = num_censored[[k]]
      Censored_Obs = X[[k]][which_censored[[k]],]
      Censor_Lower_Bound = Censored[[k]][which_censored[[k]]]
      if (n_gens != 0) {
        Mu_Survival = if (length(current_betas[[k]]) == 1) Censored_Obs * current_betas[[k]]
                      else Censored_Obs %*% current_betas[[k]]
        random_survival = rtruncnorm(n_gens, a = log(Censor_Lower_Bound),
                                     mean = Mu_Survival, sd = rep(sqrt(sigma2[i]), n_gens))
        Y[[k]][which_censored[[k]]] = random_survival
      }
    }
  }

  list(betas = betas, beta_tilde = beta_tilde, lambda2 = lambda2, sigma2 = sigma2,
       gamma = gamma, pi = pi)
}
