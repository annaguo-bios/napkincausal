#' Estimate E(Y(x)) under an intervention density for Z with covariates C
#'
#' Covariate-aware analogue of `napkin.a()`.
#'
#' The target is \eqn{E_C\{\kappa_1(C; \tilde p_z) / \kappa_2(C; \tilde p_z)\}},
#' where the inner kappas integrate over W given C and the intervention density
#' \eqn{\tilde p(z)}.
#'
#' @param x Treatment level.
#' @param z Optional intervention level for binary Z. If NULL and Z is binary,
#'   estimates are returned for z = 0, z = 1, and a combined target.
#' @param data Data frame containing the observed variables.
#' @param treatment Character string naming the treatment variable.
#' @param Z.variables Character string naming the single Z variable.
#' @param W.variables Character string naming the single binary W variable.
#' @param covariates Character vector naming baseline covariates C.
#' @param outcome Character string naming the outcome variable.
#' @param z.density Optional function for the intervention density \eqn{\tilde p(z)}
#'   when Z is continuous. If NULL, the empirical density of observed Z is used.
#' @param z_w.method Method for estimating the conditional density or probability
#'   of Z given W and C. Supported values are "dnorm", "np", or a user-specified
#'   density function.
#' @param superlearner.Y Logical; use SuperLearner for the outcome regression.
#' @param superlearner.X Logical; use SuperLearner for the treatment mechanism.
#' @param superlearner.Z Logical; use SuperLearner for the Z model.
#' @param superlearner.W Logical; use SuperLearner for the W given C model.
#' @param crossfit Logical; use cross-fitted SuperLearner nuisance estimates.
#' @param K Number of folds for cross-fitting.
#' @param lib.Y SuperLearner library for the outcome regression.
#' @param lib.X SuperLearner library for the treatment mechanism.
#' @param lib.Z SuperLearner library for the Z model.
#' @param lib.W SuperLearner library for the W given C model.
#' @param n.iter Maximum number of targeting iterations.
#' @param cvg.criteria Convergence criterion for mean EIF components.
#' @param formula.Y Regression formula for the outcome regression.
#' @param formula.X Regression formula for the treatment mechanism.
#' @param formula.Z Regression formula for the Z model.
#' @param formula.W Regression formula for the W given C model.
#' @param linkY_binary Link for binary outcome regression when using glm.
#' @param link.X Link for treatment mechanism glm.
#' @param link.Z Link for binary Z glm.
#' @param link.W Link for binary W glm.
#' @param truncate_lower.X Lower truncation bound for treatment probabilities.
#' @param truncate_upper.X Upper truncation bound for treatment probabilities.
#' @param truncate_lower.Z Lower truncation bound for Z densities/probabilities.
#' @param truncate_upper.Z Upper truncation bound for Z densities/probabilities.
#' @param truncate_lower.W Lower truncation bound for W probabilities.
#' @param truncate_upper.W Upper truncation bound for W probabilities.
#' @param minZ Lower integration/sampling bound for continuous Z.
#' @param maxZ Upper integration/sampling bound for continuous Z.
#' @param verbose Logical; print progress messages.
#' @param fast Logical retained for interface compatibility.
#' @param nMC Number of Monte Carlo draws for continuous Z integration.
#' @param boundedsubmodelY Logical retained for interface compatibility; not
#'   implemented in this covariate-aware draft.
#' @param use.kappa2.aipw Logical; use the AIPW version of kappa2(C) in
#'   denominators for the continuous Z covariate-aware estimator.
#' @return A list with `Onestep` and `TMLE` components.
#' @export
napkin.a.covariate <- function(x, z = NULL, data, treatment, Z.variables,
                               W.variables, covariates, outcome,
                               z.density = NULL, z_w.method = "dnorm",
                               superlearner.Y = FALSE,
                               superlearner.X = FALSE,
                               superlearner.Z = FALSE,
                               superlearner.W = FALSE,
                               crossfit = FALSE, K = 5,
                               lib.Y = c("SL.glm", "SL.earth", "SL.ranger", "SL.mean"),
                               lib.X = c("SL.glm", "SL.earth", "SL.ranger", "SL.mean"),
                               lib.Z = c("SL.glm", "SL.earth", "SL.ranger", "SL.mean"),
                               lib.W = c("SL.glm", "SL.earth", "SL.ranger", "SL.mean"),
                               n.iter = 500, cvg.criteria = 0.01,
                               formula.Y = "Y ~ .", formula.X = "X ~ .",
                               formula.Z = "Z ~ .", formula.W = "W ~ .",
                               linkY_binary = "logit", link.X = "logit",
                               link.Z = "logit", link.W = "logit",
                               truncate_lower.X = 0, truncate_upper.X = 1,
                               truncate_lower.Z = 0, truncate_upper.Z = 1,
                               truncate_lower.W = 0, truncate_upper.W = 1,
                               minZ = -Inf, maxZ = Inf,
                               verbose = TRUE, fast = TRUE, nMC = 5000,
                               boundedsubmodelY = FALSE,
                               use.kappa2.aipw = FALSE) {

  if (length(Z.variables) != 1) {
    stop("napkin.a.covariate currently supports a single Z variable.")
  }
  if (length(W.variables) != 1) {
    stop("napkin.a.covariate currently supports a single W variable.")
  }
  if (boundedsubmodelY) {
    warning("boundedsubmodelY is present for interface compatibility; it is not implemented in this draft.")
  }

  n <- nrow(data)

  C <- data[, covariates, drop = FALSE]
  X <- data[, treatment]
  W <- data[, W.variables]
  Z <- data[, Z.variables]
  Y <- data[, outcome]

  binaryY <- all(Y %in% c(0, 1))
  binaryZ <- all(as.vector(Z) %in% c(0, 1))

  binaryW <- all(W %in% c(0, 1))
  if (!binaryW) {
    stop("napkin.a.covariate currently supports binary numeric W only.")
  }

  formula.Y <- as.formula(formula.Y)
  formula.X <- as.formula(formula.X)
  formula.Z <- as.formula(formula.Z)
  formula.W <- as.formula(formula.W)

  dat_mpY <- data.frame(Y = Y, X = X, Z = Z, W = W, C)
  dat_mpX <- data.frame(X = X, Z = Z, W = W, C)
  dat_mpZ <- data.frame(Z = Z, W = W, C)
  dat_mpW <- data.frame(W = W, C)

  predict_cv_sl <- function(fit, newdata) {
    pred <- unlist(lapply(seq_len(K), function(k) {
      predict(fit$AllSL[[k]], newdata = newdata[fit$folds[[k]], , drop = FALSE])[[1]]
    }))
    pred[order(unlist(fit$folds))]
  }

  make_newdata_y <- function(x, z, w, c) {
    data.frame(Z = z, X = x, W = w, setNames(as.data.frame(c), covariates))
  }

  make_newdata_x <- function(z, w, c) {
    data.frame(Z = z, W = w, setNames(as.data.frame(c), covariates))
  }

  make_newdata_z <- function(w, c) {
    data.frame(W = w, setNames(as.data.frame(c), covariates))
  }

  make_newdata_w <- function(c) {
    data.frame(setNames(as.data.frame(c), covariates))
  }

  ################################################
  ############### OUTCOME REGRESSION #############
  ################################################

  fit.family.Y <- if (binaryY) binomial(linkY_binary) else gaussian()

  if (crossfit) {
    or_fit <- CV.SuperLearner(
      Y = Y, X = dat_mpY[, -1, drop = FALSE], family = fit.family.Y,
      V = K, SL.library = lib.Y, control = list(saveFitLibrary = TRUE),
      saveAll = TRUE, env = asNamespace("SuperLearner")
    )
    f.mu.xzwc <- function(x, z, w, c) {
      as.vector(predict_cv_sl(or_fit, make_newdata_y(x, z, w, c)))
    }
  } else if (superlearner.Y) {
    or_fit <- SuperLearner(
      Y = Y, X = dat_mpY[, -1, drop = FALSE], family = fit.family.Y,
      SL.library = lib.Y, env = asNamespace("SuperLearner")
    )
    f.mu.xzwc <- function(x, z, w, c) {
      as.vector(predict(or_fit, newdata = make_newdata_y(x, z, w, c))[[1]])
    }
  } else {
    or_fit <- glm(formula.Y, data = dat_mpY, family = fit.family.Y)
    f.mu.xzwc <- function(x, z, w, c) {
      as.vector(predict(or_fit, newdata = make_newdata_y(x, z, w, c), type = "response"))
    }
  }

  if (verbose) cat("Outcome regression Y|X,Z,W,C done.\n\n")

  ################################################
  ############### PROPENSITY SCORE ###############
  ################################################

  if (crossfit) {
    ps_fit <- CV.SuperLearner(
      Y = X, X = dat_mpX[, -1, drop = FALSE], family = binomial(),
      V = K, SL.library = lib.X, control = list(saveFitLibrary = TRUE),
      saveAll = TRUE, env = asNamespace("SuperLearner")
    )
    f.x.zwc <- function(x, z, w, c, truncate_lower = 0, truncate_upper = 1) {
      p.x1 <- as.vector(predict_cv_sl(ps_fit, make_newdata_x(z, w, c)))
      p.x1 <- pmin(pmax(p.x1, truncate_lower), truncate_upper)
      x * p.x1 + (1 - x) * (1 - p.x1)
    }
  } else if (superlearner.X) {
    ps_fit <- SuperLearner(
      Y = X, X = dat_mpX[, -1, drop = FALSE], family = binomial(),
      SL.library = lib.X, env = asNamespace("SuperLearner")
    )
    f.x.zwc <- function(x, z, w, c, truncate_lower = 0, truncate_upper = 1) {
      p.x1 <- as.vector(predict(ps_fit, newdata = make_newdata_x(z, w, c), type = "response")[[1]])
      p.x1 <- pmin(pmax(p.x1, truncate_lower), truncate_upper)
      x * p.x1 + (1 - x) * (1 - p.x1)
    }
  } else {
    ps_fit <- glm(formula.X, data = dat_mpX, family = binomial(link.X))
    f.x.zwc <- function(x, z, w, c, truncate_lower = 0, truncate_upper = 1) {
      p.x1 <- as.vector(predict(ps_fit, newdata = make_newdata_x(z, w, c), type = "response"))
      p.x1 <- pmin(pmax(p.x1, truncate_lower), truncate_upper)
      x * p.x1 + (1 - x) * (1 - p.x1)
    }
  }

  if (verbose) cat("Propensity score X|Z,W,C done.\n\n")

  #########################################################
  ############### f_Z(Z | W, C) DENSITY ###################
  #########################################################

  if (binaryZ) {
    if (crossfit) {
      z_fit <- CV.SuperLearner(
        Y = Z, X = dat_mpZ[, -1, drop = FALSE], family = binomial(),
        V = K, SL.library = lib.Z, control = list(saveFitLibrary = TRUE),
        saveAll = TRUE, env = asNamespace("SuperLearner")
      )
      p.z1.wc <- as.vector(z_fit$SL.predict)
    } else if (superlearner.Z) {
      z_fit <- SuperLearner(
        Y = Z, X = dat_mpZ[, -1, drop = FALSE], family = binomial(),
        SL.library = lib.Z, env = asNamespace("SuperLearner")
      )
      p.z1.wc <- as.vector(z_fit$SL.predict)
    } else {
      z_fit <- glm(formula.Z, data = dat_mpZ, family = binomial(link.Z))
      p.z1.wc <- as.vector(predict(z_fit, newdata = dat_mpZ, type = "response"))
    }
    p.z1.wc <- pmin(pmax(p.z1.wc, truncate_lower.Z), truncate_upper.Z)
    f.z.binary.obs <- function(z) {
      ifelse(z == 1, p.z1.wc, 1 - p.z1.wc)
    }
    p.z.wc <- f.z.binary.obs(Z)
  } else if (z_w.method == "np") {
    bw <- np::npcdensbw(formula = formula.Z, data = dat_mpZ)
    z_fit <- np::npcdens(bws = bw)
    p.z.wc <- as.vector(predict(z_fit))
  } else if (z_w.method == "dnorm") {
    if (crossfit) {
      z_fit <- CV.SuperLearner(
        Y = Z, X = dat_mpZ[, -1, drop = FALSE], family = gaussian(),
        V = K, SL.library = lib.Z, control = list(saveFitLibrary = TRUE),
        saveAll = TRUE, env = asNamespace("SuperLearner")
      )
      z.mean <- as.vector(z_fit$SL.predict)
    } else if (superlearner.Z) {
      z_fit <- SuperLearner(
        Y = Z, X = dat_mpZ[, -1, drop = FALSE], family = gaussian(),
        SL.library = lib.Z, env = asNamespace("SuperLearner")
      )
      z.mean <- as.vector(z_fit$SL.predict)
    } else {
      z_fit <- glm(formula.Z, data = dat_mpZ, family = gaussian())
      z.mean <- as.vector(predict(z_fit, newdata = dat_mpZ))
    }
    z.sd <- sd(Z - z.mean)
    p.z.wc <- dnorm(Z, mean = z.mean, sd = z.sd)
  } else if (is.function(z_w.method)) {
    p.z.wc <- z_w.method(Z.variables = Z, W.variables = W, covariates = C)
  } else {
    stop("z_w.method must be 'dnorm', 'np', or a density function.")
  }

  p.z.wc <- pmin(pmax(p.z.wc, truncate_lower.Z), truncate_upper.Z)

  if (verbose) cat("Conditional density Z|W,C done.\n\n")

  #########################################################
  ############### f_W(W | C) MODEL ########################
  #########################################################

  if (crossfit) {
    w_fit <- CV.SuperLearner(
      Y = W, X = C, family = binomial(), V = K, SL.library = lib.W,
      control = list(saveFitLibrary = TRUE), saveAll = TRUE, env = asNamespace("SuperLearner")
    )
    p.w1.c <- pmin(pmax(as.vector(w_fit$SL.predict), truncate_lower.W), truncate_upper.W)
  } else if (superlearner.W) {
    w_fit <- SuperLearner(Y = W, X = C, family = binomial(),
                                        SL.library = lib.W, env = asNamespace("SuperLearner"))
    p.w1.c <- pmin(pmax(as.vector(w_fit$SL.predict), truncate_lower.W), truncate_upper.W)
  } else {
    w_fit <- glm(formula.W, data = dat_mpW, family = binomial(link.W))
    p.w1.c <- pmin(pmax(as.vector(predict(w_fit, newdata = dat_mpW, type = "response")),
                        truncate_lower.W), truncate_upper.W)
  }
  p.w0.c <- 1 - p.w1.c
  w.levels <- c(0, 1)
  w.prob <- list(`0` = p.w0.c, `1` = p.w1.c)

  if (verbose) cat("Conditional model W|C done.\n\n")

  ################################################
  ############### TARGET DENSITY #################
  ################################################

  if (!binaryZ && is.null(z.density)) {
    den.z <- density(Z)
    z.density <- function(z) approx(den.z$x, den.z$y, xout = z, rule = 2)$y
    p.z <- z.density(Z)
    Zsim <- sample(Z, nMC, replace = TRUE)
  } else if (!binaryZ) {
    p.z <- z.density(Z)
    z.range <- range(Z)
    if (is.finite(minZ) && is.finite(maxZ)) z.range <- c(minZ, maxZ)
    Zsim <- RVCompare::sampleFromDensity(z.density, nMC, z.range)
  }

  if (!binaryZ && verbose) cat("Target density tilde p(Z) evaluated.\n\n")

  ################################################
  ################## ESTIMATORS ##################
  ################################################

  make_state <- function(x) {
    # Observed nuisance predictions.
    mu.obs <- f.mu.xzwc(x, Z, W, C)
    pi.obs <- f.x.zwc(x, Z, W, C,
                      truncate_lower = max(truncate_lower.X, 0.001),
                      truncate_upper = min(truncate_upper.X, 0.999))

    # Nuisance predictions for the W|C and tilde-Z plug-in integration.
    mu.sim <- replicate(length(w.levels), matrix(NA_real_, n, nMC), simplify = FALSE)
    pi.sim <- replicate(length(w.levels), matrix(NA_real_, n, nMC), simplify = FALSE)
    names(mu.sim) <- names(pi.sim) <- as.character(w.levels)
    for (w in w.levels) {
      wi <- as.character(w)
      for (j in seq_len(nMC)) {
        mu.sim[[wi]][, j] <- f.mu.xzwc(x, Zsim[j], w, C)
        pi.sim[[wi]][, j] <- f.x.zwc(
          x, Zsim[j], w, C,
          truncate_lower = max(truncate_lower.X, 0.001),
          truncate_upper = min(truncate_upper.X, 0.999)
        )
      }
    }

    # Nuisance predictions for Phi_W: integrate over tilde-Z at observed W_i.
    mu.obsW.sim <- matrix(NA_real_, n, nMC)
    pi.obsW.sim <- matrix(NA_real_, n, nMC)
    for (j in seq_len(nMC)) {
      mu.obsW.sim[, j] <- f.mu.xzwc(x, Zsim[j], W, C)
      pi.obsW.sim[, j] <- f.x.zwc(
        x, Zsim[j], W, C,
        truncate_lower = max(truncate_lower.X, 0.001),
        truncate_upper = min(truncate_upper.X, 0.999)
      )
    }

    list(mu.obs = mu.obs, pi.obs = pi.obs, mu.sim = mu.sim, pi.sim = pi.sim,
         mu.obsW.sim = mu.obsW.sim, pi.obsW.sim = pi.obsW.sim)
  }

  make_state_binaryZ <- function(x, z) {
    # Observed nuisance predictions at the fixed binary Z value.
    mu.obs <- f.mu.xzwc(x, z, W, C)
    pi.obs <- f.x.zwc(x, z, W, C,
                      truncate_lower = max(truncate_lower.X, 0.001),
                      truncate_upper = min(truncate_upper.X, 0.999))

    # Nuisance predictions for integrating W|C at the fixed z.
    mu.w <- pi.w <- vector("list", length(w.levels))
    names(mu.w) <- names(pi.w) <- as.character(w.levels)
    for (w in w.levels) {
      wi <- as.character(w)
      mu.w[[wi]] <- f.mu.xzwc(x, z, w, C)
      pi.w[[wi]] <- f.x.zwc(
        x, z, w, C,
        truncate_lower = max(truncate_lower.X, 0.001),
        truncate_upper = min(truncate_upper.X, 0.999)
      )
    }
    mu.obsW <- mu.obs
    pi.obsW <- pi.obs

    out <- list(mu.obs = mu.obs, pi.obs = pi.obs, mu.obsW = mu.obsW, pi.obsW = pi.obsW)
    out$mu.w <- mu.w
    out$pi.w <- pi.w
    out
  }

  compute_kappas <- function(state) {
    # kappas contains C-specific plug-in pieces under W|C and tilde p(z).
    kappa1.mat <- matrix(0, n, nMC)
    kappa2.mat <- matrix(0, n, nMC)
    for (w in w.levels) {
      wi <- as.character(w)
      kappa1.mat <- kappa1.mat +
        sweep(state$mu.sim[[wi]] * state$pi.sim[[wi]], 1, w.prob[[wi]], "*")
      kappa2.mat <- kappa2.mat +
        sweep(state$pi.sim[[wi]], 1, w.prob[[wi]], "*")
    }

    kappa1 <- rowMeans(kappa1.mat)
    kappa2 <- pmax(rowMeans(kappa2.mat), 0.001)
    psi.c <- kappa1 / kappa2

    list(kappa1 = kappa1, kappa2 = kappa2, psi.c = psi.c,
         kappa1.mat = kappa1.mat, kappa2.mat = kappa2.mat)
  }

  compute_kappas_binaryZ <- function(state) {
    # kappas contains C-specific plug-in pieces at the fixed z, integrating only W|C.
    kappa1 <- rep(0, n)
    kappa2 <- rep(0, n)
    for (w in w.levels) {
      wi <- as.character(w)
      kappa1 <- kappa1 + w.prob[[wi]] * state$mu.w[[wi]] * state$pi.w[[wi]]
      kappa2 <- kappa2 + w.prob[[wi]] * state$pi.w[[wi]]
    }

    kappa2 <- pmax(kappa2, 0.001)
    psi.c <- kappa1 / kappa2

    list(kappa1 = kappa1, kappa2 = kappa2, psi.c = psi.c)
  }

  regress_on_c <- function(target) {
    if (superlearner.W) {
      fit <- SuperLearner(Y = target, X = C, family = gaussian(),
                                        SL.library = lib.W, env = asNamespace("SuperLearner"))
      as.vector(predict(fit, newdata = C)[[1]])
    } else {
      fit <- lm(target ~ ., data = data.frame(target = target, C))
      as.vector(predict(fit, newdata = C))
    }
  }

  get_pi_tilde_observed_w <- function(state) {
    rowMeans(state$pi.obsW.sim)
  }

  get_kappa2_aipw <- function(x, state) {
    target <- (p.z / p.z.wc) * ((X == x) - state$pi.obs) + get_pi_tilde_observed_w(state)
    pmax(regress_on_c(target), 0.001)
  }

  get_kappa2_denom <- function(x, state, kappas) {
    if (use.kappa2.aipw) return(get_kappa2_aipw(x, state))
    kappas$kappa2
  }

  get_kappa2_aipw_binaryZ <- function(x, z, state) {
    fz.wc <- f.z.binary.obs(z)
    target <- ((Z == z) / fz.wc) * ((X == x) - state$pi.obs) + state$pi.obsW
    pmax(regress_on_c(target), 0.001)
  }

  compute_eif <- function(x, state, kappas, kappa2.denom = NULL) {
    if (is.null(kappa2.denom)) {
      kappa2.denom <- get_kappa2_denom(x, state, kappas)
    }

    EIF.Y <- (X == x) * p.z / (kappa2.denom * p.z.wc) * (Y - state$mu.obs)
    EIF.X <- p.z / (kappa2.denom * p.z.wc) *
      (state$mu.obs - kappas$psi.c) * ((X == x) - state$pi.obs)
    EIF.W <- rowMeans(state$pi.obsW.sim *
                        sweep(state$mu.obsW.sim, 1, kappas$psi.c, "-")) /
      kappa2.denom
    EIF.C <- kappas$psi.c - mean(kappas$psi.c)
    EIF <- EIF.Y + EIF.X + EIF.W + EIF.C
    estimated <- mean(kappas$psi.c)

    list(estimated = estimated, EIF = EIF, EIF.Y = EIF.Y, EIF.X = EIF.X,
         EIF.W = EIF.W, EIF.C = EIF.C, kappa2.denom = kappa2.denom)
  }

  compute_eif_binaryZ <- function(x, z, state, kappas, kappa2.denom = NULL) {
    if (is.null(kappa2.denom)) {
      kappa2.denom <- get_kappa2_aipw_binaryZ(x, z, state)
    }

    fz.wc <- f.z.binary.obs(z)
    EIF.Y <- (X == x) * (Z == z) / (kappa2.denom * fz.wc) * (Y - state$mu.obs)
    EIF.X <- (Z == z) / (kappa2.denom * fz.wc) *
      (state$mu.obs - kappas$psi.c) * ((X == x) - state$pi.obs)
    EIF.W <- state$pi.obsW * (state$mu.obsW - kappas$psi.c) / kappa2.denom
    EIF.C <- kappas$psi.c - mean(kappas$psi.c)
    EIF <- EIF.Y + EIF.X + EIF.W + EIF.C
    estimated <- mean(kappas$psi.c)

    list(estimated = estimated, EIF = EIF, EIF.Y = EIF.Y, EIF.X = EIF.X,
         EIF.W = EIF.W, EIF.C = EIF.C, kappa2.denom = kappa2.denom)
  }

  f.onestep.x_nonbinaryZ <- function(x) {
    state <- make_state(x)
    kappas <- compute_kappas(state)
    eif.out <- compute_eif(x, state, kappas)

    estimated <- eif.out$estimated + mean(eif.out$EIF)
    lower.ci <- estimated - 1.96 * sqrt(mean(eif.out$EIF^2) / n)
    upper.ci <- estimated + 1.96 * sqrt(mean(eif.out$EIF^2) / n)

    list(estimated = estimated, EIF = eif.out$EIF, EIF.Y = eif.out$EIF.Y,
         EIF.X = eif.out$EIF.X, EIF.W = eif.out$EIF.W, EIF.C = eif.out$EIF.C,
         lower.ci = lower.ci, upper.ci = upper.ci)
  }

  f.onestep.xz_binaryZ <- function(x, z) {
    state <- make_state_binaryZ(x, z)
    kappas <- compute_kappas_binaryZ(state)
    eif.out <- compute_eif_binaryZ(x, z, state, kappas)

    estimated <- eif.out$estimated + mean(eif.out$EIF)
    lower.ci <- estimated - 1.96 * sqrt(mean(eif.out$EIF^2) / n)
    upper.ci <- estimated + 1.96 * sqrt(mean(eif.out$EIF^2) / n)

    list(estimated = estimated, EIF = eif.out$EIF, EIF.Y = eif.out$EIF.Y,
         EIF.X = eif.out$EIF.X, EIF.W = eif.out$EIF.W, EIF.C = eif.out$EIF.C,
         lower.ci = lower.ci, upper.ci = upper.ci)
  }

  f.tmle.x_nonbinaryZ <- function(x) {
    state <- make_state(x)
    kappas <- compute_kappas(state)
    kappa2.denom <- get_kappa2_denom(x, state, kappas)
    eif.out <- compute_eif(x, state, kappas, kappa2.denom)

    iter <- 0
    converged <- max(abs(mean(eif.out$EIF.Y)), abs(mean(eif.out$EIF.X)),
                     abs(mean(eif.out$EIF.W))) <= cvg.criteria

    while (!converged && iter < n.iter) {
      iter <- iter + 1

      #############################
      # update pi: p(X=x|Z,W,C)
      #############################

      iter.X <- 0
      while ((iter.X == 0 || abs(mean(eif.out$EIF.X)) > cvg.criteria) &&
             iter.X < n.iter) {
        covariate.X <- (state$mu.obs - kappas$psi.c) / kappa2.denom
        ps_model <- glm((X == x) ~ offset(qlogis(state$pi.obs)) + covariate.X - 1,
                        family = binomial(), start = 0, weights = p.z / p.z.wc)
        eps.X <- as.numeric(coef(ps_model))

        state$pi.obs <- plogis(qlogis(state$pi.obs) + eps.X * covariate.X)

        for (w in w.levels) {
          wi <- as.character(w)
          for (j in seq_len(nMC)) {
            covariate.Xj <- (state$mu.sim[[wi]][, j] - kappas$psi.c) / kappa2.denom
            state$pi.sim[[wi]][, j] <- plogis(qlogis(state$pi.sim[[wi]][, j]) +
                                                 eps.X * covariate.Xj)
          }
        }

        for (j in seq_len(nMC)) {
          covariate.Xj <- (state$mu.obsW.sim[, j] - kappas$psi.c) / kappa2.denom
          state$pi.obsW.sim[, j] <- plogis(qlogis(state$pi.obsW.sim[, j]) +
                                             eps.X * covariate.Xj)
        }

        kappas <- compute_kappas(state)
        kappa2.denom <- get_kappa2_denom(x, state, kappas)
        eif.out <- compute_eif(x, state, kappas, kappa2.denom)
        iter.X <- iter.X + 1
      }

      #############################
      # update mu: E(Y|X=x,Z,W,C)
      #############################

      weight.Y <- (X == x) * p.z / (kappa2.denom * p.z.wc)
      or_model <- glm(Y ~ offset(state$mu.obs) + 1, weights = weight.Y)
      eps.Y <- as.numeric(coef(or_model))

      state$mu.obs <- state$mu.obs + eps.Y
      state$mu.obsW.sim <- state$mu.obsW.sim + eps.Y
      for (w in w.levels) {
        wi <- as.character(w)
        state$mu.sim[[wi]] <- state$mu.sim[[wi]] + eps.Y
      }

      kappas <- compute_kappas(state)
      kappa2.denom <- get_kappa2_denom(x, state, kappas)
      eif.out <- compute_eif(x, state, kappas, kappa2.denom)

      #############################
      # update psi(C)
      #############################

      pseudo.W <- kappas$kappa1 / kappas$kappa2
      H.W <- kappas$kappa2 / kappa2.denom
      psi_model <- glm(pseudo.W ~ offset(kappas$psi.c) + 1, weights = H.W)
      eps.psi <- as.numeric(coef(psi_model))
      kappas$psi.c <- kappas$psi.c + eps.psi

      eif.out <- compute_eif(x, state, kappas, kappa2.denom)
      converged <- max(abs(mean(eif.out$EIF.Y)), abs(mean(eif.out$EIF.X)),
                       abs(mean(eif.out$EIF.W))) <= cvg.criteria

      if (verbose) {
        cat("TMLE iter:", iter,
            "EIF.Y:", mean(eif.out$EIF.Y),
            "EIF.X:", mean(eif.out$EIF.X),
            "EIF.W:", mean(eif.out$EIF.W), "\n")
      }
    }

    estimated <- eif.out$estimated
    lower.ci <- estimated - 1.96 * sqrt(mean(eif.out$EIF^2) / n)
    upper.ci <- estimated + 1.96 * sqrt(mean(eif.out$EIF^2) / n)

    list(estimated = estimated, EIF = eif.out$EIF, EIF.Y = eif.out$EIF.Y,
         EIF.X = eif.out$EIF.X, EIF.W = eif.out$EIF.W, EIF.C = eif.out$EIF.C,
         lower.ci = lower.ci, upper.ci = upper.ci, iter = iter)
  }

  f.tmle.xz_binaryZ <- function(x, z) {
    state <- make_state_binaryZ(x, z)
    kappas <- compute_kappas_binaryZ(state)
    kappa2.denom <- get_kappa2_aipw_binaryZ(x, z, state)
    eif.out <- compute_eif_binaryZ(x, z, state, kappas, kappa2.denom)

    iter <- 0
    converged <- max(abs(mean(eif.out$EIF.Y)), abs(mean(eif.out$EIF.X)),
                     abs(mean(eif.out$EIF.W))) <= cvg.criteria

    while (!converged && iter < n.iter) {
      iter <- iter + 1

      #############################
      # update pi: p(X=x|Z=z,W,C)
      #############################

      iter.X <- 0
      while ((iter.X == 0 || abs(mean(eif.out$EIF.X)) > cvg.criteria) &&
             iter.X < n.iter) {
        covariate.X <- (state$mu.obs - kappas$psi.c) / kappa2.denom
        ps_model <- glm((X == x) ~ offset(qlogis(state$pi.obs)) + covariate.X - 1,
                        family = binomial(), start = 0,
                        weights = (Z == z) / f.z.binary.obs(z))
        eps.X <- as.numeric(coef(ps_model))

        state$pi.obs <- plogis(qlogis(state$pi.obs) + eps.X * covariate.X)
        state$pi.obsW <- plogis(qlogis(state$pi.obsW) + eps.X * covariate.X)

        for (w in w.levels) {
          wi <- as.character(w)
          covariate.Xw <- (state$mu.w[[wi]] - kappas$psi.c) / kappa2.denom
          state$pi.w[[wi]] <- plogis(qlogis(state$pi.w[[wi]]) + eps.X * covariate.Xw)
        }

        kappas <- compute_kappas_binaryZ(state)
        kappa2.denom <- get_kappa2_aipw_binaryZ(x, z, state)
        eif.out <- compute_eif_binaryZ(x, z, state, kappas, kappa2.denom)
        iter.X <- iter.X + 1
      }

      #############################
      # update mu: E(Y|X=x,Z=z,W,C)
      #############################

      weight.Y <- (X == x) * (Z == z) / (kappa2.denom * f.z.binary.obs(z))
      if (binaryY) {
        state$mu.obs <- pmin(pmax(state$mu.obs, 0.001), 0.999)
        or_model <- glm(Y ~ offset(qlogis(state$mu.obs)) + 1,
                        family = binomial(), start = 0, weights = weight.Y)
        eps.Y <- as.numeric(coef(or_model))
        update_mu <- function(mu) plogis(qlogis(pmin(pmax(mu, 0.001), 0.999)) + eps.Y)
      } else {
        or_model <- glm(Y ~ offset(state$mu.obs) + 1, weights = weight.Y)
        eps.Y <- as.numeric(coef(or_model))
        update_mu <- function(mu) mu + eps.Y
      }

      state$mu.obs <- update_mu(state$mu.obs)
      state$mu.obsW <- update_mu(state$mu.obsW)
      for (w in w.levels) {
        wi <- as.character(w)
        state$mu.w[[wi]] <- update_mu(state$mu.w[[wi]])
      }

      kappas <- compute_kappas_binaryZ(state)
      kappa2.denom <- get_kappa2_aipw_binaryZ(x, z, state)
      eif.out <- compute_eif_binaryZ(x, z, state, kappas, kappa2.denom)

      #############################
      # update psi(C)
      #############################

      pseudo.W <- kappas$kappa1 / kappas$kappa2
      H.W <- kappas$kappa2 / kappa2.denom
      psi_model <- glm(pseudo.W ~ offset(kappas$psi.c) + 1, weights = H.W)
      eps.psi <- as.numeric(coef(psi_model))
      kappas$psi.c <- kappas$psi.c + eps.psi

      eif.out <- compute_eif_binaryZ(x, z, state, kappas, kappa2.denom)
      converged <- max(abs(mean(eif.out$EIF.Y)), abs(mean(eif.out$EIF.X)),
                       abs(mean(eif.out$EIF.W))) <= cvg.criteria

      if (verbose) {
        cat("TMLE binary-Z z =", z, "iter:", iter,
            "EIF.Y:", mean(eif.out$EIF.Y),
            "EIF.X:", mean(eif.out$EIF.X),
            "EIF.W:", mean(eif.out$EIF.W), "\n")
      }
    }

    estimated <- eif.out$estimated
    lower.ci <- estimated - 1.96 * sqrt(mean(eif.out$EIF^2) / n)
    upper.ci <- estimated + 1.96 * sqrt(mean(eif.out$EIF^2) / n)

    list(estimated = estimated, EIF = eif.out$EIF, EIF.Y = eif.out$EIF.Y,
         EIF.X = eif.out$EIF.X, EIF.W = eif.out$EIF.W, EIF.C = eif.out$EIF.C,
         lower.ci = lower.ci, upper.ci = upper.ci, iter = iter)
  }

  combine_binary_z <- function(out.z1, out.z0) {
    alpha.opt <- mean(out.z0$EIF * (out.z0$EIF - out.z1$EIF)) /
      mean((out.z0$EIF - out.z1$EIF)^2)
    estimated <- out.z1$estimated * alpha.opt + out.z0$estimated * (1 - alpha.opt)
    EIF <- out.z1$EIF * alpha.opt + out.z0$EIF * (1 - alpha.opt)
    lower.ci <- estimated - 1.96 * sqrt(mean(EIF^2) / n)
    upper.ci <- estimated + 1.96 * sqrt(mean(EIF^2) / n)
    list(estimated = estimated, EIF = EIF, lower.ci = lower.ci, upper.ci = upper.ci,
         alpha.opt = alpha.opt)
  }

  if (binaryZ) {
    if (is.null(z)) {
      out.z1 <- f.onestep.xz_binaryZ(x, 1)
      out.z0 <- f.onestep.xz_binaryZ(x, 0)
      out.all.z <- combine_binary_z(out.z1, out.z0)
      onestep.out <- list(out.z1 = out.z1, out.z0 = out.z0, out.all.z = out.all.z)

      out.z1 <- f.tmle.xz_binaryZ(x, 1)
      out.z0 <- f.tmle.xz_binaryZ(x, 0)
      out.all.z <- combine_binary_z(out.z1, out.z0)
      tmle.out <- list(out.z1 = out.z1, out.z0 = out.z0, out.all.z = out.all.z)
    } else {
      out.z <- f.onestep.xz_binaryZ(x, z)
      onestep.out <- setNames(list(out.z), paste0("out.z", z))

      out.z <- f.tmle.xz_binaryZ(x, z)
      tmle.out <- setNames(list(out.z), paste0("out.z", z))
    }
  } else {
    out.all.z <- f.onestep.x_nonbinaryZ(x)
    onestep.out <- list(out.all.z = out.all.z)

    out.all.z <- f.tmle.x_nonbinaryZ(x)
    tmle.out <- list(out.all.z = out.all.z)
  }

  list(Onestep = onestep.out, TMLE = tmle.out)
}
