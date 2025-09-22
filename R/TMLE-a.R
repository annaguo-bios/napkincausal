#' Estimate average counterfactual outcome E(Y(x)).
#'
#' Function for estimating the average counterfactual outcome E(Y(a)) for napkin graph.
#' @param x Treatment level.
#'          The average counterfactual outcome will be computed at the specified treatment level x.
#'          The function returns `E(Y(x))`, which is the
#'          average counterfactual outcome under the treatment level specified by `x`.
#' @param data A Dataframe contains all the variables listed in vertices parameter
#' @param treatment A character string indicating the treatment variable
#' @param Z.variables A character vector indicating the variables, the level of which are fixed in the estimation of the counterfactual outcome.
#' @param z A numeric value indicating the level of Z.variables that we want to intervene on when estimating the average counterfactual outcome.
#' The default is NULL. If Z is discrete, it will return the average counterfactual outcome under all observed level of Z and the combined estimate over all levels of Z.
#' If Z is continuous, it will return the combined estimate over all levels of Z.
#' @param W.variables A character vector indicating the parents of Z.variables.
#' @param outcome A character string indicating the outcome variable
#' @param z.density A argument that takes the user specified function for estimating the density of Z.
#' @param z.method A character string indicating the method used to estimate the conditional density of p(Z|W). There are two options: 'dnorm' and 'np'. The default is 'dnorm'.
#' If Z is binary, regression based method will be adopted. The formula of the regression can be specified by `formula.Z`.
#' If Z is continuous, the `z.method` method will be adopted. The `dnorm` method assumes that the conditional density of Z given W is normal.
#' The `np` method estimates the conditional density of Z given W via non-parametric method via the \link[np]{np} package.
#' @details In assuming normal distribution. The `dnorm` method estimates the mean of the normal distribution by fitting a linear regression model with only linear terms and no interaction terms and the variance of the normal distribution sample variance of the error term resulted from the linear regression model.
#' In assuming logistic regression, the `dnorm` method further assume the logistic regression contains only linear terms and no interaction terms.
#' @param superlearner.Y A logical indicator determines whether SuperLearner via the \link[SuperLearner]{SuperLearner} function is adopted for estimating the outcome regression.
#' @param superlearner.X A logical indicator determines whether SuperLearner via the \link[SuperLearner]{SuperLearner} function is adopted for estimating the propensity score.
#' @param superlearner.Z A logical indicator determines whether SuperLearner via the \link[SuperLearner]{SuperLearner} function is adopted for estimating the conditional density p(Z|W) when Z is binary.
#' @param crossfit A logical indicator determines whether crossfitting is adopted for SuperLearner. If crossfit is set to TRUE, the data is split into K folds as specified by the `K` parameter.
#' @param K An integer specifying the number of folds for crossfitting.
#' @param lib.Y A character vector specifying the library of algorithms to be used in the SuperLearner for outcome regression.
#' @param lib.X A character vector specifying the library of algorithms to be used in the SuperLearner for propensity score estimation.
#' @param lib.Z A character vector specifying the library of algorithms to be used in the SuperLearner for estimating the conditional density of Z given W.
#' @param formula.Y Regression formula for the outcome regression of Y on it's Markov pillow. The default is 'Y ~ .'.
#' @param formula.X Regression formula for the propensity score regression of A on it's Markov pillow. The default is 'X ~ .'.
#' @param formula.Z Regression formula for the conditional density of Z given W. The default is 'Z ~ .'.
#' @param linkY_binary The link function used for outcome regression of Y on it's Markov pillow when Y is binary and superlearner is not sued. The default is the 'logit' link.
#' @param link.X The link function used for propensity score regression of X on it's Markov pillow if superlearner is not used. The default is the 'logit' link.
#' @param link.Z The link function used for the regression of Z on it's Markov pillow if superlearner is not used. The default is the 'logit' link.
#' @param cvg.criteria A numerical value representing the convergence criteria for the iterative update of the nusiances in TMLE.
#' If the absolute of the mean of efficient influence function is less than cvg.criteria, the iterative update stops. The default value is 0.01.
#' @param n.iter The maximum number of iterations for the iterative update of the nuisances in TMLE. The default value is 500.
#' @param truncate_lower.X The lower bound for truncation of the propensity score. The default is 0, which means no truncation.
#' @param truncate_upper.X The upper bound for truncation of the propensity score. The default is 1, which means no truncation.
#' @param truncate_lower.Z The lower bound for truncation of the conditional density of Z given W. The default is 0, which means no truncation.
#' @param truncate_upper.Z The upper bound for truncation of the conditional density of Z given W. The default is 1, which means no truncation.
#' @param minZ The lower bound used for performing integration of Z when Z is continuous. The default is -Inf.
#' @param maxZ The upper bound used for performing integration of Z when Z is continuous. The default is Inf.
#' @param verbose A logical indicator determines whether the function prints out detailed progress of the estimation. The default is TRUE.
#' @param fast A logical indicator determines whether the integration involved in the estimation is performed via `integrate()` function or via Monte Carlo integration. The former is lower while the later is faster. The default is TRUE.
#' @param nMC The number of Monte Carlo samples used for integration when `fast` is set to TRUE The default is 5000.
#' @param boundedsubmodelY An indicator for whether the bounded submodel is used for targeting the outcome regression when Z is discrete. The default is FALSE.
#' @return Function outputs a list containing TMLE results and Onestep results. Access TMLE result via `$TMLE` and access onestep result via `$onestep`.
#' \describe{
#'       \item{\code{estimated}}{The estimated parameter of interest: \eqn{E(Y^x)}}
#'       \item{\code{lower.ci}}{Lower bound of the 95% confidence interval for \code{estimated}}
#'       \item{\code{upper.ci}}{Upper bound of the 95 pct confidence interval for \code{estimated}}
#'       \item{\code{EIF}}{The estimated efficient influence function evaluated at the observed data}
#'       \item{\code{EIF.Y}}{The part of the EIF corresponding to the outcome regression}
#'       \item{\code{EIF.X}}{The part of EIF corresponding to the propensity score regression}
#'       \item{\code{EDstar.record}}{The sample mean of the estimated efficient influence function mapped into the tangent space of Y, Z, W over iterations. If convergence is achieved, this should be close to 0 for TMLE estimators.}
#'       \item{\code{iter}}{Number of iterations where convergence is achieved for the iterative update of the mediator density and propensity score.}}
#' @examples
#' # E(Y(1)) estimation.
#' res <- napkin.a(x=1, z = NULL, data=data_Zbinary_Ycontinuous,
#' treatment="X", Z.variables="Z", W.variables="W", outcome="Y",
#' formula.Y="Y ~ X*Z*W", formula.X="X ~ Z*W", formula.Z="Z~.",
#' link.X="identity", link.Z="logit")
#' @importFrom dplyr %>% mutate select
#' @importFrom MASS mvrnorm
#' @importFrom SuperLearner CV.SuperLearner SuperLearner
#' @importFrom mvtnorm dmvnorm
#' @importFrom densratio densratio
#' @importFrom utils combn
#' @importFrom stats rnorm runif rbinom dnorm dbinom binomial gaussian predict glm as.formula qlogis plogis lm coef cov sd density approx integrate setNames
#' @importFrom np npcdensbw npcdens
#' @importFrom RVCompare sampleFromDensity
#' @export

napkin.a <- function(x, z = NULL,data,treatment, Z.variables, W.variables, outcome,
                     z.density=NULL, z.method="dnorm", superlearner.Y=F, superlearner.X=F,superlearner.Z=F,
                     crossfit=F,K=5,
                     lib.Y = c("SL.glm","SL.earth","SL.ranger","SL.mean"),
                     lib.X = c("SL.glm","SL.earth","SL.ranger","SL.mean"),
                     lib.Z = c("SL.glm","SL.earth","SL.ranger","SL.mean"),
                     n.iter=500, cvg.criteria=0.01,
                     formula.Y="Y ~ .", formula.X="X ~ .", formula.Z="Z~.",
                     linkY_binary="logit", link.X="logit", link.Z="logit",
                     truncate_lower.X=0, truncate_upper.X=1,
                     truncate_lower.Z=0, truncate_upper.Z=1,
                     minZ=-Inf,maxZ=Inf,
                     verbose=T,fast=T,nMC=5000,boundedsubmodelY=F){



  n <- nrow(data)

  # Variables
  X <- data[,treatment]
  W <- data[,W.variables, drop=F]
  Z <- data[,Z.variables]
  Y <- data[,outcome]

  minY <- min(Y)
  maxY <- max(Y)

  #formulas
  formula.Z <- as.formula(formula.Z)
  formula.X <- as.formula(formula.X)
  formula.Y <- as.formula(formula.Y)

  # new data sets
  dat_mpY = data.frame(X,Z,W)
  dat_mpX = data.frame(Z,W)
  dat_mpZ = data.frame(W)
  dat_ZmpZ = data.frame(Z,W)

  binaryY <- all(Y %in% c(0,1))
  binaryZ <- all(as.vector(Z) %in% c(0,1))

  z.density.provided <- is.null(z.density)

  ################################################
  ############### OUTCOME REGRESSION #############
  ################################################


  if (crossfit==T){ #### cross fitting + super learner #####

    fit.family <- if(binaryY){binomial(linkY_binary)}else{gaussian()} # family for super learner depending on whether Y is binary or continuous

    or_fit <- CV.SuperLearner(Y=Y, X=dat_mpY, family = fit.family, V = K, SL.library = lib.Y, control = list(saveFitLibrary=T), saveAll = T)

    f.mu.xzw<- function(x,z,w){

      if(is.vector(w)){w <- as.data.frame(t(w))}

      dat_mpY.xz <- data.frame(Z = z, X = x, setNames(as.data.frame(w), W.variables))

      mu.Y.xz <- unlist(lapply(1:K, function(x) predict(or_fit$AllSL[[x]], newdata=dat_mpY.xz[or_fit$folds[[x]],])[[1]] %>% as.vector()))[order(unlist(lapply(1:K, function(x) or_fit$folds[[x]])))]


      return(mu.Y.xz)
    } # end of function f.mu.xz


  } else if (superlearner.Y==T){ #### super learner #####

    fit.family <- if(binaryY){binomial(linkY_binary)}else{gaussian()} # family for super learner depending on whether Y is binary or continuous

    or_fit <- SuperLearner(Y=Y, X=dat_mpY, family = fit.family, SL.library = lib.Y)

    f.mu.xzw<- function(x,z,w){

      if(is.vector(w)){w <- as.data.frame(t(w))}

      dat_mpY.xz <- data.frame(Z = z, X = x, setNames(as.data.frame(w), W.variables))

      mu.Y.xz <- predict(or_fit, newdata=dat_mpY.xz)[[1]] %>% as.vector()

      return(mu.Y.xz)
    } # end of function f.mu.xz


  } else { #### simple linear regression with user input regression formula: default="Y ~ ." ####

    fit.family <- if(binaryY){binomial(linkY_binary)}else{gaussian()} # family for super learner depending on whether Y is binary or continuous

    or_fit <- glm(formula.Y, data=dat_mpY, family = fit.family)

    f.mu.xzw <- function(x,z,w){

      if(is.vector(w)){w <- as.data.frame(t(w))}

      dat_mpY.xz <- data.frame(Z = z, X = x, setNames(as.data.frame(w), W.variables))

      mu.Y.xz <- predict(or_fit, newdata=dat_mpY.xz, type="response")


      return(mu.Y.xz)
    }



  } # end of if-else for outcome regression fitting method


  if(boundedsubmodelY){

    if (crossfit==T){ #### cross fitting + super learner #####

      fit.family <- if(binaryY){binomial(linkY_binary)}else{gaussian()} # family for super learner depending on whether Y is binary or continuous

      or_fit <- CV.SuperLearner(Y=(Y-min(Y))/(max(Y)-min(Y)), X=dat_mpY, family = fit.family, V = K, SL.library = lib.Y, control = list(saveFitLibrary=T), saveAll = T)

      f.mu.xzw.bounded<- function(x,z,w){


        if(is.vector(w)){w <- as.data.frame(t(w))}

        dat_mpY.xz <- data.frame(Z = z, X = x, setNames(as.data.frame(w), W.variables))

        mu.Y.xz <- unlist(lapply(1:K, function(x) predict(or_fit$AllSL[[x]], newdata=dat_mpY.xz[or_fit$folds[[x]],])[[1]] %>% as.vector()))[order(unlist(lapply(1:K, function(x) or_fit$folds[[x]])))]


        return(mu.Y.xz)
      } # end of function f.mu.xz


    } else if (superlearner.Y==T){ #### super learner #####

      fit.family <- if(binaryY){binomial(linkY_binary)}else{gaussian()} # family for super learner depending on whether Y is binary or continuous

      or_fit <- SuperLearner(Y=(Y-min(Y))/(max(Y)-min(Y)), X=dat_mpY, family = fit.family, SL.library = lib.Y)

      f.mu.xzw.bounded<- function(x,z,w){


        if(is.vector(w)){w <- as.data.frame(t(w))}

        dat_mpY.xz <- data.frame(Z = z, X = x, setNames(as.data.frame(w), W.variables))

        mu.Y.xz <- predict(or_fit, newdata=dat_mpY.xz)[[1]] %>% as.vector()

        return(mu.Y.xz)
      } # end of function f.mu.xz


    } else { #### simple linear regression with user input regression formula: default="Y ~ ." ####

      fit.family <- if(binaryY){binomial(linkY_binary)}else{gaussian()} # family for super learner depending on whether Y is binary or continuous

      dat_mpY.bounded <- dat_mpY %>% mutate(Y = (Y-min(Y))/(max(Y)-min(Y)))

      or_fit <- glm(formula.Y, data=dat_mpY.bounded, family = fit.family)

      f.mu.xzw.bounded <- function(x,z,w){

        if(is.vector(w)){w <- as.data.frame(t(w))}

        dat_mpY.xz <- data.frame(Z = z, X = x, setNames(as.data.frame(w), W.variables))

        mu.Y.xz <- predict(or_fit, newdata=dat_mpY.xz, type="response")


        return(mu.Y.xz)
      }



    } # end of if-else for outcome regression fitting method

  }

  cat("outcome regression done.")

  ################################################
  ############### PROPENSITY SCORE ###############
  ################################################


  if (crossfit==T){ #### cross fitting + super learner #####

    ps_fit <- CV.SuperLearner(Y=X, X=dat_mpX, family = binomial(), V = K, SL.library = lib.X, control = list(saveFitLibrary=T),saveAll = T)

    f.x.zw <- function(x, z, w, truncate_lower, truncate_upper){


      if(is.vector(w)){w <- as.data.frame(t(w))}

      dat_mpX.z <- data.frame(Z=z, setNames(as.data.frame(w), W.variables))

      p.x1.z <- unlist(lapply(1:K, function(x) predict(ps_fit$AllSL[[x]], newdata=dat_mpX.z[ps_fit$folds[[x]],])[[1]] %>% as.vector()))[order(unlist(lapply(1:K, function(x) ps_fit$folds[[x]])))]

      # truncation
      p.x1.z[p.x1.z < truncate_lower] <- truncate_lower
      p.x1.z[p.x1.z > truncate_upper] <- truncate_upper

      p.x.z <- x*p.x1.z + (1-x)*(1-p.x1.z)

      return(p.x.z)

    } # end of function f.pi.z


  } else if (superlearner.X==T){ #### super learner #####

    ps_fit <- SuperLearner(Y=X, X=dat_mpX, family = binomial(), SL.library = lib.X)

    # p(X=1|Z=z,w)
    f.x.zw <- function(x, z, w, truncate_lower, truncate_upper){


      if(is.vector(w)){w <- as.data.frame(t(w))}

      dat_mpX.z <- data.frame(Z=z, setNames(as.data.frame(w), W.variables))

      p.x1.z <- predict(ps_fit, newdata=dat_mpX.z, type="response")[[1]] %>% as.vector()

      # truncation
      p.x1.z[p.x1.z < truncate_lower] <- truncate_lower
      p.x1.z[p.x1.z > truncate_upper] <- truncate_upper

      p.x.z <- x*p.x1.z + (1-x)*(1-p.x1.z)

      return(p.x.z)

    } # end of function f.pi.z



  } else { #### simple linear regression with user input regression formula: default="A ~ ." ####

    ps_fit <- glm(formula.X, data=dat_mpX,  family = binomial(link.X))


    f.x.zw <- function(x, z,w, truncate_lower, truncate_upper){


      if(is.vector(w)){w <- as.data.frame(t(w))}

      dat_mpX.z <- data.frame(Z=z, setNames(as.data.frame(w), W.variables))
      p.x1.z <- predict(ps_fit, newdata=dat_mpX.z, type="response")

      # truncation
      p.x1.z[p.x1.z < truncate_lower] <- truncate_lower
      p.x1.z[p.x1.z > truncate_upper] <- truncate_upper

      p.x.z <- x*p.x1.z + (1-x)*(1-p.x1.z)

      return(p.x.z)

    } # end of function f.pi.z

  } # end of if-else for propensity score fitting


  cat("propensity score regression done.")


  ################################################################
  ############### f_Z(Z|W) Conditional density ###############
  ################################################################


  ############ If Z is binary ############
  if (all(as.vector(Z) %in% c(0,1))){

    if (crossfit==T){ #### cross fitting + super learner #####

      z_fit <- CV.SuperLearner(Y=Z, X=dat_mpZ, family = binomial(), V = K, SL.library = lib.Z, control = list(saveFitLibrary=T),saveAll = T)

      f.z <- function(z, truncate_lower=0, truncate_upper=1){

        p.z1 <- z_fit$SL.predict

        # truncation
        p.z1[p.z1 < truncate_lower] <- truncate_lower
        p.z1[p.z1 > truncate_upper] <- truncate_upper

        p.z <- z*p.z1 + (1-z)*(1-p.z1)

        return(p.z)

      } # end of function f.z


    } else if (superlearner.Z==T){ #### super learner #####

      z_fit <- SuperLearner(Y=Z, X=dat_mpZ, family = binomial(), SL.library = lib.Z)

      # p(X=1|Z=z,W)
      f.z <- function(z, truncate_lower=0, truncate_upper=1){

        p.z1 <- predict(z_fit, type = "response")[[1]] %>% as.vector()

        # truncation
        p.z1[p.z1 < truncate_lower] <- truncate_lower
        p.z1[p.z1 > truncate_upper] <- truncate_upper

        p.z <- z*p.z1 + (1-z)*(1-p.z1)

        return(p.z)

      } # end of function f.pi.z



    } else { #### simple linear regression with user input regression formula: default="A ~ ." ####

      z_fit <- glm(formula.Z, data=dat_mpZ,  family = binomial(link.Z))

      f.z <- function(z, truncate_lower=0, truncate_upper=1){

        p.z1 <- predict(z_fit, type="response")

        # truncation
        p.z1[p.z1 < truncate_lower] <- truncate_lower
        p.z1[p.z1 > truncate_upper] <- truncate_upper

        p.z <- z*p.z1 + (1-z)*(1-p.z1)

        return(p.z)

      } # end of function f.pi.z

    } # end of if-else for Z|W regression fitting




  ########## If Z is NOT binary ############
  }else{


    if (is.function(z.method)){

      # make prediction for p(Z_i|W_i)
      p.z <- sapply(1:n, function(i) z.method(Z.variables = Z[i], W.variables = as.vector(W[i,])) )

    }else if (z.method=="np"){

      ## Z|W
      # Methods on density estimation:
      # np: https://cran.r-project.org/web/packages/np/np.pdf

      environment(formula.Z) <- environment()

      bw <- eval(substitute(npcdensbw(formula = FORM, data = DATA),
                            list(FORM = formula.Z, DATA = dat_ZmpZ)))

      # bw <- npcdensbw(formula=formula.Z, data=dat_ZmpZ)
      Z_fit <- npcdens(bws=bw)

      # make prediction for p(Z_i|W_i)
      p.z <- predict(Z_fit)

    }else if (z.method=="dnorm"){

      dnorm.density <- calculate_density_dnorm(Z.variables=Z.variables, W.variables=W.variables, data=data, formula.Z=formula.Z, superlearner.Z=superlearner.Z, crossfit=crossfit, K=K)

      # make prediction for p(Z_i|W_i)
      p.z <- dnorm.density[[1]]

    }else{


      stop("Invalid input. z.method must be a user-specified density function or either 'np' or 'dnorm'")

    }



  } # end of if-else for binary Z


  cat("Estimation for conditional density of Z completed.")



  #############################################
  ############### p(Z) density ###############
  #############################################

  if (z.density.provided){

    den.zi <- density(Z)

    p.zi <- approx(den.zi$x, den.zi$y, xout = Z)$y

  }else{

    p.zi <- sapply(Z, z.density)

  }


  cat("Estimation/Evaluation for p(Z) completed.")

  ##################################################################
  #################### One-step estimator ##########################
  ##################################################################


  ### Binary Z ###

  f.onestep.xz_binaryZ <- function(x,z, truncate_lower.X, truncate_upper.X ,truncate_lower.Z, truncate_upper.Z){ # one-step estimate and its EIF

    # nuisance parameters
    p.z <- f.z(z, truncate_lower.Z, truncate_upper.Z) # z-density
    p.x.z <- f.x.zw(x, z,W, truncate_lower.X, truncate_upper.X) # propensity score
    mu.xz <- f.mu.xzw(x, z,W) # outcome regression

    phi1 <- mean(mu.xz*p.x.z) # numerator of the plugin estimator
    phi2 <- mean(p.x.z) # denominator of the plugin estimator

    # uncetered EIF of phi1
    kappa1_plus_phi1 <- {(Z==z)/p.z}*{(X==x)*Y - mu.xz*p.x.z} + mu.xz*p.x.z

    # uncentered EIF of phi2
    kappa2_plus_phi2 <- {(Z==z)/p.z}*( (X==x) - p.x.z ) + p.x.z

    plugin.est <- phi1/phi2 # plug-in estimate

    # EIF
    EIF.Y <- {(X==x)*(Z==z)}/{mean(kappa2_plus_phi2)*p.z}*(Y-mu.xz) # EIF for Y
    EIF.X <- (Z==z)/{mean(kappa2_plus_phi2)*p.z}*(mu.xz-plugin.est)*( (X==x) - p.x.z ) # EIF for X
    EIF.W <- 1/mean(kappa2_plus_phi2)*{mu.xz*p.x.z - plugin.est*p.x.z}

    EIF <- EIF.Y + EIF.X + EIF.W

    # the one-step estimator
    estimated <- plugin.est + mean(EIF)

    # confidence interval
    lower.ci <- estimated-1.96*sqrt(mean(EIF^2)/n)
    upper.ci <- estimated+1.96*sqrt(mean(EIF^2)/n)

    return(list(estimated=estimated, EIF=EIF, EIF.Y=EIF.Y, EIF.X=EIF.X, EIF.W=EIF.W, lower.ci=lower.ci, upper.ci=upper.ci, kappa1_plus_phi1=kappa1_plus_phi1, kappa2_plus_phi2=kappa2_plus_phi2))


  }

  ### Non-binary Z ###
  f.onestep.x_nonbinaryZ <- function(x, truncate_lower.X, truncate_upper.X ,truncate_lower.Z, truncate_upper.Z, p.z){ # one-step estimate and its EIF

    # nuisance parameters
    p.x.z <- f.x.zw(x, Z,W, truncate_lower.X, truncate_upper.X) # propensity score
    mu.xz <- f.mu.xzw(x, Z,W) # outcome regression

    phi1 <- sapply(1:n, function(i) mean(f.x.zw(x, Z[i],W, truncate_lower.X, truncate_upper.X)*f.mu.xzw(x, Z[i],W) )) # numerator of the plugin estimator
    phi2 <- sapply(1:n, function(i) mean(f.x.zw(x, Z[i],W, truncate_lower.X, truncate_upper.X) )) # denominator of the plugin estimator

    plugin.est <- phi1/phi2 # plug-in estimate

    # EIF
    EIF.Y <- {(X==x)*(p.zi)}/{phi2*p.z}*(Y-mu.xz) # EIF for Y
    EIF.X <- (p.zi)/{phi2*p.z}*(mu.xz-plugin.est)*( (X==x) - p.x.z ) # EIF for X

    # functions for computing EIF.W and plugin psi(Z_j) ----------------- #

    # integrand for calculating the EIF for W
    integrand.EIF.W <- function(z,w){

      phi1.z <- sapply(z, function(zi) mean(f.x.zw(x, zi,W, truncate_lower.X, truncate_upper.X)*f.mu.xzw(x, zi,W) ))
      phi2.z <- sapply(z, function(zi) mean(f.x.zw(x, zi,W, truncate_lower.X, truncate_upper.X)))

      plugin.est.z <- phi1.z/phi2.z

      int.w <- z.density(z)*(f.x.zw(x, z, w, truncate_lower.X, truncate_upper.X)/phi2.z*(f.mu.xzw(x, z,w,c) - plugin.est.z) )

      return(int.w)

    }

    # used to calculate the EIF for W and plugin estimator FASTLY
    fast.EIF.W <- function(Zsim){

      mu.matrix.sim <- matrix(NA, n, nMC)

      for (i in 1:nMC) { # the columns are the predictions for p(Y_i|X=x,Z_j,W_i) for i from 1 to n, Z_j in Zsim
        mu.matrix.sim[, i] <- f.mu.xzw(x, Zsim[i],W)
      }

      p.x.matrix.sim <- matrix(NA, n, nMC)

      for (i in 1:nMC) { # the columns are the predictions for p(X=x|Z_j,W_i) for i from 1 to n, Z_j in Zsim
        p.x.matrix.sim[, i] <- f.x.zw(x, Zsim[i],W, truncate_lower.X, truncate_upper.X)
      }

      phi1.sim <- colMeans(mu.matrix.sim*p.x.matrix.sim) # numerator of the plugin estimator
      phi2.sim <- colMeans(p.x.matrix.sim) # denominator of the plugin estimator

      plugin.est.sim <- phi1.sim/phi2.sim # plug-in estimate

      # EIF.W <- matrix(phi1.Zsim.long, nrow=n, byrow = T)/phi2.Zsim - matrix(phi2.Zsim.long, nrow=n, byrow = T)*(plugin.est.z/phi2.Zsim)
      EIF.W <- rowMeans(sweep(p.x.matrix.sim*sweep(mu.matrix.sim, 2, plugin.est.sim, "-"),2, phi2.sim,"/"))

      return(list(EIF.W=EIF.W, plugin.est.z=plugin.est.sim))

    }


    # integrand for calculating the plugin estimator
    integrand.plugin <- function(z){

      phi1.z <- sapply(z, function(zi) mean(f.x.zw(x, zi,W, truncate_lower.X, truncate_upper.X)*f.mu.xzw(x, zi,W) ))
      phi2.z <- sapply(z, function(zi) mean(f.x.zw(x, zi,W, truncate_lower.X, truncate_upper.X)))

      plugin.est.z <- phi1.z/phi2.z

      int.w1 <- z.density(z)*plugin.est.z

      return(int.w1)

    }

    if (!z.density.provided){ # if the z.density function is given

      if (fast){

        if(verbose){cat("Fast integration used. Integration computed with Monte Carlo integration method.")}

        Zsim <- sampleFromDensity(z.density, nMC, range(Z))

        fast.object <- fast.EIF.W(Zsim)

        EIF.W <- fast.object$EIF.W

        plugin.est <- mean(fast.object$plugin.est.z)

      }else{

        if(verbose){cat("Integration computed with integrate() function. This step may take a while.\n If you want to speed up the process, set fast=TRUE.")}

        EIF.W <- unlist(lapply(1:n, function(i) integrate(integrand.EIF.W, lower = minZ, upper = maxZ, w=W[i,])$value))

        plugin.est <- integrate(integrand.plugin, lower = minZ, upper = maxZ)$value

      }

      EIF <- EIF.Y + EIF.X + EIF.W

      if(verbose){cat("Density function of Z is provided. \n Influence function and plug-in estimator computed using the given density function of Z.")}

    }else{ # if the z.density function is NOT given

      EIF.W <- Reduce(`+`, lapply(1:n, function(i) f.x.zw(x, Z[i],W, truncate_lower.X, truncate_upper.X)/phi2[i]*(f.mu.xzw(x, Z[i],W) - plugin.est[i])    ))/n

      EIF <- EIF.Y + EIF.X + EIF.W + plugin.est - mean(plugin.est)

      if(verbose){cat("Density function of Z is not provided. Influence function and plug-in estimator computed via density estimation.")}

    }


    estimated <- mean(EIF.Y+EIF.X+EIF.W) + mean(plugin.est)

    # confidence interval
    lower.ci <- estimated-1.96*sqrt(mean(EIF^2)/n)
    upper.ci <- estimated+1.96*sqrt(mean(EIF^2)/n)

    return(list(estimated=estimated, EIF=EIF, EIF.Y=EIF.Y, EIF.X=EIF.X, EIF.W=EIF.W, lower.ci=lower.ci, upper.ci=upper.ci))


  }

  # if Z is binary, return three one-step estimators:
  # 1. at Z=1
  # 2. at Z=0
  # 3. average of the two

  if (binaryZ ){

    if(verbose){print("Z is binary. Computing one-step estimators at Z=1, Z=0, and the average of the two.")}

    out.z1=f.onestep.xz_binaryZ(x,z=1, truncate_lower.X, truncate_upper.X ,truncate_lower.Z, truncate_upper.Z)
    out.z0=f.onestep.xz_binaryZ(x,z=0, truncate_lower.X, truncate_upper.X ,truncate_lower.Z, truncate_upper.Z)

    ## using all levels of Z

    # point estimate
    # estimated= {mean(out.z1$kappa1_plus_phi1)*mean(Z==1) + mean(out.z0$kappa1_plus_phi1)*mean(Z==0)}/{mean(out.z1$kappa2_plus_phi2)*mean(Z==1) + mean(out.z0$kappa2_plus_phi2)*mean(Z==0)}
    # estimated = mean(out.z1$EIF*mean(Z==1)+out.z0$EIF*mean(Z==0)+Z*out.z1$estimated + (1-Z)*out.z0$estimated)
    estimated = out.z1$estimated*mean(Z==1) + out.z0$estimated*mean(Z==0)

    # EIF
    EIF= out.z1$estimated*{(Z==1)-mean(Z==1)} + out.z0$estimated*{(Z==0)-mean(Z==0)} + out.z1$EIF*mean(Z==1) + out.z0$EIF*mean(Z==0)


    # confidence interval
    lower.ci <- estimated-1.96*sqrt(mean(EIF^2)/n)
    upper.ci <- estimated+1.96*sqrt(mean(EIF^2)/n)

    out.all.z = list(estimated=estimated, EIF=EIF, lower.ci=lower.ci, upper.ci=upper.ci)

    onestep.out <- list(out.z1=out.z1, out.z0=out.z0, out.all.z=out.all.z)

  }else{

    if(verbose){cat(paste0('Z is continuous. Computing one-step estimator at the ',ifelse(!z.density.provided,'given','estimated'),' density function of Z.'))}

    out.all.z <- f.onestep.x_nonbinaryZ(x, truncate_lower.X, truncate_upper.X ,truncate_lower.Z, truncate_upper.Z, p.z)

    onestep.out <- list(out.all.z=out.all.z)

  }



  ##################################################################
  #################### Estimating Equation Estimator ###############
  ##################################################################


    f.equation.xz_binaryZ <- function(x,z, truncate_lower.X, truncate_upper.X ,truncate_lower.Z, truncate_upper.Z){ # one-step estimate and its EIF

      # nuisance parameters
      p.z <- f.z(z, truncate_lower.Z, truncate_upper.Z) # z-density
      p.x.z <- f.x.zw(x, z,W,truncate_lower.X, truncate_upper.X) # propensity score
      mu.xz <- f.mu.xzw(x, z,W) # outcome regression

      plugin.est <- mean(mu.xz*p.x.z) # plug-in estimate

      # uncetered EIF of phi1
      kappa1_plus_phi1 <- {(Z==z)/p.z}*{(X==x)*Y - mu.xz*p.x.z} + mu.xz*p.x.z

      # uncentered EIF of phi2
      kappa2_plus_phi2 <- {(Z==z)/p.z}*( (X==x) - p.x.z ) + p.x.z

      # point estimate
      estimated <- mean(kappa1_plus_phi1)/mean(kappa2_plus_phi2)

      # EIF
      EIF.Y <- {(X==x)*(Z==z)}/{mean(kappa2_plus_phi2)*p.z}*(Y-mu.xz) # EIF for Y
      EIF.X <- (Z==z)/{mean(kappa2_plus_phi2)*p.z}*(mu.xz-plugin.est)*( (X==x) - p.x.z ) # EIF for X
      EIF.W <- 1/mean(kappa2_plus_phi2)*{mu.xz*p.x.z - plugin.est*p.x.z}

      EIF <- EIF.Y + EIF.X + EIF.W

      # confidence interval
      lower.ci <- estimated-1.96*sqrt(mean(EIF^2)/n)
      upper.ci <- estimated+1.96*sqrt(mean(EIF^2)/n)

      return(list(estimated=estimated, EIF=EIF, EIF.Y=EIF.Y, EIF.X=EIF.X, EIF.W=EIF.W, lower.ci=lower.ci, upper.ci=upper.ci, kappa1_plus_phi1=kappa1_plus_phi1, kappa2_plus_phi2=kappa2_plus_phi2))


    }


  f.equation.xz_nonbinaryZ <- function(x, truncate_lower.X, truncate_upper.X ,truncate_lower.Z, truncate_upper.Z, p.z){ # one-step estimate and its EIF

    # integrand for calculating the plugin estimator
    integrand.plugin <- function(z){

      phi1.z <- sapply(z, function(zi) mean(f.x.zw(x, zi,W, truncate_lower.X, truncate_upper.X)*f.mu.xzw(x, zi,W) ))
      phi2.z <- sapply(z, function(zi) mean(f.x.zw(x, zi,W, truncate_lower.X, truncate_upper.X)))

      plugin.est.z <- phi1.z/phi2.z

      int.w1 <- z.density(z)*plugin.est.z

      return(int.w1)

    }

    integrated.plugin.est <- integrate(integrand.plugin, lower = minZ, upper = maxZ)$value # integrated plugin estimate


    # nuisance parameters
    p.x.z <- f.x.zw(x, Z,W, truncate_lower.X, truncate_upper.X) # propensity score
    mu.xz <- f.mu.xzw(x, Z,W) # outcome regression

    phi1 <- sapply(1:n, function(i) mean(f.x.zw(x, Z[i],W, truncate_lower.X, truncate_upper.X)*f.mu.xzw(x, Z[i],W) )) # numerator of the plugin estimator
    phi2 <- sapply(1:n, function(i) mean(f.x.zw(x, Z[i],W, truncate_lower.X, truncate_upper.X) )) # denominator of the plugin estimator

    plugin.est <- phi1/phi2

    # estimating equation estimator
    phi1.aipw <- mean(p.zi/{phi2*p.z}*( (X==x)*Y - mu.xz*p.x.z )) + integrated.plugin.est
    phi2.aipw <- mean(p.zi/{phi2*p.z}*( (X==x) - p.x.z )) + 1

    estimated <- phi1.aipw/phi2.aipw


    ## calculate EIF for inference
    # EIF
    EIF.Y <- {(X==x)*(p.zi)}/{phi2*p.z}*(Y-mu.xz) # EIF for Y
    EIF.X <- (p.zi)/{phi2*p.z}*(mu.xz-plugin.est)*( (X==x) - p.x.z ) # EIF for X

    # functions for computing EIF.W and plugin psi(Z_j) ----------------- #

    # integrand for calculating the EIF for W
    integrand.EIF.W <- function(z,w){

      phi1.z <- sapply(z, function(zi) mean(f.x.zw(x, zi,W, truncate_lower.X, truncate_upper.X)*f.mu.xzw(x, zi,W) ))
      phi2.z <- sapply(z, function(zi) mean(f.x.zw(x, zi,W, truncate_lower.X, truncate_upper.X)))

      plugin.est.z <- phi1.z/phi2.z

      int.w <- z.density(z)*(f.x.zw(x, z, w, truncate_lower.X, truncate_upper.X)/phi2.z*(f.mu.xzw(x, z,w,c) - plugin.est.z) )

      return(int.w)

    }

    # used to calculate the EIF for W and plugin estimator FASTLY
    fast.EIF.W <- function(Zsim){

      mu.matrix.sim <- matrix(NA, n, nMC)

      for (i in 1:nMC) { # the columns are the predictions for p(Y_i|X=x,Z_j,W_i) for i from 1 to n, Z_j in Zsim
        mu.matrix.sim[, i] <- f.mu.xzw(x, Zsim[i],W)
      }

      p.x.matrix.sim <- matrix(NA, n, nMC)

      for (i in 1:nMC) { # the columns are the predictions for p(X=x|Z_j,W_i) for i from 1 to n, Z_j in Zsim
        p.x.matrix.sim[, i] <- f.x.zw(x, Zsim[i],W, truncate_lower.X, truncate_upper.X)
      }

      phi1.sim <- colMeans(mu.matrix.sim*p.x.matrix.sim) # numerator of the plugin estimator
      phi2.sim <- colMeans(p.x.matrix.sim) # denominator of the plugin estimator

      plugin.est.sim <- phi1.sim/phi2.sim # plug-in estimate

      # EIF.W <- matrix(phi1.Zsim.long, nrow=n, byrow = T)/phi2.Zsim - matrix(phi2.Zsim.long, nrow=n, byrow = T)*(plugin.est.z/phi2.Zsim)
      EIF.W <- rowMeans(sweep(p.x.matrix.sim*sweep(mu.matrix.sim, 2, plugin.est.sim, "-"),2, phi2.sim,"/"))

      return(list(EIF.W=EIF.W, plugin.est.z=plugin.est.sim))

    }


    if (!z.density.provided){ # if the z.density function is given

      if (fast){

        if(verbose){cat("Fast integration used. Integration computed with Monte Carlo integration method.")}

        Zsim <- sampleFromDensity(z.density, nMC, range(Z))

        fast.object <- fast.EIF.W(Zsim)

        EIF.W <- fast.object$EIF.W


      }else{

        if(verbose){cat("Integration computed with integrate() function. This step may take a while.\n If you want to speed up the process, set fast=TRUE.")}

        EIF.W <- unlist(lapply(1:n, function(i) integrate(integrand.EIF.W, lower = minZ, upper = maxZ, w=W[i,])$value))

      }

      EIF <- EIF.Y + EIF.X + EIF.W

      if(verbose){cat("Density function of Z is provided. \n Influence function and plug-in estimator computed using the given density function of Z.")}

    }else{ # if the z.density function is NOT given

      EIF.W <- Reduce(`+`, lapply(1:n, function(i) f.x.zw(x, Z[i],W, truncate_lower.X, truncate_upper.X)/phi2[i]*(f.mu.xzw(x, Z[i],W) - plugin.est[i])    ))/n

      EIF <- EIF.Y + EIF.X + EIF.W

      if(verbose){cat("Density function of Z is not provided. Influence function and plug-in estimator computed via density estimation.")}

    }

    # confidence interval
    lower.ci <- estimated-1.96*sqrt(mean(EIF^2)/n)
    upper.ci <- estimated+1.96*sqrt(mean(EIF^2)/n)

    return(list(estimated=estimated, EIF=EIF, EIF.Y=EIF.Y, EIF.X=EIF.X, EIF.W=EIF.W, lower.ci=lower.ci, upper.ci=upper.ci))


  }


    # if Z is binary, return three estimating equation based estimators:
    # 1. at Z=1
    # 2. at Z=0
    # 3. average of the two
    if (binaryZ ){

      out.z1=f.equation.xz_binaryZ(x,z=1, truncate_lower.X, truncate_upper.X ,truncate_lower.Z, truncate_upper.Z)
      out.z0=f.equation.xz_binaryZ(x,z=0, truncate_lower.X, truncate_upper.X ,truncate_lower.Z, truncate_upper.Z)

      ## using all levels of Z

      # point estimate
      # estimated= {mean(out.z1$kappa1_plus_phi1)*mean(Z==1) + mean(out.z0$kappa1_plus_phi1)*mean(Z==0)}/{mean(out.z1$kappa2_plus_phi2)*mean(Z==1) + mean(out.z0$kappa2_plus_phi2)*mean(Z==0)}
      # estimated = mean(out.z1$EIF*mean(Z==1)+out.z0$EIF*mean(Z==0)+Z*out.z1$estimated + (1-Z)*out.z0$estimated)
      estimated = out.z1$estimated*mean(Z==1) + out.z0$estimated*mean(Z==0)

      # EIF
      EIF= out.z1$estimated*{(Z==1)-mean(Z==1)} + out.z0$estimated*{(Z==0)-mean(Z==0)} + out.z1$EIF*mean(Z==1) + out.z0$EIF*mean(Z==0)


      # confidence interval
      lower.ci <- estimated-1.96*sqrt(mean(EIF^2)/n)
      upper.ci <- estimated+1.96*sqrt(mean(EIF^2)/n)

      out.all.z = list(estimated=estimated, EIF=EIF, lower.ci=lower.ci, upper.ci=upper.ci)

      equation.out <- list(out.z1=out.z1, out.z0=out.z0, out.all.z=out.all.z)

    }else{

      if(verbose){cat(paste0('Z is continuous. Computing estimating equation estimator at the ',ifelse(!z.density.provided,'given','estimated'),' density function of Z.'))}

      out.all.z <- f.equation.xz_nonbinaryZ(x, truncate_lower.X, truncate_upper.X ,truncate_lower.Z, truncate_upper.Z, p.z)

      equation.out <- list(out.all.z=out.all.z)

    }

  # when Z is continuous, the estimating equation based estimator is the same as the one-step estimator




  ####################################################
  #################### TMLE ##########################
  ####################################################

    # tmle estimate and its EIF at a given level of z
    f.tmle.xz_binaryZ <- function(x,z, truncate_lower.X, truncate_upper.X ,truncate_lower.Z, truncate_upper.Z){

      # initialize estimates
      p.x.z <- f.x.zw(x, z,W, truncate_lower.X, truncate_upper.X) # propensity score
      p.z <- f.z(z, truncate_lower.Z, truncate_upper.Z) # z-density
      mu.xz <- f.mu.xzw(x, z,W) # outcome regression

      #############################
      # update mu: E(Y|X=x,Z=z,W)
      #############################

      # weight.Y <- {(X==x)*(Z==z)}/{mean(p.x.z)*p.z} # weight for Y
      weight.Y <- {(X==x)*(Z==z)}/{p.z} # weight for Y

      if (binaryY|boundedsubmodelY){ # binary Y

        if(boundedsubmodelY){

          if(verbose){
            cat("boundedsubmodelY is TRUE, rescaling Y to [0,1] \n")
          }

          minY <- min(Y)
          maxY <- max(Y)

          Y <- (Y-minY)/(maxY-minY) # rescale Y to [0,1]
          mu.xz <- f.mu.xzw.bounded(x, z,W) # outcome regression

          mu.xz[mu.xz<0] <- 0.001
          mu.xz[mu.xz>1] <- 0.999

        }

        # one iteration
        or_model <- glm(
          Y ~ offset(qlogis(mu.xz))+weight.Y-1, family=binomial(), start=0
        )

        eps.Y = coef(or_model)

        # updated outcome regression
        mu.xz = plogis(qlogis(mu.xz)+eps.Y*p.z)

        # if(boundedsubmodelY){
        #
        #   Y <- Y*(maxY-minY)+minY # rescale Y back
        #
        #   mu.xz <- mu.xz*(max.mu-min.mu)+min.mu # rescale mu back
        #
        # }

      } else { # continuous Y

        # one iteration
        or_model <- glm(

          Y ~ offset(mu.xz)+1, weights = weight.Y
        )

        eps.Y <- coef(or_model)

        mu.xz <- mu.xz+eps.Y

      } # end of if-else for updating mu

      # update nuisance that depend on outcome regression
      # EIF of phi1 plus phi1
      kappa1_plus_phi1 <- {(Z==z)/p.z}*{(X==x)*Y - mu.xz*p.x.z} + mu.xz*p.x.z

      # EIF of phi2 plus phi2
      kappa2_plus_phi2 <- {(Z==z)/p.z}*( (X==x) - p.x.z ) + p.x.z

      phi1 <- mean(kappa1_plus_phi1) # numerator of the plugin estimator
      phi2 <- mean(kappa2_plus_phi2) # denominator of the plugin estimator

      estimated <- mean(mu.xz*p.x.z)/mean(p.x.z) # point estimate

      EIF.X <- {(Z==z)*(mu.xz-estimated)}/{p.z}*( (X==x) - p.x.z )

      EIF.Y <- {(X==x)*(Z==z)}/{phi2*p.z}*(Y-mu.xz) # EIF for Y


      ######################
      # update p(X=1|Z=z,W,C)
      ######################

      # iter <- 0

      # while(abs(mean(EIF.X)) > cvg.criteria & iter < n.iter){

      # derive eps3
      ps_model <- glm(
        (X==x)  ~ offset(qlogis(p.x.z))+mu.xz+1, family=binomial(), start=c(0,0), weights={(Z==z)/p.z}
      )

      eps.X <- coef(ps_model)


      # updated propensity score
      p.x.z <- plogis(qlogis(p.x.z)+eps.X[2]*mu.xz+eps.X[1])


      # update nuisance that depend on the propensity score

      # EIF of phi1 plus phi1
      kappa1_plus_phi1 <- {(Z==z)/p.z}*{(X==x)*Y - mu.xz*p.x.z} + mu.xz*p.x.z

      # EIF of phi2 plus phi2
      kappa2_plus_phi2 <- {(Z==z)/p.z}*( (X==x) - p.x.z ) + p.x.z

      estimated <- mean(mu.xz*p.x.z)/mean(p.x.z) # point estimate

      EIF.X <- {(Z==z)*(mu.xz-estimated)}/{p.z}*( (X==x) - p.x.z )

      # iter <- iter + 1

      # }

      # update nuisance that depend on the propensity score

      # EIF of phi1 plus phi1
      kappa1_plus_phi1 <- {(Z==z)/p.z}*{(X==x)*Y - mu.xz*p.x.z} + mu.xz*p.x.z

      # EIF of phi2 plus phi2
      kappa2_plus_phi2 <- {(Z==z)/p.z}*( (X==x) - p.x.z ) + p.x.z

      phi2 <- mean(p.x.z) # denominator of the plugin estimator

      # EIF
      EIF.Y <- {(X==x)*(Z==z)}/{mean(kappa2_plus_phi2)*p.z}*(Y-mu.xz) # EIF for Y
      EIF.X <- (Z==z)/{mean(kappa2_plus_phi2)*p.z}*(mu.xz-estimated)*( (X==x) - p.x.z ) # EIF for X
      EIF.W <- 1/mean(kappa2_plus_phi2)*{mu.xz*p.x.z - estimated*p.x.z}

      EIF <- EIF.Y + EIF.X + EIF.W

      EDstar <- mean(EIF)


      # confidence interval
      lower.ci <- estimated-1.96*sqrt(mean(EIF^2)/n)
      upper.ci <- estimated+1.96*sqrt(mean(EIF^2)/n)

      return(list(estimated=estimated, EIF=EIF, EIF.Y=EIF.Y, EIF.X=EIF.X, EIF.W=EIF.W, lower.ci=lower.ci, upper.ci=upper.ci, kappa1_plus_phi1=kappa1_plus_phi1, kappa2_plus_phi2=kappa2_plus_phi2))

    }



  # tmle estimate and its EIF at a given level of z
  f.tmle.x_nonbinaryZ <- function(x, truncate_lower.X, truncate_upper.X ,truncate_lower.Z, truncate_upper.Z,p.z, p.zi){

    # initialize estimates
    p.x.z <- f.x.zw(x, Z,W, truncate_lower.X, truncate_upper.X) # propensity score
    mu.xz <- f.mu.xzw(x, Z,W) # outcome regression

    # make prediction for p(Y_i|X=x,Z_j,W_i) for both i and j from 1 to n
    mu.matrix <- matrix(NA, n, n)

    for (i in 1:n) { # the columns are the predictions for p(Y_i|X=x,Z_j,W_i) for i from 1 to n
      mu.matrix[, i] <- f.mu.xzw(x, Z[i],W)
    }

    # make prediction for p(X=x|Z_j,W_i) for both i and j from 1 to n
    p.x.matrix <- matrix(NA, n, n)

    for (i in 1:n) { # the columns are the predictions for p(X=x|Z_j,W_i)for i from 1 to n
      p.x.matrix[, i] <- f.x.zw(x, Z[i],W, truncate_lower.X, truncate_upper.X)
    }

    phi1 <- colMeans(mu.matrix*p.x.matrix) # numerator of the plugin estimator
    phi2 <- colMeans(p.x.matrix) # denominator of the plugin estimator

    plugin.est <- phi1/phi2 # plug-in estimate

    ## create mu.matrx and p.x.matrix for Zsim
    if(z.density.provided){

      Zsim <- sample(Z, nMC, replace=TRUE) # sample Zsim from the observed data

    }else{

      Zsim <- sampleFromDensity(z.density, nMC, range(Z))

    }


    mu.matrix.sim <- matrix(NA, n, nMC)

    for (i in 1:nMC) { # the columns are the predictions for p(Y_i|X=x,Z_j,W_i) for i from 1 to n, Z_j in Zsim
      mu.matrix.sim[, i] <- f.mu.xzw(x, Zsim[i],W)
    }

    p.x.matrix.sim <- matrix(NA, n, nMC)

    for (i in 1:nMC) { # the columns are the predictions for p(X=x|Z_j,W_i) for i from 1 to n, Z_j in Zsim
      p.x.matrix.sim[, i] <- f.x.zw(x, Zsim[i],W, truncate_lower.X, truncate_upper.X)
    }

    phi1.sim <- colMeans(mu.matrix.sim*p.x.matrix.sim) # numerator of the plugin estimator
    phi2.sim <- colMeans(p.x.matrix.sim) # denominator of the plugin estimator

    plugin.est.sim <- phi1.sim/phi2.sim # plug-in estimate


    if (binaryY){ ########################## binary Y ###########################

      EIF.Y <- {(X==x)*(p.zi)}/{phi2*p.z}*(Y-mu.xz) # EIF for Y

      iter.Y <- 0 # initialize iteration counter

      while(abs(mean(EIF.Y)) > cvg.criteria & iter.Y < n.iter){


        ######################
        # update p(X=1|Z=z,W)
        ######################

        EIF.X <- {(p.zi)*(mu.xz-plugin.est)}/{p.z*phi2}*( (X==x) - p.x.z )

        iter.X <- 0

        while(abs(mean(EIF.X)) > cvg.criteria & iter.X < n.iter){

          weight.X <- p.zi/p.z

          covariate.X <- (mu.xz-plugin.est)/phi2

          # derive eps3
          ps_model <- glm(
            (X==x)  ~ offset(qlogis(p.x.z))+covariate.X-1, family=binomial(), start = 0,weights=weight.X
          )

          eps.X <- coef(ps_model)


          # updated propensity score for pi matrix
          p.x.z <- plogis(qlogis(p.x.z)+eps.X*covariate.X)

          ## update nuisance that depend on the propensity score

          # update p(x|Z_j,W_i)
          for(j in 1:n){

            covariate.Xj <- (mu.matrix[,j]-plugin.est[j])/phi2[j]

            p.x.matrix[,j] <- plogis(qlogis(p.x.matrix[,j])+eps.X*covariate.Xj)

          }

          phi1 <- colMeans(mu.matrix*p.x.matrix) # numerator of the plugin estimator
          phi2 <- colMeans(p.x.matrix) # denominator of the plugin estimator

          plugin.est <- phi1/phi2 # plug-in estimate

          # updated propensity score for Zsim  pi matrix
          for (j in 1:nMC){

            covariate.Xj <- (mu.matrix.sim[,j]-plugin.est.sim[j])/phi2.sim[j]

            p.x.matrix.sim[,j] <- plogis(qlogis(p.x.matrix.sim[,j])+eps.X*covariate.Xj)

          }

          phi1.sim <- colMeans(mu.matrix.sim*p.x.matrix.sim) # numerator of the plugin estimator
          phi2.sim <- colMeans(p.x.matrix.sim) # denominator of the plugin estimator

          plugin.est.sim <- phi1.sim/phi2.sim # plug-in estimate


          EIF.X <- {(p.zi)*(mu.xz-plugin.est)}/{p.z*phi2}*( (X==x) - p.x.z )

          iter.X <- iter.X + 1

          if(verbose){cat("TMLE updating p(X|Z,W), iter:", iter.X, "EIF.X:", mean(EIF.X), "\n")}

        }


        #############################
        # update mu: E(Y|X=x,Z=z,W)
        #############################

        weight.Y <- {(X==x)*(p.zi)}/{phi2*p.z} # weight for Y

        # one iteration
        or_model <- glm(
          Y ~ offset(qlogis(mu.xz))+1, family=binomial(), start=0,weights=weight.Y
        )

        eps.Y = coef(or_model)

        # updated outcome regression
        mu.xz = plogis(qlogis(mu.xz)+eps.Y)

        # update nuisances that depend on the updated mu
        mu.matrix[,j] <- plogis(qlogis(mu.matrix[,j])+eps.Y)

        phi1 <- colMeans(mu.matrix*p.x.matrix) # numerator of the plugin estimator
        phi2 <- colMeans(p.x.matrix) # denominator of the plugin estimator

        plugin.est <- phi1/phi2 # plug-in estimate

        # updated outcome regression for Zsim
        mu.matrix.sim[,j] <- plogis(qlogis(mu.matrix.sim[,j])+eps.Y)

        phi1.sim <- colMeans(mu.matrix.sim*p.x.matrix.sim) # numerator of the plugin estimator
        phi2.sim <- colMeans(p.x.matrix.sim) # denominator of the plugin estimator

        plugin.est.sim <- phi1.sim/phi2.sim # plug-in estimate

        # EIF
        EIF.Y <- {(X==x)*(p.zi)}/{phi2*p.z}*(Y-mu.xz) # EIF for Y

        iter.Y <- iter.Y + 1

        if(verbose){cat("TMLE updating E(Y|X,Z,W), iter:", iter.Y, "EIF.X:", mean(EIF.Y), "\n")}

      } # end of while loop over Y


    }else{ ########################## continuous Y ###########################


      ######################
      # update p(X=1|Z=z,W)
      ######################

      EIF.X <- {(p.zi)*(mu.xz-plugin.est)}/{p.z*phi2}*( (X==x) - p.x.z )

      iter.X <- 0

      while(abs(mean(EIF.X)) > cvg.criteria & iter.X < n.iter){

        weight.X <- p.zi/p.z

        covariate.X <- (mu.xz-plugin.est)/phi2

        # derive eps3
        ps_model <- glm(
          (X==x)  ~ offset(qlogis(p.x.z))+covariate.X-1, family=binomial(), start = 0,weights=weight.X
        )

        eps.X <- coef(ps_model)


        # updated propensity score for pi matrix
        p.x.z <- plogis(qlogis(p.x.z)+eps.X*covariate.X)

        ## update nuisance that depend on the propensity score

        # update p(x|Z_j,W_i)
        for(j in 1:n){

          covariate.Xj <- (mu.matrix[,j]-plugin.est[j])/phi2[j]

          p.x.matrix[,j] <- plogis(qlogis(p.x.matrix[,j])+eps.X*covariate.Xj)

        }

        phi1 <- colMeans(mu.matrix*p.x.matrix) # numerator of the plugin estimator
        phi2 <- colMeans(p.x.matrix) # denominator of the plugin estimator

        plugin.est <- phi1/phi2 # plug-in estimate

        # updated propensity score for Zsim  pi matrix
        for (j in 1:nMC){

          covariate.Xj <- (mu.matrix.sim[,j]-plugin.est.sim[j])/phi2.sim[j]

          p.x.matrix.sim[,j] <- plogis(qlogis(p.x.matrix.sim[,j])+eps.X*covariate.Xj)

        }

        phi1.sim <- colMeans(mu.matrix.sim*p.x.matrix.sim) # numerator of the plugin estimator
        phi2.sim <- colMeans(p.x.matrix.sim) # denominator of the plugin estimator

        plugin.est.sim <- phi1.sim/phi2.sim # plug-in estimate


        EIF.X <- {(p.zi)*(mu.xz-plugin.est)}/{p.z*phi2}*( (X==x) - p.x.z )

        iter.X <- iter.X + 1

        if(verbose){cat("TMLE updating p(X|Z,W), iter:", iter.X, "EIF.X:", mean(EIF.X), "\n")}

      }

      #############################
      # update mu: E(Y|X=x,Z=z,W)
      #############################

      weight.Y <- {(X==x)*(p.zi)}/{phi2*p.z} # weight for Y

      # one iteration
      or_model <- glm(

        Y ~ offset(mu.xz)+1, weights = weight.Y
      )

      eps.Y <- coef(or_model)

      mu.xz <- mu.xz+eps.Y

      # updated outcome regression score for mu matrix
      mu.matrix <- mu.matrix + eps.Y

      phi1 <- colMeans(mu.matrix*p.x.matrix) # numerator of the plugin estimator
      phi2 <- colMeans(p.x.matrix) # denominator of the plugin estimator

      plugin.est <- phi1/phi2 # plug-in estimate


      # updated outcome regression score for Zsim mu matrix
      mu.matrix.sim <- mu.matrix.sim + eps.Y

      phi1.sim <- colMeans(mu.matrix.sim*p.x.matrix.sim) # numerator of the plugin estimator
      phi2.sim <- colMeans(p.x.matrix.sim) # denominator of the plugin estimator

      plugin.est.sim <- phi1.sim/phi2.sim # plug-in estimate


    }


    # EIF
    EIF.Y <- {(X==x)*(p.zi)}/{phi2*p.z}*(Y-mu.xz) # EIF for Y
    EIF.X <- (p.zi)/{phi2*p.z}*(mu.xz-plugin.est)*( (X==x) - p.x.z ) # EIF for X

    if (!z.density.provided){

      EIF.W <- rowMeans(sweep(p.x.matrix.sim*sweep(mu.matrix.sim, 2, plugin.est.sim, "-"),2, phi2.sim,"/"))

      EIF <- EIF.Y + EIF.X + EIF.W

      estimated <- mean(EIF.Y+EIF.X+EIF.W) + mean(plugin.est.sim)

      if(verbose){print("Density function of Z is provided.")}

    }else{

      EIF.W <- rowMeans(sweep(p.x.matrix*sweep(mu.matrix, 2, plugin.est, "-"),2, phi2,"/"))

      EIF <- EIF.Y + EIF.X + EIF.W + plugin.est - mean(plugin.est)

      estimated <- mean(EIF.Y+EIF.X+EIF.W) + mean(plugin.est)

      if(verbose){print("Density function of Z is not provided. Influence function equals plug-in estimator at all observed Z minus average of the plug-in")}


    }



    # confidence interval
    lower.ci <- estimated-1.96*sqrt(mean(EIF^2)/n)
    upper.ci <- estimated+1.96*sqrt(mean(EIF^2)/n)


    # return(list(estimated=estimated, EIF=EIF, EIF.Y=EIF.Y, EIF.X=EIF.X, EIF.W=EIF.W, lower.ci=lower.ci, upper.ci=upper.ci, EDstar.record=EDstar.record, kappa1_plus_phi1=kappa1_plus_phi1, kappa2_plus_phi2=kappa2_plus_phi2))
    return(list(estimated=estimated, EIF=EIF, EIF.Y=EIF.Y, EIF.X=EIF.X, EIF.W=EIF.W, lower.ci=lower.ci, upper.ci=upper.ci))

  }




    # if Z is binary, return three TMLE estimators:
    # 1. at Z=1
    # 2. at Z=0
    # 3. average of the two
    if (binaryZ){

      out.z1=f.tmle.xz_binaryZ(x,z=1, truncate_lower.X, truncate_upper.X ,truncate_lower.Z, truncate_upper.Z)
      out.z0=f.tmle.xz_binaryZ(x,z=0, truncate_lower.X, truncate_upper.X ,truncate_lower.Z, truncate_upper.Z)

      if(boundedsubmodelY){

        out.z1$estimated <- out.z1$estimated*(maxY-minY)+minY
        out.z0$estimated <- out.z0$estimated*(maxY-minY)+minY

        out.z1$EIF <- out.z1$EIF*(maxY-minY)
        out.z0$EIF <- out.z0$EIF*(maxY-minY)

        out.z1$EIF.Y <- out.z1$EIF.Y*(maxY-minY)
        out.z0$EIF.Y <- out.z0$EIF.Y*(maxY-minY)

        out.z1$EIF.X <- out.z1$EIF.X*(maxY-minY)
        out.z0$EIF.X <- out.z0$EIF.X*(maxY-minY)

        out.z1$EIF.W <- out.z1$EIF.W*(maxY-minY)
        out.z0$EIF.W <- out.z0$EIF.W*(maxY-minY)

        out.z1$lower.ci <- out.z1$estimated-1.96*sqrt(mean(out.z1$EIF^2)/n)
        out.z0$lower.ci <- out.z0$estimated-1.96*sqrt(mean(out.z0$EIF^2)/n)

        out.z1$upper.ci <- out.z1$estimated+1.96*sqrt(mean(out.z1$EIF^2)/n)
        out.z0$upper.ci <- out.z0$estimated+1.96*sqrt(mean(out.z0$EIF^2)/n)

      }

      ## using all levels of Z

      # point estimate
      estimated = out.z1$estimated*mean(Z==1) + out.z0$estimated*mean(Z==0)


      # EIF
      EIF= out.z1$estimated*{(Z==1)-mean(Z==1)} + out.z0$estimated*{(Z==0)-mean(Z==0)} + out.z1$EIF*mean(Z==1) + out.z0$EIF*mean(Z==0)

      # confidence interval
      lower.ci <- estimated-1.96*sqrt(mean(EIF^2)/n)
      upper.ci <- estimated+1.96*sqrt(mean(EIF^2)/n)

      out.all.z = list(estimated=estimated, EIF=EIF, lower.ci=lower.ci, upper.ci=upper.ci)

      tmle.out <- list(out.z1=out.z1, out.z0=out.z0, out.all.z=out.all.z)

    }else{

      out.all.z <- f.tmle.x_nonbinaryZ(x, truncate_lower.X, truncate_upper.X ,truncate_lower.Z, truncate_upper.Z, p.z, p.zi)

      tmle.out <- list(out.all.z=out.all.z)
    }



  ## if Z is univariate binary, return
  # 1. one-step estimator
  # 2. TMLE estimator
  # 3. estimated equation

  ## if Z is not univariate binary, return 1 and 2 only because the estimating equation aligns with the one-step estimator
  out <- list(Onestep=onestep.out, TMLE=tmle.out, EstEquation=equation.out)

  return(out)




}

