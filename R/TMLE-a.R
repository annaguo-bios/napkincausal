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
#' @param z_w.method A character string indicating the method used to estimate the conditional density of p(Z|W). There are two options: 'dnorm' and 'np'. The default is 'dnorm'.
#' If Z is binary, regression based method will be adopted. The formula of the regression can be specified by `formula.Z`.
#' If Z is continuous, the `z_w.method` method will be adopted. The `dnorm` method assumes that the conditional density of Z given W is normal.
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
                     z.density=NULL, z_w.method="dnorm", superlearner.Y=F, superlearner.X=F,superlearner.Z=F,
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

  z.density.provided <- !is.null(z.density)

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

  if(verbose){cat("Outcome regression Y|X,Z,W done.\n\n")}

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


  if(verbose){cat("Propensity score X|Z,W done.\n\n")}


  ################################################################
  ############### f_Z(Z|W) Conditional density ###############
  ################################################################


  ############ If Z is binary ############
  if (all(as.vector(Z) %in% c(0,1))){

    if (crossfit==T){ #### cross fitting + super learner #####

      z_fit <- CV.SuperLearner(Y=Z, X=dat_mpZ, family = binomial(), V = K, SL.library = lib.Z, control = list(saveFitLibrary=T),saveAll = T)

      f.z <- function(z, truncate_lower=0, truncate_upper=1){

        fz_w1 <- z_fit$SL.predict

        # truncation
        fz_w1[fz_w1 < truncate_lower] <- truncate_lower
        fz_w1[fz_w1 > truncate_upper] <- truncate_upper

        fz_w <- z*fz_w1 + (1-z)*(1-fz_w1)

        return(fz_w)

      } # end of function f.z


    } else if (superlearner.Z==T){ #### super learner #####

      z_fit <- SuperLearner(Y=Z, X=dat_mpZ, family = binomial(), SL.library = lib.Z)

      # p(X=1|Z=z,W)
      f.z <- function(z, truncate_lower=0, truncate_upper=1){

        fz_w1 <- predict(z_fit, type = "response")[[1]] %>% as.vector()

        # truncation
        fz_w1[fz_w1 < truncate_lower] <- truncate_lower
        fz_w1[fz_w1 > truncate_upper] <- truncate_upper

        fz_w <- z*fz_w1 + (1-z)*(1-fz_w1)

        return(fz_w)

      } # end of function f.pi.z



    } else { #### simple linear regression with user input regression formula: default="A ~ ." ####

      z_fit <- glm(formula.Z, data=dat_mpZ,  family = binomial(link.Z))

      f.z <- function(z, truncate_lower=0, truncate_upper=1){

        fz_w1 <- predict(z_fit, type="response")

        # truncation
        fz_w1[fz_w1 < truncate_lower] <- truncate_lower
        fz_w1[fz_w1 > truncate_upper] <- truncate_upper

        fz_w <- z*fz_w1 + (1-z)*(1-fz_w1)

        return(fz_w)

      } # end of function f.pi.z

    } # end of if-else for Z|W regression fitting




  ########## If Z is NOT binary ############
  }else{


    if (is.function(z_w.method)){

      # make prediction for p(Z_i|W_i)
      fz_w <- sapply(1:n, function(i) z_w.method(Z.variables = Z[i], W.variables = as.vector(W[i,])) )

    }else if (z_w.method=="np"){

      ## Z|W
      # Methods on density estimation:
      # np: https://cran.r-project.org/web/packages/np/np.pdf

      environment(formula.Z) <- environment()

      bw <- eval(substitute(npcdensbw(formula = FORM, data = DATA),
                            list(FORM = formula.Z, DATA = dat_ZmpZ)))

      # bw <- npcdensbw(formula=formula.Z, data=dat_ZmpZ)
      Z_fit <- npcdens(bws=bw)

      # make prediction for p(Z_i|W_i)
      fz_w <- predict(Z_fit)

    }else if (z_w.method=="dnorm"){

      dnorm.density <- calculate_density_dnorm(Z.variables=Z.variables, W.variables=W.variables, data=data, formula.Z=formula.Z, superlearner.Z=superlearner.Z, crossfit=crossfit, K=K)

      # make prediction for p(Z_i|W_i)
      fz_w <- dnorm.density[[1]]

    }else{


      stop("Invalid input. z_w.method must be a user-specified density function or either 'np' or 'dnorm'.\n\n")

    }



  } # end of if-else for binary Z


  if(verbose){cat("Estimation for conditional density of Z completed.\n\n")}



  #############################################
  ############### p(Z) density ###############
  #############################################

  if (z.density.provided){

    p.zi <- sapply(Z, z.density)

  }else{

    den.zi <- density(Z)

    z.density <- function(z){
      approx(den.zi$x, den.zi$y, xout = z, yleft = 0, yright = 0)$y
    }

    p.zi <- z.density(Z)

  }


  if(verbose){cat("Estimation/Evaluation for p(Z) completed.\n\n")}

  ##################################################################
  #################### One-step estimator ##########################
  ##################################################################


  ### Binary Z ###

  f.onestep.xz_binaryZ <- function(x,z, truncate_lower.X, truncate_upper.X ,truncate_lower.Z, truncate_upper.Z){ # one-step estimate and its EIF

    # nuisance parameters
    fz_w <- f.z(z, truncate_lower.Z, truncate_upper.Z) # z-density
    p.x.z <- f.x.zw(x, z,W, truncate_lower.X, truncate_upper.X) # propensity score
    mu.xz <- f.mu.xzw(x, z,W) # outcome regression

    kappa1 <- mean(mu.xz*p.x.z) # numerator of the plugin estimator
    kappa2 <- mean(p.x.z) # denominator of the plugin estimator

    plugin.est <- kappa1/kappa2 # plug-in estimate
    kappa2_aipw <- mean( {(Z==z)/fz_w}*( (X==x) - p.x.z ) + p.x.z )

    # EIF
    EIF.Y <- {(X==x)*(Z==z)}/{kappa2_aipw*fz_w}*(Y-mu.xz) # EIF for Y
    EIF.X <- (Z==z)/{kappa2_aipw*fz_w}*(mu.xz-plugin.est)*( (X==x) - p.x.z ) # EIF for X
    EIF.W <- 1/kappa2_aipw*{mu.xz*p.x.z - plugin.est*p.x.z}

    EIF <- EIF.Y + EIF.X + EIF.W

    # the one-step estimator
    estimated <- plugin.est + mean(EIF)

    # confidence interval
    lower.ci <- estimated-1.96*sqrt(mean(EIF^2)/n)
    upper.ci <- estimated+1.96*sqrt(mean(EIF^2)/n)

    return(list(estimated=estimated, EIF=EIF, EIF.Y=EIF.Y, EIF.X=EIF.X, EIF.W=EIF.W, lower.ci=lower.ci, upper.ci=upper.ci))


  }

  ### Non-binary Z ###
  f.onestep.x_nonbinaryZ <- function(x, truncate_lower.X, truncate_upper.X ,truncate_lower.Z, truncate_upper.Z, fz_w, p.zi){ # one-step estimate and its EIF

    # nuisance parameters
    p.x.z <- f.x.zw(x, Z,W, truncate_lower.X, truncate_upper.X) # propensity score
    mu.xz <- f.mu.xzw(x, Z,W) # outcome regression

    kappa.x <- function(z){
      kappa1 <- sapply(z, function(zi) mean(f.x.zw(x, zi,W, truncate_lower.X, truncate_upper.X)*f.mu.xzw(x, zi,W) ))
      kappa2 <- sapply(z, function(zi) mean(f.x.zw(x, zi,W, truncate_lower.X, truncate_upper.X) ))

      return(list(kappa1=kappa1, kappa2=kappa2))
    }

    # the plug-in: \psi_{x_0}(Q;\tilde{p}_z) = \int kappa1(z) p(z) dz / \int kappa2(z) p(z) dz, where kappa1(z)=E{p(X=x|Z=z,W)*E(Y|X=x,Z=z,W)} and kappa2(z)=E{p(X=x|Z=z,W)}
    if (fast){

      Zsim <- sampleFromDensity(z.density, nMC, range(Z))

      kappa.sim <- kappa.x(Zsim)

      kappa1.tilde <- mean(kappa.sim$kappa1)
      kappa2.tilde <- mean(kappa.sim$kappa2)

    }else{

      integrand.kappa1 <- function(z){
        z.density(z)*kappa.x(z)$kappa1
      }

      integrand.kappa2 <- function(z){
        z.density(z)*kappa.x(z)$kappa2
      }

      kappa1.tilde <- integrate(integrand.kappa1, lower = minZ, upper = maxZ)$value
      kappa2.tilde <- integrate(integrand.kappa2, lower = minZ, upper = maxZ)$value

    }


    plugin.est <- kappa1.tilde/kappa2.tilde # plug-in estimate

    # kappa_2 aipw
    kappa2_aipw <- mean( (p.zi/fz_w)*{(X==x) - p.x.z} + kappa2.tilde)

    # EIF
    EIF.Y <- {(X==x)*(p.zi)}/{kappa2_aipw*fz_w}*(Y-mu.xz) # EIF for Y
    EIF.X <- (p.zi)/{kappa2_aipw*fz_w}*(mu.xz-plugin.est)*( (X==x) - p.x.z ) # EIF for X

    ## functions for computing EIF.W ##

    # integrand for calculating the EIF for W
    integrand.EIF.W <- function(z,w){

      int.w <- z.density(z)*(f.x.zw(x, z, w, truncate_lower.X, truncate_upper.X)/kappa2_aipw*(f.mu.xzw(x, z,w) - plugin.est) )

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

      EIF.W <- rowMeans(p.x.matrix.sim*(mu.matrix.sim - plugin.est))/kappa2.tilde

      return(list(EIF.W=EIF.W, kappa1=phi1.sim, kappa2=phi2.sim))

    }

    # EIF.W
    if (fast){

      if(verbose){cat("Compute EIF.W: fast integration with Monte Carlo method.\n")}

      Zsim <- sampleFromDensity(z.density, nMC, range(Z))

      fast.object <- fast.EIF.W(Zsim)

      EIF.W <- fast.object$EIF.W

    }else{

      if(verbose){cat("Compute EIF.W: integration with integrate().\nThis step may take a while. Set fast=TRUE for fast integration.\n")}

      EIF.W <- unlist(lapply(1:n, function(i) integrate(integrand.EIF.W, lower = minZ, upper = maxZ, w=W[i,])$value))

    }

    EIF <- EIF.Y + EIF.X + EIF.W

    estimated <- mean(EIF.Y+EIF.X+EIF.W) + plugin.est

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

    alpha_opt = mean(out.z0$EIF*(out.z0$EIF-out.z1$EIF))/mean((out.z0$EIF-out.z1$EIF)^2)

    # point estimate
    estimated = out.z1$estimated*alpha_opt + out.z0$estimated*(1-alpha_opt)

    # EIF
    # EIF= out.z1$estimated*{(Z==1)-mean(Z==1)} + out.z0$estimated*{(Z==0)-mean(Z==0)} + out.z1$EIF*mean(Z==1) + out.z0$EIF*mean(Z==0)
    EIF= out.z1$EIF*alpha_opt + out.z0$EIF*(1-alpha_opt)


    # confidence interval
    lower.ci <- estimated-1.96*sqrt(mean(EIF^2)/n)
    upper.ci <- estimated+1.96*sqrt(mean(EIF^2)/n)

    out.all.z = list(estimated=estimated, EIF=EIF, lower.ci=lower.ci, upper.ci=upper.ci)

    onestep.out <- list(out.z1=out.z1, out.z0=out.z0, out.all.z=out.all.z)

  }else{

    if(verbose){cat(paste0('Computing one-step estimator at the ',ifelse(z.density.provided,'given','estimated'),' density function of Z|W.\n'))}

    out.all.z <- f.onestep.x_nonbinaryZ(x, truncate_lower.X, truncate_upper.X ,truncate_lower.Z, truncate_upper.Z, fz_w, p.zi)

    onestep.out <- list(out.all.z=out.all.z)

  }



  ##################################################################
  #################### Estimating Equation Estimator ###############
  ##################################################################


    f.equation.xz_binaryZ <- function(x,z, truncate_lower.X, truncate_upper.X ,truncate_lower.Z, truncate_upper.Z){ # one-step estimate and its EIF

      # nuisance parameters
      fz_w <- f.z(z, truncate_lower.Z, truncate_upper.Z) # z-density
      p.x.z <- f.x.zw(x, z,W,truncate_lower.X, truncate_upper.X) # propensity score
      mu.xz <- f.mu.xzw(x, z,W) # outcome regression

      kappa1.pi <- mean(mu.xz*p.x.z)
      kappa2.pi <- mean(p.x.z)

      plugin.est <- kappa1.pi/kappa2.pi # plug-in estimate

      # estimating equation estimator
      kappa1.aipw <- mean((Z==z)/fz_w*((X==x)*Y - mu.xz*p.x.z) + kappa1.pi)
      kappa2.aipw <- mean((Z==z)/fz_w*((X==x) - p.x.z) + kappa2.pi)

      # point estimate
      estimated <- kappa1.aipw/kappa2.aipw

      # EIF
      EIF.Y <- {(X==x)*(Z==z)}/{kappa2.aipw*fz_w}*(Y-mu.xz) # EIF for Y
      EIF.X <- (Z==z)/{kappa2.aipw*fz_w}*(mu.xz-plugin.est)*( (X==x) - p.x.z ) # EIF for X
      EIF.W <- 1/kappa2.aipw*{mu.xz*p.x.z - plugin.est*p.x.z}

      EIF <- EIF.Y + EIF.X + EIF.W

      # confidence interval
      lower.ci <- estimated-1.96*sqrt(mean(EIF^2)/n)
      upper.ci <- estimated+1.96*sqrt(mean(EIF^2)/n)

      return(list(estimated=estimated, EIF=EIF, EIF.Y=EIF.Y, EIF.X=EIF.X, EIF.W=EIF.W, lower.ci=lower.ci, upper.ci=upper.ci))


    }


  f.equation.xz_nonbinaryZ <- function(x, truncate_lower.X, truncate_upper.X ,truncate_lower.Z, truncate_upper.Z, fz_w){ # one-step estimate and its EIF

    # nuisance parameters
    p.x.z <- f.x.zw(x, Z,W, truncate_lower.X, truncate_upper.X) # propensity score
    mu.xz <- f.mu.xzw(x, Z,W) # outcome regression

    kappa.x <- function(z){
      kappa1 <- sapply(z, function(zi) mean(f.x.zw(x, zi,W, truncate_lower.X, truncate_upper.X)*f.mu.xzw(x, zi,W) ))
      kappa2 <- sapply(z, function(zi) mean(f.x.zw(x, zi,W, truncate_lower.X, truncate_upper.X) ))

      return(list(kappa1=kappa1, kappa2=kappa2))
    }

    plugin.est <- mean(kappa.x(Z)$kappa1)/mean(kappa.x(Z)$kappa2) # plug-in estimate

    if (fast){

      Zsim <- sampleFromDensity(z.density, nMC, range(Z))

      kappa.sim <- kappa.x(Zsim)

      kappa1.pi <- mean(kappa.sim$kappa1)
      kappa2.pi <- mean(kappa.sim$kappa2)

    }else{

      integrand.kappa1 <- function(z){
        z.density(z)*kappa.x(z)$kappa1
      }

      integrand.kappa2 <- function(z){
        z.density(z)*kappa.x(z)$kappa2
      }

      kappa1.pi <- integrate(integrand.kappa1, lower = minZ, upper = maxZ)$value
      kappa2.pi <- integrate(integrand.kappa2, lower = minZ, upper = maxZ)$value

    }

    # estimating equation estimator
    kappa1.aipw <- mean(p.zi/fz_w*((X==x)*Y - mu.xz*p.x.z)) + kappa1.pi
    kappa2.aipw <- mean(p.zi/fz_w*((X==x) - p.x.z)) + kappa2.pi

    estimated <- kappa1.aipw/kappa2.aipw


    ## calculate EIF for inference
    # EIF
    EIF.Y <- {(X==x)*(p.zi)}/{kappa2.aipw*fz_w}*(Y-mu.xz) # EIF for Y
    EIF.X <- (p.zi)/{kappa2.aipw*fz_w}*(mu.xz-plugin.est)*( (X==x) - p.x.z ) # EIF for X

    ## functions for computing EIF.W ##

    # integrand for calculating the EIF for W
    integrand.EIF.W <- function(z,w){

      int.w <- z.density(z)*(f.x.zw(x, z, w, truncate_lower.X, truncate_upper.X)/kappa2.aipw*(f.mu.xzw(x, z,w) - plugin.est) )

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

      EIF.W <- rowMeans(p.x.matrix.sim*(mu.matrix.sim - plugin.est))/kappa2.aipw

      return(list(EIF.W=EIF.W))

    }

    # EIF.W
    if (fast){

      if(verbose){cat("Compute EIF.W: fast integration with Monte Carlo method.\n")}

      Zsim <- sampleFromDensity(z.density, nMC, range(Z))

      fast.object <- fast.EIF.W(Zsim)

      EIF.W <- fast.object$EIF.W

    }else{

      if(verbose){cat("Compute EIF.W: integration with integrate().\nThis step may take a while. Set fast=TRUE for fast integration.")}

      EIF.W <- unlist(lapply(1:n, function(i) integrate(integrand.EIF.W, lower = minZ, upper = maxZ, w=W[i,])$value))

    }

    EIF <- EIF.Y + EIF.X + EIF.W

    estimated <- mean(EIF.Y+EIF.X+EIF.W) + plugin.est

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

      # optimal weight
      alpha_opt = mean(out.z0$EIF*(out.z0$EIF-out.z1$EIF))/mean((out.z0$EIF-out.z1$EIF)^2)

      # point estimate
      estimated = out.z1$estimated*alpha_opt + out.z0$estimated*(1-alpha_opt)

      # EIF
      # EIF= out.z1$estimated*{(Z==1)-mean(Z==1)} + out.z0$estimated*{(Z==0)-mean(Z==0)} + out.z1$EIF*mean(Z==1) + out.z0$EIF*mean(Z==0)
      EIF= out.z1$EIF*alpha_opt + out.z0$EIF*(1-alpha_opt)

      # confidence interval
      lower.ci <- estimated-1.96*sqrt(mean(EIF^2)/n)
      upper.ci <- estimated+1.96*sqrt(mean(EIF^2)/n)

      out.all.z = list(estimated=estimated, EIF=EIF, lower.ci=lower.ci, upper.ci=upper.ci)

      equation.out <- list(out.z1=out.z1, out.z0=out.z0, out.all.z=out.all.z)

    }else{

      if(verbose){cat(paste0('Computing estimating equation estimator at the ',ifelse(z.density.provided,'given','estimated'),' density function of Z|W.\n'))}

      out.all.z <- f.equation.xz_nonbinaryZ(x, truncate_lower.X, truncate_upper.X ,truncate_lower.Z, truncate_upper.Z, fz_w)

      equation.out <- list(out.all.z=out.all.z)

    }

  # when Z is continuous, the estimating equation based estimator is the same as the one-step estimator




  ####################################################
  #################### TMLE ##########################
  ####################################################

    f.tmle.xZ <- function(x, z=NULL, truncate_lower.X, truncate_upper.X ,truncate_lower.Z, truncate_upper.Z, fz_w=NULL, p.zi=NULL){

      binaryZ.target <- !is.null(z)

      weight.z <- if(binaryZ.target){
        fz_w <- f.z(z, truncate_lower.Z, truncate_upper.Z)
        as.numeric(Z==z)/fz_w
      }else{
        p.zi/fz_w
      }

      p.x.z <- f.x.zw(x, if(binaryZ.target){z}else{Z}, W, truncate_lower.X, truncate_upper.X)
      mu.xz <- f.mu.xzw(x, if(binaryZ.target){z}else{Z}, W)

      if(!binaryZ.target){
        Zsim <- sampleFromDensity(z.density, nMC, range(Z))

        mu.matrix <- matrix(NA, n, nMC)
        p.x.matrix <- matrix(NA, n, nMC)

        for (j in 1:nMC) {
          mu.matrix[,j] <- f.mu.xzw(x, Zsim[j], W)
          p.x.matrix[,j] <- f.x.zw(x, Zsim[j], W, truncate_lower.X, truncate_upper.X)
        }
      }

      get.plugin <- function(){
        if(binaryZ.target){
          return(mean(mu.xz*p.x.z)/mean(p.x.z))
        }

        kappa1 <- colMeans(mu.matrix*p.x.matrix)
        kappa2 <- colMeans(p.x.matrix)

        return(mean(kappa1)/mean(kappa2))
      }

      get.kappa2 <- function(){
        if(binaryZ.target){
          return(mean(p.x.z))
        }

        return(mean(colMeans(p.x.matrix)))
      }

      get.kappa2_aipw <- function(){
        if(binaryZ.target){
          return(mean(weight.z*{(X==x) - p.x.z}) + mean(p.x.z))
        }

        return(mean(weight.z*{(X==x) - p.x.z})+mean(colMeans(p.x.matrix)))
      }

      update.mu.matrix <- function(eps.Y){
        if(!binaryZ.target){
          if(binaryY|boundedsubmodelY){
            mu.matrix <<- plogis(qlogis(mu.matrix)+eps.Y)
          }else{
            mu.matrix <<- mu.matrix+eps.Y
          }
        }
      }

      plugin.est <- get.plugin()

      #############################
      # update mu: E(Y|X=x,Z,W)
      #############################

      weight.Y <- (X==x)*weight.z

      if (binaryY|boundedsubmodelY){

        if(boundedsubmodelY){

          if(verbose){cat("boundedsubmodelY=TRUE, rescaling Y to [0,1] \n")}

          Y <- (Y-minY)/(maxY-minY)
          mu.xz <- f.mu.xzw.bounded(x, if(binaryZ.target){z}else{Z}, W)
          mu.xz[mu.xz<0.001] <- 0.001
          mu.xz[mu.xz>0.999] <- 0.999
        }

        or_model <- glm(
          Y ~ offset(qlogis(mu.xz))+1, family=binomial(), start=0, weights=weight.Y
        )

        eps.Y <- coef(or_model)
        mu.xz <- plogis(qlogis(mu.xz)+eps.Y)
        update.mu.matrix(eps.Y)

      }else{

        or_model <- glm(
          Y ~ offset(mu.xz)+1, weights=weight.Y
        )

        eps.Y <- coef(or_model)
        mu.xz <- mu.xz+eps.Y
        update.mu.matrix(eps.Y)

      }

      plugin.est <- get.plugin()

      #############################
      # update pi: p(X=x|Z,W)
      #############################

      EIF.X <- 10 # just to enter the while loop
      iter.X <- 0

      while(abs(mean(EIF.X)) > cvg.criteria & iter.X < n.iter){

        covariate.X <- mu.xz-plugin.est

        ps_model <- glm(
          (X==x) ~ offset(qlogis(p.x.z))+covariate.X-1,
          family=binomial(), start=0, weights=weight.z
        )

        eps.X <- coef(ps_model)
        p.x.z <- plogis(qlogis(p.x.z)+eps.X*covariate.X)

        if(!binaryZ.target){
          for(j in 1:nMC){
            covariate.Xj <- mu.matrix[,j]-plugin.est
            p.x.matrix[,j] <- plogis(qlogis(p.x.matrix[,j])+eps.X*covariate.Xj)
          }
        }

        plugin.est <- get.plugin()

        EIF.X <- weight.z*(mu.xz-plugin.est)*((X==x)-p.x.z)

        iter.X <- iter.X + 1

        if(verbose){cat("TMLE updating p(X|Z,W), iter:", iter.X, "EIF.X:", mean(EIF.X), "\n")}

      }

      kappa2_aipw <- get.kappa2_aipw()

      EIF.Y <- (X==x)*weight.z/kappa2_aipw*(Y-mu.xz)
      EIF.X <- weight.z/kappa2_aipw*(mu.xz-plugin.est)*((X==x)-p.x.z)

      if(binaryZ.target){
        EIF.W <- 1/kappa2_aipw*(mu.xz*p.x.z-plugin.est*p.x.z)
        EIF <- EIF.Y + EIF.X + EIF.W
      }else{
        EIF.W <- rowMeans(p.x.matrix*(mu.matrix-plugin.est))/kappa2_aipw
        EIF <- EIF.Y + EIF.X + EIF.W
      }

      estimated <- plugin.est

      lower.ci <- estimated-1.96*sqrt(mean(EIF^2)/n)
      upper.ci <- estimated+1.96*sqrt(mean(EIF^2)/n)

      return(list(estimated=estimated, EIF=EIF, EIF.Y=EIF.Y, EIF.X=EIF.X, EIF.W=EIF.W, lower.ci=lower.ci, upper.ci=upper.ci))

    }




    # if Z is binary, return three TMLE estimators:
    # 1. at Z=1
    # 2. at Z=0
    # 3. average of the two
    if (binaryZ){

      out.z1=f.tmle.xZ(x,z=1, truncate_lower.X, truncate_upper.X ,truncate_lower.Z, truncate_upper.Z)
      out.z0=f.tmle.xZ(x,z=0, truncate_lower.X, truncate_upper.X ,truncate_lower.Z, truncate_upper.Z)

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

      # optimal weight
      alpha_opt = mean(out.z0$EIF*(out.z0$EIF-out.z1$EIF))/mean((out.z0$EIF-out.z1$EIF)^2)

      # point estimate
      estimated = out.z1$estimated*alpha_opt + out.z0$estimated*(1-alpha_opt)

      # EIF
      EIF= out.z1$EIF*alpha_opt + out.z0$EIF*(1-alpha_opt)

      # confidence interval
      lower.ci <- estimated-1.96*sqrt(mean(EIF^2)/n)
      upper.ci <- estimated+1.96*sqrt(mean(EIF^2)/n)

      out.all.z = list(estimated=estimated, EIF=EIF, lower.ci=lower.ci, upper.ci=upper.ci)

      tmle.out <- list(out.z1=out.z1, out.z0=out.z0, out.all.z=out.all.z)

    }else{

      out.all.z <- f.tmle.xZ(x, truncate_lower.X=truncate_lower.X, truncate_upper.X=truncate_upper.X ,truncate_lower.Z=truncate_lower.Z, truncate_upper.Z=truncate_upper.Z, fz_w=fz_w, p.zi=p.zi)

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
