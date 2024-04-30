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
#' @param covariates A character vector indicating the pre-treatment covariates in the napkin graph.
#' @param z.method A character string indicating the method used to estimate the conditional density of p(Z|W,C). There are two options: 'dnorm' and 'np'. The default is 'dnorm'.
#' If Z is binary, regression based method will be adopted. The formula of the regression can be specified by `formula.Z`.
#' If Z is continuous, the `z.method` method will be adopted. The `dnorm` method assumes that the conditional density of Z given W and C is normal.
#' The `np` method estimates the conditional density of Z given W and C via non-parametric method via the \link[np]{np} package.
#' @details In assuming normal distribution. The `dnorm` method estimates the mean of the normal distribution by fitting a linear regression model with only linear terms and no interaction terms and the variance of the normal distribution sample variance of the error term resulted from the linear regression model.
#' In assuming logistic regression, the `dnorm` method further assume the logistic regression contains only linear terms and no interaction terms.
#' @param superlearner.Y A logical indicator determines whether SuperLearner via the \link[SuperLearner]{SuperLearner} function is adopted for estimating the outcome regression.
#' @param superlearner.X A logical indicator determines whether SuperLearner via the \link[SuperLearner]{SuperLearner} function is adopted for estimating the propensity score.
#' @param superlearner.Z A logical indicator determines whether SuperLearner via the \link[SuperLearner]{SuperLearner} function is adopted for estimating the conditional density p(Z|W,C) when Z is binary.
#' @param crossfit A logical indicator determines whether crossfitting is adopted for SuperLearner. If crossfit is set to TRUE, the data is split into K folds as specified by the `K` parameter.
#' @param K An integer specifying the number of folds for crossfitting.
#' @param lib.Y A character vector specifying the library of algorithms to be used in the SuperLearner for outcome regression.
#' @param lib.X A character vector specifying the library of algorithms to be used in the SuperLearner for propensity score estimation.
#' @param lib.Z A character vector specifying the library of algorithms to be used in the SuperLearner for estimating the conditional density of Z given W and C.
#' @param formula.Y Regression formula for the outcome regression of Y on it's Markov pillow. The default is 'Y ~ .'.
#' @param formula.X Regression formula for the propensity score regression of A on it's Markov pillow. The default is 'X ~ .'.
#' @param formula.Z Regression formula for the conditional density of Z given W and C. The default is 'Z ~ .'.
#' @param linkY_binary The link function used for outcome regression of Y on it's Markov pillow when Y is binary and superlearner is not sued. The default is the 'logit' link.
#' @param link.X The link function used for propensity score regression of X on it's Markov pillow if superlearner is not used. The default is the 'logit' link.
#' @param link.Z The link function used for the regression of Z on it's Markov pillow if superlearner is not used. The default is the 'logit' link.
#' @param cvg.criteria A numerical value representing the convergence criteria for the iterative update of the nusiances in TMLE.
#' If the absolute of the mean of efficient influence function is less than cvg.criteria, the iterative update stops. The default value is 0.01.
#' @param n.iter The maximum number of iterations for the iterative update of the nuisances in TMLE. The default value is 500.
#' @param truncate_lower.X The lower bound for truncation of the propensity score. The default is 0, which means no truncation.
#' @param truncate_upper.X The upper bound for truncation of the propensity score. The default is 1, which means no truncation.
#' @param truncate_lower.Z The lower bound for truncation of the conditional density of Z given W and C. The default is 0, which means no truncation.
#' @param truncate_upper.Z The upper bound for truncation of the conditional density of Z given W and C. The default is 1, which means no truncation.
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
#' res <- TMLE.a(x=1, z = NULL, data=data_Zbinary_Ycontinuous,
#' treatment="X", Z.variables="Z", W.variables="W", outcome="Y", covariates="C",
#' formula.Y="Y ~ C + X*Z*W", formula.X="X ~ Z*W", formula.Z="Z~.",
#' link.X="identity", link.Z="logit")
#' @importFrom dplyr %>% mutate select
#' @importFrom MASS mvrnorm
#' @importFrom SuperLearner CV.SuperLearner SuperLearner
#' @importFrom mvtnorm dmvnorm
#' @importFrom densratio densratio
#' @importFrom utils combn
#' @importFrom stats rnorm runif rbinom dnorm dbinom binomial gaussian predict glm as.formula qlogis plogis lm coef cov sd
#' @export

TMLE.a <- function(x, z = NULL,data,treatment, Z.variables, W.variables, outcome, covariates,
                     z.method="dnorm", superlearner.Y=F, superlearner.X=F,superlearner.Z=F,
                     crossfit=F,K=5,
                     lib.Y = c("SL.glm","SL.earth","SL.ranger","SL.mean"),
                     lib.X = c("SL.glm","SL.earth","SL.ranger","SL.mean"),
                     lib.Z = c("SL.glm","SL.earth","SL.ranger","SL.mean"),
                     n.iter=500, cvg.criteria=0.01,
                     formula.Y="Y ~ .", formula.X="X ~ .", formula.Z="Z~.",
                     linkY_binary="logit", link.X="logit", link.Z="logit",
                     truncate_lower.X=0, truncate_upper.X=1,
                     truncate_lower.Z=0, truncate_upper.Z=1){
#
# library(dplyr)
# library(SuperLearner)
  # data = dt
  # treatment = "X"
  # Z.variables = "Z"
  # W.variables = "W"
  # outcome = "Y"
  # covariates = "C"
  # x=1
  # z=1
  # z.method="dnorm"
  # superlearner=F
  # crossfit=F
  # K=5
  # lib.Y = c("SL.glm","SL.earth","SL.ranger","SL.mean")
  # lib.X = c("SL.glm","SL.earth","SL.ranger","SL.mean")
  # lib.Z = c("SL.glm","SL.earth","SL.ranger","SL.mean")
  # n.iter=500
  # cvg.criteria=0.01
  # formula.Y="Y ~ ."
  # formula.X="X ~ ."
  # formula.Z="Z~."
  # linkY_binary="logit"
  # link.X="identity"
  # link.Z="logit"
  # truncate_lower.X=0
  # truncate_upper.X=1
  # truncate_lower.Z=0
  # truncate_upper.Z=1

  n <- nrow(data)

  # Variables
  C <- data[,covariates, drop = F]
  X <- data[,treatment]
  W <- data[,W.variables, drop = F]
  Z <- data[,Z.variables,drop = F]
  Y <- data[,outcome]


  # new data sets
  dat_mpY = data.frame(X,Z,W, C)
  dat_mpX = data.frame(Z,W, C)
  dat_mpZ = data.frame(W, C)

  # pre function
  I.z <- function(i) {
    vec <- rep(0, n)  # Create a vector of n zeros
    vec[i] <- 1  # Set the ith element to 1
    return(vec)
  }


  ################################################
  ############### OUTCOME REGRESSION #############
  ################################################


  if (crossfit==T){ #### cross fitting + super learner #####

    fit.family <- if(all(Y %in% c(0,1))){binomial(linkY_binary)}else{gaussian()} # family for super learner depending on whether Y is binary or continuous

    or_fit <- CV.SuperLearner(Y=Y, X=dat_mpY, family = fit.family, V = K, SL.library = lib.Y, control = list(saveFitLibrary=T),saveAll = T)

    f.mu.xz<- function(x,z){

      dat_mpY.xz <- dat_mpY %>% mutate(Z=z, X=x)

      mu.Y.xz <- unlist(lapply(1:K, function(x) predict(or_fit$AllSL[[x]], newdata=dat_mpY.xz[or_fit$folds[[x]],])[[1]] %>% as.vector()))[order(unlist(lapply(1:K, function(x) or_fit$folds[[x]])))]


      return(mu.Y.xz)
    } # end of function f.mu.xz


  } else if (superlearner.Y==T){ #### super learner #####

    fit.family <- if(all(Y %in% c(0,1))){binomial(linkY_binary)}else{gaussian()} # family for super learner depending on whether Y is binary or continuous

    or_fit <- SuperLearner(Y=Y, X=dat_mpY, family = fit.family, SL.library = lib.Y)

    f.mu.xz<- function(x,z){

      dat_mpY.xz <- dat_mpY %>% mutate(Z=z, X=x)

      mu.Y.xz <- predict(or_fit, newdata=dat_mpY.xz)[[1]] %>% as.vector()

      return(mu.Y.xz)
    } # end of function f.mu.xz


  } else { #### simple linear regression with user input regression formula: default="Y ~ ." ####

    fit.family <- if(all(Y %in% c(0,1))){binomial(linkY_binary)}else{gaussian()} # family for super learner depending on whether Y is binary or continuous

    or_fit <- glm(as.formula(formula.Y), data=dat_mpY, family = fit.family)

    f.mu.xz <- function(x,z){

      dat_mpY.xz <- dat_mpY %>% mutate(Z=z, X=x)

      mu.Y.xz <- predict(or_fit, newdata=dat_mpY.xz)


      return(mu.Y.xz)
    }



  } # end of if-else for outcome regression fitting method
  print("outcome regression done.")

  ################################################
  ############### PROPENSITY SCORE ###############
  ################################################


  if (crossfit==T){ #### cross fitting + super learner #####

    ps_fit <- CV.SuperLearner(Y=X, X=dat_mpX, family = binomial(), V = K, SL.library = lib.X, control = list(saveFitLibrary=T),saveAll = T)

    f.x.z <- function(x, z, truncate_lower, truncate_upper){

      dat_mpX.z <- dat_mpX %>% mutate(Z=z)

      p.x1.z <- unlist(lapply(1:K, function(x) predict(ps_fit$AllSL[[x]], newdata=dat_mpX.z[ps_fit$folds[[x]],])[[1]] %>% as.vector()))[order(unlist(lapply(1:K, function(x) ps_fit$folds[[x]])))]

      # truncation
      p.x1.z[p.x1.z < truncate_lower] <- truncate_lower
      p.x1.z[p.x1.z > truncate_upper] <- truncate_upper

      p.x.z <- x*p.x1.z + (1-x)*(1-p.x1.z)

      return(p.x.z)

    } # end of function f.pi.z


  } else if (superlearner.X==T){ #### super learner #####

    ps_fit <- SuperLearner(Y=X, X=dat_mpX, family = binomial(), SL.library = lib.X)

    # p(X=1|Z=z,W,C)
    f.x.z <- function(x, z, truncate_lower, truncate_upper){

      dat_mpX.z <- dat_mpX %>% mutate(Z=z)

      p.x1.z <- predict(ps_fit, newdata=dat_mpX.z, type="response")[[1]] %>% as.vector()

      # truncation
      p.x1.z[p.x1.z < truncate_lower] <- truncate_lower
      p.x1.z[p.x1.z > truncate_upper] <- truncate_upper

      p.x.z <- x*p.x1.z + (1-x)*(1-p.x1.z)

      return(p.x.z)

    } # end of function f.pi.z



  } else { #### simple linear regression with user input regression formula: default="A ~ ." ####

    ps_fit <- glm(as.formula(formula.X), data=dat_mpX,  family = binomial(link.X))

    f.x.z <- function(x, z, truncate_lower, truncate_upper){

      dat_mpX.z <- dat_mpX %>% mutate(Z=z)

      p.x1.z <- predict(ps_fit, newdata=dat_mpX.z, type="response")

      # truncation
      p.x1.z[p.x1.z < truncate_lower] <- truncate_lower
      p.x1.z[p.x1.z > truncate_upper] <- truncate_upper

      p.x.z <- x*p.x1.z + (1-x)*(1-p.x1.z)

      return(p.x.z)

    } # end of function f.pi.z

  } # end of if-else for propensity score fitting


  print("propensity score regression done.")


  ################################################
  ############### Z-DENSITY ###############
  ################################################


  ## If Z is binary
  if (length(Z.variables)==1 & all(as.vector(Z[[1]]) %in% c(0,1))){

    Z <- as.vector(Z[[1]])

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

      # p(X=1|Z=z,W,C)
      f.z <- function(z, truncate_lower=0, truncate_upper=1){

        p.z1 <- predict(z_fit, type = "response")[[1]] %>% as.vector()

        # truncation
        p.z1[p.z1 < truncate_lower] <- truncate_lower
        p.z1[p.z1 > truncate_upper] <- truncate_upper

        p.z <- z*p.z1 + (1-z)*(1-p.z1)

        return(p.z)

      } # end of function f.pi.z



    } else { #### simple linear regression with user input regression formula: default="A ~ ." ####

      z_fit <- glm(as.formula(formula.Z), data=dat_mpZ,  family = binomial(link.Z))

      f.z <- function(z, truncate_lower=0, truncate_upper=1){

        p.z1 <- predict(z_fit, type="response")

        # truncation
        p.z1[p.z1 < truncate_lower] <- truncate_lower
        p.z1[p.z1 > truncate_upper] <- truncate_upper

        p.z <- z*p.z1 + (1-z)*(1-p.z1)

        return(p.z)

      } # end of function f.pi.z

    } # end of if-else for propensity score fitting




  ## If Z is NOT binary
  }else{


    if (z.method=="np"){

      print('to be filled')

    }else if (z.method=="dnorm"){

      print('to be filled')

    }else{


      stop("Invalid input. z.method must be either 'np' or 'dnorm'")

    }



  } # end of if-else for binary Z




  ##################################################################
  #################### One-step estimator ##########################
  ##################################################################


    f.onestep.xz <- function(x,z, truncate_lower.X, truncate_upper.X ,truncate_lower.Z, truncate_upper.Z){ # one-step estimate and its EIF

      # nuisance parameters
      p.z <- f.z(z, truncate_lower.Z, truncate_upper.Z) # z-density
      p.x.z <- f.x.z(x, z, truncate_lower.X, truncate_upper.X) # propensity score
      mu.xz <- f.mu.xz(x, z) # outcome regression

      phi1 <- mean(mu.xz*p.x.z) # numerator of the plugin estimator
      phi2 <- mean(p.x.z) # denominator of the plugin estimator

      # EIF of phi1 + phi1
      kappa1_plus_phi1 <- {(Z==z)/p.z}*{(X==x)*Y - mu.xz*p.x.z} + mu.xz*p.x.z

      # EIF of phi2 + phi2
      kappa2_plus_phi2 <- {(Z==z)/p.z}*( (X==x) - p.x.z ) + p.x.z

      # point estimate
      estimated <- mean(kappa1_plus_phi1)/mean(kappa2_plus_phi2)

      # EIF
      EIF.Y <- {(X==x)*(Z==z)}/{phi2*p.z}*(Y-mu.xz) # EIF for Y
      EIF.X <- (Z==z)/{phi2*p.z}*(mu.xz-estimated)*( (X==x) - p.x.z ) # EIF for X
      EIF.W <- 1/phi2*{mu.xz*p.x.z - estimated*p.x.z}

      EIF <- EIF.Y + EIF.X + EIF.W

      # confidence interval
      lower.ci <- estimated-1.96*sqrt(mean(EIF^2)/n)
      upper.ci <- estimated+1.96*sqrt(mean(EIF^2)/n)

      return(list(estimated=estimated, EIF=EIF, EIF.Y=EIF.Y, EIF.X=EIF.X, EIF.W=EIF.W, lower.ci=lower.ci, upper.ci=upper.ci, kappa1_plus_phi1=kappa1_plus_phi1, kappa2_plus_phi2=kappa2_plus_phi2))


    }


    # if Z is binary, return three one-step estimators:
    # 1. at Z=1
    # 2. at Z=0
    # 3. average of the two
    if (length(Z.variables)==1 & all(as.vector(Z[[1]]) %in% c(0,1)) ){

      out.z1=f.onestep.xz(x,z=1, truncate_lower.X, truncate_upper.X ,truncate_lower.Z, truncate_upper.Z)
      out.z0=f.onestep.xz(x,z=0, truncate_lower.X, truncate_upper.X ,truncate_lower.Z, truncate_upper.Z)

      ## using all levels of Z

      # point estimate
      estimated= {mean(out.z1$kappa1_plus_phi1)*mean(Z==1) + mean(out.z0$kappa1_plus_phi1)*mean(Z==0)}/{mean(out.z1$kappa2_plus_phi2)*mean(Z==1) + mean(out.z0$kappa2_plus_phi2)*mean(Z==0)}

      # EIF
      EIF= out.z1$estimated*{(Z==1)-mean(Z==1)} + out.z0$estimated*{(Z==0)-mean(Z==0)} + out.z1$EIF*mean(Z==1) + out.z0$EIF*mean(Z==0)


      # confidence interval
      lower.ci <- estimated-1.96*sqrt(mean(EIF^2)/n)
      upper.ci <- estimated+1.96*sqrt(mean(EIF^2)/n)

      out.all.z = list(estimated=estimated, EIF=EIF, lower.ci=lower.ci, upper.ci=upper.ci)

      onestep.out <- list(out.z1=out.z1, out.z0=out.z0, out.all.z=out.all.z)

    }else{

      ## using all levels of Z

      # point estimate at all different levels of z
      onestep.z <- lapply(Z, function(i) f.onestep.xz(x, z=i, truncate_lower.Z, truncate_upper.Z))

      estimated.z <- unlist(lapply(onestep.z, `[[`, "estimated"))
      EIF.z <- do.call(cbind, lapply(onestep.z, `[[`, "EIF"))

      # combined the above point estimates
      estimated = 1/n*{sum(estimated.z)}

      EIF <-  colSums(  do.call(cbind, lapply(1:n, function(i) estimated.z[i]*{I.z(i)-1/n} + EIF.z[,i]/n)) )

      # confidence interval
      lower.ci <- estimated-1.96*sqrt(mean(EIF^2)/n)
      upper.ci <- estimated+1.96*sqrt(mean(EIF^2)/n)

      out.all.z = list(estimated=estimated, EIF=EIF, lower.ci=lower.ci, upper.ci=upper.ci)


      # if the level of z is provided, the function returns an estimate at the given level of z apart from the combined estimate
      if (!is.null(z)){

        out.z <- f.onestep.xz(x,z=z, truncate_lower.Z, truncate_upper.Z)

        onestep.out <- list(out.z=out.z, out.all.z=out.all.z)

      }


      onestep.out <- out.all.z

    }



  ####################################################
  #################### TMLE ##########################
  ####################################################

    # tmle estimate and its EIF at a given level of z
    f.tmle.xz <- function(x,z, truncate_lower.X, truncate_upper.X ,truncate_lower.Z, truncate_upper.Z){

      # initialize estimates
      p.x.z <- f.x.z(x, z, truncate_lower.X, truncate_upper.X) # propensity score
      p.z <- f.z(z, truncate_lower.Z, truncate_upper.Z) # z-density
      mu.xz <- f.mu.xz(x, z) # outcome regression

      # initialize EDstar, the mean of EIF
      EDstar <- 10 # random large number

      # record EDstar over iterations
      EDstar.record <- data.frame(EDstar.Y=numeric(), EDstar.X=numeric(), EDstar.W=numeric(), EDstar=numeric())

      # initialize iteration counter
      iter <- 0


      while(abs(EDstar) > cvg.criteria & iter < n.iter){

        #############################
        # update mu: E(Y|X=x,Z=z,W,C)
        #############################

        weight.Y <- {(X==x)*(Z==z)}/{mean(p.x.z)*p.z} # weight for Y

        if (all(Y %in% c(0,1))){ # binary Y

          # one iteration
          or_model <- glm(
            Y ~ offset(qlogis(mu.xz))+weight.Y-1, family=binomial(), start=0
          )

          eps.Y = coef(or_model)

          # updated outcome regression
          mu.xz = plogis(qlogis(mu.xz)+eps.Y*weight.Y)

        } else { # continuous Y

          # one iteration
          or_model <- glm(

            Y ~ offset(mu.xz)+1, weights = weight.Y
          )

          eps.Y <- coef(or_model)

          mu.xz <- mu.xz+eps.Y

        } # end of if-else for updating mu

        # update nuisance that depend on outcome regression
        estimated <- mean(mu.xz*p.x.z)/mean(p.x.z)
        weight.X <- {(Z==z)*(mu.xz-estimated)}/{mean(p.x.z)*p.z}



        ######################
        # update p(X=1|Z=z,W,C)
        ######################

        # derive eps3
        ps_model <- glm(
          (X==x)  ~ offset(qlogis(p.x.z))+weight.X-1, family=binomial(), start=0
        )

        eps.X <- coef(ps_model)


        # updated propensity score
        p.x.z <- plogis(qlogis(p.x.z)+eps.X*weight.X)

        # update nuisance that depend on the propensity score

        phi1 <- mean(mu.xz*p.x.z) # numerator of the plugin estimator
        phi2 <- mean(p.x.z) # denominator of the plugin estimator

        # EIF of phi1 plus phi1
        kappa1_plus_phi1 <- {(Z==z)/p.z}*{(X==x)*Y - mu.xz*p.x.z} + mu.xz*p.x.z

        # EIF of phi2 plus phi2
        kappa2_plus_phi2 <- {(Z==z)/p.z}*( (X==x) - p.x.z ) + p.x.z

        estimated <- mean(kappa1_plus_phi1)/mean(kappa2_plus_phi2) # point estimate
        weight.Y <- {(X==x)*(Z==z)}/{mean(p.x.z)*p.z} # weight for Y

        # EIF
        EIF.Y <- weight.Y*(Y-mu.xz) # EIF for Y
        EIF.X <- {(Z==z)*(mu.xz-estimated)}/{mean(p.x.z)*p.z}*( (X==x) - p.x.z ) # EIF for X
        EIF.W <- 1/mean(p.x.z)*{mu.xz*p.x.z - estimated*p.x.z}

        EIF <- EIF.Y + EIF.X + EIF.W

        EDstar <- mean(EIF)

        EDstar.record[(iter+1),] <- c(mean(EIF.Y), mean(EIF.X), mean(EIF.W), EDstar)

        iter <- iter + 1

      }


      # confidence interval
      lower.ci <- estimated-1.96*sqrt(mean(EIF^2)/n)
      upper.ci <- estimated+1.96*sqrt(mean(EIF^2)/n)

      return(list(estimated=estimated, EIF=EIF, EIF.Y=EIF.Y, EIF.X=EIF.X, EIF.W=EIF.W, lower.ci=lower.ci, upper.ci=upper.ci, EDstar.record=EDstar.record, kappa1_plus_phi1=kappa1_plus_phi1, kappa2_plus_phi2=kappa2_plus_phi2))


    }




    # if Z is binary, return three TMLE estimators:
    # 1. at Z=1
    # 2. at Z=0
    # 3. average of the two
    if (length(Z.variables)==1 & all(as.vector(Z[[1]]) %in% c(0,1)) ){

      out.z1=f.tmle.xz(x,z=1, truncate_lower.X, truncate_upper.X ,truncate_lower.Z, truncate_upper.Z)
      out.z0=f.tmle.xz(x,z=0, truncate_lower.X, truncate_upper.X ,truncate_lower.Z, truncate_upper.Z)

      ## using all levels of Z

      # point estimate
      estimated={mean(out.z1$kappa1_plus_phi1)*mean(Z==1) + mean(out.z0$kappa1_plus_phi1)*mean(Z==0)}/{mean(out.z1$kappa2_plus_phi2)*mean(Z==1) + mean(out.z0$kappa2_plus_phi2)*mean(Z==0)}

      # EIF
      EIF= out.z1$estimated*{(Z==1)-mean(Z==1)} + out.z0$estimated*{(Z==0)-mean(Z==0)} + out.z1$EIF*mean(Z==1) + out.z0$EIF*mean(Z==0)

      # confidence interval
      lower.ci <- estimated-1.96*sqrt(mean(EIF^2)/n)
      upper.ci <- estimated+1.96*sqrt(mean(EIF^2)/n)

      out.all.z = list(estimated=estimated, EIF=EIF, lower.ci=lower.ci, upper.ci=upper.ci)

      tmle.out <- list(out.z1=out.z1, out.z0=out.z0, out.all.z=out.all.z)

    }else{

      ## using all levels of Z

      # point estimate at all different levels of z
      tmle.z <- lapply(Z, function(i) f.tmle.xz(x, z=i, truncate_lower.X, truncate_upper.X ,truncate_lower.Z, truncate_upper.Z))

      estimated.z <- unlist(lapply(tmle.z, `[[`, "estimated"))
      EIF.z <- do.call(cbind, lapply(tmle.z, `[[`, "EIF"))

      # combined the above point estimates
      estimated = 1/n*{sum(estimated.z)}

      # EIF
      EIF <-  colSums(  do.call(cbind, lapply(1:n, function(i) estimated.z[i]*{I.z(i)-1/n} + EIF.z[,i]/n)) )

      # confidence interval
      lower.ci <- estimated-1.96*sqrt(mean(EIF^2)/n)
      upper.ci <- estimated+1.96*sqrt(mean(EIF^2)/n)

      out.all.z = list(estimated=estimated, EIF=EIF, lower.ci=lower.ci, upper.ci=upper.ci)


      # if the level of z is provided, the function returns an estimate at the given level of z apart from the combined estimate
      if (!is.null(z)){

        out.z <- f.tmle.xz(x,z=z, truncate_lower.X, truncate_upper.X ,truncate_lower.Z, truncate_upper.Z)

        tmle.out <- list(out.z=out.z, out.all.z=out.all.z)

      }


      tmle.out <- out.all.z

    }



    return(list(Onestep=onestep.out, TMLE=tmle.out))




}

