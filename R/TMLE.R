#' Function for average counterfactual outcome E(Y(x)) and average causal effect (ACE) estimation. Most of the arguments are passed directly to TMLE.a.
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
#' @param covariates Optional character vector indicating baseline covariates C. If non-NULL, the covariate-aware estimator `napkin.a.covariate()` is used.
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
#' @param superlearner.W A logical indicator determines whether SuperLearner is adopted for estimating the conditional probability p(W|C) when covariates are supplied.
#' @param crossfit A logical indicator determines whether crossfitting is adopted for SuperLearner. If crossfit is set to TRUE, the data is split into K folds as specified by the `K` parameter.
#' @param K An integer specifying the number of folds for crossfitting.
#' @param lib.Y A character vector specifying the library of algorithms to be used in the SuperLearner for outcome regression.
#' @param lib.X A character vector specifying the library of algorithms to be used in the SuperLearner for propensity score estimation.
#' @param lib.Z A character vector specifying the library of algorithms to be used in the SuperLearner for estimating the conditional density of Z given W.
#' @param lib.W A character vector specifying the library of algorithms to be used in the SuperLearner for estimating W given covariates C.
#' @param formula.Y Regression formula for the outcome regression of Y on it's Markov pillow. The default is 'Y ~ .'.
#' @param formula.X Regression formula for the propensity score regression of A on it's Markov pillow. The default is 'X ~ .'.
#' @param formula.Z Regression formula for the conditional density of Z given W. The default is 'Z ~ .'.
#' @param formula.W Regression formula for W given covariates C when covariates are supplied. The default is 'W ~ .'.
#' @param linkY_binary The link function used for outcome regression of Y on it's Markov pillow when Y is binary and superlearner is not sued. The default is the 'logit' link.
#' @param link.X The link function used for propensity score regression of X on it's Markov pillow if superlearner is not used. The default is the 'logit' link.
#' @param link.Z The link function used for the regression of Z on it's Markov pillow if superlearner is not used. The default is the 'logit' link.
#' @param link.W The link function used for the regression of W on covariates C if superlearner is not used. The default is the 'logit' link.
#' @param cvg.criteria A numerical value representing the convergence criteria for the iterative update of the nusiances in TMLE.
#' If the absolute of the mean of efficient influence function is less than cvg.criteria, the iterative update stops. The default value is 0.01.
#' @param n.iter The maximum number of iterations for the iterative update of the nuisances in TMLE. The default value is 500.
#' @param truncate_lower.X The lower bound for truncation of the propensity score. The default is 0, which means no truncation.
#' @param truncate_upper.X The upper bound for truncation of the propensity score. The default is 1, which means no truncation.
#' @param truncate_lower.Z The lower bound for truncation of the conditional density of Z given W. The default is 0, which means no truncation.
#' @param truncate_upper.Z The upper bound for truncation of the conditional density of Z given W. The default is 1, which means no truncation.
#' @param truncate_lower.W The lower bound for truncation of p(W=1|C) when covariates are supplied. The default is 0.
#' @param truncate_upper.W The upper bound for truncation of p(W=1|C) when covariates are supplied. The default is 1.
#' @param minZ The lower bound used for performing integration of Z when Z is continuous. The default is -Inf.
#' @param maxZ The upper bound used for performing integration of Z when Z is continuous. The default is Inf.
#' @param verbose A logical indicator determines whether the function prints out detailed progress of the estimation. The default is TRUE.
#' @param fast A logical indicator determines whether the integration involved in the estimation is performed via `integrate()` function or via Monte Carlo integration. The former is lower while the later is faster. The default is TRUE.
#' @param nMC The number of Monte Carlo samples used for integration when `fast` is set to TRUE The default is 5000.
#' @param boundedsubmodelY An indicator for whether the bounded submodel is used for targeting the outcome regression when Z is discrete. The default is FALSE.
#' @param use.kappa2.aipw A logical indicator for whether the covariate-aware estimator uses the AIPW version of kappa2(C) in denominators. The default is FALSE.
#' @return Function outputs a list containing TMLE results (and Onestep results if 'onestep=T' is specified). When 'a=c(1,0)', function also outputs corresponding results on \eqn{E(Y^1)} and \eqn{E(Y^0)}:
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
#' res <- napkin_est(x=1, z = NULL, data=data_Zbinary_Ycontinuous,
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
#' @export

napkin_est <- function(x, z = NULL,data,treatment, Z.variables, W.variables, outcome,
                       covariates = NULL,
                       z.density=NULL, z_w.method="dnorm", superlearner.Y=F, superlearner.X=F,superlearner.Z=F,
                       superlearner.W=F,
                       crossfit=F,K=5,
                       lib.Y = c("SL.glm","SL.earth","SL.ranger","SL.mean"),
                       lib.X = c("SL.glm","SL.earth","SL.ranger","SL.mean"),
                       lib.Z = c("SL.glm","SL.earth","SL.ranger","SL.mean"),
                       lib.W = c("SL.glm","SL.earth","SL.ranger","SL.mean"),
                       n.iter=500, cvg.criteria=0.01,
                       formula.Y="Y ~ .", formula.X="X ~ .", formula.Z="Z~.", formula.W="W ~ .",
                       linkY_binary="logit", link.X="logit", link.Z="logit", link.W="logit",
                       truncate_lower.X=0, truncate_upper.X=1,
                       truncate_lower.Z=0, truncate_upper.Z=1,
                       truncate_lower.W=0, truncate_upper.W=1,
                       minZ=-Inf,maxZ=Inf,
                       verbose=T,fast=T,nMC=5000,boundedsubmodelY=F,
                       use.kappa2.aipw=F){

  # sample size

  n <- nrow(data)

  estimate_a <- function(x.value) {
    if (is.null(covariates)) {
      return(napkin.a(x.value, z = z, data, treatment, Z.variables, W.variables, outcome,
                      z.density = z.density, z_w.method = z_w.method,
                      superlearner.Y = superlearner.Y, superlearner.X = superlearner.X,
                      superlearner.Z = superlearner.Z,
                      crossfit = crossfit, K = K,
                      lib.Y = lib.Y, lib.X = lib.X, lib.Z = lib.Z,
                      n.iter = n.iter, cvg.criteria = cvg.criteria,
                      formula.Y = formula.Y, formula.X = formula.X, formula.Z = formula.Z,
                      linkY_binary = linkY_binary, link.X = link.X, link.Z = link.Z,
                      truncate_lower.X = truncate_lower.X, truncate_upper.X = truncate_upper.X,
                      truncate_lower.Z = truncate_lower.Z, truncate_upper.Z = truncate_upper.Z,
                      minZ = minZ, maxZ = maxZ, verbose = verbose, fast = fast,
                      nMC = nMC, boundedsubmodelY = boundedsubmodelY))
    }

    napkin.a.covariate(x.value, z = z, data = data, treatment = treatment,
                       Z.variables = Z.variables, W.variables = W.variables,
                       covariates = covariates, outcome = outcome,
                       z.density = z.density, z_w.method = z_w.method,
                       superlearner.Y = superlearner.Y, superlearner.X = superlearner.X,
                       superlearner.Z = superlearner.Z, superlearner.W = superlearner.W,
                       crossfit = crossfit, K = K,
                       lib.Y = lib.Y, lib.X = lib.X, lib.Z = lib.Z, lib.W = lib.W,
                       n.iter = n.iter, cvg.criteria = cvg.criteria,
                       formula.Y = formula.Y, formula.X = formula.X,
                       formula.Z = formula.Z, formula.W = formula.W,
                       linkY_binary = linkY_binary, link.X = link.X,
                       link.Z = link.Z, link.W = link.W,
                       truncate_lower.X = truncate_lower.X,
                       truncate_upper.X = truncate_upper.X,
                       truncate_lower.Z = truncate_lower.Z,
                       truncate_upper.Z = truncate_upper.Z,
                       truncate_lower.W = truncate_lower.W,
                       truncate_upper.W = truncate_upper.W,
                       minZ = minZ, maxZ = maxZ, verbose = verbose,
                       fast = fast, nMC = nMC, boundedsubmodelY = boundedsubmodelY,
                       use.kappa2.aipw = use.kappa2.aipw)
  }


  ################################## ATE

  if (is.vector(x) & length(x)>2){ ## Invalid input ==

    stop("Invalid input. Enter x=c(1,0) for Average Causal Effect estimation. Enter x=1 or x=0 for average counterfactual outcome estimation at the specified treatment level.")

  }else if (is.vector(x) & length(x)==2){ ## ATE estimate ==

    ## TMLE estimator

    out.a1 <- estimate_a(x[1])

    out.a0 <- estimate_a(x[2])


    # levels of z used for estimation
    level.z <- sub("^out", "", names(out.a1$TMLE))

    level.z.number <- sub("^out\\.z", "", names(out.a1$TMLE))
    level.z.number[level.z.number=="out.all.z"] <- "all z"

    estimators <- names(out.a1)[names(out.a1) %in% c("TMLE", "Onestep", "EstEquation")]
    output <- vector("list", length(estimators)*length(level.z)+2)

    output[[1]] <- out.a1
    output[[2]] <- out.a0
    names(output)[1:2] <- c("est.Y1","est.Y0")

    # count of method
    count <- 3

    for (m in estimators){ # loop over TMLE and onestep

        for (est in level.z){ # loop over levels of z

          out.est1 <- out.a1[[m]][[paste0('out',est)]] # get output for Y(1)
          out.est0 <- out.a0[[m]][[paste0('out',est)]] # get output for Y(0)

          hat.ATE <- out.est1$estimated - out.est0$estimated
          hat.EIF.ATE <- out.est1$EIF - out.est0$EIF
          lower.ci.ATE <- hat.ATE - 1.96*sqrt(mean(hat.EIF.ATE^2)/n)
          upper.ci.ATE <- hat.ATE + 1.96*sqrt(mean(hat.EIF.ATE^2)/n)

          out.ate <- list(ATE=hat.ATE,
                          lower.ci=lower.ci.ATE,
                          upper.ci=upper.ci.ATE,
                          EIF=hat.EIF.ATE)

          output[[count]] <- out.ate
          names(output)[count] <- paste0(m,est)
          count <- count + 1

          indicator <- which(level.z==est)

          # print estimates
          cat(paste0(m," estimated ACE for z=",level.z.number[indicator],": ",round(out.ate$ATE,2),"; 95% CI: (",round(out.ate$lower.ci,2),", ",round(out.ate$upper.ci,2),") \n"))




        }

    }

    return(output)

  } # end of ATE estimate




 ################################## E(Y(x))


  if (length(x)==1) { ## E(Y^1) estimate ==

    out.a <- estimate_a(x)

    # levels of z used for estimation
    level.z <- sub("^out", "", names(out.a$TMLE))


    level.z.number <- sub("^out\\.z", "", names(out.a$TMLE))
    level.z.number[level.z.number=="out.all.z"] <- "all z"


    estimators <- names(out.a)[names(out.a) %in% c("TMLE", "Onestep", "EstEquation")]
    output <- vector("list", length(estimators)*length(level.z))

    # count of method
    count <- 1

    for (m in estimators){ # loop over TMLE and onestep

      for (est in level.z){ # loop over levels of z

        out.est <- out.a[[m]][[paste0('out',est)]] # get output for Y(x)
        output[[count]] <- out.est
        names(output)[count] <- paste0(m,est)
        count <- count + 1

        indicator <- which(level.z==est)
        # print estimates
        cat(paste0(m," estimated E(Y(",x,")) for z=",level.z.number[indicator],": ",round(out.est$estimated,2),";95% CI: (",round(out.est$lower.ci,2),",",round(out.est$upper.ci,2),") \n"))




      }

    }

    return(output)

  } # end of E(Y(x)) estimate



} # end of TMLE function
