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
#' @param outcome A character string indicating the outcome variable
#' @param covariates A character vector indicating the pre-treatment covariates in the napkin graph.
#' @param z.density A argument that takes the user specified function for estimating the density of Z.
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
#' @param minZ The lower bound used for performing integration of Z when Z is continuous. The default is -Inf.
#' @param maxZ The upper bound used for performing integration of Z when Z is continuous. The default is Inf.
#' @param verbose A logical indicator determines whether the function prints out detailed progress of the estimation. The default is TRUE.
#' @param fast A logical indicator determines whether the integration involved in the estimation is performed via `integrate()` function or via Monte Carlo integration. The former is lower while the later is faster. The default is TRUE.
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
#' treatment="X", Z.variables="Z", W.variables="W", outcome="Y", covariates="C",
#' formula.Y="Y ~ C + X*Z*W", formula.X="X ~ Z*W", formula.Z="Z~.",
#' link.X="identity", link.Z="logit")
#' @importFrom dplyr %>% mutate select
#' @importFrom MASS mvrnorm
#' @importFrom SuperLearner CV.SuperLearner SuperLearner
#' @importFrom mvtnorm dmvnorm
#' @importFrom densratio densratio
#' @importFrom utils combn
#' @importFrom stats rnorm runif rbinom dnorm dbinom binomial gaussian predict glm as.formula qlogis plogis lm coef cov sd density approx integrate
#' @importFrom np npcdensbw npcdens
#' @export

napkin_est <- function(x, z = NULL,data,treatment, Z.variables, W.variables, outcome, covariates,
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
                       verbose=T,fast=T){

  # sample size

  n <- nrow(data)

  univariate.binaryZ <- ( length(Z.variables)==1 & all(as.vector(data[,Z.variables[1]]) %in% c(0,1)) )


  ################################## ATE

  if (is.vector(x) & length(x)>2){ ## Invalid input ==

    stop("Invalid input. Enter x=c(1,0) for Average Causal Effect estimation. Enter x=1 or x=0 for average counterfactual outcome estimation at the specified treatment level.")

  }else if (is.vector(x) & length(x)==2){ ## ATE estimate ==

    ## TMLE estimator

    out.a1 <- napkin.a(x[1], z =z, data, treatment, Z.variables, W.variables, outcome, covariates,
                       z.density=z.density,z.method = z.method, superlearner.Y = superlearner.Y, superlearner.X = superlearner.X, superlearner.Z = superlearner.Z,
                       crossfit = crossfit, K = K,
                       lib.Y = lib.Y,
                       lib.X = lib.X,
                       lib.Z = lib.Z,
                       n.iter = n.iter, cvg.criteria = cvg.criteria,
                       formula.Y = formula.Y, formula.X = formula.X, formula.Z = formula.Z,
                       linkY_binary = linkY_binary, link.X = link.X, link.Z = link.Z,
                       truncate_lower.X = truncate_lower.X, truncate_upper.X = truncate_upper.X,
                       truncate_lower.Z = truncate_lower.Z, truncate_upper.Z = truncate_upper.Z, minZ=minZ,maxZ=maxZ,verbose=verbose,fast=fast)

    out.a0 <- napkin.a(x[2], z =z, data, treatment, Z.variables, W.variables, outcome, covariates,
                       z.density=z.density,z.method = z.method, superlearner.Y = superlearner.Y, superlearner.X = superlearner.X, superlearner.Z = superlearner.Z,
                       crossfit = crossfit, K = K,
                       lib.Y = lib.Y,
                       lib.X = lib.X,
                       lib.Z = lib.Z,
                       n.iter = n.iter, cvg.criteria = cvg.criteria,
                       formula.Y = formula.Y, formula.X = formula.X, formula.Z = formula.Z,
                       linkY_binary = linkY_binary, link.X = link.X, link.Z = link.Z,
                       truncate_lower.X = truncate_lower.X, truncate_upper.X = truncate_upper.X,
                       truncate_lower.Z = truncate_lower.Z, truncate_upper.Z = truncate_upper.Z, minZ=minZ,maxZ=maxZ,verbose=verbose,fast=fast)


    # run TMLE
    TMLE_output_Y1 <- out.a1$TMLE
    TMLE_output_Y0 <- out.a0$TMLE

    # run onestep
    Onestep_output_Y1 <- out.a1$Onestep
    Onestep_output_Y0 <- out.a0$Onestep

    if (univariate.binaryZ){

      # run estimating equation
      EstEquation_output_Y1 <- out.a1$EstEquation
      EstEquation_output_Y0 <- out.a0$EstEquation

    }

    # levels of z used for estimation
    level.z <- sub("^out", "", names(TMLE_output_Y1))

    level.z.number <- sub("^out\\.z", "", names(TMLE_output_Y1))
    level.z.number[level.z.number=="out.all.z"] <- "all z"

    if (univariate.binaryZ){output <- vector("list", 3*length(level.z)+2)}else{output <- vector("list", 2*length(level.z)+2)}
    if (univariate.binaryZ){estimators <- c('TMLE','Onestep','EstEquation')}else{estimators <- c('TMLE','Onestep')}

    output[[1]] <- out.a1
    output[[2]] <- out.a0
    names(output)[1:2] <- c("est.Y1","est.Y0")

    # count of method
    count <- 3

    for (m in estimators){ # loop over TMLE and onestep

        for (est in level.z){ # loop over levels of z

          out.est1 <- get(paste0(m,'_output_Y1'))[[paste0('out',est)]] # get output for Y(1)
          out.est0 <- get(paste0(m,'_output_Y0'))[[paste0('out',est)]] # get output for Y(0)

          # assign(paste0(m,'_output_Y',i,est), out.est)

          # estimate E[Y(1)], E[Y(0)], and ATE
          assign(paste0('hat_E.Y1',est), out.est1$estimated)
          assign(paste0('hat_E.Y0',est), out.est0$estimated)
          assign(paste0('hat_ATE',est), get(paste0('hat_E.Y1',est)) - get(paste0('hat_E.Y0',est)))


          # estimated EIF
          assign(paste0('hat_EIF.Y1',est), out.est1$EIF)
          assign(paste0('hat_EIF.Y0',est), out.est0$EIF)
          assign(paste0('hat_EIF.ATE',est), get(paste0('hat_EIF.Y1',est)) - get(paste0('hat_EIF.Y0',est)))

          ## CI
          assign(paste0('lower.ci_ATE',est), get(paste0('hat_ATE',est)) - 1.96*sqrt(mean(get(paste0('hat_EIF.ATE',est))^2)/n)) # lower CI

          assign(paste0('upper.ci_ATE',est), get(paste0('hat_ATE',est)) + 1.96*sqrt(mean(get(paste0('hat_EIF.ATE',est))^2)/n)) # upper CI


          # TMLE and onestep output
          # ATE
          assign(paste0(m,est),list(ATE=get(paste0('hat_ATE',est)), # estimated parameter
                                        lower.ci=get(paste0('lower.ci_ATE',est)), # lower bound of 95% CI
                                        upper.ci=get(paste0('upper.ci_ATE',est)), # upper bound of 95% CI
                                        EIF=get(paste0('hat_EIF.ATE',est)))) # EIF

          output[[count]] <- get(paste0(m,est))
          names(output)[count] <- paste0(m,est)
          count <- count + 1

          indicator <- which(level.z==est)

          # print estimates
          cat(paste0(m," estimated ACE for z=",level.z.number[indicator],": ",round(get(paste0(m,est))$ATE,2),"; 95% CI: (",round(get(paste0(m,est))$lower.ci,2),", ",round(get(paste0(m,est))$upper.ci,2),") \n"))




        }

    }

    return(output)

  } # end of ATE estimate




 ################################## E(Y(x))


  if (length(x)==1) { ## E(Y^1) estimate ==

    out.a <- napkin.a(x, z =z, data, treatment, Z.variables, W.variables, outcome, covariates,
                      z.density=z.density,z.method = z.method, superlearner.Y = superlearner.Y, superlearner.X = superlearner.X, superlearner.Z = superlearner.Z,
                      crossfit = crossfit, K = K,
                      lib.Y = lib.Y,
                      lib.X = lib.X,
                      lib.Z = lib.Z,
                      n.iter = n.iter, cvg.criteria = cvg.criteria,
                      formula.Y = formula.Y, formula.X = formula.X, formula.Z = formula.Z,
                      linkY_binary = linkY_binary, link.X = link.X, link.Z = link.Z,
                      truncate_lower.X = truncate_lower.X, truncate_upper.X = truncate_upper.X,
                      truncate_lower.Z = truncate_lower.Z, truncate_upper.Z = truncate_upper.Z, minZ=minZ,maxZ=maxZ,verbose=verbose,fast=fast)



    # run TMLE
    TMLE_output <- out.a$TMLE


    # run onestep
    Onestep_output <- out.a$Onestep

    if(univariate.binaryZ){

      # run estimating equation
      EstEquation_output <- out.a$EstEquation

      }


    # levels of z used for estimation
    level.z <- sub("^out", "", names(TMLE_output))


    level.z.number <- sub("^out\\.z", "", names(TMLE_output))
    level.z.number[level.z.number=="out.all.z"] <- "all z"


    if (univariate.binaryZ){output <- vector("list", 3*length(level.z))}else{output <- vector("list", 2*length(level.z))}
    if (univariate.binaryZ){estimators <- c('TMLE','Onestep','EstEquation')}else{estimators <- c('TMLE','Onestep')}

    # count of method
    count <- 1

    for (m in estimators){ # loop over TMLE and onestep

      for (est in level.z){ # loop over levels of z

        output[[count]] <- get(paste0(m,'_output'))[[paste0('out',est)]] # get output for Y(x)
        names(output)[count] <- paste0(m,est)
        count <- count + 1

        indicator <- which(level.z==est)
        # print estimates
        cat(paste0(m," estimated E(Y(",x,")) for z=",level.z.number[indicator],": ",round(get(paste0(m,'_output'))[[paste0('out',est)]]$estimated,2),";95% CI: (",round(get(paste0(m,'_output'))[[paste0('out',est)]]$lower.ci,2),",",round(get(paste0(m,'_output'))[[paste0('out',est)]]$upper.ci,2),") \n"))




      }

    }

    return(output)

  } # end of E(Y(x)) estimate





} # end of TMLE function


