#' @title Calculate the density of Z given its W and covariates
#' @description This function calculates the density associated with Z: p(Z|W).
#' @param data A Dataframe contains all the variables listed in vertices parameter
#' @param Z.variables A character vector indicating the variables, the level of which are fixed in the estimation of the counterfactual outcome.
#' @param W.variables A character vector indicating the parents of Z.variables.
#' @param formula.Z A regression formula specifying the relationship between Z,W, and covariates. For example formula = "Z ~ W".
#' @param superlearner.Z A logical value indicating whether to use super learner to estimate the propensity score. Default is FALSE.
#' @param lib.Z A character vector indicating the library of the super learner. Default is c("SL.glm","SL.earth","SL.ranger","SL.mean").
#' @param crossfit A logical value indicating whether to use crossfitting. Default is FALSE.
#' @param K An integer indicating the number of folds in crossfitting. Default is 5.
#' @keywords density estimation
#' @return A vector of density for each row of the data.
#' @export
#' @importFrom mvtnorm dmvnorm
#' @examples
#' z.density <- calculate_density_dnorm(Z.variables = "Z", W.variables = "W",
#' data = data_Zcontinuous_Ycontinuous, formula = "Z ~ W")
#' head(z.density)
#'
calculate_density_dnorm <- function(Z.variables, W.variables, data, formula.Z="Z~.", superlearner.Z=F,lib.Z=c("SL.glm","SL.earth","SL.ranger","SL.mean"), crossfit=F, K=5){ # A is a vector, M and X are data frame

  # This function only allow M to be either univariate binary / continuous
  # or multivariate continuous.
  # It does not allow M to be a multivariate variable and has binary elements.

  n <- nrow(data)

  # Variables
  W <- data[,W.variables, drop = F]
  Z <- data[,Z.variables]


  # new data sets
  dat_mpZ = data.frame(W)
  dat_ZmpZ = data.frame(Z,W)


  if (crossfit==T){ #### cross fitting + super learner #####

    z_fit <- CV.SuperLearner(Y=Z, X=dat_mpZ, family = gaussian(), V = K, SL.library = lib.Z, control = list(saveFitLibrary=T),saveAll = T)

    ## make prediction for p(Z_i|W_i)
    # model prediction
    predict.Z <- z_fit$SL.predict


  } else if (superlearner.Z==T){ #### super learner #####

    z_fit <- SuperLearner(Y=Z, X=dat_mpZ, family = gaussian(), SL.library = lib.Z)

    ## make prediction for p(Z_i|W_i)
    # model prediction
    predict.Z <- z_fit$SL.predict


  }else { #### simple linear regression with user input regression formula: default="A ~ ." ####

    z_fit <- glm(as.formula(formula.Z), data=dat_mpZ,  family = gaussian())

    ## make prediction for p(Z_i|W_i)
    # model prediction
    predict.Z <- predict(z_fit)


  } # end of if-else for propensity score fitting

    # Store model errors
    model_errors <- Z - predict.Z

    z.density <- dnorm(Z, mean = predict.Z, sd = sd(model_errors))

    return(list(z.density))

}

