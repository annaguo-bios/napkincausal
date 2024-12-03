#' @title Calculate the density of Z given its W and covariates
#' @description This function calculates the density associated with Z: p(Z|W,C).
#' @param data A Dataframe contains all the variables listed in vertices parameter
#' @param Z.variables A character vector indicating the variables, the level of which are fixed in the estimation of the counterfactual outcome.
#' @param W.variables A character vector indicating the parents of Z.variables.
#' @param covariates A character vector indicating the pre-treatment covariates in the napkin graph.
#' @param formula.Z A regression formula specifying the relationship between Z,W, and covariates. For example formula = "Z ~ W + C".
#' @param superlearner.Z A logical value indicating whether to use super learner to estimate the propensity score. Default is FALSE.
#' @param lib.Z A character vector indicating the library of the super learner. Default is c("SL.glm","SL.earth","SL.ranger","SL.mean").
#' @param crossfit A logical value indicating whether to use crossfitting. Default is FALSE.
#' @param K An integer indicating the number of folds in crossfitting. Default is 5.
#' @keywords density estimation
#' @return A vector of density for each row of the data.
#' @export
#' @importFrom mvtnorm dmvnorm
#' @examples
#' z.density <- calculate_density_dnorm(Z.variables = "Z", W.variables = "W", covariates = "C",
#' data = data_Zcontinuous_Ycontinuous, formula = "Z ~ W + C")
#' head(z.density)
#'
calculate_density_dnorm <- function(Z.variables, W.variables, covariates, data, formula.Z="Z~.", superlearner.Z=F,lib.Z=c("SL.glm","SL.earth","SL.ranger","SL.mean"), crossfit=F, K=5){ # A is a vector, M and X are data frame

  # This function only allow M to be either univariate binary / continuous
  # or multivariate continuous.
  # It does not allow M to be a multivariate variable and has binary elements.

  n <- nrow(data)

  # Variables
  C <- data[,covariates, drop = F]
  W <- data[,W.variables, drop = F]
  Z <- data[,Z.variables]


  # new data sets
  dat_mpZ = data.frame(W, C)
  dat_ZmpZ = data.frame(Z,W, C)

  # for which variables, user specified formula for regression
  #
  #
  # if (length(Z.variables)>1){ ## if M is multivariate ##
  #
  #   num_columns <- length(Z.variables) # number of variables under vertex Z
  #
  #   for (var in Z.variables){
  #
  #     if (all(data[, var] %in% c(0,1))) {
  #
  #       stop("This function doesn't support multivariate variables with binary elements.")
  #
  #     }
  #
  #   } # end of for loop over variables
  #
  #   ## Fit regression model to each component of M ##
  #
  #   # model prediction
  #   predict.Z <- list() # store the prediction results for each variable in Z
  #
  #   # currently only support linear models
  #
  #   # Initialize vectors to store errors
  #   model_errors <- matrix(NA, nrow = n, ncol = num_columns)
  #
  #   for (i in 1:num_columns){ # loop over each variable in Z
  #
  #     if (Z.variables[i] %in% formula.variables){ # if the user specified formula for this variable
  #
  #       model <- lm(as.formula(formula[[Z.variables[i]]]), data=data[, c(Z.variables[i],W.variables,covariates)])
  #
  #     }else{ # if the user didn't specify formula for this variable
  #
  #       model <- lm(data[, Z.variables[i]] ~ . , data=data[, c(W.variables,covariates)])
  #
  #     }
  #
  #     model_errors[, i] <- data[, Z.variables[i]] - predict(model) # errors of the model
  #
  #     predict.Z[[Z.variables[i]]] <- predict(model) # store the prediction result
  #
  #   } # end of for loop over variables
  #
  #   # compute the variance-covariance matrix of the errors
  #   varcov <- cov(data.frame(model_errors))
  #
  #   # Define a function for the ratio calculation
  #   f.z.density <- function(j) {
  #
  #     mean.pred <- sapply(predict.Z, `[`, j) # mean of the normal distribution
  #
  #     dmvnorm(
  #       x = data[j, Z.variables],
  #       mean = mean.pred, # coeff*mp
  #       sigma = varcov
  #     )
  #   }
  #
  #   # Apply the function to each row of M
  #   z.density <- sapply(1:n, f.z.density)
  #
  #
  #
  #
  # }else{ ## if M is univariate ##

    Z <- data[, Z.variables]


    if (crossfit==T){ #### cross fitting + super learner #####

      z_fit <- CV.SuperLearner(Y=Z, X=dat_mpZ, family = gaussian(), V = K, SL.library = lib.Z, control = list(saveFitLibrary=T),saveAll = T)

      ## make prediction for p(Z_i|W_i,C_i)
      # model prediction
      predict.Z <- z_fit$SL.predict


    } else if (superlearner.Z==T){ #### super learner #####

      z_fit <- SuperLearner(Y=Z, X=dat_mpZ, family = gaussian(), SL.library = lib.Z)

      ## make prediction for p(Z_i|W_i,C_i)
      # model prediction
      predict.Z <- z_fit$SL.predict


    }else { #### simple linear regression with user input regression formula: default="A ~ ." ####

      z_fit <- glm(as.formula(formula.Z), data=dat_mpZ,  family = gaussian())

      ## make prediction for p(Z_i|W_i,C_i)
      # model prediction
      predict.Z <- predict(z_fit)


    } # end of if-else for propensity score fitting

    # Store model errors
    model_errors <- Z - predict.Z

    z.density <- dnorm(Z, mean = predict.Z, sd = sd(model_errors))

    # make prediction for p(Z_i|W_j,C_j) for both i and j from 1 to n
    p.z.matrix <- matrix(NA, n, n)

    for (i in 1:n) {

      p.z.matrix[, i] <- dnorm(Z[i], mean = predict.Z, sd = sd(model_errors))

    }

    # ## make prediction for p(Z_i|W_i,C_i)
    # # model prediction
    # predict.Z <- predict(model) # store the prediction result
    #
    # # Store model errors
    # model_errors <- Z - predict.Z
    #
    # z.density <- dnorm(Z, mean = predict.Z, sd = sd(model_errors))
    #
    # # make prediction for p(Z_i|W_j,C_j) for both i and j from 1 to n
    # p.z.matrix <- matrix(NA, n, n)
    #
    # for (i in 1:n) {
    #
    #   p.z.matrix[, i] <- dnorm(Z[i], mean = predict.Z, sd = sd(model_errors))
    #
    # }


    # if (all(Z %in% c(0,1))) { # binary variable
    #
    #   # For binary columns, use glm
    #
    #   if (Z.variables %in% formula.variables){ # if the user specified formula for this variable
    #
    #     model <- glm(formula[[Z.variables]], data=dat_ZmpZ)
    #
    #   }else{
    #
    #     model <- glm(Z ~ . , data=dat_mpZ)
    #
    #   }
    #
    #   ## make prediction for p(Z_i|W_i,C_i)
    #   # calculate the density ratio
    #   z1.density <- predict(model, type="response") # p(Z=1|W,C)
    #
    #   z.density <- Z*(z1.density) + (1-Z)*(1-z1.density) # p(Z|W,C)
    #
    #   ## make prediction for p(Z_i|W_j,C_j) for both i and j from 1 to n
    #   p.z.matrix <- matrix(NA, n, n)
    #
    #   for (i in 1:n) {
    #
    #     p.z.matrix[, i] <- Z[i]*(z1.density) + (1-Z[i])*(1-z1.density)
    #
    #   }
    #
    # } else { # continuous variable
    #
    #
    #
    # } # end of if-else

  # }

  return(list(z.density,p.z.matrix))

}

