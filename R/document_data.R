#' This is dataset with binary Z and binary outcome Y generated based on the Napkin graph
#'
#'
#' @format A data frame with 2000 rows and 5 variables: C, W, Z, X, Y.
#' \describe{
#'  \item{C}{the pre-treatment variable that follows uniform distribution in the range of 0.1 to 1}
#'  \item{W}{a variable that follows from bernoulli distribution with probability of 0.5}
#'  \item{Z}{a binary variable. Z|W,C follows bernoulli distribution with probability p(Z=1|W,C) = expit(-1 + W + C)}
#'  \item{X}{the treatment variable. X|Z,W follows bernoulli distribution with probability p(X=1|Z,W,C) = (1 -W + Z*W)/4}
#'  \item{Y}{the outcome variable. Y|X,Z,W,C follows bernoulli distribution with probability p(Y=1|X,Z,W,C) = (4+X+Z-ZW - W+(4/3)(1-W)(1-X)(1-Z)+C)/8}
#'  }

#' @name data_Zbinary_Ybinary
#' @docType data
#' @keywords napkin graph data
"data_Zbinary_Ybinary"


#' This is dataset with binary Z and continuous outcome Y generated based on the Napkin graph
#'
#' @format A data frame with 2000 rows and 5 variables: C, W, Z, X, Y.
#' \describe{
#' \item{C}{the pre-treatment variable that follows uniform distribution in the range of 0.1 to 1}
#' \item{W}{a variable that follows from bernoulli distribution with probability of 0.5}
#' \item{Z}{a binary variable. Z|W,C follows bernoulli distribution with probability p(Z=1|W,C) = expit(-1 + W + C)}
#' \item{X}{the treatment variable. X|Z,W follows bernoulli distribution with probability p(X=1|Z,W,C) = (1 -W + Z*W)/4}
#' \item{Y}{the outcome variable. Y|X,Z,W,C follows normal distribution with mean (4+X+Z-ZW - W+(4/3)(1-W)(1-X)(1-Z)+C) and standard deviation 1}
#' }

#' @name data_Zbinary_Ycontinuous
#' @docType data
#' @keywords napkin graph data
"data_Zbinary_Ycontinuous"


#' This is dataset with continuous Z and binary outcome Y generated based on the Napkin graph
#'
#' @format A data frame with 2000 rows and 5 variables: C, W, Z, X, Y.
#' \describe{
#' \item{C}{the pre-treatment variable that follows uniform distribution in the range of 0.1 to 1}
#' \item{W}{a variable that follows from bernoulli distribution with probability of 0.5}
#' \item{Z}{a continuous variable. Z|W,C follows uniform distribution in the range of 0.1 to (1+W+C)/4}
#' \item{X}{the treatment variable. X|Z,W follows bernoulli distribution with probability p(X=1|Z,W,C) = (1 -W + Z*W)/4}
#' \item{Y}{the outcome variable. Y|X,Z,W,C follows bernoulli distribution with probability p(Y=1|X,Z,W,C) = (4+X+Z-ZW - W+(4/3)(1-W)(1-X)(1-Z)+C)/8}
#' }

#' @name data_Zcontinuous_Ybinary
#' @docType data
#' @keywords napkin graph data
"data_Zcontinuous_Ybinary"


#' This is dataset with continuous Z and continuous outcome Y generated based on the Napkin graph
#'
#' @format A data frame with 2000 rows and 5 variables: C, W, Z, X, Y.
#' \describe{
#' \item{C}{the pre-treatment variable that follows uniform distribution in the range of 0.1 to 1}
#' \item{W}{a variable that follows from bernoulli distribution with probability of 0.5}
#' \item{Z}{a continuous variable. Z|W,C follows uniform distribution in the range of 0.1 to (1+W+C)/4}
#' \item{X}{the treatment variable. X|Z,W follows bernoulli distribution with probability p(X=1|Z,W,C) = (1 -W + Z*W)/4}
#' \item{Y}{the outcome variable. Y|X,Z,W,C follows normal distribution with mean (4+X+Z-ZW - W+(4/3)(1-W)(1-X)(1-Z)+C) and standard deviation 1}
#' }

#' @name data_Zcontinuous_Ycontinuous
#' @docType data
#' @keywords napkin graph data
"data_Zcontinuous_Ycontinuous"
