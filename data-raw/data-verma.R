set.seed(7)

generate_data <- function(n, 
                          Z.binary=Z.binary,
                          Y.binary=Y.binary,
                          # when Z is binary
                          parZ.binary = c(-1,1,1), 
                          # when Z is continuous
                          parZ.continuous = c(1,1,1)/4,
                          parX=c(1,-1,1)/4, 
                          # when Y continuous
                          parY.continuous=c(4,1,1,-1,-1,4/3,1),
                          # when Y binary
                          parY.binary=c(4,1,1,-1,-1,4/3,1)/8
){
  
  
  if (Z.binary){parZ <- parZ.binary} else {parZ <- parZ.continuous}
  
  if (Y.binary){parY <- parY.binary} else {parY <- parY.continuous}
  
  
  C <- runif(n,0.1,1) # p(C)
  
  # W <- rbinom(n, 1, C/(1+C)) # p(W|C)
  W <- rbinom(n,1,0.5)
  
  # p(Z|W,C)
  if (Z.binary){ Z <- rbinom(n, 1, plogis(parZ[1]+parZ[2]*W+parZ[3]*C))} else {Z <- runif(n,0.1,(parZ[1]+parZ[2]*W+parZ[3]*C)) } 
  
  X <- rbinom(n, 1, parX[1] + parX[2]*W + parX[3]*Z*W) # p(X|Z, W) it's linear link!
  
  # p(Y|Z,W,X,C)
  # if (Y.binary){ Y <- rbinom(n, 1, parY[1]+parY[2]*X+parY[3]*Z*C+parY[4]*Z*W*C+parY[5]*W+parY[6]*C*(1-W)*(1-X)*(1-Z))} else {Y <- parY[1]+parY[2]*X+parY[3]*Z*C+parY[4]*Z*W*C+parY[5]*W+parY[6]*C*(1-W)*(1-X)*(1-Z) +runif(n,0,1)} 
  if (Y.binary){ Y <- rbinom(n, 1, parY[1]+parY[2]*X+parY[3]*Z+parY[4]*Z*W +parY[5]*W+parY[6]*(1-W)*(1-X)*(1-Z)+parY[7]*C)} else {Y <- parY[1]+parY[2]*X+parY[3]*Z+parY[4]*Z*W+parY[5]*W + parY[6]*(1-W)*(1-X)*(1-Z) + parY[7]*C +runif(n,0,1)} 
  data <- data.frame(C=C, W=W, Z=Z, X=X, Y=Y)
  
  # propensity score
  ps <- X*(parX[1] + parX[2]*W + parX[3]*Z*W)+
    (1-X)*(1-(parX[1] + parX[2]*W + parX[3]*Z*W))
  
  return(list(data = data,
              Z.binary=Z.binary,
              parZ = parZ, 
              parX=parX, 
              parY=parY,
              Y.binary=Y.binary,
              ps=ps))
}


data_Zbinary_Ybinary <- generate_data(2000, Z.binary = T, Y.binary = T)$data
data_Zbinary_Ycontinuous <- generate_data(2000, Z.binary = T, Y.binary = F)$data
data_Zcontinuous_Ybinary <- generate_data(2000, Z.binary = F, Y.binary = T)$data
data_Zcontinuous_Ycontinuous <- generate_data(2000, Z.binary = F, Y.binary = F)$data

usethis::use_data(data_Zbinary_Ybinary, overwrite = T)
usethis::use_data(data_Zbinary_Ycontinuous, overwrite = T)
usethis::use_data(data_Zcontinuous_Ybinary, overwrite = T)
usethis::use_data(data_Zcontinuous_Ycontinuous, overwrite = T)
