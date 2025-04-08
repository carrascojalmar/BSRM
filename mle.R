#------------------------------------------------------------------------------#
rm(list=ls())
#------------------------------------------------------------------------------#
# Packages
library(gamlss)
library(rmutil)
#------------------------------------------------------------------------------#
# Functions
# Bivariate Simplex Regression Model
dsimplex_fgm <- function(y1,y2,x,b01,b11,b02,b12,g01,g11,g02,g12,
                         lambda,log=FALSE){
  
  mu1 <- exp(b01+b11*x)/(1+exp(b01+b11*x))
  sigma1 <- exp(g01+g11*x)
  
  mu2 <- exp(b02+b12*x)/(1+exp(b02+b12*x))
  sigma2 <- exp(g02+g12*x)
  
  dsim1 <- dsimplex(y1,mu1, sigma1)
  dsim2 <- dsimplex(y2, mu2, sigma2)
  
  u <- pSIMPLEX(y1, mu1, sigma1)
  v <- pSIMPLEX(y2, mu2, sigma2)
    
  cop_fgm <- 1+lambda*(1-2*u)*(1-2*v)
  dsim_f <- dsim1*dsim2*cop_fgm
  
  if(log==FALSE){return(dsim_f)}else{
    return(log(dsim_f))
  }
  
}

# Likelihood function

logLikFun_sfgm <- function(param,y1,y2,x){
  b01 <- param[1]
  b11 <- param[2]
  b02 <- param[3]
  b12 <- param[4]
  g01 <- param[5]
  g11 <- param[6]
  g02 <- param[7]
  g12 <- param[8]
  lambda <- param[9]
  func <- dsimplex_fgm(y1,y2,x,b01,b11,b02,b12,g01,g11,g02,g12,
                       lambda,log=TRUE)
  lfunc <- sum(func)
  return(lfunc)
}    

#------------------------------------------------------------------------------#
# data
setwd("/Users/jalmarcarrasco/Library/CloudStorage/OneDrive-Pessoal/Jalmar/Artigos/Papers_Desenvolvimento/The Simplex bivariate regression model/codes-data-MRSB")
df.alagoas <- read.table("alagoas.txt", header = T)
head(df.alagoas)
summary(df.alagoas)
attach(df.alagoas)
#------------------------------------------------------------------------------#
y1 <- idhm       
y2 <- ivs     
x <- scale(razdep,center = TRUE, scale = TRUE)
#------------------------------------------------------------------------------#
# Optimization

chute <- c(b01 = .2, b11 = .2, b02 = .2, b12 = .2, 
           g01 = .2, g11 = .2, g02 = .2, 
           g12 = .2, lambda = -.5)

op <- optim(par=chute,logLikFun_sfgm,method="L-BFGS-B",y1=y1,y2=y2,x=x, 
      lower=c(rep(-Inf,8),-1),upper=c(rep(Inf,8),1), 
      hessian=T,control=list(maxit=10000,fnscale=-1,trace=1))

ep <- sqrt(diag(-solve(op$hessian)))
zvalue <- op$par/ep
pvalue <- round(2*(1-pt(abs(zvalue), length(x)-length(op$par))),3)
results <- cbind(op$par,ep,zvalue,pvalue)
colnames(results) <- c("Estimation","Std. error","z-Stat.","p-Value")
rownames(results) <- c("beta01","beta11","beta02","beta12","gamma01",
                       "gamma11","gamma02","gamma12","lambda")
print(results, digits = 4)
#------------------------------------------------------------------------------#
