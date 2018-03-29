
library(quantreg)
library(SparseM)
library(survival)

#rankaft
rankAFT=function (x, y, delta) 
{
  ynew <- 1000 * (length(y))^2
  data1 <- data.frame(y, x)
  options(contrasts = c("contr.treatment", "contr.poly"))
  tempfit <- lm(y ~ ., x = TRUE, y = TRUE, data = data1)
  x <- as.matrix(tempfit$x[, -1])
  xn <- dimnames(x)[[2]]
  y <- tempfit$y
  dimnum <- dim(x)
  n1 <- dimnum[1]
  n2 <- dimnum[2]
  residuals <- matrix(0, nrow = n1, ncol = 3 + 1)
  yy0 <- rep(y, rep(n1, n1))
  delta1 <- rep(delta, rep(n1, n1))
  yy1 <- rep(y, n1)
  yy <- delta1 * (yy0 - yy1)
  xx0 <- matrix(rep(as.vector(x), rep(n1, n1 * n2)), nrow = n1 * n1)
  xx1 <- t(matrix(rep(as.vector(t(x)), n1), nrow = n2))
  xx <- xx0 - xx1
  xxdif <- xx * delta1
  xnew <- apply(xxdif, 2, sum)
  xnew <- rbind(xxdif, -xnew)
  yynew <- c(yy, ynew)
  dimnames(xnew) <- list(NULL, xn)
  fit <- rq(yynew ~ xnew - 1, method = "fn")
  betag <- fit$coef
  residn <- fit$resid
  residn <- (!(residn > 0))
  residn <- residn[-(length(residn))]
  betal <- betag
  for (i in 1:3) {
    fitted <- x %*% betal
    eb <- y - fitted
    ss0b <- (n1 + 1 - rank(eb))/n1
    ss0b1 <- rep(ss0b, rep(n1, n1))
    xxdifl <- xxdif/ss0b1
    xnewl <- apply(xxdifl, 2, sum)
    xnewl <- rbind(xxdifl, -xnewl)
    yyl <- c(yy/ss0b1, ynew)
    fitl <- rq(yyl ~ xnewl - 1, method = "fn")
    betal <- fitl$coef
  }
  predmatrix <- x - t(matrix(rep(apply(x, 2, mean), n1), ncol = n1))
  residuals[, 1] <- y - predmatrix %*% as.matrix(betag)
  residuals[, 2] <- y - predmatrix %*% as.matrix(betal)
  object <- list(beta = rbind(betag, betal), residuals = residuals, 
                 message = fit$message)
  class(object) <- "AFT"
  object
}

# est.F returns the estimated/empirical (unconditional) CDF of zeta
# where zeta=T-X*alpha is the residual survival time from the AFT model
est.F <- function(y.aft,x.aft,delta.aft){
  
  # y.aft: input of y outcome (observed time in the AFT model)
  # x.aft: input of X covariates in the AFT model
  # delta.aft: input of the censoring indicator in the AFT model
  
  beta.rank <- rankAFT(x=x.aft, y=y.aft, delta=delta.aft)$beta # slope estimators from "rankaft"
  beta.GH <- beta.rank[1,] # Gehan-weighted
  beta.LR <- beta.rank[2,] # logrank weighted
  beta <- beta.GH # use Gehan-weighted
  flag=0  # to flag whether the beta.GH and beta.LR differ a lot or not
  if (sum(abs(beta.GH-beta.LR))>=1) { flag=1 }
  # calculate the residual time at the estimated beta
  res <- y.aft-as.matrix(x.aft)%*%beta
  # K-M estimate of the survival function of the error term (at the estimated beta)
  KM.fit <- survfit(Surv(res,delta.aft)~1)
  # KM.fit$time is ordered
  time <- KM.fit$time
  surv <- KM.fit$surv # estimated survival function for the residual error
  q <- length(time)
  margin <- surv[q]
  jump <- c(1,surv[-q])-surv # drop the last one if it is not a failure (in that case, surv[q-1]=surv[q])
  #jump <- c(1,surv[-q])-c(surv[-q],0) #force last one to be failure if it is not a failure
  alpha <- time%*%jump # intercept estimator of the AFT model
  
  rst.1 <- cbind(time,jump,1-surv)
  rst.2 <- rst.1[rst.1[,2]>0,] # only keep these rows with jump>0
  rst.3 <- rbind(c(min(time),0,0),rst.2)
  rst.3 <- as.data.frame(rst.3) # rst.3 is ordered by time (small to large)
  colnames(rst.3) <- c("time","jump","eF")
  
  return(list(aft.coef=c(alpha,beta),F.em=rst.3,flag=flag))
  
  # F.em is the empirical "unconditional" CDF for the residual time zeta (zeta=T-aft.cov*slope with intercept absorbed)
  # for imputing these censored T, "conditional" CDF F(t|zeta>C-aft.cov*slope) is needed
  # F(t|zeta>C-aft.cov*slope) = 0 if zeta<=C-aft.cov*slope; = (F(t)-F(C-aft.cov*slope))/Pr(zeta>C-aft.cov*slope) if zeta>C-aft.cov*slope
  # where F(t) is the unconditional CDF for zeta
  # note that the empirical conditional CDF depends on aft.cov
  # empical conditional CDF is calculated in a separate function below
}
