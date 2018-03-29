
source("Rcode/rankAFT.R")

# Use STISI method to impute the covariate Z (subject to censoring)
# where the outcome Y is binary with the substantive model being the logistic regression
STISI.yLogis <- function(data.Orig,data.Boot,Y,Y.cov,Z,V,Delta.Z,aft.cov,imp,TOL,trans){
  # data.Orig: the original dataset with Z censored 
  # data.Boot: one bootstrapped dataset (from data.Orig)
  # Y: primary outcome of interest, binary (1/0)
  # Y.cov: fully observed covariate(s) that are used in the substantive model of Y
  # Z: covariate (used in the substantive model of Y) subject to censoring
  # V=G(min(Z,C)): observed data (after transformation) of Z 
  # where G(.) is a (pre-specified) monotone transformation function, T=G(Z) 
  # C: censoring variable (in LOD situation, C is the DL constant)
  # user need to apply G(.) to obtain V ahead of time
  # Delta.Z: censoring indicator for Z (1=observed, 0=censored)
  # imp: two options for the imputed residual time
  # imp="obs" impute with an observed residual time; 
  # imp="unif": impute uniformly between a (small) time interval  
  # TOL: a very small positive constant to control for the imputation bound
  # trans: transformation function G(.)
  # Four options are provided for trans: "iden", "neg.iden", "log", "neg.log"
  
  # Step 1: Obtain bootstrap data (already done, input of the function)
  # Step 2: use complete bootstrap data to fit the substantive model and obtain the regression parameter estimates
  
  form <- "data.Boot[,Y]~"
  for(k in 1:length(Y.cov)) form <- paste0(form,"data.Boot[,Y.cov[",k,"]]+")
  Z.cov="Z"
  formula.cc <- as.formula(paste0(form,"data.Boot[,Z.cov]"))
  fit.cc <- glm(formula.cc,family=binomial,subset=(data.Boot[,Delta.Z]==1))
  beta.cc <- summary(fit.cc)$coeff[,1] # include the intercept
  se.beta.cc <- summary(fit.cc)$coeff[,2]

  # Step 3: use bootstrap model to fit AFT model V~aft.cov to get F.em and aft.coef
  F.est.Boot <- est.F(y.aft=data.Boot[,V],x.aft=data.Boot[,aft.cov],delta.aft=data.Boot[,Delta.Z])
  aft.coef <- F.est.Boot$aft.coef
  F.em <- F.est.Boot$F.em
  
  # Multiple Imputation Starts 
  data.Orig.cen <- data.Orig[data.Orig[,Delta.Z]==0,]
  n.cen <- dim(data.Orig.cen)[1]
  T.imp <- rep(NA,n.cen) 
  count.draw <- rep(0,n.cen)
  
  for (i in 1:n.cen){
    # For a given covariate cov_i, generate F(t|zeta>C-cov_i*slope)
    time.lbd <- data.Orig.cen[i,V]-as.vector(as.numeric(data.Orig.cen[i,aft.cov])%*%aft.coef[-1])
    id <- sum(F.em$time<=time.lbd)
    if(id<dim(F.em)[1]) {
      # Prob is the estimated probabiliby Pr(zeta<=C-cov_i*slope)
      if(id>0) Prob <- F.em$eF[id] 
      if(id==0) Prob <- 0
      
      time.cond <- F.em$time[F.em$time>time.lbd]
      eF.new <- F.em$eF[F.em$time>time.lbd]
      eF.cond <- (eF.new-Prob)/(1-Prob)
      
      # F.cond.em is the empirical estimate of the conditional CDF F(t|zeta>C-cov_i*slope)
      F.cond.em <- as.data.frame(rbind(c(time.lbd,0),cbind(time.cond,eF.cond)))
      colnames(F.cond.em) <- c("time","eF.cond")
      
      # Step 2b: draw a value of T from the conditional empirical distribution
      # Step 3: compare with the bound to decide whether accept the draw value
      U2 <- 1 # initial value of U2 
      bdd <- 0 # initial value of bound
      # count = 0
      while(U2>bdd){
        U1 <- runif(n=1) # U1~Unif(0,1), used for draw the T value
        # find the maximum time (t) that has F.em(t)<=U1. maxT.id >=1 (since F.em starts from 0)
        maxT.id <- sum(F.cond.em$eF.cond<=U1)
        
        # find the range of the maxT (since F.em is a step function, an interval of t corresponds to the same F.em(t))
        if(maxT.id<dim(F.cond.em)[1]){
          maxT.min <- F.cond.em$time[maxT.id]
          maxT.max <- F.cond.em$time[maxT.id+1] 
          res.time <- runif(n=1,min=maxT.min,max=maxT.max) # impute the time within this interval
          #res.time <- F.cond.em$time[maxT.id] # impute the time only using the observed failure time
          T.draw <- res.time + as.vector(as.numeric(data.Orig.cen[i,aft.cov])%*%aft.coef[-1])
          
          if(trans == "neg.log") Z.draw <- exp(-T.draw)
          if(trans == "neg.iden") Z.draw <- -T.draw
          
          if(trans == "log") Z.draw <- exp(T.draw)
          if(trans == "iden") Z.draw <- T.draw
          data.Orig.cen[i,Z] <- Z.draw  
          
          # compute the bound for binary case
          bdd <- exp(c(1,as.numeric(data.Orig.cen[i,c(Y.cov,Z)])) %*% beta.cc)^data.Orig.cen[i,Y]/
            (1+exp(c(1,as.numeric(data.Orig.cen[i,c(Y.cov,Z)])) %*% beta.cc))
          
          if(bdd>TOL) U2 <- runif(n=1) # U2~Unif(0,1), used for comparison with the bound
          if(bdd<=TOL) {U2=0; T.draw=NA} # force to break out the while loop}
        }
    
        if(maxT.id==dim(F.cond.em)[1]) { 
          # this can only happen if eF.cond does not go to 1
          # which means the unconditional eF does not go to 1, 
          # i.e., the last observed time in the residual scale is a censored value
          # we don't impute in this case
          T.draw = NA # don't impute 
          bdd <- 1
          U2 <- 0 # force the inequality 'U2>bdd' do not hold and break out the while loop
        }
        
        count.draw[i] = count.draw[i] + 1
        
      }
      T.imp[i] = T.draw
    }     
        
    if(id==dim(F.em)[1]) {
      # this means "C-cov_i*slope" is beyond the last time point in F.em
      T.imp[i] = NA # don't impute       
    }    
  }
  
  if(trans=="neg.iden") Z.imp = -T.imp
  if(trans=="neg.log") Z.imp = exp(-T.imp)
  if(trans=="iden") Z.imp = T.imp
  if(trans=="log") Z.imp = exp(T.imp)
  
  data.Orig$Z.imp <- data.Orig[,Z]
  data.Orig$Z.imp[data.Orig[,Delta.Z]==0] <- Z.imp
  
  return(list(data.imp=data.Orig, count.draw=count.draw))
}