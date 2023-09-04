########################################################
# R code for Table 4 in the main text
########################################################

source('source_func.R')
library(MASS)
library(glmnet)
library(ncvreg)
library(smoothmest)

#############################################################################
#                            main functions of CVC                          #
#############################################################################

CleanErrMat <- function(err.mat) {
  # This is a subroutine called in the main function "CVPerm"
  # It removes identical columns in the cross-validated loss matrix
  n2 <- nrow(err.mat)
  M <- ncol(err.mat)
  err.mean <- apply(err.mat, 2, mean)
  err.mat.center <- err.mat - matrix(err.mean, nrow = n2, ncol = M, byrow = T)
  
  res <- list()
  list.length <- 0
  remaining.index <- 1:M
  
  while (length(remaining.index) > 0) {
    m <- remaining.index[1]
    err.diff.center <- err.mat.center[,m] - err.mat.center
    err.mean.diff <- err.mean[m] - err.mean
    sd.vec <- apply(err.diff.center, 2, sd)
    J.1 <- which(sd.vec==0 & err.mean.diff > 0)
    J.2 <- which(sd.vec==0 & err.mean.diff == 0)
    J.3 <- which(sd.vec==0 & err.mean.diff < 0)
    if (length(J.1) == 0) {
      list.length <- list.length+1
      res[[list.length]] <- J.2
    }
    remaining.index <- setdiff(remaining.index, union(J.2,J.3))
  }
  return(list(err.mat.c = err.mat[,sapply(res,function(x){x[1]})], ind=res))
}

CVPerm <- function(err.mat, B = 200, screen=T, alpha.s=0.005) {
  # This is the main function of cross-validation with confidence.
  # Requires subroutine "CleanErrMat"
  # Input: "err.mat" is a matrix where each row corresponds to a data point
  #            and each column corresponds to a candidate model.
  #            err.mat[i,j] records the cross-validated loss evaluated at
  #            the i'th data point and the j'th candidate model.
  #            For example, in least square linear regression this is the
  #            cross-validated squared residual 
  #        "B" is the number of bootstrap samples
  #        "screen" is an indicator if pre-screening is used for the test
  #        "alpha.s" is the threshold used in pre-screening
  clean.res <- CleanErrMat(err.mat)
  err.mat.c <- clean.res$err.mat.c
  err.mat.c.ind <- clean.res$ind
  n2 <- nrow(err.mat.c)
  M <- ncol(err.mat.c)
  err.mean <- apply(err.mat.c, 2, mean)
  err.mat.center <- err.mat.c - matrix(err.mean, nrow = n2, ncol = M, byrow = T)
  sign.mat <- matrix(rnorm(n2*B),ncol = B)
  norm.quantile <- qnorm(1-alpha.s/(M-1))
  screen.th <- ifelse(norm.quantile^2>=n2,Inf,norm.quantile / sqrt(1-norm.quantile^2/n2))
  sgmb.p.val <- sapply(1:M, function(m){
    err.diff.center <- err.mat.center[,m] - err.mat.center[, setdiff(1:M, m),drop=F]
    err.mean.diff <- err.mean[m] - err.mean[setdiff(1:M, m)]
    sd.vec <- apply(err.diff.center, 2, sd)
    err.mean.diff.scale <- err.mean.diff / sd.vec
    test.stat <- sqrt(n2)*max(err.mean.diff.scale)
    if (test.stat >= screen.th) {
      return(alpha.s)
    }
    if (screen) {        
      J.screen <- which(sqrt(n2) * err.mean.diff.scale > -2*screen.th)
      if (length(J.screen)==0) {return(1-alpha.s)}
    } else {
      J.screen <- 1:(M-1)
    }
    err.diff.center.scale <- err.diff.center[,J.screen] / 
      matrix(sd.vec[J.screen], nrow = n2, ncol = length(J.screen), byrow=T)
    test.stat.vec <- sapply(1:B,
                            function(ib){max(apply(err.diff.center.scale*sign.mat[,ib],2,mean))})
    return(mean(test.stat.vec>max(err.mean.diff.scale[J.screen])))
  })
  res <- rep(0,ncol(err.mat))
  for (i in 1:length(err.mat.c.ind)) {
    res[err.mat.c.ind[[i]]] <- sgmb.p.val[i]
  }
  return(res)
}

CVC <- function(X, Y, type = "stepwise", max.steps=NULL, lambda=NULL, nlambda=NULL, ols.models=NULL, n.fold = 5, B = 200, screen=T) {
  # This is the main wrapper function that applies CVC to various linear regression methods
  # Input: "X", "Y" are the covariate and response, in matrix form
  #        "type": "OLS", "stepwise", "lasso", "lars" 
  #            if type=="OLS", then a list of candidate models "ols.models" shall be provided,
  #            see "low_dim_example.R" for example
  #        "max.steps" is the maximum number of steps used in lars package, which is used to implement forward stepwise
  #        "lambda" is the vector of lambda values used in lasso, implemented by "glmnet" package
  #        "nlambda": instead of providing a sequence of lambda, one can also just provide a total number of lambda values,
  #            the sequence of lambda will then be chosen by glmnet.
  #        "ols.models": the list of candidate ols models when type=="OLS", see "low_dim_example.R" for example
  #        "n.fold" is the number of folds
  #        "B", "screen" are parameters passed to the main CVC function "CVPerm"
  n <- nrow(X)
  p <- ncol(X)
  fold.size <- floor(n/n.fold)
  n <- n.fold * fold.size
  fold.ind <- matrix(sample.int(n),nrow = n.fold)
  M <- length(ols.models)
  cv.res <- list(test.err=matrix(0,ncol=M,nrow=0), beta=array(0,c(p,M,n.fold)))
  for (i.fold in 1:n.fold) {
    test.ind <- fold.ind[i.fold,]
    train.ind <- setdiff(1:n, test.ind)
    X.train <- X[train.ind,]
    Y.train <- Y[train.ind]
    X.test <- X[test.ind,]
    Y.test <- Y[test.ind]
    beta.est <- sapply(ols.models,
                       function(x){est.coef <- rep(0,p)
                       if (length(x)>0) {
                         est.coef[x] <- lm(Y.train ~ X.train[,x]-1)$coef
                         est.coef[is.na(est.coef)] <- 0
                       }
                       return(est.coef)})
    predict.vals <- X.test %*% beta.est
    test.err <- predict.vals - Y.test
    cv.res$test.err <- rbind(cv.res$test.err, test.err)
    cv.res$beta[,,i.fold] <- beta.est
  }
  
  test.err.1 <- cv.res$test.err[1:fold.size,]
  p.vals.1 <- CVPerm(test.err.1^2, B = B)
  test.err.mean.1 <- apply(test.err.1, 2, mean)
  p.vals.c <- CVPerm(cv.res$test.err^2, B = B)
  return(list(test.err.c = cv.res$test.err, p.vals.c = p.vals.c, test.err.1 = test.err.1, p.vals.1 = p.vals.1,
              beta.1 = cv.res$beta[,,1], beta = cv.res$beta))
}

GenBeta <- function(s=3,p=200,mix=F) {
  b <- c(rep(1,s),rep(0,p-s))
  if (mix) {b[(s+1):(2*s)]<-c(0.7,0.5,0.3)}
  return(b)
}

GenData <- function(b,n,p,rho,err.sigma,standardize=F) {
  Sigma <- matrix(rho, ncol = p, nrow = p) + diag(rep(1-rho, p))
  X <- mvrnorm(n, rep(0, p), Sigma)
  if (standardize) {
    marg.s <- sqrt(apply(X^2, 2, mean))
    marg.s.inv <- 1 / marg.s
    X <- X %*% diag(marg.s.inv)
  }
  Y <- X %*% b + rnorm(n,0,err.sigma)
  return(list(X=X,Y=Y))
}

type <- "OLS" 
n.cand <- c(40, 80, 160, 320, 640)
p <- 5 
n.fold <- 5 
err.sigma.cand <- c(1,4)
b <- c(2,9,0,4,8) 
alpha.cand <- c(0.05,0.1)
ols.models <- vector(mode = 'list', length = 2^p)
s0 <- matrix(c(1,rep(0,p-1)),1,p)
s1 <- matrix(c(rep(0,p)),1,p)
for (i in 2:2^p){
  s1 <- s0 + s1
  for (j in 1:p){
    if (s1[1,j] == 2){
      s1[1,j+1] <- s1[1,j+1]+1
      s1[1,j] <- 0
    }
  }           
  ols.models[[i]] <- which(s1!=0)
}

# setting up the simulation
n.rep <- 100
r <- 200
full.var <- paste("x",1:p,sep="")
rho <- 0
sim.res <- array(0, dim=c(length(alpha.cand),length(err.sigma.cand),length(n.cand),n.rep,6))

for(h0 in 1:length(alpha.cand)){
  alpha <- alpha.cand[h0]
  for(h1 in 1:length(err.sigma.cand)){
    err.sigma <- err.sigma.cand[h1]
    for(h2 in 1:length(n.cand)){
      n <- n.cand[h2]
      for (i.rep in 1:n.rep) {
        
        train.data <- GenData(b,n,p,rho,err.sigma)
        X <- train.data$X
        Y <- train.data$Y
        
        ## MCB
        time.mcb <- system.time({
          var_mcp <- matrix(0,nrow=r,ncol=p)
          fit <- cv.ncvreg(X=X,y=Y,penalty='MCP')
          fit2 <- fit$fit
          beta <- fit2$beta[,fit$min]
          res_original <- Y - X %*% beta[2:(p+1)] - beta[1]
          res_after_center <- res_original - mean(res_original)
          constant <- X %*% beta[2:(p+1)] + beta[1]
          for(j in 1:r) { 
            ind=sample(1:nrow(X),nrow(X),replace=T)
            new_response <- constant + res_after_center[ind]
            tem <- cv.ncvreg(X=X,y=new_response,penalty='MCP')
            fit2_tem<-tem$fit
            beta_tem<-fit2_tem$beta[,tem$min]
            beta_tem <- beta_tem[2:(p+1)]
            var_mcp[j,]<-full.var%in%full.var[which(beta_tem!=0)]
          }
          CI_res <- get_opt_mcb_binary_search(var.matrix = var_mcp, confidence.level=1-alpha, p=p, r=r)
          opt.lower <- CI_res$optimal.mcbs.lower
          opt.upper <- CI_res$optimal.mcbs.upper
          opt.width <- CI_res$optimal.width
          opt.cr <- CI_res$optimal.cr
          truecapture <- cap_ture(low_list=opt.lower,upper_list=opt.upper,true.model=abs(sign(b)))
        })[3]
        
        ## CVC
        time.cvc <- system.time({
          cvc.res <- CVC(X,Y,type=type,ols.models=ols.models)
        })[3]
        k.cvc.vec <- which(cvc.res$p.vals.c > alpha)
        n.sel.cvc <- length(k.cvc.vec)
        cover.cnt <- ifelse(28%in%k.cvc.vec,1,0)
        sim.res[h0,h1,h2,i.rep,] <- c(time.mcb, time.cvc, truecapture, cover.cnt, opt.width, n.sel.cvc)
      }
      cat('---------n = ', n.cand[h2], '---done!!!---------\n')
    }
    cat('---------err.sigma = ', err.sigma.cand[h1], '---done!!!---------\n')
  }
  cat('---------alpha = ', alpha.cand[h0], '---done!!!---------\n')
}


file_name <- paste("sim_cvc_mcb", "p", p, "rho", rho, "rep", n.rep, sep = "_")
save.image(file = paste("F:/xnhu/MoST/",file_name, ".RData", sep = ""))
