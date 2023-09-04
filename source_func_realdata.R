########################################################
# Main functions for real data analysis
########################################################

library(MASS)
library(ncvreg)
library(glmnet)

get_all_methods_mcb <- function(data_input, r, p, confidence.level, method = 'adalasso', threshold=0.05){
  
  if(method == "adalasso"){
    boot_res <- boot_var_adalasso(x = data_input, r = r, p = p, confidence.level = confidence.level)
    opt.cr <- boot_res$freq
    opt.width <- boot_res$width
    var.set <- boot_res$var_adalasso_set
    samp.set <- boot_res$samp_adalasso_set
  }
  else if(method == "scad"){
    boot_res <- boot_var_scad(x = data_input, r = r, p = p, pnlt="SCAD", confidence.level = confidence.level)
    opt.cr <- boot_res$freq
    opt.width <- boot_res$width
    var.set <- boot_res$var_set
    samp.set <- boot_res$samp_set
  }
  else if(method == "mcp"){
    boot_res <- boot_var_scad(x = data_input, r = r, p = p, pnlt="MCP", confidence.level = confidence.level)
    opt.cr <- boot_res$freq
    opt.width <- boot_res$width
    var.set <- boot_res$var_set
    samp.set <- boot_res$samp_set
  }
  else if(method == "lasso"){
    boot_res <- boot_var_lasso(x = data_input, r = r, p = p, threshold = threshold, confidence.level = confidence.level)
    opt.cr <- boot_res$freq
    opt.width <- boot_res$width
    var.set <- boot_res$var_lasso_set
    samp.set <- boot_res$samp_lasso_set
  }
  else if(method == "scaledlasso"){
    boot_res <- bcr_tcr_res_scalreg(x = data_input, r = r, p = p, confidence.level = confidence.level)
    opt.cr <- boot_res$freq
    opt.width <- boot_res$width
    var.set <- boot_res$var_scal_set
    samp.set <- boot_res$samp_scal_set
  }
  
  return(list(boot_res=boot_res, 
              opt.cr=opt.cr,
              opt.width=opt.width,
              var.set=var.set, 
              samp.set=samp.set))
  
}

bcr_tcr_res_scalreg <- function(x, r, p, confidence.level) 
{
  var_scal <- matrix(0,nrow=r,ncol=p)
  var_scal_set <- vector(mode = 'list', length = r)
  fit <- scalreg(as.matrix(x[,1:p]), as.vector(x[,(p+1)]), lam0 = 'univ')
  res_original <- fit$residuals
  res_after_center <- res_original - mean(res_original)
  constant <- fit$fitted.values
  
  for(j in 1:r) { 
    ind <- sample(1:nrow(x),nrow(x),replace=T)
    new_response <- constant + res_after_center[ind]
    boot.data <- cbind(x[,1:p], new_response)
    tem <- scalreg(as.matrix(boot.data[,1:p]), as.vector(boot.data[,(p+1)]), lam0 = 'univ')
    var_scal[j,] <- full.var%in%full.var[which(tem$coefficients!=0)]
    var_scal_set[[j]] <- full.var[which(tem$coefficients!=0)]
  }
  
  CI_res <- get_opt_mcb_binary_search(var_scal, confidence.level, p, r)
  opt.lower <- CI_res$optimal.mcbs.lower
  opt.upper <- CI_res$optimal.mcbs.upper
  opt.width <- CI_res$optimal.width
  opt.cr <- CI_res$optimal.cr
  
  results <- list(var_scal_set=var_scal_set,
                  samp_scal_set=full.var[which(fit$coefficients!=0)],
                  freq=opt.cr,
                  lower=opt.lower,
                  upper=opt.upper,
                  width=opt.width)
  return(results)
}

boot_var_adalasso <- function(x, r, p, confidence.level)
{
  var_adalasso <- matrix(0,nrow=r,ncol=p)
  var_adalasso_set <- vector(mode = 'list', length = r)
  adalasso_res <- adalasso(X=as.matrix(x[,1:p]),y=x[,(p+1)],k=10)
  res_original <- x[,(p+1)] - as.matrix(x[,1:p]) %*% adalasso_res$coefficients.adalasso - adalasso_res$intercept.adalasso
  res_after_center <- res_original - mean(res_original)
  constant <- as.matrix(x[,1:p]) %*% adalasso_res$coefficients.adalasso + adalasso_res$intercept.adalasso
  for(j in 1:r) { 
    ind <- sample(1:nrow(x),nrow(x),replace=T)
    new_response <- constant + res_after_center[ind]
    boot.data <- cbind(x[,1:p], new_response)
    adalasso_tem <- adalasso(X=as.matrix(boot.data[,1:p]),y=boot.data[,(p+1)],k=10)
    var_adalasso[j,] <- full.var%in%full.var[which(adalasso_tem$coefficients.adalasso!=0)]
    var_adalasso_set[[j]] <- full.var[which(adalasso_tem$coefficients.adalasso!=0)]
  }
  
  CI_res <- get_opt_mcb_binary_search(var_adalasso, confidence.level, p, r)
  opt.lower <- CI_res$optimal.mcbs.lower
  opt.upper <- CI_res$optimal.mcbs.upper
  opt.width <- CI_res$optimal.width
  opt.cr <- CI_res$optimal.cr
  
  results <- list(var_adalasso_set=var_adalasso_set,
                  samp_adalasso_set=full.var[which(adalasso_res$coefficients.adalasso!=0)],
                  freq=opt.cr,
                  lower=opt.lower,
                  upper=opt.upper,
                  width=opt.width)
  return(results)
}


boot_var_scad <- function(x, r, p, pnlt="SCAD", confidence.level)
{
  var_scad <- matrix(0,nrow=r,ncol=p)
  var_scad_set <- vector(mode = 'list', length = r)
  fit <- cv.ncvreg(X=as.matrix(x[,1:p]),y=x[,(p+1)],penalty=pnlt)
  fit2 <- fit$fit
  beta <- fit2$beta[,fit$min]
  res_original <- x[,(p+1)] - as.matrix(x[,1:p]) %*% beta[2:(p+1)] - beta[1]
  res_after_center <- res_original - mean(res_original)
  constant <- as.matrix(x[,1:p]) %*% beta[2:(p+1)] + beta[1]
  for(j in 1:r) { 
    ind=sample(1:nrow(x),nrow(x),replace=T)
    new_response <- constant + res_after_center[ind]
    boot.data <- cbind(x[,1:p], new_response)
    tem <- cv.ncvreg(X=as.matrix(boot.data[,1:p]),y=boot.data[,(p+1)],penalty=pnlt)
    fit2_tem<-tem$fit
    beta_tem<-fit2_tem$beta[,tem$min]
    beta_tem <- beta_tem[2:(p+1)]
    var_scad[j,]<-full.var%in%names(beta_tem)[which(beta_tem!=0)]
    var_scad_set[[j]]<-names(beta_tem)[which(beta_tem!=0)]
  }
  
  CI_res <- get_opt_mcb_binary_search(var_scad, confidence.level, p, r)
  opt.lower <- CI_res$optimal.mcbs.lower
  opt.upper <- CI_res$optimal.mcbs.upper
  opt.width <- CI_res$optimal.width
  opt.cr <- CI_res$optimal.cr
  
  results <- list(var_set=var_scad_set,
                  samp_set=names(beta[2:(p+1)])[which(beta[2:(p+1)]!=0)],
                  freq=opt.cr,
                  lower=opt.lower,
                  upper=opt.upper,
                  width=opt.width)
  return(results)
}

boot_var_lasso <- function(x, r, p, threshold, confidence.level)
{
  var_lasso <- matrix(0,nrow=r,ncol=p)
  var_lasso_set <- vector(mode = 'list', length = r)
  opt.lambda <- cv.glmnet(x=as.matrix(x[,1:p]),y=x[,(p+1)],alpha=1)$lambda.min
  lasso.fit <- glmnet(x=as.matrix(x[,1:p]),y=x[,(p+1)],family='gaussian',alpha=1)
  beta <- coef(lasso.fit,s=opt.lambda)[,1]; beta[2:(p+1)] <- ifelse(beta[2:(p+1)] > threshold, beta[2:(p+1)], 0)
  res_original <- x[,(p+1)] - as.matrix(x[,1:p]) %*% beta[2:(p+1)] - beta[1]
  res_after_center <- res_original - mean(res_original)
  constant <- as.matrix(x[,1:p]) %*% beta[2:(p+1)] + beta[1]
  
  for (i in 1:r){
    ind=sample(1:nrow(x),nrow(x),replace=T)
    new_response <- constant + res_after_center[ind]
    boot.data <- cbind(x[,1:p], new_response)
    opt.lambda.boot<-cv.glmnet(x=as.matrix(boot.data[,1:p]),y=boot.data[,(p+1)],alpha=1)$lambda.min
    lasso.fit.boot<-glmnet(x=as.matrix(boot.data[,1:p]),y=boot.data[,(p+1)],family='gaussian',alpha=1)
    beta.boot <- coef(lasso.fit.boot,s=opt.lambda.boot)[-1,1][which(abs(coef(lasso.fit.boot,s=opt.lambda.boot)[-1,1]) > threshold)]
    var_lasso[i,]<-full.var%in%names(beta.boot)
    var_lasso_set[[i]] <- names(beta.boot)
  }
  
  CI_res <- get_opt_mcb_binary_search(var_lasso, confidence.level, p, r)
  opt.lower <- CI_res$optimal.mcbs.lower
  opt.upper <- CI_res$optimal.mcbs.upper
  opt.width <- CI_res$optimal.width
  opt.cr <- CI_res$optimal.cr

  results <- list(var_lasso_set=var_lasso_set,
                  samp_lasso_set=names(beta[2:(p+1)])[which(beta[2:(p+1)]!=0)],
                  freq=opt.cr,
                  lower=opt.lower,
                  upper=opt.upper,
                  width=opt.width)
  return(results)
}


alter_model <- function(datasets, p, method='scaledlasso', threshold=0.05) 
{
  if(method=='scaledlasso'){
    tem <- scalreg(as.matrix(datasets[,1:p]), as.vector(datasets[,(p+1)]), lam0 = 'univ')$coefficients
    alter_model_est <- names(tem)[which(tem!=0)]
  }else if(method=='lasso'){
    opt.lambda<-cv.glmnet(x=as.matrix(datasets[,-(p+1)]),y=datasets[,(p+1)],alpha=1)$lambda.min
    lasso.fit<-glmnet(x=as.matrix(datasets[,-(p+1)]),y=datasets[,(p+1)],family='gaussian',alpha=1)
    tem <- coef(lasso.fit,s=opt.lambda)[,1]; tem[2:(p+1)] <- ifelse(tem[2:(p+1)] > threshold, tem[2:(p+1)], 0)
    alter_model_est <- names(tem[2:(p+1)])[which(tem[2:(p+1)]!=0)]
  }else if(method=='adalasso'){
    tem<-adalasso(X=as.matrix(datasets[,1:p]),y=datasets[,(p+1)],k=10)$coefficients.adalasso
    alter_model_est <- colnames(datasets)[1:p][which(tem!=0)]
  }else if(method=='scad'){
    fit<-cv.ncvreg(X=as.matrix(datasets[,1:p]),y=datasets[,(p+1)],penalty="SCAD")
    tem<-fit$fit$beta[,fit$min]
    alter_model_est <- names(tem[2:(p+1)])[which(tem[2:(p+1)]!=0)]
  }else if(method=='mcp'){
    fit<-cv.ncvreg(X=as.matrix(datasets[,1:p]),y=datasets[,(p+1)],penalty="MCP")
    tem<-fit$fit$beta[,fit$min]
    alter_model_est <- names(tem[2:(p+1)])[which(tem[2:(p+1)]!=0)]
  }else if(method=='screening'){
    pvalue.screen<-c()
    datasets <- as.data.frame(datasets)
    for(i in 1:p){
      screen.olm<-lm(as.numeric(datasets[,(p+1)])~datasets[,i],data=datasets)
      pvalue.screen[i]<-summary(screen.olm)$coefficients[2,4]
    }
    alter_model_est <- colnames(datasets)[1:p][which(pvalue.screen<=0.05)]
    data.new <- datasets[,alter_model_est]
    olm.new <- lm(as.numeric(data.new[,ncol(data.new)])~.,data=data.new)
    tem <- olm.new$coefficients
  }
  return(list(alter_model_est=alter_model_est,
              beta.est=tem))
}
