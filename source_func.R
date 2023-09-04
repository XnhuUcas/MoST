########################################################
# Main functions
########################################################

library(MASS)
library(ncvreg)
library(glmnet)

Generate_Dataset <- function(size, p, sig.num, sd.epsi, rho, decay.factor = NULL, sig.val = 1, var.dis = 'normal', sig.val.def = NULL) # generating datasets
{
  omiga<-matrix(0,nrow=p,ncol=p)
  for(i in 1:p)
  {
    for(j in 1:p)
    {
      omiga[i,j]<-rho^(abs(i-j))
    }
  }
  data<-mvrnorm(n=size,mu=rep(0,p),omiga)
  colnames(data)=paste("x",1:p,sep="")
  
  if (var.dis == 't') u<-rt(size, df)
  if (var.dis == 'normal') u<-rnorm(size,mean=0,sd = sd.epsi)
  if (var.dis == 'cauchy') u<-rcauchy(size,location=0,scale=1)
  if (var.dis == 'Laplase') u<-rdoublex(size,mu=0,lambda =sqrt(sd.epsi^2/2))
  
  if(is.null(sig.val.def)){
    if(is.null(decay.factor)){
      b<-matrix(sig.val,nrow=1,ncol=sig.num)
    }else{
      b<-matrix(NA,nrow=1,ncol=sig.num)
      if(sig.num == 1){
        b[1] = sig.val
      }else{
        b[1] = sig.val
        for(i in 1:(sig.num-1)){
          b[,i+1] = b[,i]*decay.factor
        }
      }
    }
  }else{
    b <- matrix(sig.val.def,1,sig.num)
  }
  
  y<-t(b%*%t(data[,1:sig.num])+u)
  colnames(y)=c('y')
  return(data.frame(data,y))
}


########### PI-MCB algorithm for computing MCBs candidate
mcb_candidate <- function(var.matrix, p, r) 
{
  #--------------------
  # var.list: selected variables set
  # var.matrix: 0/1 matrix for selected variables
  #-------------------

  colsum<-apply(var.matrix,2,sum)
  order<-order(colsum,decreasing = T)
  freq<-vector(length=p+1)
  lower<-matrix(0,nrow=(p+1),ncol=p)
  upper<-matrix(0,nrow=(p+1),ncol=p)
  oper_times<-0
  
  freq[1]<-0
  for(i in 0:p)
  { 
    cap<-vector(length=p-i+1);cap[1]<-0 
    for(j in 0:(p-i))
    {
      if(j==0 && i==0){
        lowtest<-rep(0,p)
        upptest<-rep(0,p)
      }else{
        if(j==0 && i!=0){
          lowtest<-rep(0,p)
          upptest<-rep(0,p)
          upptest[order[1:i]]<-1
        }else{
          lowtest<-rep(0,p)
          upptest<-rep(0,p)
          lowtest[order[1:j]]<-1
          upptest[order[1:(i+j)]]<-1
        }
      }
      for(m in 1:r){
        if(all(all(lowtest<=var.matrix[m,]),all(var.matrix[m,]<=upptest))) cap[j+1]<-cap[j+1]+1
      }
    }
    freq[i+1]<-max(cap)/r
    maxlocation<-which.max(cap)
    oper_times<-oper_times+(p-i+1)
    
    if(i==0 && maxlocation==1)
    {
      lower[i+1,]<-rep(0,p)
      upper[i+1,]<-rep(0,p)
    }else{
      if(i!=0 & maxlocation==1){
        lower[i+1,]<-rep(0,p)
        upper[i+1,order[1:i]]<-1
      }else{
        lower[i+1,order[1:(maxlocation-1)]]<-1
        upper[i+1,order[1:(maxlocation-1+i)]]<-1
      }
    }

  }   

  result<-list(freq=freq,lower=lower,upper=upper, oper_times=oper_times)
  return(result)
}


get_opt_mcb_pi<-function(result, confidence.level){
  fit.index<-which(result$freq >= confidence.level)
  optimal.mcbs.lower<-result$lower[min(fit.index),]
  optimal.mcbs.upper<-result$upper[min(fit.index),]
  optimal.width<-min(fit.index)-1
  optimal.cr<-result$freq[min(fit.index)]
  result<-list(width=optimal.width, 
               mcb.lower=optimal.mcbs.lower, 
               mcb.upper=optimal.mcbs.upper,
               mcb.bcr=optimal.cr,
               oper_times=result$oper_times)
  return(result)
}


######### BIS-MCB for computing optimal MCB
get_opt_mcb_binary_search <- function(var.matrix, confidence.level, p, r) 
{

  colsum<-apply(var.matrix,2,sum)
  order<-order(colsum,decreasing = T)
  freq<-c()
  maxlocation<-c()
  iter_width<-c()
  i<-1
  width.lower<-c();width.lower[1]<-0 
  width.iter<-c();width.iter[1]<-floor(p/2) 
  width.upper<-c();width.upper[1]<-p 
  oper_times<-0
  
  while((width.upper[i]-width.lower[i])>0){
    cap<-vector(length=p-width.iter[i]+1);cap[1]<-0  
    for(j in 0:(p-width.iter[i]))
    {
      if(j==0 && width.iter[i]==0){
        lowtest<-rep(0,p)
        upptest<-rep(0,p)
      }else{
        if(j==0 && width.iter[i]!=0){
          lowtest<-rep(0,p)
          upptest<-rep(0,p)
          upptest[order[1:width.iter[i]]]<-1
        }else{
          lowtest<-rep(0,p)
          upptest<-rep(0,p)
          lowtest[order[1:j]]<-1
          upptest[order[1:(width.iter[i]+j)]]<-1
        }
      }
      for(m in 1:r){
        if(all(all(lowtest<=var.matrix[m,]),all(var.matrix[m,]<=upptest))) cap[j+1]<-cap[j+1]+1
      }
    }

    freq[i]<-max(cap)/r
    maxlocation[i]<-which.max(cap)
    iter_width[i]<-width.iter[i]
    oper_times<-oper_times+(p-width.iter[i]+1)

    if((width.upper[i]-width.lower[i])==1)
    { 
      if(freq[i]>=confidence.level)
      {
        if(width.iter[i]==0 && maxlocation[i]==1)
        {
          opt_lower<-rep(0,p)
          opt_upper<-rep(0,p)
        }else{
          if(maxlocation[i]==1 && width.iter[i]!=0){
            opt_lower<-rep(0,p)
            opt_upper<-rep(0,p)
            opt_upper[order[1:width.iter[i]]]<-1
          }else{
            opt_lower<-rep(0,p)
            opt_upper<-rep(0,p)
            opt_lower[order[1:(maxlocation[i]-1)]]<-1
            opt_upper[order[1:(maxlocation[i]-1+width.iter[i])]]<-1
          }
        }
        return(list(optimal.width=width.iter[i], 
                    optimal.mcbs.lower=opt_lower, 
                    optimal.mcbs.upper=opt_upper,
                    optimal.cr=freq[i],
                    oper_times=oper_times))
      }else{ 
        if(width.upper[i]==p){
          return(list(optimal.width=p, 
                      optimal.mcbs.lower=rep(0,p), 
                      optimal.mcbs.upper=rep(1,p),
                      optimal.cr=1,
                      oper_times=oper_times))
          }
          else{
            maxlocation = maxlocation[which(iter_width==width.upper[i])]
            if(width.upper[i]==0 && maxlocation==1)
            {
              opt_lower<-rep(0,p)
              opt_upper<-rep(0,p)
            }else{
              if(maxlocation==1 && width.upper[i]!=0){
                opt_lower<-rep(0,p)
                opt_upper<-rep(0,p)
                opt_upper[order[1:width.upper[i]]]<-1
              }else{
                opt_lower<-rep(0,p)
                opt_upper<-rep(0,p)
                opt_lower[order[1:(maxlocation-1)]]<-1
                opt_upper[order[1:(maxlocation-1+width.upper[i])]]<-1
              }
            }
            return(list(optimal.width=width.upper[i], 
                        optimal.mcbs.lower=opt_lower, 
                        optimal.mcbs.upper=opt_upper,
                        optimal.cr=freq[which(iter_width==width.upper[i])],
                        oper_times=oper_times))
          }
      }
    }else{
      if(freq[i]>=confidence.level){
        width.lower[i+1]<-width.lower[i]
        width.upper[i+1]<-width.iter[i]
        width.iter[i+1]<-floor((width.lower[i+1]+width.upper[i+1])/2)
      }else{
        width.lower[i+1]<-width.iter[i]
        width.upper[i+1]<-width.upper[i]
        width.iter[i+1]<-floor((width.lower[i+1]+width.upper[i+1])/2)
      }
    }

    i<-i+1
  }
}


cap_ture<-function(low_list,upper_list,true.model){
  if(all(all(low_list<=true.model),all(true.model<=upper_list))) truecapture<-1 else truecapture<-0  
  return(truecapture)
}


boot_var_adalasso <- function(x, r, p, confidence.level, true.model)
{
  var_adalasso <- matrix(0,nrow=r,ncol=p)
  var_adalasso_set <- vector(mode = 'list', length = r)
  adalasso_res <- adalasso(X=as.matrix(x[,1:p]),y=x$y,k=10)
  res_original <- x$y - as.matrix(x[,1:p]) %*% adalasso_res$coefficients.adalasso - adalasso_res$intercept.adalasso
  res_after_center <- res_original - mean(res_original)
  constant <- as.matrix(x[,1:p]) %*% adalasso_res$coefficients.adalasso + adalasso_res$intercept.adalasso
  for(j in 1:r) { 
    ind <- sample(1:nrow(x),nrow(x),replace=T)
    new_response <- constant + res_after_center[ind]
    boot.data <- cbind(x[,1:p], new_response)
    colnames(boot.data)[p+1] <- "y"
    adalasso_tem <- adalasso(X=as.matrix(boot.data[,1:p]),y=boot.data$y,k=10)
    var_adalasso[j,] <- full.var%in%full.var[which(adalasso_tem$coefficients.adalasso!=0)]
    var_adalasso_set[[j]] <- full.var[which(adalasso_tem$coefficients.adalasso!=0)]
  }

  CI_res <- get_opt_mcb_binary_search(var_adalasso, confidence.level, p, r)
  opt.lower <- CI_res$optimal.mcbs.lower
  opt.upper <- CI_res$optimal.mcbs.upper
  opt.cr <- CI_res$optimal.cr
  opt.width <- CI_res$optimal.width
  truecapture <- cap_ture(opt.lower,opt.upper,true.model)
  
  results <- list(var_adalasso_set=var_adalasso_set,
                  samp_adalasso_set=full.var[which(adalasso_res$coefficients.adalasso!=0)],
                  freq=opt.cr,
                  width=opt.width,
                  truecapture=truecapture,
                  lower=opt.lower,
                  upper=opt.upper)
  return(results)
}


boot_var_scad <- function(x, r, p, pnlt="SCAD", confidence.level, true.model)
{
  var_scad <- matrix(0,nrow=r,ncol=p)
  var_scad_set <- vector(mode = 'list', length = r)
  fit <- cv.ncvreg(X=as.matrix(x[,1:p]),y=x$y,penalty=pnlt)
  fit2 <- fit$fit
  beta <- fit2$beta[,fit$min]
  res_original <- x$y - as.matrix(x[,1:p]) %*% beta[2:(p+1)] - beta[1]
  res_after_center <- res_original - mean(res_original)
  constant <- as.matrix(x[,1:p]) %*% beta[2:(p+1)] + beta[1]
  for(j in 1:r) { 
    ind=sample(1:nrow(x),nrow(x),replace=T)
    new_response <- constant + res_after_center[ind]
    boot.data <- cbind(x[,1:p], new_response)
    tem <- cv.ncvreg(X=as.matrix(boot.data[,1:p]),y=boot.data$new_response,penalty=pnlt)
    fit2_tem<-tem$fit
    beta_tem<-fit2_tem$beta[,tem$min]
    beta_tem <- beta_tem[2:(p+1)]
    var_scad[j,]<-full.var%in%names(beta_tem)[which(beta_tem!=0)]
    var_scad_set[[j]]<-names(beta_tem)[which(beta_tem!=0)]
  }
  
  CI_res <- get_opt_mcb_binary_search(var_scad, confidence.level, p, r)
  opt.lower <- CI_res$optimal.mcbs.lower
  opt.upper <- CI_res$optimal.mcbs.upper
  opt.cr <- CI_res$optimal.cr
  opt.width <- CI_res$optimal.width
  truecapture <- cap_ture(opt.lower,opt.upper,true.model)
  
  results <- list(var_set=var_scad_set,
                samp_set=names(beta[2:(p+1)])[which(beta[2:(p+1)]!=0)],
                freq=opt.cr,
                width=opt.width,
                truecapture=truecapture,
                lower=opt.lower,
                upper=opt.upper)
  return(results)
}


boot_var_lasso <- function(x, r, p, threshold, confidence.level, true.model)
{
  var_lasso <- matrix(0,nrow=r,ncol=p)
  var_lasso_set <- vector(mode = 'list', length = r)
  opt.lambda <- cv.glmnet(x=as.matrix(x[,1:p]),y=x$y,alpha=1)$lambda.min
  lasso.fit <- glmnet(x=as.matrix(x[,1:p]),y=x$y,family='gaussian',alpha=1)
  beta <- coef(lasso.fit,s=opt.lambda)[,1]; beta[2:(p+1)] <- ifelse(beta[2:(p+1)] > threshold, beta[2:(p+1)], 0)
  res_original <- x$y - as.matrix(x[,1:p]) %*% beta[2:(p+1)] - beta[1]
  res_after_center <- res_original - mean(res_original)
  constant <- as.matrix(x[,1:p]) %*% beta[2:(p+1)] + beta[1]

  for (i in 1:r){
    ind=sample(1:nrow(x),nrow(x),replace=T)
    new_response <- constant + res_after_center[ind]
    boot.data <- cbind(x[,1:p], new_response)
    colnames(boot.data)[p+1] <- "y"
    opt.lambda.boot<-cv.glmnet(x=as.matrix(boot.data[,1:p]),y=boot.data$y,alpha=1)$lambda.min
    lasso.fit.boot<-glmnet(x=as.matrix(boot.data[,1:p]),y=boot.data$y,family='gaussian',alpha=1)
    beta.boot <- coef(lasso.fit.boot,s=opt.lambda.boot)[-1,1][which(abs(coef(lasso.fit.boot,s=opt.lambda.boot)[-1,1]) > threshold)]
    var_lasso[i,]<-full.var%in%names(beta.boot)
    var_lasso_set[[i]] <- names(beta.boot)
  }

  CI_res <- get_opt_mcb_binary_search(var_lasso, confidence.level, p, r)
  opt.lower <- CI_res$optimal.mcbs.lower
  opt.upper <- CI_res$optimal.mcbs.upper
  opt.cr <- CI_res$optimal.cr
  opt.width <- CI_res$optimal.width
  truecapture <- cap_ture(opt.lower,opt.upper,true.model)
  
  results <- list(var_lasso_set=var_lasso_set,
                  samp_lasso_set=names(beta[2:(p+1)])[which(beta[2:(p+1)]!=0)],
                  freq=opt.cr,
                  width=opt.width,
                  truecapture=truecapture,
                  lower=opt.lower,
                  upper=opt.upper)
  return(results)
}


bcr_tcr_res_scalreg <- function(x, r, p, confidence.level, true.model) 
{
  var_scal <- matrix(0,nrow=r,ncol=p)
  var_scal_set <- vector(mode = 'list', length = r)
  fit <- scalreg(as.matrix(x[,1:p]), as.vector(x$y), lam0 = 'univ')
  res_original <- fit$residuals
  res_after_center <- res_original - mean(res_original)
  constant <- fit$fitted.values
  
  for(j in 1:r) { 
    ind <- sample(1:nrow(x),nrow(x),replace=T)
    new_response <- constant + res_after_center[ind]
    boot.data <- cbind(x[,1:p], new_response)
    tem <- scalreg(as.matrix(boot.data[,1:p]), as.vector(boot.data$new_response), lam0 = 'univ')
    var_scal[j,] <- full.var%in%full.var[which(tem$coefficients!=0)]
    var_scal_set[[j]] <- full.var[which(tem$coefficients!=0)]
  }
  
  CI_res <- get_opt_mcb_binary_search(var_scal, confidence.level, p, r)
  opt.lower <- CI_res$optimal.mcbs.lower
  opt.upper <- CI_res$optimal.mcbs.upper
  opt.cr <- CI_res$optimal.cr
  opt.width <- CI_res$optimal.width
  truecapture <- cap_ture(opt.lower,opt.upper,true.model)
  
  results <- list(var_scal_set=var_scal_set,
                samp_scal_set=full.var[which(fit$coefficients!=0)],
                freq=opt.cr,
                width=opt.width,
                truecapture=truecapture,
                lower=opt.lower,
                upper=opt.upper)
  return(results)
}


get_all_methods_mcb <- function(data_input, r, p, confidence.level, method = 'adalasso', true.model, threshold=0.05){
  
  if(method == "adalasso"){
    boot_res <- boot_var_adalasso(x = data_input, r = r, p = p, confidence.level = confidence.level, true.model = true.model)
    opt.cr <- boot_res$freq
    capture <- boot_res$truecapture
    var.set <- boot_res$var_adalasso_set
    samp.set <- boot_res$samp_adalasso_set
  }
  else if(method == "scad"){
    boot_res <- boot_var_scad(x = data_input, r = r, p = p, pnlt="SCAD", confidence.level = confidence.level, true.model = true.model)
    opt.cr <- boot_res$freq
    capture <- boot_res$truecapture
    var.set <- boot_res$var_set
    samp.set <- boot_res$samp_set
  }
  else if(method == "mcp"){
    boot_res <- boot_var_scad(x = data_input, r = r, p = p, pnlt="MCP", confidence.level = confidence.level, true.model = true.model)
    opt.cr <- boot_res$freq
    capture <- boot_res$truecapture
    var.set <- boot_res$var_set
    samp.set <- boot_res$samp_set
  }
  else if(method == "lasso"){
    boot_res <- boot_var_lasso(x = data_input, r = r, p = p, threshold = threshold, confidence.level = confidence.level, true.model = true.model)
    opt.cr <- boot_res$freq
    capture <- boot_res$truecapture
    var.set <- boot_res$var_lasso_set
    samp.set <- boot_res$samp_lasso_set
  }
  else if(method == "scaledlasso"){
    boot_res <- bcr_tcr_res_scalreg(x = data_input, r = r, p = p, confidence.level = confidence.level, true.model = true.model)
    opt.cr <- boot_res$freq
    capture <- boot_res$truecapture
    var.set <- boot_res$var_scal_set
    samp.set <- boot_res$samp_scal_set
  }
  
  return(list(boot_res=boot_res, 
              opt.cr=opt.cr, 
              capture=capture, 
              var.set=var.set, 
              samp.set=samp.set))
}



alter_model <- function(datasets, p, method='scaledlasso') 
{
  full.var<-paste("x",1:p,sep="")
  if(method=='scaledlasso'){
    tem <- scalreg(as.matrix(datasets[,1:p]), as.vector(datasets$y), lam0 = 'univ')
    alter_model_est <- full.var[which(tem$coefficients!=0)]
  }else if(method=='lasso'){
    opt.lambda<-cv.glmnet(x=as.matrix(datasets[,-(p+1)]),y=datasets[,(p+1)],alpha=1)$lambda.min
    tem<-glmnet(x=as.matrix(datasets[,-(p+1)]),y=datasets[,(p+1)],family='gaussian',alpha=1)
    alter_model_est <- full.var[which(coef(tem,s=opt.lambda)[,1]!=0)]
  }else if(method=='adalasso'){
    tem<-adalasso(X=as.matrix(datasets[,1:p]),y=datasets$y,k=10)
    alter_model_est <- full.var[which(tem$coefficients.adalasso!=0)]
  }else if(method=='scad'){
    fit<-cv.ncvreg(X=as.matrix(datasets[,1:p]),y=datasets$y,penalty="SCAD")
    fit2<-fit$fit
    beta<-fit2$beta[,fit$min]
    alter_model_est <- full.var[which(beta[2:(p+1)]!=0)]
  }else if(method=='mcp'){
    fit<-cv.ncvreg(X=as.matrix(datasets[,1:p]),y=datasets$y,penalty="MCP")
    fit2<-fit$fit
    beta<-fit2$beta[,fit$min]
    alter_model_est <- full.var[which(beta[2:(p+1)]!=0)]
  }
  return(alter_model_est)
}

