rm(list=ls())
gc()
set.seed(123456)
path="F:/xnhu/MoST/MoST_codes/"
setwd(path)
source('source_func.R')
source('mylars.R')
source('adalasso.R')
library(parallel)
library(scalreg)
library(MASS)
library(ggplot2)
library(RColorBrewer)
library(reshape2)
library(simpleboot)
library(fANCOVA)
library(glmnet)
library(methods)
library(leaps)
library(ncvreg)
library(smoothmest)

# -------------- Run ----------------
size=100
p=20
sig.num=5
sd.epsi=1
rho=0.3
r=200
confidence.level <- c(0.9,0.95)
s=100
full.var<-paste("x",1:p,sep="")
threshold=0.05

# True Model
true.model<-rep(0,p)
true.model[1:sig.num]<-1

comp_methods <- c("adalasso","scad","mcp","lasso","scaledlasso")
timetotal <- system.time({
  
  opt.width <- array(0,dim=c(length(comp_methods), s, length(confidence.level)))
  opt.cr <- array(0,dim=c(length(comp_methods), s, length(confidence.level)))
  truecapture <- array(0,dim=c(length(comp_methods), s, length(confidence.level)))
  
  for(mm in 1:length(comp_methods)){
    
    cat('//////////////////////////////////////////////////////////////////\n')
    cat('----------Now is the', comp_methods[mm], '(Laplase+sd1)... -----------\n')
    method =  comp_methods[mm]
    
    for(k in 1:s){
      
      cat('----------Now replication is the', k, '-th rounds -----------\n')
      
      x <- Generate_Dataset(size,p,sig.num,sd.epsi,rho, sig.val = 3, var.dis = 'Laplase')
      if(method == "adalasso"){
        var_adalasso <- matrix(0,nrow=r,ncol=p)
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
        }
        for(i in 1:length(confidence.level)){
          cat('----------Now confidence.level=', confidence.level[i], '.. -----------\n')
          CI_res <- get_opt_mcb_binary_search(var.matrix = var_adalasso, confidence.level=confidence.level[i], p=p, r=r)
          opt.lower <- CI_res$optimal.mcbs.lower
          opt.upper <- CI_res$optimal.mcbs.upper
          opt.width[mm,k,i] <- CI_res$optimal.width
          opt.cr[mm,k,i] <- CI_res$optimal.cr
          truecapture[mm,k,i]<-cap_ture(low_list=opt.lower,upper_list=opt.upper,true.model=true.model)
        }
      }
      else if(method == "scad"){
        # residual bootstrap
        var_scad <- matrix(0,nrow=r,ncol=p)
        fit <- cv.ncvreg(X=as.matrix(x[,1:p]),y=x$y,penalty='SCAD')
        fit2 <- fit$fit
        beta <- fit2$beta[,fit$min]
        res_original <- x$y - as.matrix(x[,1:p]) %*% beta[2:(p+1)] - beta[1]
        res_after_center <- res_original - mean(res_original)
        constant <- as.matrix(x[,1:p]) %*% beta[2:(p+1)] + beta[1]
        for(j in 1:r) { 
          ind=sample(1:nrow(x),nrow(x),replace=T)
          new_response <- constant + res_after_center[ind]
          boot.data <- cbind(x[,1:p], new_response)
          tem <- cv.ncvreg(X=as.matrix(boot.data[,1:p]),y=boot.data$new_response,penalty='SCAD')
          fit2_tem<-tem$fit
          beta_tem<-fit2_tem$beta[,tem$min]
          beta_tem <- beta_tem[2:(p+1)]
          var_scad[j,]<-full.var%in%names(beta_tem)[which(beta_tem!=0)]
        }
        for(i in 1:length(confidence.level)){
          cat('----------Now confidence.level=', confidence.level[i], '.. -----------\n')
          CI_res <- get_opt_mcb_binary_search(var.matrix = var_scad, confidence.level=confidence.level[i], p=p, r=r)
          opt.lower <- CI_res$optimal.mcbs.lower
          opt.upper <- CI_res$optimal.mcbs.upper
          opt.width[mm,k,i] <- CI_res$optimal.width
          opt.cr[mm,k,i] <- CI_res$optimal.cr
          truecapture[mm,k,i]<-cap_ture(low_list=opt.lower,upper_list=opt.upper,true.model=true.model)
        }
      }
      else if(method == "mcp"){
        # residual bootstrap
        var_mcp <- matrix(0,nrow=r,ncol=p)
        fit <- cv.ncvreg(X=as.matrix(x[,1:p]),y=x$y,penalty='MCP')
        fit2 <- fit$fit
        beta <- fit2$beta[,fit$min]
        res_original <- x$y - as.matrix(x[,1:p]) %*% beta[2:(p+1)] - beta[1]
        res_after_center <- res_original - mean(res_original)
        constant <- as.matrix(x[,1:p]) %*% beta[2:(p+1)] + beta[1]
        for(j in 1:r) { 
          ind=sample(1:nrow(x),nrow(x),replace=T)
          new_response <- constant + res_after_center[ind]
          boot.data <- cbind(x[,1:p], new_response)
          tem <- cv.ncvreg(X=as.matrix(boot.data[,1:p]),y=boot.data$new_response,penalty='MCP')
          fit2_tem<-tem$fit
          beta_tem<-fit2_tem$beta[,tem$min]
          beta_tem <- beta_tem[2:(p+1)]
          var_mcp[j,]<-full.var%in%names(beta_tem)[which(beta_tem!=0)]
        }
        for(i in 1:length(confidence.level)){
          cat('----------Now confidence.level=', confidence.level[i], '.. -----------\n')
          CI_res <- get_opt_mcb_binary_search(var.matrix = var_mcp, confidence.level=confidence.level[i], p=p, r=r)
          opt.lower <- CI_res$optimal.mcbs.lower
          opt.upper <- CI_res$optimal.mcbs.upper
          opt.width[mm,k,i] <- CI_res$optimal.width
          opt.cr[mm,k,i] <- CI_res$optimal.cr
          truecapture[mm,k,i]<-cap_ture(low_list=opt.lower,upper_list=opt.upper,true.model=true.model)
        }
      }
      else if(method == "lasso"){
        # residual bootstrap
        var_lasso <- matrix(0,nrow=r,ncol=p)
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
        }
        
        for(i in 1:length(confidence.level)){
          cat('----------Now confidence.level=', confidence.level[i], '.. -----------\n')
          CI_res <- get_opt_mcb_binary_search(var.matrix = var_lasso, confidence.level=confidence.level[i], p=p, r=r)
          opt.lower <- CI_res$optimal.mcbs.lower
          opt.upper <- CI_res$optimal.mcbs.upper
          opt.width[mm,k,i] <- CI_res$optimal.width
          opt.cr[mm,k,i] <- CI_res$optimal.cr
          truecapture[mm,k,i]<-cap_ture(low_list=opt.lower,upper_list=opt.upper,true.model=true.model)
        }
      }
      else if(method == "scaledlasso"){
        # residual bootstrap
        var_scal <- matrix(0,nrow=r,ncol=p)
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
        }
        for(i in 1:length(confidence.level)){
          cat('----------Now confidence.level=', confidence.level[i], '.. -----------\n')
          CI_res <- get_opt_mcb_binary_search(var.matrix = var_scal, confidence.level=confidence.level[i], p=p, r=r)
          opt.lower <- CI_res$optimal.mcbs.lower
          opt.upper <- CI_res$optimal.mcbs.upper
          opt.width[mm,k,i] <- CI_res$optimal.width
          opt.cr[mm,k,i] <- CI_res$optimal.cr
          truecapture[mm,k,i]<-cap_ture(low_list=opt.lower,upper_list=opt.upper,true.model=true.model)
        }
      }
      else{
        warnings("You must input the argment of method!!")
      }
      
    }
  }
  
})


file_name <- paste("sim_muc","size",size, "p", p, "sd", sd.epsi, "b", r, "rep",s, 
                   "sig", sig.num, "laplase", "penall", sep = "_")
save.image(file = paste("F:/xnhu/MoST/",file_name, ".RData", sep = ""))


