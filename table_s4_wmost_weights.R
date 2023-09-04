########################################################
# R code to generate table s4 in the Supplementary Materials
########################################################

source('source_func.R')
source('mylars.R')
source('adalasso.R')
library(scalreg)
library(MASS)
library(glmnet)
library(ncvreg)
library(smoothmest)

size = 100 
p = 20
sig.num = 5
alter.beta.num = 3
sd.epsi = 1
rho = 0.3
boot_num = 200 
alpha = 0.05
theta1 = c(0.7, 0.3, 0.5, 0.8)
theta2 = c(0.3, 0.7, 0.5, 0.6)
full.var <- paste("x",1:p,sep="")
comp_methods <- c("adalasso","scad","mcp","lasso","scaledlasso")
var.dis = 'normal'
B=200
s=100
true.model <- rep(0,p)
true.model[1:sig.num] <- 1
var.ind <- paste('x',1:p, sep = '')


################# Type I error rates ################
most_type1err <- array(0, dim=c(5, length(comp_methods)))
timetotal <- system.time({
  for(h in 1:5){
    for(k in 1:length(comp_methods)){
      most_pval <- c()
      for(rep in 1:s){
        whole.data.null <- Generate_Dataset(size=size, p=p, sig.num=sig.num, sd.epsi = sd.epsi, rho = rho, sig.val=3, var.dis = var.dis)
        res.boot.null <- get_all_methods_mcb(data_input=whole.data.null, r=boot_num, p=p, confidence.level=1-alpha, method = comp_methods[k], true.model=true.model)
        lbm.null <- res.boot.null$boot_res$lower; lbm.null <- var.ind[lbm.null==1]
        ubm.null <- res.boot.null$boot_res$upper; ubm.null <- var.ind[ubm.null==1]
        var.set.null <- res.boot.null$var.set
        var.samp.alter <- res.boot.null$samp.set
        var_null_mat <- t(sapply(var.set.null, function(x) as.numeric(full.var%in%x)))
        weight.var<-apply(var_null_mat,2,sum)/boot_num
        w_in_null <- sapply(var.set.null, function(x) mean(weight.var[which(full.var%in%x==T)], na.rm=TRUE))
        w_out_null <- sapply(var.set.null, function(x) (length(union(setdiff(lbm.null,x), setdiff(x,ubm.null)))+1)/(p-length(ubm.null)+length(lbm.null)+1))
        w_in_alter <- mean(weight.var[which(full.var%in%var.samp.alter==T)], na.rm=TRUE)
        w_out_alter <- (length(union(setdiff(lbm.null,var.samp.alter), setdiff(var.samp.alter,ubm.null)))+1)/(p-length(ubm.null)+length(lbm.null)+1)
        most_pval[rep] <- sum((theta1[h]*abs(w_in_null-1)+theta2[h]*abs(w_out_null-0))>=(theta1[h]*abs(w_in_alter-1)+theta2[h]*abs(w_out_alter-0)))/boot_num
      }
      most_type1err[h,k] <- sum(most_pval<alpha[kk])/length(most_pval)
      rm(whole.data.null)
      gc()
      cat('---------comp_methods=', comp_methods[k], '---done!!!---------\n')
    }
    cat('---------(theta1,theta2)=', c(theta1[h], theta2[h]), '---done!!!---------\n')
  }
})


################# Power ################
most_power <- array(0, dim=c(5, length(comp_methods)))
timetotal2 <- system.time({
  for(h in 1:5){
    for(k in 1:length(comp_methods)){
      most_pval <- c()
      whole.data.null <- Generate_Dataset(size=size, p=p, sig.num=sig.num, sd.epsi = sd.epsi, rho = rho, sig.val=3, var.dis=var.dis)
      res.boot.null <- get_all_methods_mcb(data_input=whole.data.null, r=boot_num, p=p, confidence.level=1-alpha, method = "scaledlasso", true.model=true.model)
      lbm.null <- res.boot.null$boot_res$lower; lbm.null <- var.ind[lbm.null==1]
      ubm.null <- res.boot.null$boot_res$upper; ubm.null <- var.ind[ubm.null==1]
      var.set.null <- res.boot.null$var.set
      var_null_mat <- t(sapply(var.set.null, function(x) as.numeric(full.var%in%x)))
      weight.var<-apply(var_null_mat,2,sum)/boot_num
      w_in_null <- sapply(var.set.null, function(x) mean(weight.var[which(full.var%in%x==T)], na.rm=TRUE))
      w_out_null <- sapply(var.set.null, function(x) (length(union(setdiff(lbm.null,x), setdiff(x,ubm.null)))+1)/(p-length(ubm.null)+length(lbm.null)+1))
      for(rep in 1:s){
        whole.data.alter <- Generate_Dataset(size=size, p=p, sig.num=alter.beta.num, sd.epsi = sd.epsi, rho = rho, sig.val=3, var.dis=var.dis)
        var.samp.alter <- alter_model(datasets=whole.data.alter, p=p, method = comp_methods[k])
        w_in_alter <- mean(weight.var[which(full.var%in%var.samp.alter==T)], na.rm=TRUE)
        w_out_alter <- (length(union(setdiff(lbm.null,var.samp.alter), setdiff(var.samp.alter,ubm.null)))+1)/(p-length(ubm.null)+length(lbm.null)+1)
        most_pval[rep] <- sum((theta1[h]*abs(w_in_null-1)+theta2[h]*abs(w_out_null-0))>=(theta1[h]*abs(w_in_alter-1)+theta2[h]*abs(w_out_alter-0)))/boot_num
      }
      most_power[h,k] <- sum(most_pval<alpha)/length(most_pval)
      rm(whole.data.null,whole.data.alter)
      gc()
      cat('---------comp_methods=', comp_methods[k], '---done!!!---------\n')
    }
    cat('---------(theta1,theta2)=', c(theta1[h], theta2[h]), '---done!!!---------\n')
  }
})

