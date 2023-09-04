############################################################################################################
# R code to generate table 5-7 in the main text, and table s1-s2 in the Supplementary Materials
############################################################################################################

source('source_func.R')
source('mylars.R')
source('adalasso.R')
library(scalreg)
library(MASS)
library(ggplot2)
library(RColorBrewer)
library(reshape2)
library(glmnet)
library(ncvreg)
library(smoothmest)


############# Type I error rate ###############

size = 100 
p = 20
sig.num = 5
sd.epsi = 1
rho = 0.3
boot_num = 200 
alpha = c(0.05,0.1)
full.var <- paste("x",1:p,sep="")
comp_methods <- c("adalasso","scad","mcp","lasso","scaledlasso")
var.dis = 'Laplase'
B=200
s=100
true.model <- rep(0,p)
true.model[1:sig.num] <- 1
var.ind <- paste('x',1:p, sep = '')

## fixed coeff 
most_type1err_fix <- array(0, dim=c(length(comp_methods), length(alpha)))
timetotal <- system.time({
  for(k in 1:length(comp_methods)){
    for(kk in 1:length(alpha)){
      most_pval <- c()
      for(rep in 1:s){
        whole.data.null <- Generate_Dataset(size=size, p=p, sig.num=sig.num, sd.epsi = sd.epsi, rho = rho, sig.val=3, var.dis=var.dis)
        res.boot.null <- get_all_methods_mcb(data_input=whole.data.null, r=boot_num, p=p, confidence.level=1-alpha[kk], method = comp_methods[k], true.model=true.model)
        lbm.null <- res.boot.null$boot_res$lower; lbm.null <- var.ind[lbm.null==1]
        ubm.null <- res.boot.null$boot_res$upper; ubm.null <- var.ind[ubm.null==1]
        var.set.null <- res.boot.null$var.set
        var.samp.alter <- res.boot.null$samp.set
        most_null <- sapply(var.set.null, function(x) (length(union(setdiff(lbm.null,x), setdiff(x,ubm.null)))+1)/(p-length(ubm.null)+length(lbm.null)+1))
        most_alter <- (length(union(setdiff(lbm.null,var.samp.alter), setdiff(var.samp.alter,ubm.null)))+1)/(p-length(ubm.null)+length(lbm.null)+1)
        most_pval[rep] <- sum(most_null>=most_alter)/length(most_null)
      }
      most_type1err_fix[k,kk] <- sum(most_pval<alpha[kk])/length(most_pval)
      rm(whole.data.null)
      gc()
    }
  }
})

## decaying coeff 
most_type1err_pd <- array(0, dim=c(length(comp_methods), length(alpha)))
timetotal <- system.time({
  for(k in 1:length(comp_methods)){
    for(kk in 1:length(alpha)){
      most_pval <- c()
      for(rep in 1:s){
        whole.data.null <- Generate_Dataset(size=size, p=p, sig.num=sig.num, sd.epsi = sd.epsi, rho = rho, decay.factor = 0.8, sig.val=3, var.dis=var.dis)
        res.boot.null <- get_all_methods_mcb(data_input=whole.data.null, r=boot_num, p=p, confidence.level=1-alpha[kk], method = comp_methods[k], true.model=true.model)
        lbm.null <- res.boot.null$boot_res$lower; lbm.null <- var.ind[lbm.null==1]
        ubm.null <- res.boot.null$boot_res$upper; ubm.null <- var.ind[ubm.null==1]
        var.set.null <- res.boot.null$var.set
        var.samp.alter <- res.boot.null$samp.set
        most_null <- sapply(var.set.null, function(x) (length(union(setdiff(lbm.null,x), setdiff(x,ubm.null)))+1)/(p-length(ubm.null)+length(lbm.null)+1))
        most_alter <- (length(union(setdiff(lbm.null,var.samp.alter), setdiff(var.samp.alter,ubm.null)))+1)/(p-length(ubm.null)+length(lbm.null)+1)
        most_pval[rep] <- sum(most_null>=most_alter)/length(most_null)
      }
      most_type1err_pd[k,kk] <- sum(most_pval<alpha[kk])/length(most_pval)
      rm(whole.data.null)
      gc()
    }
  }
})


############# Power ###############

size = 100 
p = 20 
sig.num = 5 
alter.beta.num = c(3,8)
sd.epsi = 1
rho=0.3
boot_num = 200 
alpha = 0.05
full.var <- paste("x",1:p,sep="")
comp_methods <- c("adalasso","scad","mcp","lasso","scaledlasso")
var.dis = 'Laplase'
B=200
s=100
true.model <- rep(0,p)
true.model[1:sig.num] <- 1
var.ind <- paste('x',1:p, sep = '')


## decaying coeff 
most_power_pd <- array(0, dim=c(length(comp_methods), length(comp_methods), length(alter.beta.num)))
timetotal3 <- system.time({
  for(k0 in 1:length(comp_methods)){
    for(k in 1:length(comp_methods)){
      for(kk in 1:length(alter.beta.num)){
        most_pval <- c()
        whole.data.null <- Generate_Dataset(size=size, p=p, sig.num=sig.num, sd.epsi = sd.epsi, rho = rho, decay.factor = 0.8, sig.val=3, var.dis=var.dis)
        res.boot.null <- get_all_methods_mcb(data_input=whole.data.null, r=boot_num, p=p, confidence.level=1-alpha, method = comp_methods[k0], true.model=true.model)
        lbm.null <- res.boot.null$boot_res$lower; lbm.null <- var.ind[lbm.null==1]
        ubm.null <- res.boot.null$boot_res$upper; ubm.null <- var.ind[ubm.null==1]
        var.set.null <- res.boot.null$var.set
        most_null <- sapply(var.set.null, function(x) (length(union(setdiff(lbm.null,x), setdiff(x,ubm.null)))+1)/(p-length(ubm.null)+length(lbm.null)+1))
        for(rep in 1:s){
          whole.data.alter <- Generate_Dataset(size=size, p=p, sig.num=alter.beta.num[kk], sd.epsi = sd.epsi, rho = rho, decay.factor = 0.8, sig.val=3, var.dis=var.dis)
          var.samp.alter <- alter_model(datasets=whole.data.alter, p=p, method = comp_methods[k])
          most_alter <- (length(union(setdiff(lbm.null,var.samp.alter), setdiff(var.samp.alter,ubm.null)))+1)/(p-length(ubm.null)+length(lbm.null)+1)
          most_pval[rep] <- sum(most_null>=most_alter)/length(most_null)
        }
        most_power_pd[k,kk] <- sum(most_pval<alpha)/length(most_pval)
        rm(whole.data.null,whole.data.alter)
        gc()
      }
    }
  }
})


## fixed coeff 
most_power_fix <- array(0, dim=c(length(comp_methods), length(comp_methods), length(alter.beta.num)))
timetotal3 <- system.time({
  for(k0 in 1:length(comp_methods)){
    for(k in 1:length(comp_methods)){
      for(kk in 1:length(alter.beta.num)){
        most_pval <- c()
        whole.data.null <- Generate_Dataset(size=size, p=p, sig.num=sig.num, sd.epsi = sd.epsi, rho = rho, sig.val=3, var.dis=var.dis)
        res.boot.null <- get_all_methods_mcb(data_input=whole.data.null, r=boot_num, p=p, confidence.level=1-alpha, method = comp_methods[k0], true.model=true.model)
        lbm.null <- res.boot.null$boot_res$lower; lbm.null <- var.ind[lbm.null==1]
        ubm.null <- res.boot.null$boot_res$upper; ubm.null <- var.ind[ubm.null==1]
        var.set.null <- res.boot.null$var.set
        most_null <- sapply(var.set.null, function(x) (length(union(setdiff(lbm.null,x), setdiff(x,ubm.null)))+1)/(p-length(ubm.null)+length(lbm.null)+1))
        for(rep in 1:s){
          whole.data.alter <- Generate_Dataset(size=size, p=p, sig.num=alter.beta.num[kk], sd.epsi = sd.epsi, rho = rho, sig.val=3, var.dis=var.dis)
          var.samp.alter <- alter_model(datasets=whole.data.alter, p=p, method = comp_methods[k])
          most_alter <- (length(union(setdiff(lbm.null,var.samp.alter), setdiff(var.samp.alter,ubm.null)))+1)/(p-length(ubm.null)+length(lbm.null)+1)
          most_pval[rep] <- sum(most_null>=most_alter)/length(most_null)
        }
        most_power_fix[k,kk] <- sum(most_pval<alpha)/length(most_pval)
        rm(whole.data.null,whole.data.alter)
        gc()
      }
    }
  }
})


