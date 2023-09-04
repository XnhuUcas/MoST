########################################################
# R code for Figure S3 in the Supplementary Materials
########################################################

source('source_func.R')
library(scalreg)
library(MASS)
library(glmnet)
library(ncvreg)
library(smoothmest)

size = 100 
p.cand = c(20, 200) 
sig.num = 5
alter.beta.num = 3
sd.epsi = 1
rho = 0.3
boot_num = 200 
alpha = 0.05
comp_methods <- c("scad","scaledlasso")
var.dis = 'normal'
B=200
s=100


#### Type I error rates

res_total_type1err <- array(0, dim=c(length(p.cand), length(comp_methods), 4))
for(k0 in 1:length(p.cand)){
  
  p <- p.cand[k0]
  full.var <- paste("x",1:p,sep="")
  true.model <- rep(0,p)
  true.model[1:sig.num] <- 1
  var.ind <- paste('x',1:p, sep = '')
  whole.data.null <- Generate_Dataset(size=size, p=p, sig.num=sig.num, sd.epsi = sd.epsi, rho = rho, sig.val=3, var.dis=var.dis)
  
  for(k in 1:length(comp_methods)){
    
    timetotal.type1.single <- system.time({
      res.boot.null <- get_all_methods_mcb(data_input=whole.data.null, r=boot_num, p=p, confidence.level=1-alpha, method = comp_methods[k], true.model=true.model)
      lbm.null <- res.boot.null$boot_res$lower; lbm.null <- var.ind[lbm.null==1]
      ubm.null <- res.boot.null$boot_res$upper; ubm.null <- var.ind[ubm.null==1]
      var.set.null <- res.boot.null$var.set
      var.samp.alter <- res.boot.null$samp.set
      most_null <- sapply(var.set.null, function(x) (length(union(setdiff(lbm.null,x), setdiff(x,ubm.null)))+1)/(p-length(ubm.null)+length(lbm.null)+1))
      most_pval <- c()
      for(rep in 1:s){
        whole.data.alter <- Generate_Dataset(size=size, p=p, sig.num=sig.num, sd.epsi = sd.epsi, rho = rho, sig.val=3, var.dis=var.dis)
        var.samp.alter <- alter_model(datasets=whole.data.alter, p=p, method = comp_methods[k])
        most_alter <- (length(union(setdiff(lbm.null,var.samp.alter), setdiff(var.samp.alter,ubm.null)))+1)/(p-length(ubm.null)+length(lbm.null)+1)
        most_pval[rep] <- sum(most_null>=most_alter)/length(most_null)
        rm(whole.data.alter)
        gc()
      }
      res_total_type1err[k0,k,1] <- sum(most_pval<alpha)/length(most_pval)
      cat('---------single bootstrap method---', comp_methods[k], '---done!!!---------\n')
    })
    
    
    timetotal.type1.double <- system.time({
      most_null_double <- c()
      for(bb in 1:boot_num){
        ind <- sample(1:nrow(whole.data.null),nrow(whole.data.null),replace=T)
        boot.data <- whole.data.null[ind,]
        res.boot.null.double <- get_all_methods_mcb(data_input=boot.data, r=boot_num, p=p, confidence.level=1-alpha, method = comp_methods[k], true.model=true.model)
        lbm.null.double <- res.boot.null.double$boot_res$lower; lbm.null.double <- var.ind[lbm.null.double==1]
        ubm.null.double <- res.boot.null.double$boot_res$upper; ubm.null.double <- var.ind[ubm.null.double==1]
        var.samp.null.double <- res.boot.null.double$samp.set
        most_null_double[bb] <- (length(union(setdiff(lbm.null.double,var.samp.null.double), setdiff(var.samp.null.double,ubm.null.double)))+1)/(p-length(ubm.null.double)+length(lbm.null.double)+1)
      }
      
      res.boot.null <- get_all_methods_mcb(data_input=whole.data.null, r=boot_num, p=p, confidence.level=1-alpha, method = comp_methods[k], true.model=true.model)
      lbm.null <- res.boot.null$boot_res$lower; lbm.null <- var.ind[lbm.null==1]
      ubm.null <- res.boot.null$boot_res$upper; ubm.null <- var.ind[ubm.null==1]
      most_pval_double <- c()
      for(rep in 1:s){
        whole.data.alter.double <- Generate_Dataset(size=size, p=p, sig.num=sig.num, sd.epsi = sd.epsi, rho = rho, sig.val=3, var.dis=var.dis)
        var.samp.alter.double <- alter_model(datasets=whole.data.alter.double, p=p, method = comp_methods[k])
        most_alter_double <- (length(union(setdiff(lbm.null,var.samp.alter.double), setdiff(var.samp.alter.double,ubm.null)))+1)/(p-length(ubm.null)+length(lbm.null)+1)
        most_pval_double[rep] <- sum(most_null_double>=most_alter_double)/length(most_null_double)
        rm(whole.data.alter.double)
        gc()
      }
      res_total_type1err[k0,k,2] <- sum(most_pval_double<alpha)/length(most_pval_double)
      cat('---------twolayer bootstrap method---', comp_methods[k], '---done!!!---------\n')
    })
    
    res_total_type1err[k0,k,3] <- as.double(timetotal.type1.single[3])
    res_total_type1err[k0,k,4] <- as.double(timetotal.type1.double[3])
    
  }
}



#### Power 

res_total_power <- array(0, dim=c(length(p.cand), length(comp_methods), 4))
for(k0 in 1:length(p.cand)){
  
  p <- p.cand[k0]
  full.var <- paste("x",1:p,sep="")
  true.model <- rep(0,p)
  true.model[1:sig.num] <- 1
  var.ind <- paste('x',1:p, sep = '')
  whole.data.null <- Generate_Dataset(size=size, p=p, sig.num=sig.num, sd.epsi = sd.epsi, rho = rho, sig.val=3, var.dis=var.dis)
  
  for(k in 1:length(comp_methods)){
    
    timetotal.power.single <- system.time({
      most_pval <- c()
      res.boot.null <- get_all_methods_mcb(data_input=whole.data.null, r=boot_num, p=p, confidence.level=1-alpha, method = "scaledlasso", true.model=true.model)
      lbm.null <- res.boot.null$boot_res$lower; lbm.null <- var.ind[lbm.null==1]
      ubm.null <- res.boot.null$boot_res$upper; ubm.null <- var.ind[ubm.null==1]
      var.set.null <- res.boot.null$var.set
      most_null <- sapply(var.set.null, function(x) (length(union(setdiff(lbm.null,x), setdiff(x,ubm.null)))+1)/(p-length(ubm.null)+length(lbm.null)+1))
      for(rep in 1:s){
        whole.data.alter <- Generate_Dataset(size=size, p=p, sig.num=alter.beta.num, sd.epsi = sd.epsi, rho = rho, sig.val=3, var.dis=var.dis)
        var.samp.alter <- alter_model(datasets=whole.data.alter, p=p, method = comp_methods[k])
        most_alter <- (length(union(setdiff(lbm.null,var.samp.alter), setdiff(var.samp.alter,ubm.null)))+1)/(p-length(ubm.null)+length(lbm.null)+1)
        most_pval[rep] <- sum(most_null>=most_alter)/length(most_null)
        rm(whole.data.alter)
        gc()
      }
      res_total_power[k0,k,1] <- sum(most_pval<alpha)/length(most_pval)
      cat('---------single bootstrap method---', comp_methods[k], '---done!!!---------\n')
    })
    
    
    timetotal.power.double <- system.time({
      most_null_double <- c()
      for(bb in 1:boot_num){
        ind <- sample(1:nrow(whole.data.null),nrow(whole.data.null),replace=T)
        boot.data <- whole.data.null[ind,]
        res.boot.null.double <- get_all_methods_mcb(data_input=boot.data, r=boot_num, p=p, confidence.level=1-alpha, method = "scaledlasso", true.model=true.model)
        lbm.null.double <- res.boot.null.double$boot_res$lower; lbm.null.double <- var.ind[lbm.null.double==1]
        ubm.null.double <- res.boot.null.double$boot_res$upper; ubm.null.double <- var.ind[ubm.null.double==1]
        var.samp.null.double <- res.boot.null.double$samp.set
        most_null_double[bb] <- (length(union(setdiff(lbm.null.double,var.samp.null.double), setdiff(var.samp.null.double,ubm.null.double)))+1)/(p-length(ubm.null.double)+length(lbm.null.double)+1)
      }
      
      res.boot.null <- get_all_methods_mcb(data_input=whole.data.null, r=boot_num, p=p, confidence.level=1-alpha, method = "scaledlasso", true.model=true.model)
      lbm.null <- res.boot.null$boot_res$lower; lbm.null <- var.ind[lbm.null==1]
      ubm.null <- res.boot.null$boot_res$upper; ubm.null <- var.ind[ubm.null==1]
      most_pval_double <- c()
      for(rep in 1:s){
        whole.data.alter.double <- Generate_Dataset(size=size, p=p, sig.num=alter.beta.num, sd.epsi = sd.epsi, rho = rho, sig.val=3, var.dis=var.dis)
        var.samp.alter.double <- alter_model(datasets=whole.data.alter.double, p=p, method = comp_methods[k])
        most_alter_double <- (length(union(setdiff(lbm.null,var.samp.alter.double), setdiff(var.samp.alter.double,ubm.null)))+1)/(p-length(ubm.null)+length(lbm.null)+1)
        most_pval_double[rep] <- sum(most_null_double>=most_alter_double)/length(most_null_double)
        rm(whole.data.alter.double)
        gc()
      }
      res_total_power[k0,k,2] <- sum(most_pval_double<alpha)/length(most_pval_double)
      cat('---------twolayer bootstrap method---', comp_methods[k], '---done!!!---------\n')
    })
    
    res_total_power[k0,k,3] <- timetotal.power.single[3]
    res_total_power[k0,k,4] <- timetotal.power.double[3]

  }
}


