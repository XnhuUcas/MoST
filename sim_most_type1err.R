rm(list=ls())
gc()
path="/Users/huxiaonan/Desktop/【CSSC-0903Due】MoST_Review/revise_sim/code/"
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
library(glmnet)
library(ncvreg)
library(smoothmest)

############################
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

# True Model
true.model <- rep(0,p)
true.model[1:sig.num] <- 1
var.ind <- paste('x',1:p, sep = '')

# Type I error rates
most_type1err <- array(NA, dim=c(length(comp_methods), length(alpha)))
timetotal <- system.time({
  for(k in 1:length(comp_methods)){
    for(kk in 1:length(alpha)){
      
      cat('----------Now is the (', c(comp_methods[k],alpha[kk]), ')running------------\n')
      most_pval <- c()
      for(rep in 1:s){
        cat('----------Now replication is the ', rep, '-th running------------\n')
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
      most_type1err[k,kk] <- sum(most_pval<alpha[kk])/length(most_pval)
      
      rm(whole.data.null)
      gc()
      
    }
  }
})


most_type1err <- array(NA, dim=c(length(comp_methods), length(alpha)))
timetotal <- system.time({
  for(k in 1:length(comp_methods)){
    for(kk in 1:length(alpha)){
      
      cat('----------Now is the (', c(comp_methods[k],alpha[kk]), ')running(def)------------\n')
      most_pval <- c()
      for(rep in 1:s){
        cat('----------Now replication is the ', rep, '-th running------------\n')
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
      most_type1err[k,kk] <- sum(most_pval<alpha[kk])/length(most_pval)
      
      rm(whole.data.null)
      gc()
      
    }
  }
})


file_name <- paste("most_type1err_def","size",size, "p", p, "s", sig.num, "dis", var.dis, sep = "_")
save.image(file = paste("F://xnhu/MoST/",file_name, ".RData", sep = ""))


