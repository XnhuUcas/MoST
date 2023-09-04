########################################################
# R code for lung squamous carcinoma data analysis in the main text
########################################################

source('source_func_realdata.R')
source('mylars.R')
source('adalasso.R')
library(scalreg)
library(MASS)
library(ggplot2)
library(RColorBrewer)
library(reshape2)
library(glmnet)
library(ncvreg)

scale.func <- function(x){
  return((x-min(x))/(max(x)-min(x)))
}

############### TCGA lung squamous carcinoma data

# Load original data
data_response<-read.table("surv.txt",col.names = c('y','delta'))
data_response<-data_response[2:nrow(data_response),]
data_gene<-read.table("Z.txt")

# screening Original Dataset
pvalue.screen<-c()
datasets <- as.data.frame(cbind(data_response['y'], data_gene)); datasets$y <- as.numeric(datasets$y)
for(i in 1:ncol(data_gene)){
  cat('--',i,'--\n')
  screen.olm<-lm(datasets$y~datasets[,(i+1)],data=datasets)
  pvalue.screen[i]<-summary(screen.olm)$coefficients[2,4]
}
# screen_set <- colnames(datasets)[-1][which(pvalue.screen<=0.1)]
# sort(colnames(datasets)[-1][order(pvalue.screen)[1:33]])

# Selected top rank genes
topvar <- c('AP1S2','BTD','C10ORF54','CA5BP1','CAPN1','CEBPB','EFNA1','FAM107B','FLRT3','KATNB1','LRRC1','LYRM5','AP1S2.1','MYO18A',
            'NOD1','NPLOC4','PLEKHO2','POLR3GL','RAB27A','SECISBP2L','SGTB','STRADB','SWSAP1','TEP1','TERF1','THOC1','TIGD5','TK2',
            'TMEM54','TMEM106A','TOMM7','TRIM34','YARS2')
y <- data_response['y']
datagene.mcb <- data_gene[,topvar]

# Randomly select some genes
rancova <- 67
rsv.ind <- sample(ncol(data_gene)-length(topvar),rancova,replace = F)
data_rangene <- data_gene[,setdiff(names(data_gene),topvar)[rsv.ind]]
data.model <- data.frame(datagene.mcb,data_rangene,y); data.model$y <- as.numeric(data.model$y)

p <- ncol(data.model)-1
boot_num = 2000
confidence.level = 0.95
alpha = 0.05
full.var <- colnames(data.model)[1:p]

#### MCB
comp_methods <- c("adalasso","scad","mcp","lasso","scaledlasso")
timetotal <- system.time({
  opt.lower <- matrix(NA, length(comp_methods), p)
  opt.upper <- matrix(NA, length(comp_methods), p)
  opt.width <- c()
  opt.cr <- c()
  samp.set <- vector(mode = 'list', length = length(comp_methods))
  for(mm in 1:length(comp_methods)){
    method =  comp_methods[mm]
    CI_res <- get_all_methods_mcb(data_input=data.model, r=boot_num, p=p, confidence.level=confidence.level, method = comp_methods[mm])
    opt.lower[mm,] <- CI_res$boot_res$lower
    opt.upper[mm,] <- CI_res$boot_res$upper
    opt.width[mm] <- CI_res$opt.width
    opt.cr[mm] <- CI_res$opt.cr
    samp.set[[mm]] <- CI_res$samp.set
  }
})

## settings
boot_num = 2000
s=1000
most_pval <- matrix(NA, nrow=s, ncol=4)
most_pval2 <- matrix(NA, nrow=s, ncol=4)
pred.mse <- matrix(NA, nrow=s, ncol=4)
for(rr in 1:s){

  samp_tr_te <- sample(1:nrow(data.model),304)
  data.model.train <- as.matrix(data.model[samp_tr_te,])
  data.model.test <- as.matrix(data.model[-samp_tr_te,])
  
  ## testing
  res.boot.null <- get_all_methods_mcb(data_input=data.model.train, r=boot_num, p=p, confidence.level=confidence.level, method = "scaledlasso")
  lbm.null <- res.boot.null$boot_res$lower; lbm.null <- full.var[lbm.null==1]
  ubm.null <- res.boot.null$boot_res$upper; ubm.null <- full.var[ubm.null==1]
  var.set.null <- res.boot.null$var.set
  var_null_mat <- t(sapply(var.set.null, function(x) as.numeric(full.var%in%x)))
  weight.var<-apply(var_null_mat,2,sum)/boot_num
  w_in_null <- sapply(var.set.null, function(x) mean(weight.var[which(full.var%in%x==T)], na.rm=TRUE))
  w_out_null <- sapply(var.set.null, function(x) (length(union(setdiff(lbm.null,x), setdiff(x,ubm.null)))+1)/(p-length(ubm.null)+length(lbm.null)+1))
  most_null <- sapply(var.set.null, function(x) (length(union(setdiff(lbm.null,x), setdiff(x,ubm.null)))+1)/(p-length(ubm.null)+length(lbm.null)+1))
  
  
  # random model
  tempvar = setdiff(full.var,ubm.null)
  var.samp.alter <- vector(mode = 'list', length = 4)
  var.samp.alter[[1]] <- ubm.null
  var.samp.alter[[2]] <- c(ubm.null,tempvar[sample(1:length(tempvar),1)])
  var.samp.alter[[3]] <- c(ubm.null,tempvar[sample(1:length(tempvar),2)])
  var.samp.alter[[4]] <- lbm.null
  for(i in 1:length(var.samp.alter)){

    w_in_alter <- mean(weight.var[which(full.var%in%var.samp.alter[[i]]==T)])
    w_out_alter <- (length(union(setdiff(lbm.null,var.samp.alter[[i]]), setdiff(var.samp.alter[[i]],ubm.null)))+1)/(p-length(ubm.null)+length(lbm.null)+1)
    most_pval[rr,i] <- sum((abs(w_in_null-1)+abs(w_out_null-0))>=(abs(w_in_alter-1)+abs(w_out_alter-0)))/boot_num
    most_pval2[rr,i] <- sum(most_null>=w_out_alter)/length(most_null)
    
    ## prediction
    datasets <- as.data.frame(data.model.train)
    data.new <- cbind(datasets[var.samp.alter[[i]]], datasets$y); colnames(data.new)[ncol(data.new)] = 'y'
    olm.new <- lm(data.new$y~.,data=data.new)
    beta.est <- olm.new$coefficients

    datasets.test <- as.data.frame(data.model.test)
    cova.test <- datasets.test[var.samp.alter[[i]]]
    if(length(var.samp.alter[[i]])==0){
      pred.mse[rr,i] <- sum((as.matrix(datasets.test$y) - beta.est[1])^2)/100
    }else{
      pred.mse[rr,i] <- sum((as.matrix(datasets.test$y) - as.matrix(cova.test)%*%as.matrix(beta.est[2:length(beta.est)])-beta.est[1])^2)/100
    }
  }
}

pred.mse.mean <- apply(pred.mse,2,mean, na.rm=TRUE)
pred.mse.med <- apply(pred.mse,2,median, na.rm=TRUE)
pred.mse.sd <- apply(pred.mse,2,sd, na.rm=TRUE)
most_pval.trans <- ifelse(most_pval<alpha,1,0)
most_pval.trans2 <- ifelse(most_pval2<alpha,1,0)
sig.times <- apply(most_pval.trans,2,sum)
sig.times2 <- apply(most_pval.trans2,2,sum)


######################## visualization MCB 
library(funprog)
library(ggplot2)

mcb.lower.adalasso <- which(opt.lower[1,]!=0)
mcb.upper.adalasso <- which(opt.upper[1,]!=0)
samp.adalasso <- which(full.var%in%samp.set[[1]]==TRUE)

mcb.lower.scad <- which(opt.lower[2,]!=0)
mcb.upper.scad <- which(opt.upper[2,]!=0)
samp.scad <- which(full.var%in%samp.set[[2]]==TRUE)

mcb.lower.mcp <- which(opt.lower[3,]!=0)
mcb.upper.mcp <- which(opt.upper[3,]!=0)
samp.mcp <- which(full.var%in%samp.set[[3]]==TRUE)

mcb.lower.lasso <- which(opt.lower[4,]!=0)
mcb.upper.lasso <- which(opt.upper[4,]!=0)
samp.lasso <- which(full.var%in%samp.set[[4]]==TRUE)

mcb.lower.scaledlasso <- which(opt.lower[5,]!=0)
mcb.upper.scaledlasso <- which(opt.upper[5,]!=0)
samp.scaledlasso <- which(full.var%in%samp.set[[5]]==TRUE)

# mcb.lower.adalasso <- which(opt.lower[1,1:50]!=0)
# mcb.upper.adalasso <- which(opt.upper[1,1:50]!=0)
# samp.adalasso <- which(full.var[1:50]%in%samp.set[[1]]==TRUE)
# 
# mcb.lower.scad <- which(opt.lower[2,1:50]!=0)
# mcb.upper.scad <- which(opt.upper[2,1:50]!=0)
# samp.scad <- which(full.var[1:50]%in%samp.set[[2]]==TRUE)
# 
# mcb.lower.mcp <- which(opt.lower[3,1:50]!=0)
# mcb.upper.mcp <- which(opt.upper[3,1:50]!=0)
# samp.mcp <- which(full.var[1:50]%in%samp.set[[3]]==TRUE)
# 
# mcb.lower.lasso <- which(opt.lower[4,1:50]!=0)
# mcb.upper.lasso <- which(opt.upper[4,1:50]!=0)
# samp.lasso <- which(full.var[1:50]%in%samp.set[[4]]==TRUE)
# 
# mcb.lower.scaledlasso <- which(opt.lower[5,1:50]!=0)
# mcb.upper.scaledlasso <- which(opt.upper[5,1:50]!=0)
# samp.scaledlasso <- which(full.var[1:50]%in%samp.set[[5]]==TRUE)

# ind.union <- Reduce(union, list(mcb.upper.adalasso,
#                                 samp.adalasso,
#                                 mcb.upper.scad,
#                                 samp.scad,
#                                 mcb.upper.mcp,
#                                 samp.mcp,
#                                 mcb.upper.lasso,
#                                 samp.lasso,
#                                 mcb.upper.scaledlasso,
#                                 samp.scaledlasso))
# full.var.comb <- colnames(data.model)[ind.union]
full.var.comb <- full.var

obs.model.adalasso <- which(full.var.comb%in%samp.set[[1]]!=0)
obs.model.scad <- which(full.var.comb%in%samp.set[[2]]!=0)
obs.model.mcp <- which(full.var.comb%in%samp.set[[3]]!=0)
obs.model.lasso <- which(full.var.comb%in%samp.set[[4]]!=0)
obs.model.scaledlasso <- which(full.var.comb%in%samp.set[[5]]!=0)

methods<-c(rep(1.2,length(mcb.lower.adalasso)),rep(1.4,length(obs.model.adalasso)),rep(1.6,length(mcb.upper.adalasso)),
           rep(2.2,length(mcb.lower.scad)),rep(2.4,length(obs.model.scad)),rep(2.6,length(mcb.upper.scad)),
           rep(3.2,length(mcb.lower.mcp)),rep(3.4,length(obs.model.mcp)),rep(3.6,length(mcb.upper.mcp)),
           rep(4.2,length(mcb.lower.lasso)),rep(4.4,length(obs.model.lasso)),rep(4.6,length(mcb.upper.lasso)),
           rep(5.2,length(mcb.lower.scaledlasso)),rep(5.4,length(obs.model.scaledlasso)),rep(5.6,length(mcb.upper.scaledlasso)))

genes<-c(mcb.lower.adalasso,obs.model.adalasso,mcb.upper.adalasso,
         mcb.lower.scad,obs.model.scad,mcb.upper.scad,
         mcb.lower.mcp,obs.model.mcp,mcb.upper.mcp,
         mcb.lower.lasso,obs.model.lasso,mcb.upper.lasso,
         mcb.lower.scaledlasso,obs.model.scaledlasso,mcb.upper.scaledlasso)

legends<-c(rep('LBM',length(mcb.lower.adalasso)),rep('Selected Model',length(obs.model.adalasso)),rep('UBM',length(mcb.upper.adalasso)),
           rep('LBM',length(mcb.lower.scad)),rep('Selected Model',length(obs.model.scad)),rep('UBM',length(mcb.upper.scad)),
           rep('LBM',length(mcb.lower.mcp)),rep('Selected Model',length(obs.model.mcp)),rep('UBM',length(mcb.upper.mcp)),
           rep('LBM',length(mcb.lower.lasso)),rep('Selected Model',length(obs.model.lasso)),rep('UBM',length(mcb.upper.lasso)),
           rep('LBM',length(mcb.lower.scaledlasso)),rep('Selected Model',length(obs.model.scaledlasso)),rep('UBM',length(mcb.upper.scaledlasso)))

dataplot<-data.frame(Genes=genes,Methods=methods,Legend=legends)

xind<-c(1.4,2.4,3.4,4.4,5.4)
xnew<-c('Adaptive Lasso', 'SCAD', 'MCP', 'Lasso', 'Scaled Lasso')
dindx<-data.frame(xind,xnew)
yind<-seq(1,length(full.var.comb))
ynew<-full.var.comb
dindy<-data.frame(yind,ynew)

pdf('fig_mcb_lungdata2.pdf', width = 8, height = 10)
ggplot(data=dataplot, aes(x=Methods, y=Genes, color = Legend, shape=Legend)) +
  geom_point(aes(color=Legend, shape=Legend), size=2) +
  labs(x="", y="") +
  scale_x_continuous(breaks=dindx$xind, labels = dindx$xnew) + 
  scale_y_continuous(breaks=dindy$yind, labels = dindy$ynew) + 
  theme_set(theme_bw()) +
  theme(panel.grid =element_blank(),
        legend.title = element_blank(),
        axis.text.y = element_text(size=7),
        axis.text.x = element_text(size=10),
        legend.text = element_text(size = 10),
        legend.position = 'bottom',
        plot.margin = margin(t=10,r=50,b=1,l=1))
dev.off()

