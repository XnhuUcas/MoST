########################################################
# R code for Table 1 and Figure 1 in the main text
########################################################

source('source_func.R')
source('mylars.R')
source('adalasso.R')

size=200;rho=0.3
p=seq(100,1000,100)
sig.num=10
sd.epsi=1
r=200
confidence.level=0.95
threshold=0.05

var.instances <- vector(length=length(p), mode = 'list')
for(i in 1:length(p)){
  full.var<-paste("x",1:p[i],sep="")
  true.model<-rep(0,p[i])
  true.model[1:sig.num]<-1
  test.data<-Generate_Dataset(size,p[i],sig.num,sd.epsi,rho)
  var.instances[[i]]<-boot_var_adalasso(x = test.data,r=r, p=p[i],confidence.level=confidence.level, true.model=true.model)
}


time_a1=vector(length = length(p),mode='list')
time_a3=vector(length = length(p),mode='list')
opt.pimcb=vector(length = length(p),mode='list')
opt.bimcb=vector(length = length(p),mode='list')
var_mat=vector(length = length(p),mode='list')

for(i in 1:length(p)){
  full.var<-paste("x",1:p[i],sep="")
  var_mat[[i]] <- matrix(0,nrow=r,ncol=p[i])
  for(j in 1:r){
    var_mat[[i]][j,] <- full.var%in%var.instances[[i]]$var_adalasso_set[[j]]
  }
}

for(k in 1:length(p)){
  cat('-------- p is',p[k],'----------\n')
  time_a1[[k]] = system.time({
    result<-mcb_candidate(var_mat[[k]],p[k],r)
    opt.pimcb[[k]]<-get_opt_mcb_pi(result, confidence.level)
  })
  
  time_a3[[k]] = system.time({
    opt.bimcb[[k]]<-get_opt_mcb_binary_search(var_mat[[k]],confidence.level,p[k],r)
  })
}


## BCR and Width

file_names <- list.files(path)
out.tab <- NULL
for(hh in 1:length(file_names)){
  cat('-------- file is',file_names[hh],'----------\n')
  load(file_names[hh])
  res.total <- matrix(NA, nrow=4, ncol=length(p))
  for(k in 1:length(p)){
    cat('-------- p is',p[k],'----------\n')
    result<-mcb_candidate(var_mat[[k]],p[k],r)
    opt.pimcb<-get_opt_mcb_pi(result, confidence.level)
    opt.bimcb<-get_opt_mcb_binary_search(var_mat[[k]],confidence.level,p[k],r)
    res.total[1,k] <- opt.pimcb$width
    res.total[2,k] <- opt.bimcb$optimal.width
    res.total[3,k] <- opt.pimcb$mcb.bcr
    res.total[4,k] <- opt.bimcb$optimal.cr
  }
  rownames(res.total) <- c('PIMCB-width', 'BIS-PIMCB-width', 'PIMCB-bcr', 'BIS-PIMCB-bcr')
  colnames(res.total) <- as.character(p)
  out.tab <- rbind(out.tab, res.total)
}

write.csv(out.tab, file='res_alg_wid_cover.csv')


###########################################

library(ggplot2)
file_names <- list.files(path)
time_consume_all <- NULL
for(kk in 1:length(file_names)){
  load(file_names[kk])
  time_pimcb <- unname(sapply(time_a1,function(x) x[1]))
  time_new <- unname(sapply(time_a3,function(x) x[1]))
  time_consume_all <- cbind(time_consume_all, c(log(1+time_pimcb),log(1+time_new)))
}

var_num<-rep(p,10)
methods<-rep(c(rep('PIMCB',length(p)), rep('BIS-PIMCB',length(p))), 5)
var_sele <- c(rep('Adaptive Lasso', 20), rep('Lasso', 20), rep('MCP', 20), rep('SCAD', 20), rep('Scaled Lasso', 20))
df.res <- data.frame(time_consume = as.vector(time_consume_all),var_num=var_num, methods=methods, var_sele=var_sele)
df.res$var_sele <- factor(df.res$var_sele,
                             levels = c("Adaptive Lasso", "SCAD", "MCP", "Lasso", "Scaled Lasso"))

pdf('fig_alg_time_all.pdf', width = 20, height = 8)
ggplot(data=df.res, aes(x=factor(var_num), y=time_consume, color = methods, group = methods)) +
  geom_line(aes(color = methods, linetype=methods), size=1.5) +
  geom_point(aes(color=methods, shape=methods), size=6) +
  labs(x='p',y='Log(1+Time Consumption(s))') +
  scale_x_discrete(breaks = c(100,300,500,700,900), labels = c(100,300,500,700,900))+
  theme_set(theme_bw()) +
  theme(panel.grid =element_blank(),
        legend.title = element_blank(),
        axis.title.y = element_text(size=25),
        axis.text.y = element_text(size=25),
        axis.title.x = element_text(size=25),
        axis.text.x = element_text(size=25),
        legend.text = element_text(size = 25),
        legend.key.size = unit(0.8, "inches"),
        legend.position = "bottom") +
  facet_wrap(~var_sele, scales='free', nrow = 1) +
  theme(strip.text.x = element_text(size = 25),
        strip.text.y = element_text(size = 25),
        plot.title = element_text(size = 25))
dev.off()
  

