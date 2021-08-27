# source("E:\\DB\\Dropbox\\Qinfei\\Simulation of SC based on splatter\\Code\\BCV plots.R")

##plot the mean variance plot using CPM
gen_mean_variance <- function(data,type){
  data <- edgeR::cpm(data,log = T,prior.count = 1)
  mean_variance <- matrix(rep(0,3*nrow(data)),nrow=nrow(data))
  mean_variance[,1] <- apply(data,1,mean)
  mean_variance[,2] <- apply(data,1,var)
  mean_variance <- as.data.frame(mean_variance)
  colnames(mean_variance) <- c("mean","variance","Type") 
  mean_variance[,3] <- type
  return(mean_variance) 
}


## Generate meanv-variance files to prepare for mean-variance plots ##
fun.mean.variance <- function(name){

  setwd(paste0("E:\\DB\\Dropbox\\Qinfei\\Simulation of SC based on splatter\\Submition\\Bioinformatics\\Simulation data\\",name))

  real.res <- get(load(file="counts.RData"))
  colnames(real.res) <- paste("cell",c(1:ncol(real.res)),sep="")
  
  GP.res <- get(load(file="SCRIP.RData"))
  colnames(GP.res) <- paste("GP-trendedBCV",c(1:ncol(GP.res)),sep="")
  
  BGP.res <- get(load(file="BGP-trendedBCV.RData"))
  colnames(BGP.res) <- paste("BGP-trendedBCV",c(1:ncol(BGP.res)),sep="")
  
  scDesign.res <- get(load(file="scDesign.RData"))
  colnames(scDesign.res) <- paste("cell",c(1:ncol(scDesign.res)),sep="")

  scDD.res <- get(load(file="scDD.RData"))
  colnames(scDD.res) <- paste("scDD",c(1:ncol(scDD.res)),sep="")

  powsimR.res <- get(load(file="powsimR.RData"))
  colnames(powsimR.res) <- paste("powsimR",c(1:ncol(powsimR.res)),sep="")
  
  SPARSim.res <- get(load(file="SPARSim.RData"))
  colnames(SPARSim.res) <- paste("SPARSim",c(1:ncol(SPARSim.res)),sep="")
  
  dyngen.res <- get(load(file="dyngen.RData"))
  colnames(dyngen.res) <- paste("dyngen",c(1:ncol(dyngen.res)),sep="")
  
  SymSim.res <- get(load(file="SymSim.RData"))
  colnames(SymSim.res) <- paste("SymSim",c(1:ncol(SymSim.res)),sep="")

  ##plot the mean variance plot using CPM
  gen_mean_variance <- function(data,type){
    data <- edgeR::cpm(data,log = T,prior.count = 1)
    mean_variance <- matrix(rep(0,3*nrow(data)),nrow=nrow(data))
    mean_variance[,1] <- apply(data,1,mean)
    mean_variance[,2] <- apply(data,1,var)
    mean_variance <- as.data.frame(mean_variance)
    colnames(mean_variance) <- c("mean","variance","Type") 
    mean_variance[,3] <- type
    return(mean_variance) 
  }
  
  
  real_mean_variance <- gen_mean_variance(data=real.res, type="real")
  scDesign_mean_variance <- gen_mean_variance(data=scDesign.res, type="scDesign")
  GP_mean_variance <- gen_mean_variance(GP.res, type="GP-trendedBCV")
  BGP_mean_variance <- gen_mean_variance(BGP.res, type="BGP-trendedBCV")
  powsimR_mean_variance <- gen_mean_variance(powsimR.res, type="powsimR")
  SPARSim_mean_variance <- gen_mean_variance(SPARSim.res, type="SPARSim")
  scDD_mean_variance <- gen_mean_variance(scDD.res, type="scDD")
  dyngen_mean_variance <- gen_mean_variance(dyngen.res+0.01, type="dyngen")
  SymSim_mean_variance <- gen_mean_variance(SymSim.res, type="SymSim")
  

  #####  plot  #####
  ##################
  Rol <- loess.smooth(real_mean_variance$mean + runif(nrow(real_mean_variance),-0.1,0.1),
                      real_mean_variance$variance + runif(nrow(real_mean_variance),-0.1,0.1),
                    degree=2,family = c("gaussian"),span = 1/2)
  scDesign_ol <- loess.smooth(scDesign_mean_variance$mean + runif(nrow(real_mean_variance),-0.1,0.1),
                              scDesign_mean_variance$variance + runif(nrow(real_mean_variance),-0.1,0.1),
                       degree=2,family = c("gaussian"),span = 1/2)
  GP_ol <- loess.smooth(GP_mean_variance$mean + runif(nrow(real_mean_variance),-0.1,0.1),
                        GP_mean_variance$variance + runif(nrow(real_mean_variance),-0.1,0.1),
                            degree=2,family = c("gaussian"),span = 1/2)
  BGP_ol <- loess.smooth(BGP_mean_variance$mean + runif(nrow(real_mean_variance),-0.1,0.1),
                         BGP_mean_variance$variance + runif(nrow(real_mean_variance),-0.1,0.1),
                           degree=2,family = c("gaussian"),span = 1/2)
  powsimR_ol <- loess.smooth(powsimR_mean_variance$mean + runif(nrow(real_mean_variance),-0.1,0.1),
                             powsimR_mean_variance$variance + runif(nrow(real_mean_variance),-0.1,0.1),
                         degree=2,family = c("gaussian"),span = 1/2)
  SPARSim_ol <- loess.smooth(SPARSim_mean_variance$mean + runif(nrow(real_mean_variance),-0.1,0.1),
                             SPARSim_mean_variance$variance + runif(nrow(real_mean_variance),-0.1,0.1),
                             degree=2,family = c("gaussian"),span = 1/2)
  scDD_ol <- loess.smooth(scDD_mean_variance$mean + runif(nrow(real_mean_variance),-0.1,0.1),
                          scDD_mean_variance$variance + runif(nrow(real_mean_variance),-0.1,0.1),
                         degree=2,family = c("gaussian"),span = 1/2)
  dyngen_ol <- loess.smooth(dyngen_mean_variance$mean + runif(nrow(real_mean_variance),-0.1,0.1),
                            dyngen_mean_variance$variance + runif(nrow(real_mean_variance),-0.1,0.1),
                             degree=2,family = c("gaussian"),span = 1/2)
  SymSim_ol <- loess.smooth(SymSim_mean_variance$mean + runif(nrow(real_mean_variance),-0.1,0.1),
                            SymSim_mean_variance$variance + runif(nrow(real_mean_variance),-0.1,0.1),
                          degree=2,family = c("gaussian"),span = 1/2)
  
  data <- list(real_mean_variance=real_mean_variance, 
               scDesign_mean_variance=scDesign_mean_variance, 
               GP_mean_variance=GP_mean_variance,
               BGP_mean_variance=BGP_mean_variance,
               powsimR_mean_variance=powsimR_mean_variance, 
               SPARSim_mean_variance=SPARSim_mean_variance, 
               scDD_mean_variance=scDD_mean_variance,
               dyngen_mean_variance=dyngen_mean_variance, 
               SymSim_mean_variance=SymSim_mean_variance,
      
               Rol=Rol, scDesign_ol=scDesign_ol, GP_ol=GP_ol, BGP_ol=BGP_ol,
               powsimR_ol=powsimR_ol, SPARSim_ol=SPARSim_ol, scDD_ol=scDD_ol,
               dyngen_ol=dyngen_ol, SymSim_ol=SymSim_ol)
  save(data, file=paste0("files.for.mean.variance.other.simultors.RData"))
}


fun.mean.variance(name='MuSiC')
fun.mean.variance(name='Klein')
fun.mean.variance(name='Camp')
fun.mean.variance(name='Tung')
fun.mean.variance(name='Tirosh')
fun.mean.variance(name='Zhou')
fun.mean.variance(name='Seale/Human')
fun.mean.variance(name='Seale/Mice')


### Plotting mean variance plots for comparing other simulators ###
###################################################################
###################################################################
dir <- c("E:\\DB\\Dropbox\\Qinfei\\Simulation of SC based on splatter\\Submition\\Bioinformatics\\Simulation data")
setwd(dir)
pdf("Comparison to other simulators2.pdf", width=7, height=9)
par(mfrow=c(4,2),mai=c(0.52,0.55,0.28,0.1),mgp=c(2.5, 1, 0))
name.list = c("MuSiC","Klein","Tung","Camp","Tirosh","Zhou","Seale/Human","Seale/Mice")
title.list =  c("Xin","Klein","Tung","Camp","Tirosh","Zhou","Seale(human)","Seale(mouse)")
for (i in 5:8) {
  name <- name.list[i]
  title <- title.list[i]
  if (name=="MuSiC") { 
    xlim=c(0,14); ylim=c(0,9)
  }
  if (name=="Tung") { 
    xlim=c(2,14); ylim=c(0,2)
  }
  if (name=="Klein") { 
    xlim=c(4,14); ylim=c(0,2)
  }
  if (name=="Camp") { 
    xlim=c(0,14); ylim=c(0,18)
  }
  if (name=="Tirosh") { 
    xlim=c(6,11); ylim=c(0,3)
  }
  if (name=="Zhou") { 
    xlim=c(6,14); ylim=c(0,2)
  }
  if (name=="Seale/Human") { 
    xlim=c(5,14); ylim=c(0,2)
  }
  if (name=="Seale/Mice") { 
    xlim=c(5,14); ylim=c(0,2)
  }
  
  data <- get(load(file=paste0(dir,"\\",name,"\\files.for.mean.variance.other.simultors.RData")))
  library(scales)
  
  plot(1,lwd=0,type="l",main=title,
       xlab="Mean log2(CPM+1)",ylab="Variance log2(CPM+1)",ylim=ylim,
       xlim=xlim,cex.main=2.0, cex.lab=1.8, cex.axis=1.8,col="red")
  points(data$real_mean_variance$mean, data$real_mean_variance$variance,pch=19,col=alpha("red",0.05),cex=0.01)
  points(data$scDesign_mean_variance$mean, data$scDesign_mean_variance$variance,pch=19,col=alpha("gray30",0.05),cex=0.01)
  points(data$GP_mean_variance$mean, data$GP_mean_variance$variance,pch=19,col=alpha("forestgreen",0.05),cex=0.01)
  points(data$BGP_mean_variance$mean, data$BGP_mean_variance$variance,pch=19,col=alpha("royalblue2",0.05),cex=0.01)
  points(data$powsimR_mean_variance$mean, data$powsimR_mean_variance$variance,pch=19,col=alpha("purple",0.05),cex=0.01)
  points(data$SPARSim_mean_variance$mean, data$SPARSim_mean_variance$variance,pch=19,col=alpha("darkorange",0.05),cex=0.01)
  points(data$scDD_mean_variance$mean, data$scDD_mean_variance$variance,pch=19,col=alpha("khaki3",0.05),cex=0.01)
  points(data$dyngen_mean_variance$mean, data$dyngen_mean_variance$variance,pch=19,col=alpha("cyan",0.05),cex=0.01)
  points(data$SymSim_mean_variance$mean, data$SymSim_mean_variance$variance,pch=19,col=alpha("brown",0.05),cex=0.01)
  
  lines(data$Rol$x, data$Rol$y,lwd=2,col="red" )
  lines(data$scDesign_ol$x, data$scDesign_ol$y,lwd=2,col="gray30")
  lines(data$GP_ol$x, data$GP_ol$y,lwd=2,col="forestgreen" )
  lines(data$BGP_ol$x, data$BGP_ol$y,lwd=2,col="royalblue2" )
  lines(data$powsimR_ol$x, data$powsimR_ol$y,lwd=2,col="purple" )
  lines(data$SPARSim_ol$x, data$SPARSim_ol$y,lwd=2,col="darkorange" )
  lines(data$scDD_ol$x, data$scDD_ol$y,lwd=2,col="khaki3")
  lines(data$dyngen_ol$x, data$dyngen_ol$y,lwd=2,col="cyan" )
  lines(data$SymSim_ol$x, data$SymSim_ol$y,lwd=2,col="brown")
  

  # legend("topright", legend=c("Real","GP","BGP","scDesign","powsimR","SPARSim","scDD","dyngen","SymSim"),
  #        col=c("red","forestgreen","royalblue2","gray30","purple","darkorange","khaki3","cyan", "brown"),
  #        lty=1, cex=1.6, lwd = 2, text.font = 2.0)
}

dev.off()
  

