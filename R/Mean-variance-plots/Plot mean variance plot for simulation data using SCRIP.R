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
  
  setwd(paste0("E:/DB/Dropbox/Qinfei/Simulation of SC based on splatter/Submition/Bioinformatics/Simulation data/",name))
  
  real.res <- get(load(file="counts.RData"))
  colnames(real.res) <- paste("cell",c(1:ncol(real.res)),sep="")
  
  GP_commonBCV.res <- get(load(file="GP-commonBCV.RData"))
  colnames(GP_commonBCV.res) <- paste("cell",c(1:ncol(GP_commonBCV.res)),sep="")
  
  GP_trendedBCV.res <- get(load(file="GP-trendedBCV.RData"))
  colnames(GP_trendedBCV.res) <- paste("cell",c(1:ncol(GP_trendedBCV.res)),sep="")
  
  BP.res <- get(load(file="BP.RData"))
  colnames(BP.res) <- paste("cell",c(1:ncol(BP.res)),sep="")  
  
  BGP_commonBCV.res <- get(load(file="BGP-commonBCV.RData"))
  colnames(BGP_commonBCV.res) <- paste("cell",c(1:ncol(BGP_commonBCV.res)),sep="")
  
  BGP_trendedBCV.res <- get(load(file="BGP-trendedBCV.RData"))
  colnames(BGP_trendedBCV.res) <- paste("cell",c(1:ncol(BGP_trendedBCV.res)),sep="")
  
  real_mean_variance <- gen_mean_variance(data=real.res, type="real")
  GP_commonBCV_mean_variance <- gen_mean_variance(data=GP_commonBCV.res, type="GP-commonBCV")
  GP_trendedBCV_mean_variance <- gen_mean_variance(GP_trendedBCV.res, type="GP-trendedBCV")
  BP_mean_variance <- gen_mean_variance(BP.res, type="BP")
  BGP_commonBCV_mean_variance <- gen_mean_variance(BGP_commonBCV.res, type="BGP-commonBCV")
  BGP_trendedBCV_mean_variance <- gen_mean_variance(BGP_trendedBCV.res, type="BGP-trendedBCV")
  
  #####  plot  #####
  ##################
  Rol <- loess.smooth(real_mean_variance$mean + runif(nrow(GP_commonBCV_mean_variance),-0.1,0.1),
                      real_mean_variance$variance,
                      degree=2,family = c("gaussian"),span = 1/2)
  GP_commonBCV_ol <- loess.smooth(GP_commonBCV_mean_variance$mean + runif(nrow(GP_commonBCV_mean_variance),-0.1,0.1),
                                  GP_commonBCV_mean_variance$variance,
                                  degree=2,family = c("gaussian"),span = 1/2)
  GP_trendedBCV_ol <- loess.smooth(GP_trendedBCV_mean_variance$mean + runif(nrow(GP_commonBCV_mean_variance),-0.1,0.1),
                                   GP_trendedBCV_mean_variance$variance,
                                    degree=2,family = c("gaussian"),span = 1/2)
  BP_ol <- loess.smooth(BP_mean_variance$mean + runif(nrow(GP_commonBCV_mean_variance),-0.1,0.1),
                        BP_mean_variance$variance,
                        degree=2,family = c("gaussian"),span = 1/2)
  BGP_commonBCV_ol <- loess.smooth(BGP_commonBCV_mean_variance$mean + runif(nrow(GP_commonBCV_mean_variance),-0.1,0.1),
                                    BGP_commonBCV_mean_variance$variance,
                                    degree=2,family = c("gaussian"),span = 1/2)
  BGP_trendedBCV_ol <- loess.smooth(BGP_trendedBCV_mean_variance$mean + runif(nrow(GP_commonBCV_mean_variance),-0.1,0.1),
                                    BGP_trendedBCV_mean_variance$variance,
                                    degree=2,family = c("gaussian"),span = 1/2)
 
  data <- list(real_mean_variance=real_mean_variance, 
               GP_commonBCV_mean_variance=GP_commonBCV_mean_variance, 
               GP_trendedBCV_mean_variance=GP_trendedBCV_mean_variance,
               BP_mean_variance=BP_mean_variance, 
               BGP_commonBCV_mean_variance=BGP_commonBCV_mean_variance, 
               BGP_trendedBCV_mean_variance=BGP_trendedBCV_mean_variance,
                       
               Rol=Rol, GP_commonBCV_ol=GP_commonBCV_ol, GP_trendedBCV_ol=GP_trendedBCV_ol, 
               BP_ol=BP_ol, BGP_commonBCV_ol=BGP_commonBCV_ol, BGP_trendedBCV_ol=BGP_trendedBCV_ol)
  save(data, file=paste0("files.for.mean.variance.RData"))
}


fun.mean.variance(name='MuSiC')
fun.mean.variance(name='Klein')
fun.mean.variance(name='Camp')
fun.mean.variance(name='Tung')
fun.mean.variance(name='Tirosh')
fun.mean.variance(name='Zhou')
fun.mean.variance(name='Seale/Human')
fun.mean.variance(name='Seale/Mice')


### Plotting mean-variance for comparing different methods in SCRIP ######
##########################################################################
##########################################################################
dir <- c("E:/DB/Dropbox/Qinfei/Simulation of SC based on splatter/Submition/Bioinformatics/Simulation data")
setwd(dir)
pdf("Comparison of different methods in SCRIP1.pdf", width=7, height=9)
par(mfrow=c(4,2),mai=c(0.52,0.55,0.28,0.1),mgp=c(2.5, 1, 0))
name.list = c("MuSiC","Klein","Tung","Camp","Tirosh","Zhou","Seale/Human","Seale/Mice")
title.list =  c("Xin","Klein","Tung","Camp","Tirosh","Zhou","Seale(human)","Seale(mouse)")
for (i in 1:8) {
  name <- name.list[i]
  title <- title.list[i]
  if (name=="MuSiC") { 
    xlim=c(1,14); ylim=c(0,12)
  }
  if (name=="Tung") { 
    xlim=c(4,14); ylim=c(0,2)
  }
  if (name=="Klein") { 
    xlim=c(5,14); ylim=c(0,2)
  }
  if (name=="Camp") { 
    xlim=c(0,14); ylim=c(0,18)
  }
  name <- name.list[i]
  if (name=="Tirosh") { 
    xlim=c(6,12); ylim=c(0,3)
  }
  if (name=="Zhou") { 
    xlim=c(8.5,14); ylim=c(0,2)
  }
  if (name=="Seale/Human") { 
    xlim=c(8,15); ylim=c(0,2)
  }
  if (name=="Seale/Mice") { 
    xlim=c(7.5,14); ylim=c(0,2)
  }
  data <- get(load(file=paste0(dir,"/",name,"/files.for.mean.variance.RData")))
  library(scales)
  
  plot(1,lwd=0,type="l",main=title,
       xlab="Mean log2(CPM+1)",ylab="Variance log2(CPM+1)",ylim=ylim,
       xlim=xlim,cex.main=2.0, cex.lab=1.6, cex.axis=1.6,col="red")
  points(data$real_mean_variance$mean, data$real_mean_variance$variance,pch=19,col=alpha("red",0.05),cex=0.01)
  points(data$GP_commonBCV_mean_variance$mean, data$GP_commonBCV_mean_variance$variance,pch=19,col=alpha("gray30",0.05),cex=0.01)
  points(data$GP_trendedBCV_mean_variance$mean, data$GP_trendedBCV_mean_variance$variance,pch=19,col=alpha("forestgreen",0.05),cex=0.01)
  points(data$BP_mean_variance$mean, data$BP_mean_variance$variance,pch=19,col=alpha("gold4",0.05),cex=0.01)
  points(data$BGP_commonBCV_mean_variance$mean, data$BGP_commonBCV_mean_variance$variance,pch=19,col=alpha("violet",0.05),cex=0.01)
  points(data$BGP_trendedBCV_mean_variance$mean, data$BGP_trendedBCV_mean_variance$variance,pch=19,col=alpha("royalblue2",0.05),cex=0.01)
 
  lines(data$Rol$x, data$Rol$y,lwd=2,col="red" )
  lines(data$GP_commonBCV_ol$x, data$GP_commonBCV_ol$y,lwd=2,col="gray30")
  lines(data$GP_trendedBCV_ol$x, data$GP_trendedBCV_ol$y,lwd=2,col="forestgreen" )
  lines(data$BP_ol$x, data$BP_ol$y,lwd=2,col="gold4" )
  lines(data$BGP_commonBCV_ol$x, data$BGP_commonBCV_ol$y,lwd=2,col="violet" )
  lines(data$BGP_trendedBCV_ol$x, data$BGP_trendedBCV_ol$y,lwd=2,col="royalblue2")
  
  
  # legend("topright", legend=c("Real","GP-trendedBCV","GP-commonBCV","BP","BGP-commonBCV","BGP-trendedBCV","dyngen","SymSim"), 
  #        col=c("red","forestgreen","gray30","purple","darkorange","blue","cyan", "brown"), 
  #        lty=1, cex=1.6, lwd = 4, text.font = 2.0)
}

dev.off()
