library("edgeR")
##plot the mean variance plot using CPM
gen_mean_variance <- function(data, type){
  data <- cpm(data,log = T,prior.count = 1)
  mean_variance <- matrix(rep(0,3*nrow(data)),nrow=nrow(data))
  mean_variance[,1] <- apply(data,1,mean)
  mean_variance[,2] <- apply(data,1,var)
  mean_variance <- as.data.frame(mean_variance)
  colnames(mean_variance) <- c("mean","variance","Type") 
  mean_variance[,3] <- type
  return(mean_variance) 
}


##plot the mean-variance plot using log(counts+1)
gen_mean_variance_logcounts <- function(data,type){
  data_CPM <- cpm(data,log=T,prior.count=1)
  mean_variance <- matrix(rep(0,3*nrow(data)),nrow=nrow(data))
  mean_variance[,1] <- apply(log2(data+1),1,mean)
  mean_variance[,2] <- apply(data_CPM,1,var)
  mean_variance <- as.data.frame(mean_variance)
  colnames(mean_variance) <- c("mean","variance","Type") 
  mean_variance[,3] <- type
  return(mean_variance) 
}


gen_sum_cell_counts <- function(data,type){
  sum_cell_counts <- matrix(rep(0,2*ncol(data)),nrow=ncol(data))
  sum_cell_counts[,1] <- apply(data,2,sum)
  sum_cell_counts <- as.data.frame(sum_cell_counts)
  colnames(sum_cell_counts) <- c("Sum_count","Type") 
  sum_cell_counts[,2] <- type
  return(sum_cell_counts) 
}


gen_meancounts_zeropercent <- function(data,type){
  mean_zero <- matrix(rep(0,4*nrow(data)),nrow=nrow(data))
  mean_zero[,1] <- apply(data,1,mean)
  mean_zero[,2] <- apply(data==0,1,sum)
  mean_zero[,3] <- 100*mean_zero[,2]/ncol(data)
  mean_zero <- as.data.frame(mean_zero)
  colnames(mean_zero) <- c("Mean_count","counts_zeros","Percentage_zeros","Type") 
  mean_zero[,4] <- type
  return(mean_zero) 
}


gen_meancounts_zeropercent_percell <- function(data,type){
  mean_zero <- matrix(rep(0,3*ncol(data)),nrow=ncol(data))
  mean_zero[,1] <- apply(data==0,2,sum)
  mean_zero[,2] <- 100*mean_zero[,1]/nrow(data)
  mean_zero <- as.data.frame(mean_zero)
  colnames(mean_zero) <- c("counts_zeros","Percentage_zeros","Type") 
  mean_zero[,3] <- type
  return(mean_zero) 
}




library("edgeR")
##plot the mean variance plot using CPM
gen_mean_variance <- function(data, type){
  data <- cpm(data,log = T,prior.count = 1)
  mean_variance <- matrix(rep(0,3*nrow(data)),nrow=nrow(data))
  mean_variance[,1] <- apply(data,1,mean)
  mean_variance[,2] <- apply(data,1,var)
  mean_variance <- as.data.frame(mean_variance)
  colnames(mean_variance) <- c("mean","variance","Type") 
  mean_variance[,3] <- type
  return(mean_variance) 
}


##plot the mean-variance plot using log(counts+1)
gen_mean_variance_logcounts <- function(data,type){
  data_CPM <- cpm(data,log=T,prior.count=1)
  mean_variance <- matrix(rep(0,3*nrow(data)),nrow=nrow(data))
  mean_variance[,1] <- apply(log2(data+1),1,mean)
  mean_variance[,2] <- apply(data_CPM,1,var)
  mean_variance <- as.data.frame(mean_variance)
  colnames(mean_variance) <- c("mean","variance","Type") 
  mean_variance[,3] <- type
  return(mean_variance) 
}


gen_sum_cell_counts <- function(data,type){
  sum_cell_counts <- matrix(rep(0,2*ncol(data)),nrow=ncol(data))
  sum_cell_counts[,1] <- apply(data,2,sum)
  sum_cell_counts <- as.data.frame(sum_cell_counts)
  colnames(sum_cell_counts) <- c("Sum_count","Type") 
  sum_cell_counts[,2] <- type
  return(sum_cell_counts) 
}


gen_meancounts_zeropercent <- function(data,type){
  mean_zero <- matrix(rep(0,4*nrow(data)),nrow=nrow(data))
  mean_zero[,1] <- apply(data,1,mean)
  mean_zero[,2] <- apply(data==0,1,sum)
  mean_zero[,3] <- 100*mean_zero[,2]/ncol(data)
  mean_zero <- as.data.frame(mean_zero)
  colnames(mean_zero) <- c("Mean_count","counts_zeros","Percentage_zeros","Type") 
  mean_zero[,4] <- type
  return(mean_zero) 
}


gen_meancounts_zeropercent_percell <- function(data,type){
  mean_zero <- matrix(rep(0,3*ncol(data)),nrow=ncol(data))
  mean_zero[,1] <- apply(data==0,2,sum)
  mean_zero[,2] <- 100*mean_zero[,1]/nrow(data)
  mean_zero <- as.data.frame(mean_zero)
  colnames(mean_zero) <- c("counts_zeros","Percentage_zeros","Type") 
  mean_zero[,3] <- type
  return(mean_zero) 
}






Gen_heatmatdata <- function(expre_data_params,simu_GP,simu_BGP, simu_scDesign,simu_scDD,simu_SPARSim,simu_powsimR,simu_dyngen,simu_SymSim){
  
  index <- 1:min(ncol(simu_GP),ncol(simu_BGP),ncol(simu_scDesign),ncol(simu_scDD),ncol(simu_SPARSim),ncol(simu_powsimR),ncol(simu_dyngen),ncol(simu_SymSim))
  expre_data_params <- expre_data_params[,index]
  simu_GP <- simu_GP[,index]
  simu_BGP <- simu_BGP[,index]
  simu_scDesign <- simu_scDesign[,index]
  simu_scDD <- simu_scDD[,index]
  simu_SPARSim <- simu_SPARSim[,index]
  simu_powsimR <- simu_powsimR[,index]
  simu_dyngen <- simu_dyngen[,index]
  simu_SymSim <- simu_SymSim[,index]
  
  
  ##plot the mean variance plot uisng CPM
  real_mean_variance <- gen_mean_variance(data=expre_data_params,type="Real")
  simu_GP_mean_variance <- gen_mean_variance(data=simu_GP,type="GP-trendedBCV")
  simu_BGP_mean_variance <- gen_mean_variance(data=simu_BGP,type="BGP-trendedBCV")
  simu_scDesign_mean_variance <- gen_mean_variance(data=simu_scDesign,type="scDesign")
  simu_scDD_mean_variance <- gen_mean_variance(data=simu_scDD,type="scDD")
  simu_SPARSim_mean_variance <- gen_mean_variance(data=simu_SPARSim,type="SPARSim")
  simu_powsimR_mean_variance <- gen_mean_variance(data=simu_powsimR, type="powsimR")
  simu_dyngen_mean_variance <- gen_mean_variance(data=simu_dyngen,type="dyngen")
  simu_SymSim_mean_variance <- gen_mean_variance(data=simu_SymSim, type="SymSim")
  
  
  final_mean_variance6 <- rbind(real_mean_variance,
                                simu_GP_mean_variance,
                                simu_BGP_mean_variance,
                                simu_scDesign_mean_variance,
                                simu_scDD_mean_variance,
                                simu_SPARSim_mean_variance,
                                simu_powsimR_mean_variance,
                                simu_dyngen_mean_variance,
                                simu_SymSim_mean_variance)
  final_mean_variance6$Type <- factor(final_mean_variance6$Type,levels=c("Real","GP-trendedBCV","BGP-trendedBCV",
                                                                         "scDesign","scDD",
                                                                         "SPARSim","powsimR",
                                                                         "dyngen","SymSim"))
  
  ###describe the mean_rank_diff to real plot####
  ###############################################
  ###############################################
  final_mean_variance <- final_mean_variance6
  final_mean_variance_diffrankmean <- NULL
  for (i in unique(final_mean_variance$Type)){
    data1 <- final_mean_variance[final_mean_variance$Type==i,]
    data1 <- data1[order(data1$mean),]
    real_data <- final_mean_variance[which(final_mean_variance$Type=="Real"),]
    real_data <- real_data[order(real_data$mean),]
    data1$Diff_mean <- abs(log2(pmax(data1$mean,0.1)/pmax(real_data$mean,0.1)))
    final_mean_variance_diffrankmean <- rbind(final_mean_variance_diffrankmean,data1)
  }
  final_mean_variance_diffrankmean <- final_mean_variance_diffrankmean[final_mean_variance_diffrankmean$Type!="Real",]
  
  mean <- aggregate(final_mean_variance_diffrankmean$Diff_mean, list(final_mean_variance_diffrankmean$Type), median)
  mean$Mean <- rank(mean$x)
  rownames(mean) <- mean$Group.1
  mean <- mean[c("GP-trendedBCV","BGP-trendedBCV",
                 "scDesign","scDD",
                 "SPARSim","powsimR","dyngen","SymSim"),c("x","Mean")]
  
  ###describe the variance_rank_diff to real plot
  final_mean_variance <- final_mean_variance6
  final_mean_variance_diffrankvariance <- NULL
  for (i in unique(final_mean_variance$Type)){
    data1 <- final_mean_variance[final_mean_variance$Type==i,]
    data1 <- data1[order(data1$variance),]
    real_data <- final_mean_variance[which(final_mean_variance$Type=="Real"),]
    real_data <- real_data[order(real_data$variance),]
    data1$Diff_variance <-  abs(log2(pmax(data1$variance,0.1)/pmax(real_data$variance,0.1)))
    final_mean_variance_diffrankvariance <- rbind(final_mean_variance_diffrankvariance,data1)
  }
  final_mean_variance_diffrankvariance <- final_mean_variance_diffrankvariance[final_mean_variance_diffrankvariance$Type!="Real",]
  
  variance <- aggregate(final_mean_variance_diffrankvariance$Diff_variance, list(final_mean_variance_diffrankvariance$Type), median)
  variance$Variance <- rank(variance$x)
  rownames(variance) <- variance$Group.1
  variance <- variance[c("GP-trendedBCV","BGP-trendedBCV",
                         "scDesign","scDD",
                         "SPARSim","powsimR","dyngen","SymSim"),c("x","Variance")]
  
  ##plot the library size (sum counts for each cell) boxplot
  real_sum_counts <- gen_sum_cell_counts(data=expre_data_params,type="Real")
  simu_powsimR_sum_counts <- gen_sum_cell_counts(data=simu_powsimR, type="powsimR")
  simu_scDesign_sum_counts <- gen_sum_cell_counts(data=simu_scDesign,type="scDesign")
  simu_SPARSim_sum_counts <- gen_sum_cell_counts(data=simu_SPARSim,type="SPARSim")
  simu_GP_sum_counts <- gen_sum_cell_counts(data=simu_GP,type="GP-trendedBCV")
  simu_BGP_sum_counts <- gen_sum_cell_counts(data=simu_BGP,type="BGP-trendedBCV")
  simu_scDD_sum_counts <- gen_sum_cell_counts(data=simu_scDD,type="scDD")
  simu_dyngen_sum_counts <- gen_sum_cell_counts(data=simu_dyngen,type="dyngen")
  simu_SymSim_sum_counts <- gen_sum_cell_counts(data=simu_SymSim,type="SymSim")
  
  final_sum_counts <- rbind(real_sum_counts,
                            simu_GP_sum_counts,
                            simu_BGP_sum_counts,
                            simu_scDesign_sum_counts,
                            simu_scDD_sum_counts,
                            simu_SPARSim_sum_counts,
                            simu_powsimR_sum_counts,
                            simu_dyngen_sum_counts,
                            simu_SymSim_sum_counts)
  final_sum_counts$Type <- factor(final_sum_counts$Type,levels=c("Real","GP-trendedBCV","BGP-trendedBCV",
                                                                 "scDesign","scDD",
                                                                 "SPARSim","powsimR","dyngen","SymSim"))
  
  
  final_librarysize_diffranklibrarysize <- NULL
  real_data <- final_sum_counts[which(final_sum_counts$Type=="Real"),]
  set.seed(101)
  real_data <- real_data[sample(1:nrow(real_data),min(400,nrow(real_data))),]
  real_data <- real_data[order(real_data$Sum_count),]
  for (i in unique(final_sum_counts$Type)){
    data1 <- final_sum_counts[final_sum_counts$Type==i,]
    data1 <- data1[sample(1:nrow(data1),min(nrow(data1),nrow(real_data))),]
    data1 <- data1[order(data1$Sum_count),]
    data1$diffranklibrarysize <- abs(log2(pmax(data1$Sum_count,1)/pmax(real_data$Sum_count,1)))
    final_librarysize_diffranklibrarysize <- rbind(final_librarysize_diffranklibrarysize,data1)
  }
  final_librarysize_diffranklibrarysize <- final_librarysize_diffranklibrarysize[final_librarysize_diffranklibrarysize$Type!="Real",]
  
  librarysize <- aggregate(final_librarysize_diffranklibrarysize$diffranklibrarysize, list(final_librarysize_diffranklibrarysize$Type), median)
  librarysize$LibrarySize <- rank(librarysize$x)
  rownames(librarysize) <- librarysize$Group.1  
  librarysize <- librarysize[c("GP-trendedBCV","BGP-trendedBCV",
                               "scDesign","scDD",
                               "SPARSim","powsimR","dyngen","SymSim"),c("x","LibrarySize")]
  
  
  ##plot the Percentage zeros per gene  box plot
  real_mean_zero <- gen_meancounts_zeropercent(data=expre_data_params,type="Real")
  simu_powsimR_mean_zero <- gen_meancounts_zeropercent(data=simu_powsimR, type="powsimR")
  simu_scDesign_mean_zero <- gen_meancounts_zeropercent(data=simu_scDesign,type="scDesign")
  simu_SPARSim_mean_zero <- gen_meancounts_zeropercent(data=simu_SPARSim,type="SPARSim")
  simu_GP_mean_zero <- gen_meancounts_zeropercent(data=simu_GP,type="GP-trendedBCV")
  simu_BGP_mean_zero <- gen_meancounts_zeropercent(data=simu_BGP,type="BGP-trendedBCV")
  simu_scDD_mean_zero <- gen_meancounts_zeropercent(data=simu_scDD,type="scDD")
  simu_dyngen_mean_zero <- gen_meancounts_zeropercent(data=simu_dyngen,type="dyngen")
  simu_SymSim_mean_zero <- gen_meancounts_zeropercent(data=simu_SymSim,type="SymSim")
  
  final_mean_zero6 <- rbind(real_mean_zero,
                            simu_GP_mean_zero,
                            simu_BGP_mean_zero,
                            simu_scDesign_mean_zero,
                            simu_scDD_mean_zero,
                            simu_SPARSim_mean_zero,
                            simu_powsimR_mean_zero,
                            simu_dyngen_mean_zero,
                            simu_SymSim_mean_zero)
  final_mean_zero6$Type <- factor(final_mean_zero6$Type,levels=c("Real","GP-trendedBCV","BGP-trendedBCV",
                                                                 "scDesign","scDD",
                                                                 "SPARSim","powsimR","dyngen","SymSim"))
  
  
  
  ##plot the Percentage zeros per cell
  real_percentage_zero_percell <- gen_meancounts_zeropercent_percell(data=expre_data_params,type="Real")
  simu_powsimR_percentage_zero_percell <- gen_meancounts_zeropercent_percell(data=simu_powsimR, type="powsimR")
  simu_scDesign_percentage_zero_percell <- gen_meancounts_zeropercent_percell(data=simu_scDesign,type="scDesign")
  simu_SPARSim_percentage_zero_percell <- gen_meancounts_zeropercent_percell(data=simu_SPARSim,type="SPARSim")
  simu_GP_percentage_zero_percell <- gen_meancounts_zeropercent_percell(data=simu_GP,type="GP-trendedBCV")
  simu_BGP_percentage_zero_percell <- gen_meancounts_zeropercent_percell(data=simu_BGP,type="BGP-trendedBCV")
  simu_scDD_percentage_zero_percell <- gen_meancounts_zeropercent_percell(data=simu_scDD,type="scDD")
  simu_dyngen_percentage_zero_percell <- gen_meancounts_zeropercent_percell(data=simu_dyngen,type="dyngen")
  simu_SymSim_percentage_zero_percell <- gen_meancounts_zeropercent_percell(data=simu_SymSim,type="SymSim")
  
  
  final_percentage_zero_percell <- rbind(real_percentage_zero_percell,
                                         simu_GP_percentage_zero_percell,
                                         simu_BGP_percentage_zero_percell,
                                         simu_scDesign_percentage_zero_percell,
                                         simu_scDD_percentage_zero_percell,
                                         simu_SPARSim_percentage_zero_percell,
                                         simu_powsimR_percentage_zero_percell,
                                         simu_dyngen_percentage_zero_percell,
                                         simu_SymSim_percentage_zero_percell)
  final_percentage_zero_percell$Type <- factor(final_percentage_zero_percell$Type,levels=c("Real","GP-trendedBCV","BGP-trendedBCV",
                                                                                           "scDesign","scDD",
                                                                                           "SPARSim","powsimR","dyngen","SymSim"))
  
  ###plot the zero_diff to real vs rank_counts plot per gene
  gen_diffzero_pergene <- function(final_mean_zero){
    final_mean_zero_rankcount <- NULL
    for (i in unique(final_mean_zero$Type)){
      data1 <- final_mean_zero[final_mean_zero$Type==i,]
      # data1$rank <- rank(data1$Mean_count)
      data1 <- data1[order(data1$Mean_count),]
      real_data <- final_mean_zero[which(final_mean_zero$Type=="Real"),]
      real_data <- real_data[order(real_data$Mean_count),]
      data1$Diff_zero <- abs(log2(pmax(data1$Percentage_zeros,0.1)/pmax(real_data$Percentage_zeros,0.1)))
      final_mean_zero_rankcount <- rbind(final_mean_zero_rankcount,data1)
    }
    final_mean_zero_rankcount <- final_mean_zero_rankcount[final_mean_zero_rankcount$Type!="Real",]
    return(final_mean_zero_rankcount)
  }
  
  final_diff_zero_pergene <- gen_diffzero_pergene(final_mean_zero=final_mean_zero6)
  Diff_zero_pergene <- aggregate(final_diff_zero_pergene$Diff_zero, list(final_diff_zero_pergene$Type), median)
  Diff_zero_pergene$Zeros_Gene <- rank(Diff_zero_pergene$x)
  rownames(Diff_zero_pergene) <- Diff_zero_pergene$Group.1
  Diff_zero_pergene <- Diff_zero_pergene[c("GP-trendedBCV","BGP-trendedBCV",
                                           "scDesign","scDD",
                                           "SPARSim","powsimR","dyngen","SymSim"),c("x","Zeros_Gene")]
  
  
  
  ###plot the zero_diff to real vs rank_counts plot per cell
  gen_diffzero_percell <- function(final_mean_zero){
    final_mean_zero_rankcount <- NULL
    for (i in unique(final_mean_zero$Type)){
      data1 <- final_mean_zero[final_mean_zero$Type==i,]
      data1 <- data1[order(data1$Percentage_zeros),]
      real_data <- final_mean_zero[which(final_mean_zero$Type=="Real"),]
      real_data <- real_data[order(real_data$Percentage_zeros),]
      data2 <- as.data.frame(abs(log2((pmax(data1$Percentage_zeros[min(nrow(data1),nrow(real_data))],0.1)/pmax(real_data$Percentage_zeros[min(nrow(data1),nrow(real_data))],0.1)))))
      colnames(data2) <- "Diff_zero"
      data2$Type <- i
      final_mean_zero_rankcount <- rbind(final_mean_zero_rankcount,data2)
    }
    final_mean_zero_rankcount <- final_mean_zero_rankcount[final_mean_zero_rankcount$Type!="Real",]
    
    
    return(final_mean_zero_rankcount)
  }
  
  final_diff_zero_percell <- gen_diffzero_percell(final_mean_zero=final_percentage_zero_percell)
  
  Diff_zero_percell <- aggregate(final_diff_zero_percell$Diff_zero, list(final_diff_zero_percell$Type), median)
  Diff_zero_percell$Zeros_Cell <- rank(Diff_zero_percell$x)
  rownames(Diff_zero_percell) <- Diff_zero_percell$Group.1
  Diff_zero_percell <- Diff_zero_percell[c("GP-trendedBCV","BGP-trendedBCV",
                                           "scDesign","scDD",
                                           "SPARSim","powsimR","dyngen","SymSim"),c("x","Zeros_Cell")]
  
  
  Heatmap_data <- cbind(mean,variance,librarysize,Diff_zero_pergene,Diff_zero_percell)
  Heatmap_data <- Heatmap_data[,c(1:5)*2-1]
  colnames(Heatmap_data) <- c("Mean","Variance","LibrarySize","Zeros_Gene","Zeros_Cell")
  library(BBmisc)
  Heatmap_data <- t(Heatmap_data)
  # Heatmap_data <- t(scale(Heatmap_data))
  # Heatmap_data <- t(Heatmap_data[,c("Mean","Variance","LibrarySize","Zeros_Gene","Zeros_Cell")])
  return(Heatmap_data)
}






Fun_heatmap_data <- function(name){
  
  setwd(paste0("E:\\DB\\Dropbox\\Qinfei\\Simulation of SC based on splatter\\Submition\\Bioinformatics\\Simulation data\\", name))
  expre_data_params <- get(load(file="counts.RData"))
  simu_GP <- get(load(file="SCRIP.RData"))
  simu_BGP <- get(load(file="BGP-trendedBCV.RData"))  
  simu_scDesign <- get(load(file="scDesign.RData"))  
  simu_scDD <- get(load(file="scDD.RData"))  
  simu_SPARSim <- get(load(file="SPARSim.RData"))  
  simu_powsimR <- get(load(file="powsimR.RData"))  
  simu_SymSim <- get(load(file="SymSim.RData"))  
  simu_dyngen <- get(load(file="dyngen.RData"))  
  
  heatmap_data <- Gen_heatmatdata(expre_data_params=expre_data_params,
                                  simu_GP=simu_GP,
                                  simu_BGP=simu_BGP,
                                  simu_scDesign=simu_scDesign,
                                  simu_scDD=simu_scDD,
                                  simu_SPARSim=simu_SPARSim,
                                  simu_powsimR=simu_powsimR, 
                                  simu_SymSim=simu_SymSim,
                                  simu_dyngen=simu_dyngen+0.01)
  save(heatmap_data, file="heatmap_data.RData")
  return(heatmap_data)
}






Heatmap_data_MuSiC <- Fun_heatmap_data(name="MuSiC")
Heatmap_data_Camp <- Fun_heatmap_data(name="Camp")
Heatmap_data_Klein <- Fun_heatmap_data(name="Klein")
Heatmap_data_Tung <- Fun_heatmap_data(name="Tung")
Heatmap_data_Tirosh <- Fun_heatmap_data(name="Tirosh")
Heatmap_data_Zhou <- Fun_heatmap_data(name="Zhou")
Heatmap_data_SealeH <- Fun_heatmap_data(name="Seale/Human")
Heatmap_data_SealeM <- Fun_heatmap_data(name="Seale/Mice")

dir <- c("E:\\DB\\Dropbox\\Qinfei\\Simulation of SC based on splatter\\Submition\\Bioinformatics\\Simulation data")

Heatmap_data_MuSiC <- get(load(paste0(dir,"\\MuSiC\\heatmap_data.RData")))
Heatmap_data_Camp <- get(load(paste0(dir,"\\Camp\\heatmap_data.RData")))
Heatmap_data_Klein <- get(load(paste0(dir,"\\Klein\\heatmap_data.RData")))
Heatmap_data_Tung <- get(load(paste0(dir,"\\Tung\\heatmap_data.RData")))
Heatmap_data_Tirosh <- get(load(paste0(dir,"\\Tirosh\\heatmap_data.RData")))
Heatmap_data_Zhou <- get(load(paste0(dir,"\\Zhou\\heatmap_data.RData")))
Heatmap_data_SealeH <- get(load(paste0(dir,"\\Seale\\Human\\heatmap_data.RData")))
Heatmap_data_SealeM <- get(load(paste0(dir,"\\Seale\\Mice\\heatmap_data.RData")))


Heatmap_data <- rbind(Heatmap_data_MuSiC,Heatmap_data_Klein,Heatmap_data_Tung,Heatmap_data_Camp,
                      Heatmap_data_Tirosh,Heatmap_data_Zhou,Heatmap_data_SealeH,Heatmap_data_SealeM)
Heatmap_data <- Heatmap_data[,-c(4,7)]
rownames(Heatmap_data) <- paste("measure",1:40,sep="")
# Heatmap_data<- t(scale(t(Heatmap_data),scale = F))
cols = colorRampPalette(c("royalblue3","white"))(30)

library(pheatmap)
data <- Heatmap_data[c((0:7)*5+1,(0:7)*5+2,(0:7)*5+3,(0:7)*5+4,(0:7)*5+5),]




features <- data.frame(Features=rep(c("Mean","Variance","library size","Zeros per cell","Zeros per gene"),each=8))
rownames(features) <- rownames(data)
features$Features <- factor(features$Features, levels=c("Mean","Variance","library size","Zeros per cell","Zeros per gene"))

setwd("E:\\DB\\Dropbox\\Qinfei\\Simulation of SC based on splatter\\Submition\\Bioinformatics\\Simulation data")
pdf("Heatmap comparion other simulators for eight datasets log2 FC3.pdf",width=15,height=17)
pheatmap(data,color=cols,cluster_rows=F,cluster_cols=F,
         fontsize = 25,angle_col=45,annotation_row=features,
         labels_row = rep(c("Xin","Klein","Tung","Camp","Tirosh","Zhou","Seale(Human)","Seale(Mice)"),5),cellheight=25, cellwidth=45)
dev.off()






data1 <- data.frame(value=c(as.vector(data[1:8,1:ncol(data)]),as.vector(data[9:16,1:ncol(data)]),
                            as.vector(data[17:24,1:ncol(data)]), as.vector(data[25:32,1:ncol(data)]),
                            as.vector(data[33:40,1:ncol(data)])), 
                    type=rep(c("Mean","Variance","library size","Zeros per cell","Zeros per gene"),each=48))
data1$type <- factor(data1$type, levels = c(c("Mean","Variance","library size","Zeros per cell","Zeros per gene")))



library("ggplot2")
setwd("E:\\DB\\Dropbox\\Qinfei\\Simulation of SC based on splatter\\Submition\\Bioinformatics\\Simulation data")
pdf("other simulators boxplot for different features.pdf",height=8,width=6)
ggplot(data = data1, aes(x=type,y=value)) +
  labs(y="MAD", x = "Type")+
  geom_boxplot(lwd=1,fatten=1)+
  scale_y_continuous(limits=c(0,3.5)) +
  geom_point(size = 1) +
  theme(legend.position="none",
        plot.title = element_text(face="bold", size = 20, hjust = 0.5),
        axis.text.x =element_text(size=25,angle=30, hjust=1,colour="black"),
        axis.text.y =element_text(size=30,colour="black"),
        axis.title.x = element_text(size=30,face="bold",margin = margin(t = 0, r = 20, b = 0, l = 15)),
        axis.title.y = element_text(size=30,face="bold",margin = margin(t = 0, r = 20, b = 0, l = 15)),
        panel.background = element_rect(fill='white',colour="black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) 

dev.off()



data2 <- data[c(9:16),]
setwd("E:\\DB\\Dropbox\\Qinfei\\Simulation of SC based on splatter\\Submition\\Bioinformatics\\Simulation data")
pdf("Heatmap comparion other simulators for variance using eight datasets log2 FC3.pdf",width=7,height=6)
pheatmap(data2,color=cols,cluster_rows=F,cluster_cols=F,
         fontsize = 15,angle_col=45,
         labels_row = rep(c("Xin","Klein","Tung","Camp","Tirosh","Zhou","Seale(Human)","Seale(Mice)"),1),cellheight=35, cellwidth=35)
dev.off()


