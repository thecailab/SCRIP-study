## Below are R functions and detailed analysis steps to evaluate simulation data from SCRIP using MAD with five characteristics 
## (i.e., gene-wise expression mean, variance, cell-wise zero proportion, gene-wise zero proportion, and library size) in
## eight real datasets (Xin, Klein, Tung, Camp, Tirosh, Zhou, Seale human, Seale mouse).


##plot the mean variance plot using CPM
gen_mean_variance <- function(data,type){
  library(edgeR)
  data <- cpm(data,log = T,prior.count = 1)
  mean_variance <- matrix(rep(0,3*nrow(data)),nrow=nrow(data))
  mean_variance[,1] <- apply(data,1,mean)
  mean_variance[,2] <- apply(data,1,var)
  mean_variance <- as.data.frame(mean_variance)
  colnames(mean_variance) <- c("mean","variance","Type") 
  mean_variance[,3] <- type
  return(mean_variance) 
}

## Calculate library size
gen_sum_cell_counts <- function(data,type){
  sum_cell_counts <- matrix(rep(0,2*ncol(data)),nrow=ncol(data))
  sum_cell_counts[,1] <- apply(data,2,sum)
  sum_cell_counts <- as.data.frame(sum_cell_counts)
  colnames(sum_cell_counts) <- c("Sum_count","Type") 
  sum_cell_counts[,2] <- type
  return(sum_cell_counts) 
}


## Calculate zero proportion per gene
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


## Calculate zero proportion per cell
gen_meancounts_zeropercent_percell <- function(data,type){
  mean_zero <- matrix(rep(0,3*ncol(data)),nrow=ncol(data))
  mean_zero[,1] <- apply(data==0,2,sum)
  mean_zero[,2] <- 100*mean_zero[,1]/nrow(data)
  mean_zero <- as.data.frame(mean_zero)
  colnames(mean_zero) <- c("counts_zeros","Percentage_zeros","Type") 
  mean_zero[,3] <- type
  return(mean_zero) 
}




Gen_heatmatdata <- function(res){
  expre_data_params <- res$counts
  nobursting <- res$GP_trendedBCV
  bursting <- res$BGP_trendedBCV
  simu_common_burst <- res$BGP_commonBCV
  simu_common_noburst <- res$GP_commonBCV
  simu_noBCV_burst <- res$BP
  ##plot the mean variance plot uisng CPM
  real_mean_variance <- gen_mean_variance(data=expre_data_params,type="Real")
  simu_common_noburstring_mean_variance <- gen_mean_variance(data=simu_common_noburst,type="GP_commonBCV")
  simu_nobursting_mean_variance <- gen_mean_variance(data=nobursting,type="GP_trendedBCV")
  simu_noBCV_burst_mean_variance <- gen_mean_variance(data=simu_noBCV_burst,type="BP")
  simu_common_burstring_mean_variance <- gen_mean_variance(data=simu_common_burst,type="BGP_commonBCV")
  simu_bursting_mean_variance <- gen_mean_variance(data=bursting, type="BGP_trendedBCV")
  
  final_mean_variance6 <- rbind(real_mean_variance,
                                simu_noBCV_burst_mean_variance,
                                simu_common_noburstring_mean_variance,
                                simu_nobursting_mean_variance,
                                simu_common_burstring_mean_variance,
                                simu_bursting_mean_variance)
  final_mean_variance6$Type <- factor(final_mean_variance6$Type,levels=c("Real","GP_commonBCV",
                                                                         "GP_trendedBCV","BP",
                                                                         "BGP_commonBCV","BGP_trendedBCV"))
  
  
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
    data1$Diff_mean <- abs(log2(pmax(data1$mean,0.0000001)/pmax(real_data$mean,0.0000001)))
    final_mean_variance_diffrankmean <- rbind(final_mean_variance_diffrankmean,data1)
  }
  final_mean_variance_diffrankmean <- final_mean_variance_diffrankmean[final_mean_variance_diffrankmean$Type!="Real",]
  
  mean <- aggregate(final_mean_variance_diffrankmean$Diff_mean, list(final_mean_variance_diffrankmean$Type), median)
  mean$Mean <- rank(mean$x)
  rownames(mean) <- mean$Group.1
  mean <- mean[c("GP_commonBCV",
                 "GP_trendedBCV","BP",
                 "BGP_commonBCV","BGP_trendedBCV"),c("x","Mean")]
  
  ###describe the variance_rank_diff to real plot
  final_mean_variance <- final_mean_variance6
  final_mean_variance_diffrankvariance <- NULL
  for (i in unique(final_mean_variance$Type)){
    data1 <- final_mean_variance[final_mean_variance$Type==i,]
    data1 <- data1[order(data1$variance),]
    real_data <- final_mean_variance[which(final_mean_variance$Type=="Real"),]
    real_data <- real_data[order(real_data$variance),]
    data1$Diff_variance <-  abs(log2(pmax(data1$variance,0.0000001)/pmax(real_data$variance,0.0000001)))
    final_mean_variance_diffrankvariance <- rbind(final_mean_variance_diffrankvariance,data1)
  }
  final_mean_variance_diffrankvariance <- final_mean_variance_diffrankvariance[final_mean_variance_diffrankvariance$Type!="Real",]
  
  variance <- aggregate(final_mean_variance_diffrankvariance$Diff_variance, list(final_mean_variance_diffrankvariance$Type), median)
  variance$Variance <- rank(variance$x)
  rownames(variance) <- variance$Group.1
  variance <- variance[c("GP_commonBCV",
                         "GP_trendedBCV","BP",
                         "BGP_commonBCV","BGP_trendedBCV"),c("x","Variance")]

  
  ##plot the library size (sum counts for each cell) boxplot
  real_sum_counts <- gen_sum_cell_counts(data=expre_data_params,type="Real")
  simu_bursting_sum_counts <- gen_sum_cell_counts(data=bursting, type="BGP_trendedBCV")
  simu_nobursting_sum_counts <- gen_sum_cell_counts(data=nobursting,type="GP_trendedBCV")
  simu_common_bursting_sum_counts <- gen_sum_cell_counts(data=simu_common_burst,type="BGP_commonBCV")
  simu_common_nobursting_sum_counts <- gen_sum_cell_counts(data=simu_common_noburst,type="GP_commonBCV")
  simu_noBCV_burst_sum_counts <- gen_sum_cell_counts(data=simu_noBCV_burst,type="BP")
  
  final_sum_counts <- rbind(real_sum_counts,
                            simu_noBCV_burst_sum_counts,
                            simu_common_nobursting_sum_counts,
                            simu_nobursting_sum_counts,
                            simu_common_bursting_sum_counts,
                            simu_bursting_sum_counts)
  final_sum_counts$Type <- factor(final_sum_counts$Type,levels=c("Real","GP_commonBCV",
                                                                 "GP_trendedBCV","BP",
                                                                 "BGP_commonBCV","BGP_trendedBCV"))
  
  
  final_librarysize_diffranklibrarysize <- NULL
  real_data <- final_sum_counts[which(final_sum_counts$Type=="Real"),]
  # set.seed(101)
  # real_data <- real_data[sample(1:nrow(real_data),min(400,nrow(real_data))),]
  real_data <- real_data[order(real_data$Sum_count),]
  for (i in unique(final_sum_counts$Type)){
    data1 <- final_sum_counts[final_sum_counts$Type==i,]
    # data1 <- data1[sample(1:nrow(data1),min(nrow(data1),nrow(real_data))),]
    data1 <- data1[order(data1$Sum_count),]
    data1$diffranklibrarysize <- abs(log2(pmax(data1$Sum_count,0.000001)/pmax(real_data$Sum_count,0.000001)))
    final_librarysize_diffranklibrarysize <- rbind(final_librarysize_diffranklibrarysize,data1)
  }
  
  final_librarysize_diffranklibrarysize <- final_librarysize_diffranklibrarysize[final_librarysize_diffranklibrarysize$Type!="Real",]
  
  librarysize <- aggregate(final_librarysize_diffranklibrarysize$diffranklibrarysize, list(final_librarysize_diffranklibrarysize$Type), median)
  librarysize$LibrarySize <- rank(librarysize$x)
  rownames(librarysize) <- librarysize$Group.1
  librarysize <- librarysize[c("GP_commonBCV",
                               "GP_trendedBCV","BP",
                               "BGP_commonBCV","BGP_trendedBCV"),c("x","LibrarySize")]
  
  ##plot the Percentage zeros per gene  box plot
  real_mean_zero <- gen_meancounts_zeropercent(data=expre_data_params,type="Real")
  simu_bursting_mean_zero <- gen_meancounts_zeropercent(data=bursting, type="BGP_trendedBCV")
  simu_nobursting_mean_zero <- gen_meancounts_zeropercent(data=nobursting,type="GP_trendedBCV")
  simu_common_burstring_mean_zero <- gen_meancounts_zeropercent(data=simu_common_burst,type="BGP_commonBCV")
  simu_common_noburstring_mean_zero <- gen_meancounts_zeropercent(data=simu_common_noburst,type="GP_commonBCV")
  simu_noBCV_burst_mean_zero <- gen_meancounts_zeropercent(data=simu_noBCV_burst,type="BP")
  
  final_mean_zero <- rbind(real_mean_zero,
                            simu_noBCV_burst_mean_zero,
                            simu_common_noburstring_mean_zero,
                            simu_nobursting_mean_zero,
                            simu_common_burstring_mean_zero,
                            simu_bursting_mean_zero)
  final_mean_zero$Type <- factor(final_mean_zero$Type,levels=c("Real","GP_commonBCV",
                                                                 "GP_trendedBCV","BP",
                                                                 "BGP_commonBCV","BGP_trendedBCV"))
  
  
  
  ##plot the Percentage zeros per cell
  real_percentage_zero_percell <- gen_meancounts_zeropercent_percell(data=expre_data_params,type="Real")
  simu_bursting_percentage_zero_percell <- gen_meancounts_zeropercent_percell(data=bursting, type="BGP_trendedBCV")
  simu_nobursting_percentage_zero_percell <- gen_meancounts_zeropercent_percell(data=nobursting,type="GP_trendedBCV")
  simu_common_burstring_percentage_zero_percell <- gen_meancounts_zeropercent_percell(data=simu_common_burst,type="BGP_commonBCV")
  simu_common_noburstring_percentage_zero_percell <- gen_meancounts_zeropercent_percell(data=simu_common_noburst,type="GP_commonBCV")
  simu_noBCV_burst_percentage_zero_percell <- gen_meancounts_zeropercent_percell(data=simu_noBCV_burst,type="BP")
  
  
  final_percentage_zero_percell <- rbind(real_percentage_zero_percell,
                                         simu_noBCV_burst_percentage_zero_percell,
                                         simu_common_noburstring_percentage_zero_percell,
                                         simu_nobursting_percentage_zero_percell,
                                         simu_common_burstring_percentage_zero_percell,
                                         simu_bursting_percentage_zero_percell)
  final_percentage_zero_percell$Type <- factor(final_percentage_zero_percell$Type,levels=c("Real","GP_commonBCV",
                                                                                           "GP_trendedBCV","BP",
                                                                                           "BGP_commonBCV","BGP_trendedBCV"))
  
  ###plot the zero_diff to real vs rank_counts plot per gene
  gen_diffzero_pergene <- function(final_mean_zero){
    final_mean_zero_rankcount <- NULL      
    real_data <- final_mean_zero[which(final_mean_zero$Type=="Real"),]
    real_data <- real_data[order(real_data$counts_zeros),]
    for (i in unique(final_mean_zero$Type)){
      data1 <- final_mean_zero[final_mean_zero$Type==i,]
      # data1$rank <- rank(data1$Mean_count)
      data1 <- data1[order(data1$counts_zeros),]

      data1$Diff_zero <- abs(log2(pmax(data1$Percentage_zeros,0.000001)/pmax(real_data$Percentage_zeros,0.000001)))
      final_mean_zero_rankcount <- rbind(final_mean_zero_rankcount,data1)
    }
    final_mean_zero_rankcount <- final_mean_zero_rankcount[final_mean_zero_rankcount$Type!="Real",]
    return(final_mean_zero_rankcount)
  }
  
  final_diff_zero_pergene <- gen_diffzero_pergene(final_mean_zero=final_mean_zero)
  # final_diff_zero_pergene <- final_diff_zero_pergene[is.na(final_diff_zero_pergene$Diff_zero)==F,]
  Diff_zero_pergene <- aggregate(final_diff_zero_pergene$Diff_zero, list(final_diff_zero_pergene$Type), median)
  Diff_zero_pergene$Zeros_Gene <- rank(Diff_zero_pergene$x)
  rownames(Diff_zero_pergene) <- Diff_zero_pergene$Group.1
  Diff_zero_pergene <- Diff_zero_pergene[c("GP_commonBCV",
                                           "GP_trendedBCV","BP",
                                           "BGP_commonBCV","BGP_trendedBCV"),c("x","Zeros_Gene")]
  
  
  
  ###plot the zero_diff to real vs rank_counts plot per cell
  gen_diffzero_percell <- function(final_mean_zero){
    final_mean_zero_rankcount <- NULL
    for (i in unique(final_mean_zero$Type)){
      data1 <- final_mean_zero[final_mean_zero$Type==i,]
      data1 <- data1[order(data1$Percentage_zeros),]
      real_data <- final_mean_zero[which(final_mean_zero$Type=="Real"),]
      real_data <- real_data[order(real_data$Percentage_zeros),]
      data1$Diff_zero_cell <- abs(log2(pmax(data1$Percentage_zeros,0.000001)/pmax(real_data$Percentage_zeros,0.000001)))
      # colnames(data2) <- "Diff_zero"
      # data2$Type <- i
      final_mean_zero_rankcount <- rbind(final_mean_zero_rankcount,data1)
    }
    final_mean_zero_rankcount <- final_mean_zero_rankcount[final_mean_zero_rankcount$Type!="Real",]
    
    
    return(final_mean_zero_rankcount)
  }
  
  final_diff_zero_percell <- gen_diffzero_percell(final_mean_zero=final_percentage_zero_percell)
  
  Diff_zero_percell <- aggregate(final_diff_zero_percell$Diff_zero_cell, list(final_diff_zero_percell$Type), median)
  Diff_zero_percell$Zeros_Cell <- rank(Diff_zero_percell$x)
  rownames(Diff_zero_percell) <- Diff_zero_percell$Group.1
  Diff_zero_percell <- Diff_zero_percell[c("GP_commonBCV",
                                           "GP_trendedBCV","BP",
                                           "BGP_commonBCV","BGP_trendedBCV"),c("x","Zeros_Cell")]
  
  
  Heatmap_data <- cbind(mean,variance,librarysize,Diff_zero_pergene,Diff_zero_percell)
  Heatmap_data <- Heatmap_data[,c(1:5)*2-1]
  colnames(Heatmap_data) <- c("Mean","Variance","LibrarySize","Zeros_Gene","Zeros_Cell")
  library(BBmisc)
  Heatmap_data <- t(Heatmap_data)
  # Heatmap_data <- t(scale(Heatmap_data))
  # Heatmap_data <- t(Heatmap_data[,c("Mean","Variance","LibrarySize","Zeros_Gene","Zeros_Cell")])
  return(Heatmap_data)
}





# Below We evaluated the performance of simulation data from SCIRP to recover the key characteristics from #
# real data. Characteristics included gene-wise expression mean, variance, cell-wise zero proportion, #
# gene-wise zero proportion, and library size. To know more details about how to use SCRIP generate these # 
# simulation data, you can check https://github.com/thecailab/SCRIP/blob/main/vignettes/SCRIPsimu.pdf # 


################################################# Klein data ##############################################
###########################################################################################################
###########################################################################################################
setwd("E:/DB/Dropbox/Qinfei/Simulation of SC based on splatter/Submition/Bioinformatics/Simulation data/Klein")
library("Biobase")
library('edgeR')
library('stats')
library(scales)

Klein <- read.csv("Klein.csv",header=T,row.names = 1)
expre_data <- as.matrix(Klein)

counts <- expre_data

GP_commonBCV <- get(load("GP-commonBCV.RData"))
colnames(GP_commonBCV) <- paste("cell",c(1:ncol(GP_commonBCV)),sep="")
GP_commonBCV <- round(GP_commonBCV)

GP_trendedBCV <- get(load("GP-trendedBCV.RData"))
colnames(GP_trendedBCV) <- paste("cell",c(1:ncol(GP_trendedBCV)),sep="")
GP_trendedBCV <- round(GP_trendedBCV)

BGP_commonBCV <- get(load("BGP-commonBCV.RData"))
colnames(BGP_commonBCV) <- paste("cell",c(1:ncol(BGP_commonBCV)),sep="")
BGP_commonBCV <- round(BGP_commonBCV)

BGP_trendedBCV <- get(load("BGP-trendedBCV.RData"))
colnames(BGP_trendedBCV) <- paste("cell",c(1:ncol(BGP_trendedBCV)),sep="")
BGP_trendedBCV <- round(BGP_trendedBCV)

BP <- get(load("BP.RData"))
colnames(BP) <- paste("cell",c(1:ncol(BP)),sep="")
BP <- round(BP)

res <- list()
res$counts <- counts
res$GP_commonBCV <- GP_commonBCV
res$GP_trendedBCV <- GP_trendedBCV
res$BGP_commonBCV <- BGP_commonBCV
res$BGP_trendedBCV <- BGP_trendedBCV
res$BP <- BP

Heatmap_data_Klein <- Gen_heatmatdata(res)




#######Tung data####
#####################
#####################
#### Tung ####
setwd("E:/DB/Dropbox/Qinfei/Simulation of SC based on splatter/Submition/Bioinformatics/Simulation data/Tung")
library("Biobase")
library('edgeR')
library('stats')
library(scales)
Tung <- read.table("Tung.txt",sep="\t",header=T,row.names=1)
expre_data <- as.matrix(Tung)
counts <- expre_data

GP_commonBCV <- get(load("GP-commonBCV.RData"))
colnames(GP_commonBCV) <- paste("cell",c(1:ncol(GP_commonBCV)),sep="")
GP_commonBCV <- round(GP_commonBCV)

GP_trendedBCV <- get(load("GP-trendedBCV.RData"))
colnames(GP_trendedBCV) <- paste("cell",c(1:ncol(GP_trendedBCV)),sep="")
GP_trendedBCV <- round(GP_trendedBCV)

BGP_commonBCV <- get(load("BGP-commonBCV.RData"))
colnames(BGP_commonBCV) <- paste("cell",c(1:ncol(BGP_commonBCV)),sep="")
BGP_commonBCV <- round(BGP_commonBCV)

BGP_trendedBCV <- get(load("BGP-trendedBCV.RData"))
colnames(BGP_trendedBCV) <- paste("cell",c(1:ncol(BGP_trendedBCV)),sep="")
BGP_trendedBCV <- round(BGP_trendedBCV)

BP <- get(load("BP.RData"))
colnames(BP) <- paste("cell",c(1:ncol(BP)),sep="")
BP <- round(BP)

res <- list()
res$counts <- counts
res$GP_commonBCV <- GP_commonBCV
res$GP_trendedBCV <- GP_trendedBCV
res$BGP_commonBCV <- BGP_commonBCV
res$BGP_trendedBCV <- BGP_trendedBCV
res$BP <- BP

Heatmap_data_Tung <- Gen_heatmatdata(res)




#######Camp data####
#####################
#####################
#### Camp ####
setwd("E:/DB/Dropbox/Qinfei/Simulation of SC based on splatter/Submition/Bioinformatics/Simulation data/Camp")
library("Biobase")
library('edgeR')
library('stats')
library(scales)
Camp <- read.table("Camp.txt",sep="\t",header=T,row.names=1)
Camp <- Camp[,c(7:ncol(Camp))]
expre_data_Camp <- Camp
colnames(expre_data_Camp) <- c(paste("Cell",1:ncol(expre_data_Camp),sep=""))
expre_data <- as.matrix(expre_data_Camp)

counts <- expre_data

GP_commonBCV <- get(load("GP-commonBCV.RData"))
colnames(GP_commonBCV) <- paste("cell",c(1:ncol(GP_commonBCV)),sep="")
GP_commonBCV <- round(GP_commonBCV)

GP_trendedBCV <- get(load("GP-trendedBCV.RData"))
colnames(GP_trendedBCV) <- paste("cell",c(1:ncol(GP_trendedBCV)),sep="")
GP_trendedBCV <- round(GP_trendedBCV)

BGP_commonBCV <- get(load("BGP-commonBCV.RData"))
colnames(BGP_commonBCV) <- paste("cell",c(1:ncol(BGP_commonBCV)),sep="")
BGP_commonBCV <- round(BGP_commonBCV)

BGP_trendedBCV <- get(load("BGP-trendedBCV.RData"))
colnames(BGP_trendedBCV) <- paste("cell",c(1:ncol(BGP_trendedBCV)),sep="")
BGP_trendedBCV <- round(BGP_trendedBCV)

BP <- get(load("BP.RData"))
colnames(BP) <- paste("cell",c(1:ncol(BP)),sep="")
BP <- round(BP)

res <- list()
res$counts <- counts
res$GP_commonBCV <- GP_commonBCV
res$GP_trendedBCV <- GP_trendedBCV
res$BGP_commonBCV <- BGP_commonBCV
res$BGP_trendedBCV <- BGP_trendedBCV
res$BP <- BP

Heatmap_data_Camp <- Gen_heatmatdata(res)








#######Xin data####
#####################
#####################
setwd("E:/DB/Dropbox/Qinfei/Simulation of SC based on splatter/Submition/Bioinformatics/Simulation data/MuSiC")
EMTAB.eset <- readRDS("EMTABesethealthy.rds")
library("Biobase")
library('edgeR')
library('stats')
library(scales)
expre_data = exprs(EMTAB.eset)
pheno_data=pData(EMTAB.eset)
expre_data=expre_data[,which(pheno_data[,3] == 2)]
counts <- expre_data

GP_commonBCV <- get(load("GP-commonBCV.RData"))
colnames(GP_commonBCV) <- paste("cell",c(1:ncol(GP_commonBCV)),sep="")

GP_trendedBCV <- get(load("GP-trendedBCV.RData"))
colnames(GP_trendedBCV) <- paste("cell",c(1:ncol(GP_trendedBCV)),sep="")

BGP_commonBCV <- get(load("BGP-commonBCV.RData"))
colnames(BGP_commonBCV) <- paste("cell",c(1:ncol(BGP_commonBCV)),sep="")

BGP_trendedBCV <- get(load("BGP-trendedBCV.RData"))
colnames(BGP_trendedBCV) <- paste("cell",c(1:ncol(BGP_trendedBCV)),sep="")

BP <- get(load("BP.RData"))
colnames(BP) <- paste("cell",c(1:ncol(BP)),sep="")
res <- list()
res$counts <- counts
res$GP_commonBCV <- GP_commonBCV
res$GP_trendedBCV <- GP_trendedBCV
res$BGP_commonBCV <- BGP_commonBCV
res$BGP_trendedBCV <- BGP_trendedBCV
res$BP <- BP

Heatmap_data_MuSiC <- Gen_heatmatdata(res)








#######Tirosh data####
#####################
#####################
#### Tirosh ####
setwd("E:/DB/Dropbox/Qinfei/Simulation of SC based on splatter/Submition/Bioinformatics/Simulation data/Tirosh")
library("Biobase")
library('edgeR')
library('stats')
library(scales)

counts <- get(load("counts.RData"))

GP_commonBCV <- get(load("GP-commonBCV.RData"))
colnames(GP_commonBCV) <- paste("cell",c(1:ncol(GP_commonBCV)),sep="")
GP_commonBCV <- round(GP_commonBCV)

GP_trendedBCV <- get(load("GP-trendedBCV.RData"))
colnames(GP_trendedBCV) <- paste("cell",c(1:ncol(GP_trendedBCV)),sep="")
GP_trendedBCV <- round(GP_trendedBCV)

BGP_commonBCV <- get(load("BGP-commonBCV.RData"))
colnames(BGP_commonBCV) <- paste("cell",c(1:ncol(BGP_commonBCV)),sep="")
BGP_commonBCV <- round(BGP_commonBCV)

BGP_trendedBCV <- get(load("BGP-trendedBCV.RData"))
colnames(BGP_trendedBCV) <- paste("cell",c(1:ncol(BGP_trendedBCV)),sep="")
BGP_trendedBCV <- round(BGP_trendedBCV)

BP <- get(load("BP.RData"))
colnames(BP) <- paste("cell",c(1:ncol(BP)),sep="")
BP <- round(BP)

res <- list()
res$counts <- counts
res$GP_commonBCV <- GP_commonBCV
res$GP_trendedBCV <- GP_trendedBCV
res$BGP_commonBCV <- BGP_commonBCV
res$BGP_trendedBCV <- BGP_trendedBCV
res$BP <- BP

Heatmap_data_Tirosh <- Gen_heatmatdata(res)










#######Zhou data####
#####################
#####################
#### Zhou  ####
setwd("E:/DB/Dropbox/Qinfei/Simulation of SC based on splatter/Submition/Bioinformatics/Simulation data/Zhou")
library("Biobase")
library('edgeR')
library('stats')
library(scales)

counts <- get(load("counts.RData"))

GP_commonBCV <- get(load("GP-commonBCV.RData"))
colnames(GP_commonBCV) <- paste("cell",c(1:ncol(GP_commonBCV)),sep="")
GP_commonBCV <- round(GP_commonBCV)

GP_trendedBCV <- get(load("GP-trendedBCV.RData"))
colnames(GP_trendedBCV) <- paste("cell",c(1:ncol(GP_trendedBCV)),sep="")
GP_trendedBCV <- round(GP_trendedBCV)

BGP_commonBCV <- get(load("BGP-commonBCV.RData"))
colnames(BGP_commonBCV) <- paste("cell",c(1:ncol(BGP_commonBCV)),sep="")
BGP_commonBCV <- round(BGP_commonBCV)

BGP_trendedBCV <- get(load("BGP-trendedBCV.RData"))
colnames(BGP_trendedBCV) <- paste("cell",c(1:ncol(BGP_trendedBCV)),sep="")
BGP_trendedBCV <- round(BGP_trendedBCV)

BP <- get(load("BP.RData"))
colnames(BP) <- paste("cell",c(1:ncol(BP)),sep="")
BP <- round(BP)

res <- list()
res$counts <- counts
res$GP_commonBCV <- GP_commonBCV
res$GP_trendedBCV <- GP_trendedBCV
res$BGP_commonBCV <- BGP_commonBCV
res$BGP_trendedBCV <- BGP_trendedBCV
res$BP <- BP

Heatmap_data_Zhou <- Gen_heatmatdata(res)









#######Seale Human data####
#####################
#####################
#### Seale Human ####
setwd("E:/DB/Dropbox/Qinfei/Simulation of SC based on splatter/Submition/Bioinformatics/Simulation data/Seale/Human")
library("Biobase")
library('edgeR')
library('stats')
library(scales)

counts <- get(load("counts.RData"))

GP_commonBCV <- get(load("GP-commonBCV.RData"))
colnames(GP_commonBCV) <- paste("cell",c(1:ncol(GP_commonBCV)),sep="")
GP_commonBCV <- round(GP_commonBCV)

GP_trendedBCV <- get(load("GP-trendedBCV.RData"))
colnames(GP_trendedBCV) <- paste("cell",c(1:ncol(GP_trendedBCV)),sep="")
GP_trendedBCV <- round(GP_trendedBCV)

BGP_commonBCV <- get(load("BGP-commonBCV.RData"))
colnames(BGP_commonBCV) <- paste("cell",c(1:ncol(BGP_commonBCV)),sep="")
BGP_commonBCV <- round(BGP_commonBCV)

BGP_trendedBCV <- get(load("BGP-trendedBCV.RData"))
colnames(BGP_trendedBCV) <- paste("cell",c(1:ncol(BGP_trendedBCV)),sep="")
BGP_trendedBCV <- round(BGP_trendedBCV)

BP <- get(load("BP.RData"))
colnames(BP) <- paste("cell",c(1:ncol(BP)),sep="")
BP <- round(BP)

res <- list()
res$counts <- counts
res$GP_commonBCV <- GP_commonBCV
res$GP_trendedBCV <- GP_trendedBCV
res$BGP_commonBCV <- BGP_commonBCV
res$BGP_trendedBCV <- BGP_trendedBCV
res$BP <- BP

Heatmap_data_Seale_Human <- Gen_heatmatdata(res)










#######Seale Mice data####
#####################
#####################
#### Seale Mice ####
setwd("E:/DB/Dropbox/Qinfei/Simulation of SC based on splatter/Submition/Bioinformatics/Simulation data/Seale/Mice")
library("Biobase")
library('edgeR')
library('stats')
library(scales)

counts <- get(load("counts.RData"))

GP_commonBCV <- get(load("GP-commonBCV.RData"))
colnames(GP_commonBCV) <- paste("cell",c(1:ncol(GP_commonBCV)),sep="")
GP_commonBCV <- round(GP_commonBCV)

GP_trendedBCV <- get(load("GP-trendedBCV.RData"))
colnames(GP_trendedBCV) <- paste("cell",c(1:ncol(GP_trendedBCV)),sep="")
GP_trendedBCV <- round(GP_trendedBCV)

BGP_commonBCV <- get(load("BGP-commonBCV.RData"))
colnames(BGP_commonBCV) <- paste("cell",c(1:ncol(BGP_commonBCV)),sep="")
BGP_commonBCV <- round(BGP_commonBCV)

BGP_trendedBCV <- get(load("BGP-trendedBCV.RData"))
colnames(BGP_trendedBCV) <- paste("cell",c(1:ncol(BGP_trendedBCV)),sep="")
BGP_trendedBCV <- round(BGP_trendedBCV)

BP <- get(load("BP.RData"))
colnames(BP) <- paste("cell",c(1:ncol(BP)),sep="")
BP <- round(BP)

res <- list()
res$counts <- counts
res$GP_commonBCV <- GP_commonBCV
res$GP_trendedBCV <- GP_trendedBCV
res$BGP_commonBCV <- BGP_commonBCV
res$BGP_trendedBCV <- BGP_trendedBCV
res$BP <- BP

Heatmap_data_Seale_Mice <- Gen_heatmatdata(res)







tran_data <- function(data){
  for (i in 1:nrow(data)){
    mean <- mean(data[i,])
    
  }
}


#### Now plot the heatmap (Figure S1) of MAD with five characteristics for differert methods in SCRIP and eight datasets ####
#############################################################################################################################
setwd("E:\\DB\\Dropbox\\Qinfei\\Simulation of SC based on splatter\\Submition\\Bioinformatics\\Simulation data")
# load("Heatmap_data.RData")
# load("Heatmap_data_mean.RData")
Heatmap_data <- rbind(Heatmap_data_MuSiC,Heatmap_data_Klein,Heatmap_data_Tung,Heatmap_data_Camp,
                      Heatmap_data_Tirosh,Heatmap_data_Zhou,Heatmap_data_Seale_Human,Heatmap_data_Seale_Mice)
rownames(Heatmap_data) <- paste0(rep(c("Xin","Klein","Tung","Camp","Tirosh","Zhou","Seale(Human)","Seale(Mice)"),each=5),
                                 rep(c("Mean","Variance","library size","Zeros per gene","Zeros per cell"),8))
cols = colorRampPalette(c("royalblue3","white"))(15)

data <- Heatmap_data[c((0:7)*5+1,(0:7)*5+2,(0:7)*5+3,(0:7)*5+4,(0:7)*5+5),]
# data <- t(scale(t(data), scale=F))

library(pheatmap)
features <- data.frame(Features=rep(c("Mean","Variance","library size","Zeros per cell","Zeros per gene"),each=8))
rownames(features) <- rownames(data)
features$Features <- factor(features$Features, levels=c("Mean","Variance","library size","Zeros per cell","Zeros per gene"))

colnames(data) <- c("GP-commonBCV(Splat)","GP-trendedBCV","BP","BP-commonBCV","BP-trendedBCV")

pdf("Heatmap mean comparion models in SCRIP for eights datasets log2 FC.pdf",width=15,height=15)
pheatmap(data,color=cols,cluster_rows=F,cluster_cols=F,
         fontsize = 20,angle_col=45,annotation_row=features,
         labels_row = rep(c("Xin","Klein","Tung","Camp","Tirosh","Zhou","Seale(Human)","Seale(Mice)"),5),
         cellheight=20, cellwidth = 60)
dev.off()





#### transfer the heatmap into boxplots (Figure 2B) and heatmap with only variance (Figure 2C) 
################################################################################################
data1 <- data.frame(value=c(as.vector(data[1:8,1:ncol(data)]),as.vector(data[9:16,1:ncol(data)]),
                            as.vector(data[17:24,1:ncol(data)]), as.vector(data[25:32,1:ncol(data)]),
                            as.vector(data[33:40,1:ncol(data)])), 
                    type=rep(c("Mean","Variance","library size","Zeros per cell","Zeros per gene"),each=40))
data1$type <- factor(data1$type, levels = c(c("Mean","Variance","library size","Zeros per cell","Zeros per gene")))


library("ggplot2")
setwd("E:\\DB\\Dropbox\\Qinfei\\Simulation of SC based on splatter\\Submition\\Bioinformatics\\Simulation data")
pdf("SCIRP boxplot for different features.pdf",height=8,width=6)
ggplot(data = data1, aes(x=type,y=value)) +
  labs(y="MAD", x = "Type")+
  geom_boxplot(lwd=1,fatten=1)+
  # scale_y_continuous(limits=c(0,4)) +
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
pdf("Heatmap comparion SCRIP for variance using eight datasets log2 FC3.pdf",width=9,height=6)
pheatmap(data2,color=cols,cluster_rows=F,cluster_cols=F,
         fontsize = 20,angle_col=45,
         labels_row = rep(c("Xin","Klein","Tung","Camp","Tirosh","Zhou","Seale(Human)","Seale(Mice)"),1),cellheight=35, cellwidth=35)
dev.off()
