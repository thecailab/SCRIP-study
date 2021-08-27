
setwd("/home/fqin/SCRIP/Cluster")
####Zhou
library(Seurat)
Zhou <- readRDS(file="GSE139827_mice.rds")
expre_data <- as.matrix(Zhou@assays$RNA@counts)
celltype  <- as.character(Idents(Zhou))
table(celltype)
CTlist <- unique(celltype)

library("ggplot2")

Group.sim <- function(data){
  library(dplyr)
  library(Seurat)
  library(cowplot)
  rownames(data) <- paste0('gene',1:nrow(data))
  colnames(data) <- paste0('cell',1:ncol(data))
  
  pbmc <- CreateSeuratObject(counts = data, project = "MuSiC" )
  pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
  pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 1000)
  all.genes <- rownames(pbmc)
  pbmc <- ScaleData(pbmc, features = all.genes)
  pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
  pbmc <- RunUMAP(pbmc, dims = 1:20)
  # pbmc <- RunTSNE(pbmc)
  return(pbmc)
}

### plot UMAP for the real dataset with predefined celltypes
real.data <-  expre_data[,which(celltype %in% CTlist)]
size=18
real.pbmc <- Group.sim(data=real.data)
real.pbmc@meta.data$Celltype <- celltype[celltype %in% CTlist]
real.umap <- DimPlot(real.pbmc, reduction = "umap",group.by = "Celltype",pt.size=1)+
  labs(title="Real")+
  scale_x_continuous(limits = c(-25, 25))+
  scale_y_continuous(limits = c(-25, 25))+
  
  theme(plot.title = element_text(face="bold",size=size,hjust=0.5),
        legend.text=element_text(size=size),
        axis.text=element_text(size=size,face="bold"),
        axis.title.x = element_text(size=size,face="bold"),
        axis.title.y = element_text(size=size,face="bold"),
        panel.background = element_rect(fill='white',colour="black"))

pdf("Zhou.UMAP.pdf")
real.umap
dev.off()



nsimu=1
combins <- combn(CTlist, 2)

library(ggplot2)
cor.matrix <- matrix(NA,11*ncol(combins),(nsimu+2))
cor.matrix[,1] <- rep(c("real","GP-commonBCV","GP-trendedBCV", "BP",
                        "BGP-commonBCV","BGP-trendedBCV",
                        "scDesign","powsimR","SPARSim","dyngen","SymSim"),ncol(combins))
colnames(cor.matrix) <- c("dataset","CTs",paste0("simu",1:nsimu))
simu=1
for (simu in 1:nsimu){
  ### plot the tSNE plot to check the clustering performance
  if (genetype=="whole"){ 
    GPcommonBCV.data <- NULL
    GPtrendedBCV.data <- NULL
    BGPtrendedBCV.data <- NULL
    BGPcommonBCV.data <- NULL
    BP.data <- NULL
    scDesign.data <- NULL
    powsimR.data <- NULL
    SPARSim.data <- NULL
    dyngen.data <- NULL
    SymSim.data <- NULL
    
    SCRIP.CT <- NULL
    scDesign.CT <- NULL
    powsimR.CT <- NULL
    SPARSim.CT <- NULL
    dyngen.CT <- NULL
    SymSim.CT <- NULL
    for (CT in CTlist){
      res <- get(load(paste0(CT,".group.simulation",simu,"Zhou.VEG1000Genes.addmean.RData")))  
      GPcommonBCV.data <- cbind(GPcommonBCV.data, res$SCRIP.GPcommonBCV)
      GPtrendedBCV.data <- cbind(GPtrendedBCV.data, res$SCRIP.GPtrendedBCV)
      BP.data <- cbind(BP.data, res$SCRIP.BP)  
      BGPcommonBCV.data <- cbind(BGPcommonBCV.data, res$SCRIP.BGPcommonBCV)
      BGPtrendedBCV.data <- cbind(BGPtrendedBCV.data, res$SCRIP.BGPtrendedBCV)
      scDesign.data <- cbind(scDesign.data, res$scdesign)
      powsimR.data <- cbind(powsimR.data, res$powsimR)
      SPARSim.data <- cbind(SPARSim.data, res$SPARSim)
      dyngen.data <- cbind(dyngen.data, res$dyngen)
      SymSim.data <- cbind(SymSim.data, res$SymSim)
      
      SCRIP.CT <- c(SCRIP.CT, rep(CT,ncol(res$SCRIP.GPcommonBCV)))
      scDesign.CT <- c(scDesign.CT, rep(CT,ncol(res$scdesign)))
      powsimR.CT <- c(powsimR.CT, rep(CT,ncol(res$powsimR)))
      SPARSim.CT <- c(SPARSim.CT, rep(CT,ncol(res$SPARSim)))
      dyngen.CT <- c(dyngen.CT, rep(CT,ncol(res$dyngen)))
      SymSim.CT <- c(SymSim.CT, rep(CT,ncol(res$SymSim)))
    }
    
    size=8
    ## Real data  
    real.data <-  expre_data[,which(celltype %in% CTlist)]
  }
  
  real.pbmc <- Group.sim(data=real.data)
  real.pbmc@meta.data$Celltype <- celltype[celltype %in% CTlist]
  real.umap <- DimPlot(real.pbmc, reduction = "umap",group.by = "Celltype",pt.size=1)+
    labs(title="Real")+
    scale_x_continuous(limits = c(-25, 25))+
    scale_y_continuous(limits = c(-25, 25))+
    theme(plot.title = element_text(face="bold",size=size,hjust=0.5),
          legend.text=element_text(size=size),
          axis.text=element_text(size=size,face="bold"),
          axis.title.x = element_text(size=size,face="bold"),
          axis.title.y = element_text(size=size,face="bold"),
          panel.background = element_rect(fill='white',colour="black"))
  
  
  
  GPcommonBCV.pbmc <- Group.sim(data=GPcommonBCV.data)
  GPcommonBCV.pbmc@meta.data$Celltype <- SCRIP.CT
  GPcommonBCV.umap <- DimPlot(GPcommonBCV.pbmc, reduction = "umap",group.by = "Celltype",pt.size=1)+
    labs(title="SCRIP GP-commonBCV")+
    scale_x_continuous(limits = c(-25, 25))+
    scale_y_continuous(limits = c(-25, 25))+
    theme(plot.title = element_text(face="bold",size=size,hjust=0.5),
          legend.text=element_text(size=size),
          axis.text=element_text(size=size,face="bold"),
          axis.title.x = element_text(size=size,face="bold"),
          axis.title.y = element_text(size=size,face="bold"),
          panel.background = element_rect(fill='white',colour="black"))
  
  
  GPtrendedBCV.pbmc <- Group.sim(data=GPtrendedBCV.data)
  GPtrendedBCV.pbmc@meta.data$Celltype <- SCRIP.CT
  GPtrendedBCV.umap <- DimPlot(GPtrendedBCV.pbmc, reduction = "umap",group.by = "Celltype",pt.size=1)+
    labs(title="SCRIP GP-trendedBCV")+
    scale_x_continuous(limits = c(-25, 25))+
    scale_y_continuous(limits = c(-25, 25))+
    theme(plot.title = element_text(face="bold",size=size,hjust=0.5),
          legend.text=element_text(size=size),
          axis.text=element_text(size=size,face="bold"),
          axis.title.x = element_text(size=size,face="bold"),
          axis.title.y = element_text(size=size,face="bold"),
          panel.background = element_rect(fill='white',colour="black"))
  
  
  
  BP.pbmc <- Group.sim(data=BP.data)
  BP.pbmc@meta.data$Celltype <- SCRIP.CT
  BP.umap <- DimPlot(BP.pbmc, reduction = "umap",group.by = "Celltype",pt.size=1)+
    labs(title="SCRIP BP")+
    scale_x_continuous(limits = c(-25, 25))+
    scale_y_continuous(limits = c(-25, 25))+
    theme(plot.title = element_text(face="bold",size=size,hjust=0.5),
          legend.text=element_text(size=size),
          axis.text=element_text(size=size,face="bold"),
          axis.title.x = element_text(size=size,face="bold"),
          axis.title.y = element_text(size=size,face="bold"),
          panel.background = element_rect(fill='white',colour="black"))
  
  
  
  ## Common BCV bursting 
  BGPcommonBCV.pbmc <- Group.sim(data=BGPcommonBCV.data)
  BGPcommonBCV.pbmc@meta.data$Celltype <- SCRIP.CT
  BGPcommonBCV.umap  <-  DimPlot(BGPcommonBCV.pbmc, reduction = "umap",group.by = "Celltype",pt.size=1)+
    labs(title="SCRIP BGP-commonBCV")+
    scale_x_continuous(limits = c(-25, 25))+
    scale_y_continuous(limits = c(-25, 25))+
    theme(plot.title = element_text(face="bold",size=size,hjust=0.5),
          legend.text=element_text(size=size),
          axis.text=element_text(size=size,face="bold"),
          axis.title.x = element_text(size=size,face="bold"),
          axis.title.y = element_text(size=size,face="bold"),
          panel.background = element_rect(fill='white',colour="black"))
  
  
  
  BGPtrendedBCV.pbmc <- Group.sim(data=BGPtrendedBCV.data)
  BGPtrendedBCV.pbmc@meta.data$Celltype <- SCRIP.CT
  BGPtrendedBCV.umap  <- DimPlot(BGPtrendedBCV.pbmc, reduction = "umap",group.by = "Celltype",pt.size=1)+
    labs(title="SCRIP BGP-trendedBCV")+
    scale_x_continuous(limits = c(-25, 25))+
    scale_y_continuous(limits = c(-25, 25))+
    theme(plot.title = element_text(face="bold",size=size,hjust=0.5),
          legend.text=element_text(size=size),
          axis.text=element_text(size=size,face="bold"),
          axis.title.x = element_text(size=size,face="bold"),
          axis.title.y = element_text(size=size,face="bold"),
          panel.background = element_rect(fill='white',colour="black"))
  
  
  
  
  ## scDesign 
  scDesign.pbmc <- Group.sim(data=scDesign.data)
  scDesign.pbmc@meta.data$Celltype <- scDesign.CT
  scDesign.umap  <- DimPlot(scDesign.pbmc, reduction = "umap",group.by = "Celltype",pt.size=1)+
    labs(title="scDesign")+
    scale_x_continuous(limits = c(-25, 25))+
    scale_y_continuous(limits = c(-25, 25))+
    theme(plot.title = element_text(face="bold",size=size,hjust=0.5),
          legend.text=element_text(size=size),
          axis.text=element_text(size=size,face="bold"),
          axis.title.x = element_text(size=size,face="bold"),
          axis.title.y = element_text(size=size,face="bold"),
          panel.background = element_rect(fill='white',colour="black"))
  
  
  
  ## powsimR 
  powsimR.pbmc <- Group.sim(data=powsimR.data)
  powsimR.pbmc@meta.data$Celltype <-  powsimR.CT
  powsimR.umap  <- DimPlot(powsimR.pbmc, reduction = "umap",group.by = "Celltype",pt.size=1)+
    labs(title="powsimR")+
    scale_x_continuous(limits = c(-25, 25))+
    scale_y_continuous(limits = c(-25, 25))+
    theme(plot.title = element_text(face="bold",size=size,hjust=0.5),
          legend.text=element_text(size=size),
          axis.text=element_text(size=size,face="bold"),
          axis.title.x = element_text(size=size,face="bold"),
          axis.title.y = element_text(size=size,face="bold"),
          panel.background = element_rect(fill='white',colour="black"))
  
  
  
  ## SPARSim 
  SAPRSim.pbmc <- Group.sim(data=SPARSim.data)
  SAPRSim.pbmc@meta.data$Celltype <- SPARSim.CT
  SAPRSim.umap   <-   DimPlot(SAPRSim.pbmc, reduction = "umap",group.by = "Celltype",pt.size=1)+
    labs(title="SPARSim")+
    scale_x_continuous(limits = c(-25, 25))+
    scale_y_continuous(limits = c(-25, 25))+
    theme(plot.title = element_text(face="bold",size=size,hjust=0.5),
          legend.text=element_text(size=size),
          axis.text=element_text(size=size,face="bold"),
          axis.title.x = element_text(size=size,face="bold"),
          axis.title.y = element_text(size=size,face="bold"),
          panel.background = element_rect(fill='white',colour="black"))
  
  
  
  
  
  
  
  ## dyngen 
  dyngen.pbmc <- Group.sim(data=as.matrix(dyngen.data))
  dyngen.pbmc@meta.data$Celltype <- dyngen.CT
  dyngen.umap   <-   DimPlot(dyngen.pbmc, reduction = "umap",group.by = "Celltype",pt.size=1)+
    labs(title="dyngen")+
    scale_x_continuous(limits = c(-25, 25))+
    scale_y_continuous(limits = c(-25, 25))+
    theme(plot.title = element_text(face="bold",size=size,hjust=0.5),
          legend.text=element_text(size=size),
          axis.text=element_text(size=size,face="bold"),
          axis.title.x = element_text(size=size,face="bold"),
          axis.title.y = element_text(size=size,face="bold"),
          panel.background = element_rect(fill='white',colour="black"))
  
  
  
  
  ## SymSim
  SymSim.pbmc <- Group.sim(data=SymSim.data)
  SymSim.pbmc@meta.data$Celltype <- SymSim.CT
  SymSim.umap   <-   DimPlot(SymSim.pbmc, reduction = "umap",group.by = "Celltype",pt.size=1)+
    labs(title="SymSim")+
    scale_x_continuous(limits = c(-25, 25))+
    scale_y_continuous(limits = c(-25, 25))+
    theme(plot.title = element_text(face="bold",size=size,hjust=0.5),
          legend.text=element_text(size=size),
          axis.text=element_text(size=size,face="bold"),
          axis.title.x = element_text(size=size,face="bold"),
          axis.title.y = element_text(size=size,face="bold"),
          panel.background = element_rect(fill='white',colour="black"))
  
  
  
  library(ggpubr)
  ## plot UMAPs for real data and simulation data together
  pdf(paste0("UMPA.ordered.Zhou.all.methods.simu.",simu,".pdf"),height=6,width=7)
  ggarrange(real.umap, GPcommonBCV.umap, GPtrendedBCV.umap, BP.umap, BGPcommonBCV.umap, BGPtrendedBCV.umap,
            scDesign.umap, powsimR.umap, SAPRSim.umap, dyngen.umap, SymSim.umap, 
            common.legend = TRUE, legend = "right", font.label = list(size = 10),
            nrow=4, ncol=3)
  dev.off()
  
  
  
  
  real.xy <- real.pbmc@reductions$umap@cell.embeddings
  group.real <- as.character(real.pbmc@meta.data$Celltype)
  
  GPcommonBCV.xy <- GPcommonBCV.pbmc@reductions$umap@cell.embeddings
  group.GPcommonBCV <- as.character(GPcommonBCV.pbmc@meta.data$Celltype)
  
  GPtrendedBCV.xy <- GPtrendedBCV.pbmc@reductions$umap@cell.embeddings
  group.GPtrendedBCV <- as.character(GPtrendedBCV.pbmc@meta.data$Celltype)
  
  BP.xy <- BP.pbmc@reductions$umap@cell.embeddings
  group.BP <- as.character(BP.pbmc@meta.data$Celltype)
  
  BGPcommonBCV.xy <- BGPcommonBCV.pbmc@reductions$umap@cell.embeddings
  group.BGPcommonBCV <- as.character(BGPcommonBCV.pbmc@meta.data$Celltype)
  
  BGPtrendedBCV.xy <- BGPtrendedBCV.pbmc@reductions$umap@cell.embeddings
  group.BGPtrendedBCV <- as.character(BGPtrendedBCV.pbmc@meta.data$Celltype)
  
  scDesign.xy <- scDesign.pbmc@reductions$umap@cell.embeddings
  group.scDesign <- as.character(scDesign.pbmc@meta.data$Celltype)
  
  
  powsimR.xy <- powsimR.pbmc@reductions$umap@cell.embeddings
  group.powsimR <- as.character(powsimR.pbmc@meta.data$Celltype)
  
  SAPRSim.xy <- SAPRSim.pbmc@reductions$umap@cell.embeddings
  group.SAPRSim <- as.character(SAPRSim.pbmc@meta.data$Celltype)
  
  dyngen.xy <- dyngen.pbmc@reductions$umap@cell.embeddings
  group.dyngen <- as.character(dyngen.pbmc@meta.data$Celltype)
  
  SymSim.xy <- SymSim.pbmc@reductions$umap@cell.embeddings
  group.SymSim <- as.character(SymSim.pbmc@meta.data$Celltype)
  
  
  
  library(pdist)
  library(RPMG)
  ## generate the distance matrix for pairwise cell types
  combins <- combn(CTlist, 2)
  for (s in 1:ncol(combins)){
    i=combins[1,s]
    j=combins[2,s]
    
    real.dist <- pdist(colMeans(real.xy[group.real==i,]), colMeans(real.xy[group.real==j,]))@dist
    GPcommonBCV.dist <- pdist(colMeans(GPcommonBCV.xy[group.GPcommonBCV==i,]), colMeans(GPcommonBCV.xy[group.GPcommonBCV==j,]))@dist
    GPtrendedBCV.dist <- pdist(colMeans(GPtrendedBCV.xy[group.GPtrendedBCV==i,]), colMeans(GPtrendedBCV.xy[group.GPtrendedBCV==j,]))@dist
    BP.dist <- pdist(colMeans(BP.xy[group.BP==i,]), colMeans(BP.xy[group.BP==j,]))@dist
    BGPcommonBCV.dist <- pdist(colMeans(BGPcommonBCV.xy[group.BGPcommonBCV==i,]), colMeans(BGPcommonBCV.xy[group.BGPcommonBCV==j,]))@dist
    BGPtrendedBCV.dist <- pdist(colMeans(BGPtrendedBCV.xy[group.BGPtrendedBCV==i,]), colMeans(BGPtrendedBCV.xy[group.BGPtrendedBCV==j,]))@dist
    scDesign.dist <- pdist(colMeans(scDesign.xy[group.scDesign==i,]), colMeans(scDesign.xy[group.scDesign==j,]))@dist
    powsimR.dist <- pdist(colMeans(powsimR.xy[group.powsimR==i,]), colMeans(powsimR.xy[group.powsimR==j,]))@dist
    SAPRSim.dist <- pdist(colMeans(SAPRSim.xy[group.SAPRSim==i,]), colMeans(SAPRSim.xy[group.SAPRSim==j,]))@dist
    dyngen.dist <- pdist(colMeans(dyngen.xy[group.dyngen==i,]), colMeans(dyngen.xy[group.dyngen==j,]))@dist
    SymSim.dist <- pdist(colMeans(SymSim.xy[group.SymSim==i,]), colMeans(SymSim.xy[group.SymSim==j,]))@dist
    
    
    cor.matrix[(s-1)*11+1,simu+2] <- real.dist
    cor.matrix[(s-1)*11+2,simu+2] <- GPcommonBCV.dist
    cor.matrix[(s-1)*11+3,simu+2] <- GPtrendedBCV.dist
    cor.matrix[(s-1)*11+4,simu+2] <- BP.dist
    cor.matrix[(s-1)*11+5,simu+2] <- BGPcommonBCV.dist
    cor.matrix[(s-1)*11+6,simu+2] <- BGPtrendedBCV.dist
    cor.matrix[(s-1)*11+7,simu+2] <- scDesign.dist
    cor.matrix[(s-1)*11+8,simu+2] <- powsimR.dist
    cor.matrix[(s-1)*11+9,simu+2] <- SAPRSim.dist
    cor.matrix[(s-1)*11+10,simu+2] <- dyngen.dist
    cor.matrix[(s-1)*11+11,simu+2] <- SymSim.dist
    
    cor.matrix[(((s-1)*11+1):((s-1)*11+11)),2] <- paste0(i," ", j)
    
  } 
}

cor.matrix <- as.data.frame(cor.matrix)
cor.matrix$simu1 <- as.numeric(cor.matrix$simu1)
save(cor.matrix,file="Zhou.cor.matrix.RData")


nsimu=1
cor.matrix$dataset <- factor(cor.matrix$dataset, levels=c("real","GP-commonBCV","GP-trendedBCV","BP",
                                                          "BGP-commonBCV","BGP-trendedBCV",
                                                          "scDesign","powsimR","SPARSim","dyngen","SymSim"))

## the distance of each pair of cell types need to scale by divided by the total distance of all the pair of cell types
for (type in c("real","GP-commonBCV","GP-trendedBCV","BP",
               "BGP-commonBCV","BGP-trendedBCV",
               "scDesign","powsimR","SPARSim","dyngen","SymSim")){
  data <- cor.matrix[which(cor.matrix$dataset==type),]
  datasum <- sum(data[,3:(nsimu+2)])
  cor.matrix[which(cor.matrix$dataset==type),3:(nsimu+2)] <-  cor.matrix[which(cor.matrix$dataset==type),3:(nsimu+2)]/datasum
}


## the scaled distance of each pair of cell types from  simulation data minus that from real data ##
## Thus, the samller the value, the better ##
####################################################################################################
for (i in c("GP-commonBCV","GP-trendedBCV","BP",
            "BGP-commonBCV","BGP-trendedBCV",
            "scDesign","powsimR","SPARSim","dyngen","SymSim")){
  
  data <- cor.matrix[which(cor.matrix$dataset==i),3:(nsimu+2)]
  real <- cor.matrix[which(cor.matrix$dataset=="real"),3:(nsimu+2)]
  cor.matrix[which(cor.matrix$dataset==i),3:(nsimu+2)] <- abs(data-real)
}


cor.matrix <- cor.matrix[which(cor.matrix$dataset!="real"),]
cor.matrix$mean <- cor.matrix[,3:(nsimu+2)]*100

cor.matrix$dataset <- factor(cor.matrix$dataset, levels=c("GP-commonBCV","GP-trendedBCV","BP",
                                                          "BGP-commonBCV","BGP-trendedBCV",
                                                          "scDesign","SPARSim","powsimR","dyngen","SymSim"))


mean <- aggregate(cor.matrix[, 4], list(cor.matrix$dataset), mean)
mean$CTs <- "Mean"
mean$simu1=0
mean <- mean[,c(1,3,4,2)]
colnames(mean) <- c("dataset","CTs","simu1","mean")
cor.matrix <- rbind(cor.matrix, mean)


## Generate the boxplots which can show the results (the number of pair of cell types) for each simulator
library("ggplot2")
pdf("Cluster Zhou.pdf",height=8,width=20)
cor.matrix1 <- cor.matrix[cor.matrix$CTs!="Mean",]
ggplot(data = cor.matrix1, aes(x=dataset,y=mean)) +
  labs(title="Clustering comparison",y="D.prop diff (%)", x = "Datasets")+
  geom_boxplot(lwd=1,fatten=1)+
  geom_point(size = 2) +
  # scale_colour_manual(values=c("red","gray30","forestgreen","gold4","violet","royalblue2"))+
  # scale_y_continuous(limits=c(0,100000)) +
  theme(legend.position="none",
        plot.title = element_text(face="bold", size = 40, hjust = 0.5),
        legend.title = element_text(size=35,face="bold",hjust=0.5),
        legend.text=element_text(size=35),
        axis.text.x =element_text(size=30,face="bold",angle=30, hjust=1,colour="black"),
        axis.text.y =element_text(size=30,face="bold",colour="black"),
        axis.title.x = element_text(size=35,face="bold",margin = margin(t = 0, r = 20, b = 0, l = 15)),
        axis.title.y = element_text(size=35,face="bold",margin = margin(t = 0, r = 20, b = 0, l = 15)),
        panel.background = element_rect(fill='white',colour="black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        plot.margin = unit(c(0.7,1.5,1.0,1.5), "inches")) 

dev.off()

