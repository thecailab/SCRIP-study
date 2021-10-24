## Below are functions about how we use simulation data to generate trajectory plots using Monocle R package and plot heatmap for cell distance evaluation.

library("monocle")
library(Rcpp)
### evaluate simulation using trajectory ####   
trajectplot <- function(counts){
  library(monocle)
  cds <- newCellDataSet(counts)
  # cds <- preprocess_cds(cds, num_dim = 50)
  cds <- estimateSizeFactors(cds)
  cds <- estimateDispersions(cds)
  
  cds <- detectGenes(cds, min_expr = 0.1)
  
  disp_table <- dispersionTable(cds)
  ordering_genes <- subset(disp_table, mean_expression >= 0.1)
  cds <- setOrderingFilter(cds, ordering_genes)
  # plot_ordering_genes(cds)
  HSMM_myo <- reduceDimension(cds, max_components = 2, num_dim = 2,
                              method = 'DDRTree')
  HSMM_myo <- orderCells(HSMM_myo)
  return(HSMM_myo)
}



setwd("E:\\DB\\Dropbox\\Qinfei\\Simulation of SC based on splatter\\Submition\\Bioinformatics\\Simulation data\\Trajectory")
setwd("/home/fqin/SCRIP/Cluster")
library(Seurat)
Tirosh <- readRDS(file="GSE72056_Melanoma.RDS")
expre_data <- as.matrix(Tirosh@assays$RNA@counts)
celltype  <- as.character(Idents(Tirosh))
CTlist <- c("B cell","CAF cell","Endo. cell","Macro cell","NK cell","T cell")
celltype <- celltype[which(celltype %in% CTlist)]

table(celltype)

expre_data <-  expre_data[,which(celltype %in% CTlist)]
real <- trajectplot(as.matrix(expre_data))
save(real, file="Tirosh.real.trajectory.RData")


GPcommonBCV <- get(load(file="GP-commonBCV.RData"))
GPcommonBCV.traj <- trajectplot(GPcommonBCV)
save(GPcommonBCV.traj,file="Tirosh.GPcommonBCV.trajectory.RData")


GPtrendedBCV <- get(load(file="GP-trendedBCV.RData"))
GPtrendedBCV.traj <- trajectplot(GPtrendedBCV)
save(GPtrendedBCV.traj ,file="Tirosh.GPtrendedBCV.trajectory.RData")


BP <- get(load(file="BP.RData"))
BP.traj <- trajectplot(BP)
save(BP.traj ,file="Tirosh.BP.trajectory.RData")


Splat <- get(load(file="Splat.RData"))
Splat.traj <- trajectplot(BP)
save(Splat.traj ,file="Tirosh.Splat.trajectory.RData")


scdesign <- get(load(file="scDesign.RData"))
scdesign.traj <- trajectplot(scdesign)
save(scdesign.traj ,file="Tirosh.scdesign.trajectory.RData")


powsimR <- get(load(file="powsimR.RData"))
powsimR.traj <- trajectplot(powsimR)
save(powsimR.traj ,file="Tirosh.powsimR.trajectory.RData")


SPARSim <- get(load(file="SPARSim.RData"))
SPARSim.traj <- trajectplot(SPARSim)
save(SPARSim.traj ,file="Tirosh.SPARSim.trajectory.RData")

SymSim <- get(load(file="SymSim.RData"))
SymSim.traj <- trajectplot(SymSim)
save(SymSim.traj ,file="Tirosh.SymSim.trajectory.RData")


dyngen <- get(load(file="dyngen.RData"))
dyngen.traj <- trajectplot(dyngen)
save(dyngen.traj ,file="Tirosh.dyngen.trajectory.RData")






library(monocle)
setwd("E:\\DB\\Dropbox\\Qinfei\\Simulation of SC based on splatter\\Submition\\Bioinformatics\\Simulation data\\Trajectory")
real <- get(load("Tirosh.real.trajectory.RData"))
real <- plot_cell_trajectory(real,cell_size = 3,cex=1.5,cell_name_size = 3)+
  ggtitle("Real")+
  theme(plot.title = element_text(face="bold",size=20,hjust=0.5),
        legend.title = element_text(size=20, face="bold"),
        legend.text=element_text(size=20),
        axis.text=element_text(size=20),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20))




GPcommonBCV <- get(load("Tirosh.GPcommonBCV.trajectory.RData"))
GPcommonBCV <- plot_cell_trajectory(GPcommonBCV,cell_size = 3,cex=1.5,cell_name_size = 3)+
  ggtitle("SCRIP GP-commonBCV nobursting")+
  theme(plot.title = element_text(face="bold",size=20,hjust=0.5),
        legend.title = element_text(size=20, face="bold"),
        legend.text=element_text(size=20),
        axis.text=element_text(size=20),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20))


GPtrendedBCV <- get(load("Tirosh.GPtrendedBCV.trajectory.RData"))
GPtrendedBCV <- plot_cell_trajectory(GPtrendedBCV,cell_size = 3,cex=1.5,cell_name_size = 3)+
  ggtitle("SCRIP GP-trendedBCV nobursting")+
  theme(plot.title = element_text(face="bold",size=20,hjust=0.5),
        legend.title = element_text(size=20, face="bold"),
        legend.text=element_text(size=20),
        axis.text=element_text(size=20),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20))


BP <- get(load("Tirosh.BP.trajectory.RData"))
BP <- plot_cell_trajectory(BP,cell_size = 3,cex=1.5,cell_name_size = 3)+
  ggtitle("SCRIP BP")+
  theme(plot.title = element_text(face="bold",size=20,hjust=0.5),
        legend.title = element_text(size=20, face="bold"),
        legend.text=element_text(size=20),
        axis.text=element_text(size=20),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20))



Splat <- get(load("Tirosh.Splat.trajectory.RData"))
Splat <- plot_cell_trajectory(Splat, cell_size = 3, cex=1.5, cell_name_size = 3)+
  ggtitle("Splat")+
  theme(plot.title = element_text(face="bold",size=20,hjust=0.5),
        legend.title = element_text(size=20, face="bold"),
        legend.text=element_text(size=20),
        axis.text=element_text(size=20),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20))




scdesign <- get(load("Tirosh.scdesign.trajectory.RData"))
scdesign <- plot_cell_trajectory(scdesign, cell_size = 3, cex=1.5, cell_name_size = 3)+
  ggtitle("scDesign")+
  theme(plot.title = element_text(face="bold",size=20,hjust=0.5),
        legend.title = element_text(size=20, face="bold"),
        legend.text=element_text(size=20),
        axis.text=element_text(size=20),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20))



powsimR <- get(load("Tirosh.powsimR.trajectory.RData"))
powsimR <- plot_cell_trajectory(powsimR, cell_size = 3, cex=1.5, cell_name_size = 3)+
  ggtitle("powsimR")+
  theme(plot.title = element_text(face="bold",size=20,hjust=0.5),
        legend.title = element_text(size=20, face="bold"),
        legend.text=element_text(size=20),
        axis.text=element_text(size=20),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20))



SPARSim <- get(load("Tirosh.SPARSim.trajectory.RData"))
SPARSim <- plot_cell_trajectory(SPARSim, cell_size = 3, cex=1.5, cell_name_size = 3)+
  ggtitle("SPARSim")+
  theme(plot.title = element_text(face="bold",size=20,hjust=0.5),
        legend.title = element_text(size=20, face="bold"),
        legend.text=element_text(size=20),
        axis.text=element_text(size=20),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20))



dyngen <- get(load("Tirosh.dyngen.trajectory.RData"))
dyngen <- plot_cell_trajectory(dyngen, cell_size = 3, cex=1.5, cell_name_size = 3)+
  ggtitle("dyngen")+
  theme(plot.title = element_text(face="bold",size=20,hjust=0.5),
        legend.title = element_text(size=20, face="bold"),
        legend.text=element_text(size=20),
        axis.text=element_text(size=20),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20))



SymSim <- get(load("Tirosh.SymSim.trajectory.RData"))
SymSim <- plot_cell_trajectory(SymSim, cell_size = 3, cex=1.5, cell_name_size = 3)+
  ggtitle("SymSim")+
  theme(plot.title = element_text(face="bold",size=20,hjust=0.5),
        legend.title = element_text(size=20, face="bold"),
        legend.text=element_text(size=20),
        axis.text=element_text(size=20),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20))



library(ggpubr)
pdf("trajectplot Tirosh.pdf",height=20,width=18)
ggarrange(real, GPcommonBCV, GPtrendedBCV, BP, BGPcommonBCV, BGPtrendedBCV,  scdesign, powsimR, SPARSim,  dyngen, SymSim, 
          font.label = list(size = 28),
          ncol = 3, nrow = 4)
dev.off()






Compon1 <- data.frame(real=real$data$data_dim_1, GPcommon=GPcommonBCV$data$data_dim_1, 
                      GPtrend=GPtrendedBCV$data$data_dim_1, BP=BP$data$data_dim_1, Splat=Splat$data$data_dim_1, 
                      scdesign=scdesign$data$data_dim_1, powsimR=powsimR$data$data_dim_1,
                      SPARSim=SPARSim$data$data_dim_1,
                      dyngen=dyngen$data$data_dim_1, SymSim=SymSim$data$data_dim_1)

Compon2 <- data.frame(real=real$data$data_dim_2, GPcommon=GPcommonBCV$data$data_dim_2, 
                      GPtrend=GPtrendedBCV$data$data_dim_2, BP=BP$data$data_dim_2, Splat=Splat$data$data_dim_2, 
                      scdesign=scdesign$data$data_dim_2, powsimR=powsimR$data$data_dim_2,
                      SPARSim=SPARSim$data$data_dim_2,
                      dyngen=dyngen$data$data_dim_2, SymSim=SymSim$data$data_dim_2)

pseudotime <- data.frame(real=real$data$Pseudotime, GPcommon=GPcommonBCV$data$Pseudotime, 
                         GPtrend=GPtrendedBCV$data$Pseudotime, BP=BP$data$Pseudotime, Splat=Splat$data$Pseudotime, 
                         scdesign=scdesign$data$Pseudotime, powsimR=powsimR$data$Pseudotime,
                         SPARSim=SPARSim$data$Pseudotime,
                         dyngen=dyngen$data$Pseudotime, SymSim=SymSim$data$Pseudotime)

data.order <- list(Compon1=Compon1, Compon2=Compon2, pseudotime=pseudotime)
save(data.order, file="data.order.Tirosh.RData")






data.order <- get(load("data.order.Tirosh.RData"))
Compon1 <- data.order$Compon1
Compon2 <- data.order$Compon2
pseudotime <- data.order$pseudotime

res.matrix <- matrix(NA,9,4)
rownames(res.matrix) <- c("GPcommon","Splat","GPtrend","BP",
                          "scdesign","powsimR","SPARSim","dyngen","SymSim")
colnames(res.matrix)=c("cp1_posi","cp1_nega","cp2_posi","cp2_nega")
for (i in 1:9){
  res.matrix[i,1] <- mean(abs(pseudotime[rank(Compon1[,i+1]),i+1]-pseudotime[rank(Compon1[,1]),1]))
  res.matrix[i,2] <- mean(abs(pseudotime[rank(-Compon1[,i+1]),i+1]-pseudotime[rank(Compon1[,1]),1]))
  res.matrix[i,3] <-mean(abs(pseudotime[rank(Compon2[,i+1]),i+1]-pseudotime[rank(Compon2[,1]),1]))
  res.matrix[i,4] <- mean(abs(pseudotime[rank(-Compon2[,i+1]),i+1]-pseudotime[rank(Compon2[,1]),1]))
}

res <- as.data.frame(res.matrix)
res$order_cp1.posi <- rank(res$cp1_posi)
res$order_cp1.nega <- rank(res$cp1_nega)
res$order_cp2.posi <- rank(res$cp2_posi)
res$order_cp2_nega <- rank(res$cp2_nega)

res$best_order <- apply(res[,5:8],1,min)
res.Tirosh <- res
res.Tirosh


## Xin
setwd("E:\\DB\\Dropbox\\Qinfei\\Simulation of SC based on splatter\\Submition\\Bioinformatics\\Simulation data\\Trajectory\\Xin")
data.order <- get(load("data.order.Xin.RData"))
Compon1 <- data.order$Compon1
Compon2 <- data.order$Compon2
pseudotime <- data.order$pseudotime

res.matrix <- matrix(NA,9,4)
rownames(res.matrix) <- c("GPcommon","Splat","GPtrend","BP",
                          "scdesign","powsimR","SPARSim","dyngen","SymSim")
colnames(res.matrix)=c("cp1_posi","cp1_nega","cp2_posi","cp2_nega")
for (i in 1:9){
  res.matrix[i,1] <- mean(abs(pseudotime[rank(Compon1[,i+1]),i+1]-pseudotime[rank(Compon1[,1]),1]))
  res.matrix[i,2] <- mean(abs(pseudotime[rank(-Compon1[,i+1]),i+1]-pseudotime[rank(Compon1[,1]),1]))
  res.matrix[i,3] <-mean(abs(pseudotime[rank(Compon2[,i+1]),i+1]-pseudotime[rank(Compon2[,1]),1]))
  res.matrix[i,4] <- mean(abs(pseudotime[rank(-Compon2[,i+1]),i+1]-pseudotime[rank(Compon2[,1]),1]))
}

res <- as.data.frame(res.matrix)
res$order_cp1.posi <- rank(res$cp1_posi)
res$order_cp1.nega <- rank(res$cp1_nega)
res$order_cp2.posi <- rank(res$cp2_posi)
res$order_cp2_nega <- rank(res$cp2_nega)

res$best_order <- apply(res[,5:8],1,min)
res.Xin <- res
res.Xin


## Seale Mice
setwd("E:\\DB\\Dropbox\\Qinfei\\Simulation of SC based on splatter\\Submition\\Bioinformatics\\Simulation data\\Trajectory\\Seale_M")
data.order <- get(load("data.order.Seale_M.RData"))
Compon1 <- data.order$Compon1
Compon2 <- data.order$Compon2
pseudotime <- data.order$pseudotime

res.matrix <- matrix(NA,9,4)
rownames(res.matrix) <- c("GPcommon","Splat","GPtrend","BP",
                          "scdesign","powsimR","SPARSim","dyngen","SymSim")
colnames(res.matrix)=c("cp1_posi","cp1_nega","cp2_posi","cp2_nega")
for (i in 1:9){
  res.matrix[i,1] <- mean(abs(pseudotime[rank(Compon1[,i+1]),i+1]-pseudotime[rank(Compon1[,1]),1]))
  res.matrix[i,2] <- mean(abs(pseudotime[rank(-Compon1[,i+1]),i+1]-pseudotime[rank(Compon1[,1]),1]))
  res.matrix[i,3] <-mean(abs(pseudotime[rank(Compon2[,i+1]),i+1]-pseudotime[rank(Compon2[,1]),1]))
  res.matrix[i,4] <- mean(abs(pseudotime[rank(-Compon2[,i+1]),i+1]-pseudotime[rank(Compon2[,1]),1]))
}

res <- as.data.frame(res.matrix)
res$order_cp1.posi <- rank(res$cp1_posi)
res$order_cp1.nega <- rank(res$cp1_nega)
res$order_cp2.posi <- rank(res$cp2_posi)
res$order_cp2_nega <- rank(res$cp2_nega)

res$best_order <- apply(res[,5:8],1,min)
res.Seale_M <- res
res.Seale_M



## Zhou
setwd("E:\\DB\\Dropbox\\Qinfei\\Simulation of SC based on splatter\\Submition\\Bioinformatics\\Simulation data\\Trajectory\\Zhou")
data.order <- get(load("data.order.Zhou.RData"))
Compon1 <- data.order$Compon1
Compon2 <- data.order$Compon2
pseudotime <- data.order$pseudotime

res.matrix <- matrix(NA,9,4)
rownames(res.matrix) <- c("GPcommon","Splat","GPtrend","BP",
                          "scdesign","powsimR","SPARSim","dyngen","SymSim")
colnames(res.matrix)=c("cp1_posi","cp1_nega","cp2_posi","cp2_nega")
for (i in 1:9){
  res.matrix[i,1] <- mean(abs(pseudotime[rank(Compon1[,i+1]),i+1]-pseudotime[rank(Compon1[,1]),1]))
  res.matrix[i,2] <- mean(abs(pseudotime[rank(-Compon1[,i+1]),i+1]-pseudotime[rank(Compon1[,1]),1]))
  res.matrix[i,3] <-mean(abs(pseudotime[rank(Compon2[,i+1]),i+1]-pseudotime[rank(Compon2[,1]),1]))
  res.matrix[i,4] <- mean(abs(pseudotime[rank(-Compon2[,i+1]),i+1]-pseudotime[rank(Compon2[,1]),1]))
}

res <- as.data.frame(res.matrix)
res$order_cp1.posi <- rank(res$cp1_posi)
res$order_cp1.nega <- rank(res$cp1_nega)
res$order_cp2.posi <- rank(res$cp2_posi)
res$order_cp2_nega <- rank(res$cp2_nega)

res$best_order <- apply(res[,5:8],1,min)
res.Zhou <- res
res.Zhou



final.res <- data.frame(Datasets= c(rep("Seale_Mice",9), rep("Tirosh",9), rep("Xin", 9), rep("Zhou", 9)),
                        Order=c(res.Seale_M$best_order, res.Tirosh$best_order, res.Xin$best_order, res.Zhou$best_order),
                        Methods=c("GPcommon","Splat","GPtrend","BP",
                                  "scdesign","powsimR","SPARSim","dyngen","SymSim"))

final.res


final.res1 <- data.frame(Seale_Mice = res.Seale_M$best_order, Tirosh=res.Tirosh$best_order, 
                         Xin= res.Xin$best_order, Zhou=res.Zhou$best_order)
reorder <- order(apply(final.res1,1,mean))
rownames(final.res1) <- c("GP-commonBCV","Splat","GP-trendedBCV","BP",
                          "scdesign","powsimR","SPARSim","dyngen","SymSim")
final.res1 <- final.res1[reorder,]
setwd("E:\\DB\\Dropbox\\Qinfei\\Simulation of SC based on splatter\\Submition\\Bioinformatics\\Simulation data\\Trajectory")
pdf("Trajecotry evaluation heatmap plots 1.2.pdf", height = 8, width=10)
library(pheatmap)
pheatmap(t(final.res1),
         main="Trajectory",
         color=colorRampPalette(c("coral1","white"))(8),
         cluster_rows=F,cluster_cols=F,
         fontsize = 20,angle_col=45,
         cellheight=40, cellwidth=32)
dev.off()




