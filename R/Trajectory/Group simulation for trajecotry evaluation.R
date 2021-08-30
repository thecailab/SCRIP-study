setwd("E:/Dropbox/scRNA-deconvolution/Splatter simulator/Data from Tung")


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
  plot <- plot_cell_trajectory(HSMM_myo,cell_size = 5,cex=3,cell_name_size = 4)+
    theme(plot.title = element_text(face="bold",size=30,hjust=0.5),
          legend.title = element_text(size=30, face="bold"),
          legend.text=element_text(size=30),
          axis.text=element_text(size=30),
          axis.title.x = element_text(size=30,face="bold"),
          axis.title.y = element_text(size=30,face="bold"))
  return(plot)
}

tiff("Tung dataset Trajectories splatter.tiff",units="in",height=10,width=10,res=300)
trajectplot(round(counts(exp1)))
dev.off()

setwd("E:\\DB\\Dropbox\\Qinfei\\Simulation of SC based on splatter\\Submition\\Bioinformatics\\Simulation data\\Trajectory")
Tung <- read.table("Tung.txt",sep="\t",header=T,row.names=1)
real <- trajectplot(as.matrix(Tung))
save(real, file="Tung.real.trajectory.RData")

GPcommonBCV <- get(load(file="GP-commonBCV.RData"))
GPcommonBCV.traj <- trajectplot(GPcommonBCV)
save(GPcommonBCV.traj,file="Tung.GPcommonBCV.trajectory.RData")

BGPtrendedBCV <- get(load(file="GP-trendedBCV.RData"))
BGPtrendedBCV.traj <- trajectplot(BGPtrendedBCV)
save(BGPtrendedBCV.traj ,file="Tung.BGPtrendedBCV.trajectory.RData")

BP <- get(load(file="BP.RData"))
BP.traj <- trajectplot(BP)
save(BP.traj ,file="Tung.BP.trajectory.RData")

BGPcommonBCV <- get(load(file="BGP-trendedBCV.RData"))
BGPcommonBCV.traj <- trajectplot(BGPcommonBCV)
save(BGPcommonBCV.traj ,file="Tung.BGPcommonBCV.trajectory.RData")

BGPtrendedBCV <- get(load(file="BGP-trendedBCV.RData"))
BGPtrendedBCV.traj <- trajectplot(BGPtrendedBCV)
save(BGPtrendedBCV.traj ,file="Tung.BGPtrendedBCV.trajectory.RData")

scdesign <- get(load(file="scDesign.RData"))
scdesign.traj <- trajectplot(scdesign)
save(scdesign.traj ,file="Tung.scdesign.trajectory.RData")

scDD <- get(load(file="scDD.RData"))
scDD.traj <- trajectplot(scDD)
save(scDD.traj ,file="Tung.scDD.trajectory.RData")

powsimR <- get(load(file="powsimR.RData"))
powsimR.traj <- trajectplot(powsimR)
save(powsimR.traj ,file="Tung.powsimR.trajectory.RData")

SPARSim <- get(load(file="SPARSim.RData"))
SPARSim.traj <- trajectplot(SPARSim)
save(SPARSim.traj ,file="Tung.SPARSim.trajectory.RData")

SymSim <- get(load(file="SymSim.RData"))
SymSim.traj <- trajectplot(SymSim)
save(SymSim.traj ,file="Tung.SymSim.trajectory.RData")

dyngen <- get(load(file="dyngen.RData"))
dyngen.traj <- trajectplot(dyngen)
save(dyngen.traj ,file="Tung.dyngen.trajectory.RData")






setwd("E:\\DB\\Dropbox\\Qinfei\\Simulation of SC based on splatter\\Tung trajectory")
data <- get(load("Tung.real.trajectory.RData"))

library(monocle)
setwd("E:\\DB\\Dropbox\\Qinfei\\Simulation of SC based on splatter\\Submition\\Bioinformatics\\Simulation data\\Trajectory")
real <- get(load("Tung.real.trajectory.RData"))
real <- plot_cell_trajectory(data,cell_size = 3,cex=1.5,cell_name_size = 3)+
  ggtitle("Real")+
  theme(plot.title = element_text(face="bold",size=20,hjust=0.5),
        legend.title = element_text(size=20, face="bold"),
        legend.text=element_text(size=20),
        axis.text=element_text(size=20),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20))
real 

GPcommonBCV <- get(load("Tung.GPcommonBCV.trajectory.RData"))
GPcommonBCV <- plot_cell_trajectory(GPcommonBCV,cell_size = 3,cex=1.5,cell_name_size = 3)+
  ggtitle("SCRIP GP-commonBCV nobursting")+
  theme(plot.title = element_text(face="bold",size=20,hjust=0.5),
        legend.title = element_text(size=20, face="bold"),
        legend.text=element_text(size=20),
        axis.text=element_text(size=20),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20))


GPtrendedBCV <- get(load("Tung.GPtrendedBCV.trajectory.RData"))
GPtrendedBCV <- plot_cell_trajectory(GPtrendedBCV,cell_size = 3,cex=1.5,cell_name_size = 3)+
  ggtitle("SCRIP GP-trendedBCV nobursting")+
  theme(plot.title = element_text(face="bold",size=20,hjust=0.5),
        legend.title = element_text(size=20, face="bold"),
        legend.text=element_text(size=20),
        axis.text=element_text(size=20),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20))


BP <- get(load("Tung.BP.trajectory.RData"))
BP <- plot_cell_trajectory(BP,cell_size = 3,cex=1.5,cell_name_size = 3)+
  ggtitle("SCRIP BP")+
  theme(plot.title = element_text(face="bold",size=20,hjust=0.5),
        legend.title = element_text(size=20, face="bold"),
        legend.text=element_text(size=20),
        axis.text=element_text(size=20),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20))


BGPcommonBCV <- get(load("Tung.BGPcommonBCV.trajectory.RData"))
BGPcommonBCV <- plot_cell_trajectory(BGPcommonBCV, cell_size = 3, cex=1.5, cell_name_size = 3)+
  ggtitle("SCRIP BGP-commonBCV bursting")+
  theme(plot.title = element_text(face="bold",size=20,hjust=0.5),
        legend.title = element_text(size=20, face="bold"),
        legend.text=element_text(size=20),
        axis.text=element_text(size=20),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20))

BGPtrendedBCV <- get(load("Tung.BGPtrendedBCV.trajectory.RData"))
BGPtrendedBCV <- plot_cell_trajectory(BGPtrendedBCV, cell_size = 3, cex=1.5, cell_name_size = 3)+
  ggtitle("SCRIP BGP-trendedBCV bursting")+
  theme(plot.title = element_text(face="bold",size=20,hjust=0.5),
        legend.title = element_text(size=20, face="bold"),
        legend.text=element_text(size=20),
        axis.text=element_text(size=20),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20))


scdesign <- get(load("Tung.scdesign.trajectory.RData"))
scdesign <- plot_cell_trajectory(scdesign, cell_size = 3, cex=1.5, cell_name_size = 3)+
  ggtitle("scDesign")+
  theme(plot.title = element_text(face="bold",size=20,hjust=0.5),
        legend.title = element_text(size=20, face="bold"),
        legend.text=element_text(size=20),
        axis.text=element_text(size=20),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20))


powsimR <- get(load("Tung.powsimR.trajectory.RData"))
powsimR <- plot_cell_trajectory(powsimR, cell_size = 3, cex=1.5, cell_name_size = 3)+
  ggtitle("powsimR")+
  theme(plot.title = element_text(face="bold",size=20,hjust=0.5),
        legend.title = element_text(size=20, face="bold"),
        legend.text=element_text(size=20),
        axis.text=element_text(size=20),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20))


SPARSim <- get(load("Tung.SPARSim.trajectory.RData"))
SPARSim <- plot_cell_trajectory(SPARSim, cell_size = 3, cex=1.5, cell_name_size = 3)+
  ggtitle("SPARSim")+
  theme(plot.title = element_text(face="bold",size=20,hjust=0.5),
        legend.title = element_text(size=20, face="bold"),
        legend.text=element_text(size=20),
        axis.text=element_text(size=20),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20))


scDD <- get(load("Tung.scDD.trajectory.RData"))
scDD <- plot_cell_trajectory(scDD, cell_size = 3, cex=1.5, cell_name_size = 3)+
  ggtitle("scDD")+
  theme(plot.title = element_text(face="bold",size=20,hjust=0.5),
        legend.title = element_text(size=20, face="bold"),
        legend.text=element_text(size=20),
        axis.text=element_text(size=20),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20))


dyngen <- get(load("Tung.dyngen.trajectory.RData"))
dyngen <- plot_cell_trajectory(dyngen, cell_size = 3, cex=1.5, cell_name_size = 3)+
  ggtitle("dyngen")+
  theme(plot.title = element_text(face="bold",size=20,hjust=0.5),
        legend.title = element_text(size=20, face="bold"),
        legend.text=element_text(size=20),
        axis.text=element_text(size=20),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20))


SymSim <- get(load("Tung.SymSim.trajectory.RData"))
SymSim <- plot_cell_trajectory(SymSim, cell_size = 3, cex=1.5, cell_name_size = 3)+
  ggtitle("SymSim")+
  theme(plot.title = element_text(face="bold",size=20,hjust=0.5),
        legend.title = element_text(size=20, face="bold"),
        legend.text=element_text(size=20),
        axis.text=element_text(size=20),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20))



library(ggpubr)
pdf("trajectplot.pdf",height=15,width=18)
ggarrange(real, GPcommonBCV, GPtrendedBCV, BP, BGPcommonBCV, BGPtrendedBCV, 
          # scdesign, powsimR, SPARSim, scDD, dyngen, SymSim,
          font.label = list(size = 28),
          ncol = 3, nrow = 3)
dev.off()



