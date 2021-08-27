source("E:\\DB\\Dropbox\\Qinfei\\Simulation of SC based on splatter\\Code\\BCV plots.R")

setwd("E:\\DB\\Dropbox\\Qinfei\\Simulation of SC based on splatter\\Code\\powsimR")
for (i in list.files()){
  source(i)
}

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


fun <- function(counts.matrix){
  rownames(counts.matrix) <- paste0("Gene",1:nrow(counts.matrix))
  colnames(counts.matrix) <- paste0("Cell",1:ncol(counts.matrix))
  
  #### SCRIP ####
  start_time <- Sys.time()  
  message("Starting simulating SCRIP")
  library(splatter)
  params <- splatEstimate(counts.matrix)
  parshape = params@mean.shape
  parrate = params@mean.rate
  parnCells = params@batchCells
  parnGene = params@nGenes  
  library(SCRIP)
  single_cellsample <- SCRIPsimu(data=counts.matrix, params=params,
                                 mode="GP-trendedBCV")
  SCRIP <- counts(single_cellsample)
  end_time <- Sys.time()
  print(paste0("the time for running SCRIP is"))
  print(end_time-start_time)
  save(SCRIP, file="SCRIP.RData")


  #### scDesign ####
  start_time <- Sys.time()
  message("Starting simulating scDesign")
  library(scDesign)
  librarysize <- sum(apply(counts.matrix,2,sum))
  scDesign <- design_data(counts.matrix, S = librarysize, ncell=parnCells, ngroup = 1, pUp = 0.05,
                          pDown = 0.05, fU = 5, fL = 1.5, ncores = 1)
  end_time <- Sys.time()
  print(paste0("the time for running scDesign is"))
  print(end_time-start_time)
  save(scDesign, file="scDesign.RData")


  #### powsimR ####
  start_time <- Sys.time()
  message("Starting simulating powsimR")
  estparam <- estimateParam(countData = counts.matrix,
                            readData = NULL,
                            MeanFragLengths = NULL,
                            RNAseq = 'singlecell', Protocol = 'Read',
                            Distribution = 'NB', Normalisation = "scran",
                            GeneFilter = 0.0, SampleFilter = 5,
                            sigma = 1.96, NCores = NULL, verbose = TRUE)
  # plotParam(estparam, Annot = FALSE)
  p.lfc <- function(x) sample(c(-1,1), size=x,replace=T)*rgamma(x, shape = parshape, rate = parrate)
  setupres <- Setup(ngenes = parnGene, nsims = 1,
                    p.DE = 0.0, pLFC = p.lfc,
                    n1 = parnCells/2, n2 = parnCells/2,
                    Thinning = NULL, LibSize = 'given',
                    estParamRes = estparam,
                    estSpikeRes = NULL,
                    DropGenes = FALSE,
                    setup.seed = 2020, verbose = TRUE)
  powsimR <- simulateDE_new(SetupRes = setupres,
                            Prefilter = NULL, Imputation = NULL,
                            Normalisation = 'sctransform', Label = 'none',
                            DEmethod = "limma-trend", DEFilter = FALSE,
                            NCores = NULL, verbose = TRUE)
  end_time <- Sys.time()
  print(paste0("the time for running powsimR is"))
  print(end_time-start_time)
  save(powsimR, file="powsimR.RData")


  ## SPARSim ##
  start_time <- Sys.time()
  message("Starting simulating SPARSim")
  library(SPARSim)
  # count_norm <- scran_normalization(counts.matrix)
  count_conditions <- list(cond_A=c(1:(ncol(counts.matrix))))
  SPARSim_sim_param <- SPARSim_estimate_parameter_from_data(raw_data = counts.matrix,
                                                            norm_data =  counts.matrix,
                                                            conditions = count_conditions)
  SPARSim <- SPARSim_simulation(dataset_parameter = SPARSim_sim_param)
  SPARSim <- SPARSim$count_matrix
  end_time <- Sys.time()
  print(paste0("the time for running SPARSim is"))
  print(end_time-start_time)
  save(SPARSim, file="SPARSim.RData")

  
  
  ## scDD ##
  library(scDD)
  start_time <- Sys.time()  
  message("Starting simulating scDD")
  # library(SingleCellExperiment)
  condition <- sample(1:2, ncol(counts.matrix), replace = TRUE)
  names(condition) <- colnames(counts.matrix)
  scDD.param <- scDDEstimate(counts.matrix,conditions=condition)
  
  scDD <- scDDSimulate(params = scDD.param)
  scDD <- as.matrix(assay(scDD))
  scDD <- round(scDD)
  end_time <- Sys.time()
  print(paste0("the time for running scDD is"))
  print(end_time-start_time)
  save(scDD, file="scDD.RData")
  
  
  ## dyngen ##
  library(dyngen)
  start_time <- Sys.time()  
  counts.matrix1 <- Matrix::Matrix(counts.matrix, sparse = TRUE)
  message("Starting simulating dyngen")
  
  backbone <- backbone_bifurcating_loop()
  
  num_cells <- nrow(counts.matrix)
  num_feats <- ncol(counts.matrix)
  num_tfs <- nrow(backbone$module_info)
  num_tar <- round((num_feats - num_tfs) / 2)
  num_hks <- num_feats - num_tfs - num_tar
  
  # the simulation is being sped up because rendering all vignettes with one core
  # for pkgdown can otherwise take a very long time
  set.seed(1)
  
  config <-
    initialise_model(
      backbone = backbone,
      num_cells = num_cells,
      num_tfs = num_tfs,
      num_targets = num_tar,
      num_hks = num_hks,
      verbose = interactive(),
      download_cache_dir = tools::R_user_dir("dyngen", "data"),
      simulation_params = simulation_default(
        total_time = 1000,
        census_interval = 2, 
        ssa_algorithm = ssa_etl(tau = 300/3600),
        experiment_params = simulation_type_wild_type(num_simulations = 100)
      ),
      experiment_params = experiment_snapshot(
        realcount = counts.matrix1
      )
    )
  
  out <- generate_dataset(config, make_plots = TRUE)
  dyngen <- out$dataset$counts
  dyngen <- as.matrix(dyngen)
  
  end_time <- Sys.time()
  print(paste0("the time for running dyngen is "))
  print(end_time-start_time)
  save(dyngen, file="dyngen.RData") 
  
  
  
  ## SymSim ## 
  library(SymSim)
  message("Starting simulating SymSim")
  start_time <- Sys.time()  
  
  best_matches_UMI <- BestMatchParams_new('UMI',counts.matrix)
  true_counts_res <- SimulateTrueCounts(ncells_total=ncol(counts.matrix), ngenes=nrow(counts.matrix), evf_type="one.population",
                                        gene_effects_sd=best_matches_UMI$gene_effects_sd[1], Sigma=best_matches_UMI$Sigma[1],
                                        scale_s=best_matches_UMI$scale_s[1], gene_effect_prob=best_matches_UMI$gene_effect_prob[1],
                                        randseed=0)
  data(gene_len_pool)
  gene_len <- sample(gene_len_pool, nrow(counts.matrix), replace = FALSE)
  observed_counts <- True2ObservedCounts(true_counts=true_counts_res[[1]], meta_cell=true_counts_res[[3]], protocol="UMI", 
                                         alpha_mean=best_matches_UMI$alpha_mean, alpha_sd=best_matches_UMI$alpha_sd, 
                                         gene_len=gene_len, depth_mean=best_matches_UMI$depth_mean, depth_sd=best_matches_UMI$depth_sd)
  
  SymSim <- observed_counts
  SymSim <- SymSim$counts
  
  end_time <- Sys.time()
  print(paste0("the time for running SymSim is"))
  print(end_time-start_time)
  save(SymSim, file="SymSim.RData")
  
  
  # res <- list(SCRIP=SCRIP, scdesign=scdesign, powsimR=powsimR, SPARSim=SPARSim, scDD=scDD, dyngen=dyngen, SymSim=SymSim)
  # return(res)
}


#### Tung ####
library(devtools)
library(scDesign)
setwd("E:/DB/Dropbox/Qinfei/Simulation of SC based on splatter/Splatter simulator/Data from Tung")
Tung <- read.table("Tung.txt",sep="\t",header=T,row.names=1)
counts.matrix <- as.matrix(Tung)
librarysize <- median(apply(counts.matrix,2,sum))
Tung.res <- fun(counts.matrix)
setwd("E:\\DB\\Dropbox\\Qinfei\\Simulation of SC based on splatter\\Other methods")
save(Tung.res,file="Tung.simu.methods.RData")


#### Xin ####
EMTAB.eset <- readRDS("E:/DB/Dropbox/Qinfei/Simulation of SC based on splatter/Splatter simulator/Data from MuSiC/EMTABesethealthy.rds")
library(Biobase)
library(edgeR)
library(stats)
expre_data = exprs(EMTAB.eset)
pheno_data=pData(EMTAB.eset)
expre_data_params <- expre_data[,which((pheno_data$sampleID %in% c(5)) & (pheno_data$cellTypeID == 1))]
counts.matrix <- as.matrix(expre_data_params)
rownames(counts.matrix) <- paste0("gene",1:nrow(counts.matrix))
librarysize <- median(apply(counts.matrix,2,sum))
Xin.res <- fun(counts.matrix)
setwd("E:\\DB\\Dropbox\\Qinfei\\Simulation of SC based on splatter\\Other methods")
save(Xin.res,file="Xin.simu.methods.RData")


#### Klein ####
setwd("E:/DB/Dropbox/Qinfei/Simulation of SC based on splatter/Splatter simulator/Data from Klein")
Klein <- read.csv("Klein.csv",header=T,row.names = 1)
counts.matrix <- as.matrix(Klein)
rownames(counts.matrix) <- paste0("gene",1:nrow(counts.matrix))
librarysize <- median(apply(counts.matrix,2,sum))
Klein.res <- fun(counts.matrix)
setwd("E:\\DB\\Dropbox\\Qinfei\\Simulation of SC based on splatter\\Other methods")
save(Klein.res,file="Klein.simu.methods.RData")



#### Camp ####
setwd("E:/DB/Dropbox/Qinfei/Simulation of SC based on splatter/Splatter simulator/Data from Camp")
Camp <- read.table("Camp.txt",sep="\t",header=T,row.names=1)
Camp <- Camp[,c(7:ncol(Camp))]
expre_data_Camp <- Camp
colnames(expre_data_Camp) <- c(paste("Cell",1:ncol(expre_data_Camp),sep=""))
counts.matrix <- as.matrix(Camp)
librarysize <- median(apply(counts.matrix,2,sum))
Camp.res <- fun(counts.matrix)
setwd("E:\\DB\\Dropbox\\Qinfei\\Simulation of SC based on splatter\\Other methods")
save(Camp.res,file="Camp.simu.methods.RData")




#### Seale human
setwd("E:\\DB\\Dropbox\\Qinfei\\Simulation of SC based on splatter\\Submition\\Bioinformatics\\Simulation data\\Seale\\Human")
library(Seurat)

Seale_H <- readRDS(file="GSE128890_human.rds")
expre_data <- as.matrix(Seale_H@assays$RNA@counts)
celltype  <- as.character(Idents(Seale_H))
table(celltype)
expre_data <- expre_data[,celltype=="0"]
dim(expre_data)
expre_data1 <- expre_data[,1:3000]
counts.matrix1 <- as.matrix(expre_data1)
rownames(counts.matrix1) <- paste0("gene",1:nrow(counts.matrix))
librarysize <- median(apply(counts.matrix,2,sum))
fun(counts.matrix)



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
  

