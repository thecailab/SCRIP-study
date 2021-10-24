## The functions below show how we generated simulation data for each celltype using SCRIP with VEGs and other simulators.
## Then simulation data for each cell type were combined as the final simulation data for each dataset. Tirosh dataset was used as example.


setwd("/home/fqin/SCRIP/powsimR")
for (i in list.files()){
  source(i)
}


fun <- function(counts.matrix, base_allcellmeans){
  message("Starting simulating SCRIP")
  rownames(counts.matrix) <- paste0("Gene",1:nrow(counts.matrix))
  colnames(counts.matrix) <- paste0("Cell",1:ncol(counts.matrix))
  
  library(Seurat)
  pbmc <- CreateSeuratObject(counts = counts.matrix, project = "pbmc3k", min.cells = 3, min.features = 200)
  seurat_stim <- NormalizeData(pbmc,
                               normalization.method = "LogNormalize",
                               scale.factor = 10000)
  stim <- FindVariableFeatures(object = seurat_stim,
                               selection.method = "vst",
                               nfeatures = 1000)
  top1000 <- head(VariableFeatures(stim),1000)
  EGcounts <- counts.matrix[top1000,]
  
  com_base_cellmeans <- base_allcellmeans
  names(com_base_cellmeans) <- rownames(counts.matrix)
  
  com_base_cellmeans[top1000] <- apply(EGcounts,1,mean)
  
  #### SCRIP ####
  library(splatter)
  params <- splatEstimate(counts.matrix)
  
  library(SCRIP)
  message("Starting simulating SCRIP GP-commonBCV")
  SCRIP.GPcommonBCV <- SCRIPsimu(data=counts.matrix, params=params, base_allcellmeans_SC=com_base_cellmeans,
                                 mode="GP-commonBCV") 
  SCRIP.GPcommonBCV <- counts(SCRIP.GPcommonBCV)
  
  message("Starting simulating SCRIP GP-trendedBCV")
  SCRIP.GPtrendedBCV <-  SCRIPsimu(data=counts.matrix, params=params, base_allcellmeans_SC=com_base_cellmeans,
                                   mode="GP-trendedBCV") 
  SCRIP.GPtrendedBCV <- counts(SCRIP.GPtrendedBCV)
  
  message("Starting simulating SCRIP BP")
  SCRIP.BP <- SCRIPsimu(data=counts.matrix, params=params, base_allcellmeans_SC=com_base_cellmeans,
                        mode="BP") 
  SCRIP.BP <- counts(SCRIP.BP)
  
  message("Starting simulating Splat")
  Splat <- SCRIPsimu(data=counts.matrix, params=params, mode="GP-commonBCV") 
  Splat <- counts(Splat)
  
  #### scDesign ####  
  message("Starting simulating scDesign")
  library(scDesign)
  librarysize <- sum(apply(counts.matrix,2,sum))
  scdesign <- design_data(counts.matrix, S = librarysize, ncell=parnCells, ngroup = 1, pUp = 0.05,
                          pDown = 0.05, fU = 5, fL = 1.5, ncores = 1)
  
  #### powsimR ####
  message("Starting simulating powsimR")
  estparam <- estimateParam(countData = counts.matrix,
                            readData = NULL,
                            MeanFragLengths = NULL,
                            RNAseq = 'singlecell', Protocol = 'Read',
                            Distribution = 'NB', Normalisation = "scran",
                            GeneFilter = 0.0, SampleFilter = 5,
                            sigma = 1.96, NCores = NULL, verbose = TRUE)
  plotParam(estparam, Annot = FALSE)
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
  
  ## SPARSim ##
  message("Starting simulating SPARSim")
  library(SPARSim)
  # count_norm <- scran_normalization(counts.matrix)
  count_conditions <- list(cond_A=c(1:(ncol(counts.matrix))))
  SPARSim_sim_param <- SPARSim_estimate_parameter_from_data(raw_data = counts.matrix,
                                                            norm_data =  counts.matrix,
                                                            conditions = count_conditions)
  SPARSim <- SPARSim_simulation(dataset_parameter = SPARSim_sim_param)
  SPARSim <- SPARSim$count_matrix
  
  
  ## scDD ##
  # library(scDD)
  # message("Starting simulating scDD")
  # # library(SingleCellExperiment)
  # condition <- sample(1:2, ncol(counts.matrix), replace = TRUE)
  # names(condition) <- colnames(counts.matrix)
  # scDD.param <- scDDEstimate(counts.matrix,conditions=condition)
  # 
  # scDD <- scDDSimulate(params = scDD.param)
  # scDD <- assay(scDD)
  
  
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
  
  
  res <- list(SCRIP.GPcommonBCV=SCRIP.GPcommonBCV, SCRIP.GPtrendedBCV=SCRIP.GPtrendedBCV,
              SCRIP.BP=SCRIP.BP, Splat=Splat, scdesign=scdesign, powsimR=powsimR, 
              SPARSim=SPARSim, dyngen=dyngen, SymSim=SymSim)
  return(res)
}



####Tirosh dataset Malanoma tumor
setwd("/home/fqin/SCRIP/Tirosh")
library(Seurat)
Tirosh <- readRDS(file="GSE72056_Melanoma.RDS")
expre_data <- as.matrix(Tirosh@assays$RNA@counts)
celltype  <- as.character(Idents(Tirosh))
table(celltype)

library(Biobase)
library(splatter)
CTlist <- unique(celltype)[-c(1,2,6)]

expre_data <- expre_data[,celltype %in% CTlist]
celltype <- celltype[celltype %in% CTlist]

library(splatter)
params <- splatEstimate(expre_data)
parnGene <- params@nGenes
parshape <- params@mean.shape
parrate <- params@mean.rate
parnCells <- params@batchCells
base_allcellmeans=rgamma(parnGene, shape=parshape, rate=parrate)

setwd("/home/fqin/SCRIP/Cluster")
GPcommon.final <- NULL
GPtrend.final <- NULL
BP.final <- NULL
Splat.final <- NULL
scdesign.final <- NULL
powsimR.final <- NULL
SAPRSim.final <- NULL
dyngen.final <- NULL
SymSim.final <- NULL
for (CT in CTlist){
  message(paste0("Starting simulating CT", CT))
  counts <- expre_data[,which(celltype==CT)]
  res <- fun(counts.matrix = counts, base_allcellmeans=base_allcellmeans)
  GPcommon.final <- cbind(GPcommon.final, res$SCRIP.GPcommonBCV)
  GPtrend.final <- cbind(GPtrend.final, res$SCRIP.GPtrendedBCV)
  BP.final <- cbind(BP.final, res$SCRIP.BP)
  Splat.final <- cbind(Splat.final, res$Splat)
  scdesign.final <- cbind(scdesign.final, res$scdesign)
  powsimR.final <- cbind(powsimR, res$powsimR)
  SAPRSim.final <- cbind(SAPRSim.final, res$SAPRSim)
  dyngen.final <- cbind(dyngen.final, res$dyngen)
  SymSim.final <- cbind(SymSim.final, res$SymSim)
}

save(GPcommon.final, file="GP-commonBCV.RData")
save(GPtrend.final, file="GP-trendedBCV.RData")
save(BP.final, file="BP.RData")
save(Splat.final, file="Splat.RData")
save(scdesign.final, file="scDesign.RData")
save(powsimR.final, file="powsimR.RData")
save(SAPRSim.final, file="SAPRSim.RData")
save(dyngen.final, file="dyngen.RData")
save(SymSim.final, file="SymSim.RData")

