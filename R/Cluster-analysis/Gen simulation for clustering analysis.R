setwd("/home/fqin/SCRIP/powsimR")
for (i in list.files()){
  source(i)
}

fun <- function(counts.matrix){
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
  parnGene <- params@nGenes
  parshape <- params@mean.shape
  parrate <- params@mean.rate
  parnCells <- params@batchCells
  
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
  
  message("Starting simulating SCRIP BGP-commonBCV")
  SCRIP.BGPcommonBCV <-  SCRIPsimu(data=counts.matrix, params=params, base_allcellmeans_SC=com_base_cellmeans,
                                   mode="BGP-commonBCV") 
  SCRIP.BGPcommonBCV <- counts(SCRIP.BGPcommonBCV)
  
  message("Starting simulating SCRIP BGP-trendedBCV")
  SCRIP.BGPtrendedBCV <-  SCRIPsimu(data=counts.matrix, params=params, base_allcellmeans_SC=com_base_cellmeans,
                                    mode="BGP-trendedBCV") 
  SCRIP.BGPtrendedBCV <- counts(SCRIP.BGPtrendedBCV)
  
  counts.matrix <- counts.matrix  + matrix(rgamma(ncol(counts.matrix)*nrow(counts.matrix),0.1, 1),nrow(counts.matrix),ncol(counts.matrix))
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
              SCRIP.BP=SCRIP.BP, SCRIP.BGPcommonBCV=SCRIP.BGPcommonBCV, SCRIP.BGPtrendedBCV=SCRIP.BGPtrendedBCV,
              scdesign=scdesign, powsimR=powsimR, SPARSim=SPARSim, dyngen=dyngen, SymSim=SymSim)
  return(res)
}







####Zhou
setwd("/home/fqin/SCRIP/Zhou")
library(Seurat)
Zhou <- readRDS(file="GSE139827_mice.rds")
expre_data <- as.matrix(Zhou@assays$RNA@counts)
celltype  <- as.character(Idents(Zhou))
table(celltype)

library(Biobase)
library(splatter)
CTlist <- unique(celltype)[5:6]
params <- splatEstimate(expre_data)
parnGene <- params@nGenes
parshape <- params@mean.shape
parrate <- params@mean.rate
parnCells <- params@batchCells
base_allcellmeans=rgamma(parnGene, shape=parshape, rate=parrate)

### all the simulators should be applied to each celltype from real dataset
### For cell types not clearly clustered in UMAP, can consider delete ##
setwd("/home/fqin/SCRIP/Cluster")
for (i in 1){
  for (CT in CTlist){
    set.seed(i)
    message(paste0("Starting simulating CT", CT))
    counts <- expre_data[,which(celltype==CT)]
    counts <- counts[,-c(1)]
    res <- fun(counts.matrix = counts)
    save(res, file=paste0(CT,".group.simulation",i,"Zhou.VEG1000Genes.addmean.RData"))
  }
}
