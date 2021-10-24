## the functions below show how simulation data were generated for DE analysis methods. Also, after DE analyses were conducted with 
## methods including edgeR, DESeq2. Limma-voom, MAST, ZINB-Wave.edgeR, ZINB-Wave.DESeq2 and ZINB-Wave.limma_voom, how we describe the evaluation results

library("Biobase")
library(splatter)
library(gtools)
library(lcmix)


setwd("E:\\DB\\Dropbox\\Qinfei\\Research\\deconvolution\\deconvolution")
EMTAB.eset <- readRDS("EMTABesethealthy.rds")
expre_data <- exprs(EMTAB.eset)
pheno_data <- pData(EMTAB.eset)

expre_data_sub <- expre_data[,which((pheno_data$sampleID %in% c(5)) & (pheno_data$cellTypeID == 1))][1:1000,1:50]
pheno_data_sub <- pheno_data[which((pheno_data$sampleID %in% c(5)) & (pheno_data$cellTypeID == 1)),]

params1 <- splatEstimate(expre_data_sub)

nGenes <- params1@nGenes


base_allcellmeans=rgamma(nGenes, shape=params1@mean.shape, rate = params1@mean.rate)
res <- SCRIPsimu(data=expre_data_sub,params=params1)



#### Tung ####
library(devtools)
library(scDesign)
setwd("E:/DB/Dropbox/scRNA-deconvolution/Splatter simulator/Data from Tung")
Tung <- read.table("Tung.txt",sep="\t",header=T,row.names=1)
counts <- as.matrix(Tung)
library("Biobase")
library(splatter)
library(gtools)
library(lcmix)
expre_data_sub <- counts[1:2000,1:100]
params1 <- splatEstimate(expre_data_sub)
nGenes <- params1@nGenes





fun <- function(ngene,nDE,bulk_cell2,fold,Dropout_rate=NULL,libsize=NULL,pre.bcv.df=NULL,bcv.shrink=1){

  params2=setParams(params1,update=list("batchCells"=bulk_cell2))
  params2=setParams(params1,update=list("nGenes"=ngene))

  a=getParams(params2, c("mean.rate", "mean.shape"))
  rate=a[[1]]
  shape=a[[2]]
  ngene=getParams(params2,"nGenes")[[1]]
  library("TruncatedDistributions")
  set.seed(2020)
  base_allcellmeans=rtgamma(ngene, shape=shape, scale=1/rate, a=1, b=3)
  DEgene <- sample(1:length(base_allcellmeans),nDE,replace=F)


  base_allcellmeansDE <- base_allcellmeans
  base_allcellmeansDE[DEgene[1:nDE/2]] <- base_allcellmeans[DEgene[1:nDE/2]]*fold
  base_allcellmeansDE[DEgene[(nDE/2+1):nDE]] <- base_allcellmeans[DEgene[(nDE/2+1):nDE]]*1/fold

  params <- setParams(params2,  nGenes=nGenes, batchCells=bulk_cell2,  seed=100)
  
  sim <- SCRIPsimu(data=expre_data_sub, params=params, batchCells=bulk_cell2, CT.index=1, nc=1, libsize=libsize, method="single", bcv.shrink=bcv.shrink,
                                 base_allcellmeans_SC=base_allcellmeans,Dropout_rate=Dropout_rate,pre.bcv.df=pre.bcv.df,
                                 trend.m="locfit", mode="Trend nobursting") 
  exps <- counts(sim)
  
  
  simDE <-  SCRIPsimu(data=expre_data_sub, params=params, batchCells=bulk_cell2, CT.index=1, nc=1, libsize=libsize, method="single", bcv.shrink=bcv.shrink,
                                   base_allcellmeans_SC=base_allcellmeansDE,Dropout_rate=Dropout_rate,pre.bcv.df=pre.bcv.df,
                                   trend.m="locfit", mode="Trend nobursting") 
  expsDE <- counts(simDE)
  

  counts <- cbind(exps,expsDE)
  colnames(counts) <- paste0("cell",1:ncol(counts))
  rownames(counts) <- paste0("gene",1:nrow(counts))

  genenameDE <- paste0("gene",DEgene)
  res.data <- list(sim,simDE,genenameDE)
}



setwd(paste0("E:\\DB\\Dropbox\\Qinfei\\Simulation of SC based on splatter\\DE\\DEdata\\Tung"))
res1 <- fun(ngene=nGenes,nDE=200,bulk_cell2=100,fold=1.2)
save(res1,file="simu data2000 fold1.2.RData")
res2 <- fun(ngene=nGenes,nDE=200,bulk_cell2=100,fold=1.4)
save(res2,file="simu data2000 fold1.4.RData")
res3 <- fun(ngene=nGenes,nDE=200,bulk_cell2=100,fold=1.6)
save(res3,file="simu data2000 fold1.6.RData")
res4 <- fun(ngene=nGenes,nDE=200,bulk_cell2=100,fold=1.8)
save(res4,file="simu data2000 fold1.8.RData")


setwd(paste0("E:\\DB\\Dropbox\\Qinfei\\Simulation of SC based on splatter\\DE\\DEdata\\Tung"))
res1 <- fun(ngene=nGenes,nDE=200,bulk_cell2=100,fold=1.4,Dropout_rate=0.1)
save(res1,file="simu data2000 fold1.4 Drop0.1.RData")
res2 <- fun(ngene=nGenes,nDE=200,bulk_cell2=100,fold=1.4,Dropout_rate=0.2)
save(res2,file="simu data2000 fold1.4 Drop0.2.RData")
res3 <- fun(ngene=nGenes,nDE=200,bulk_cell2=100,fold=1.4,Dropout_rate=0.3)
save(res3,file="simu data2000 fold1.4 Drop0.3.RData")
res4 <- fun(ngene=nGenes,nDE=200,bulk_cell2=100,fold=1.4,Dropout_rate=0.4)
save(res4,file="simu data2000 fold1.4 Drop0.4.RData")


setwd(paste0("E:\\DB\\Dropbox\\Qinfei\\Simulation of SC based on splatter\\DE\\DEdata\\Tung"))
res1 <- fun(ngene=nGenes,nDE=200,bulk_cell2=100,fold=1.4,libsize=5000)
save(res1,file="simu data2000 fold1.4 lib5k.RData")
res2 <- fun(ngene=nGenes,nDE=200,bulk_cell2=100,fold=1.4,libsize=10000)
save(res2,file="simu data2000 fold1.4 lib10k.RData")
res3 <- fun(ngene=nGenes,nDE=200,bulk_cell2=100,fold=1.4,libsize=20000)
save(res3,file="simu data2000 fold1.4 lib20k.RData")
res4 <- fun(ngene=nGenes,nDE=200,bulk_cell2=100,fold=1.4,libsize=40000)
save(res4,file="simu data2000 fold1.4 lib40k.RData")


setwd(paste0("E:\\DB\\Dropbox\\Qinfei\\Simulation of SC based on splatter\\DE\\DEdata\\Tung"))
res2 <- fun(ngene=nGenes,nDE=200,bulk_cell2=100,fold=1.4,pre.bcv.df=2)
save(res2,file="simu data2000 fold1.4 bcv.df2.RData")
res1 <- fun(ngene=nGenes,nDE=200,bulk_cell2=100,fold=1.4,pre.bcv.df=5)
save(res1,file="simu data2000 fold1.4 bcv.df5.RData")
res1 <- fun(ngene=nGenes,nDE=200,bulk_cell2=100,fold=1.4,pre.bcv.df=10)
save(res1,file="simu data2000 fold1.4 bcv.df10.RData")
res2 <- fun(ngene=nGenes,nDE=200,bulk_cell2=100,fold=1.4,pre.bcv.df=20)
save(res2,file="simu data2000 fold1.4 bcv.df20.RData")



setwd(paste0("E:\\DB\\Dropbox\\Qinfei\\Simulation of SC based on splatter\\DE\\DEdata\\Tung"))
res2 <- fun(ngene=nGenes,nDE=200,bulk_cell2=100,fold=1.4,bcv.shrink = 1.0)
save(res2,file="simu data2000 fold1.4 bcv.shrink1.0.RData")
res7 <- fun(ngene=nGenes,nDE=200,bulk_cell2=100,fold=1.4,bcv.shrink = 0.5)
save(res7,file="simu data2000 fold1.4 bcv.shrink0.5.RData")
res8 <- fun(ngene=nGenes,nDE=200,bulk_cell2=100,fold=1.4,bcv.shrink = 1.5)
save(res8,file="simu data2000 fold1.4 bcv.shrink1.5.RData")
res9 <- fun(ngene=nGenes,nDE=200,bulk_cell2=100,fold=1.4,bcv.shrink = 2.0)
save(res9,file="simu data2000 fold1.4 bcv.shrink2.0.RData")








fun_DE <- function(sim, simDE, genenameDE) {
  library("Biobase")
  library(splatter)
  library(gtools)
  library(lcmix)  
  exps <- counts(sim)
  expsDE <- counts(simDE)
  counts <- cbind(exps,expsDE)
  sum=apply(counts,1,sum)
  counts=counts[which(sum>0),]

  colnames(counts) <- paste0("cell",1:ncol(counts))
  rownames(counts) <- paste0("gene",1:nrow(counts))
  
  truth <- rownames(counts)
  truth[which(rownames(counts) %in% genenameDE)] <- "DE"
  truth[-which(rownames(counts) %in% genenameDE)] <- "common"
  truth <- as.factor(truth)
  
##################### DE methods ##############################
#################### ZINB WAVE##################################
  library(zinbwave)
  library(DESeq2)
  group <- factor(c(rep(1,ncol(exps)),rep(2,ncol(expsDE))))
  coldata <- data.frame(condition=group)
  rownames(coldata) <- colnames(counts)
  fluidigm <- DESeqDataSetFromMatrix(countData = counts,
                                     colData = coldata,
                                     design = ~ condition)
  zinb <- zinbFit(fluidigm, K=2, epsilon=1000)
  fluidigm_zinb <- zinbwave(fluidigm, fitted_model = zinb, K = 2, epsilon=1000,
                          observationalWeights = TRUE)
  weights <- assay(fluidigm_zinb, "weights")

  
  
  ### DESeq2 ##
  #############
  library(DESeq2)
  dds <- DESeqDataSet(fluidigm, design = ~ condition)
  dds <- DESeq(dds, sfType="poscounts", useT=TRUE, minmu=1e-6)
  
  res <- results(dds)
  res <- as.data.frame(res)
  DESeq2=res
  DESeq2 <- DESeq2[rownames(counts),]
  
  
  
  ### DESeq2 zinbwave ##
  #############
  library(DESeq2)
  dds <- DESeqDataSet(fluidigm_zinb, design = ~ condition)
  dds <- DESeq(dds, sfType="poscounts", useT=TRUE, minmu=1e-6)
  
  res <- results(dds)
  res <- as.data.frame(res)
  DESeq2.weighted=res
  DESeq2.weighted <- DESeq2.weighted[rownames(counts),]
  
  
  
  ### edgeR ###
  library(edgeR)
  dgList <- DGEList(assay(fluidigm))
  countsPerMillion <- cpm(dgList)
  countCheck <- countsPerMillion > 1
  keep <- which(rowSums(countCheck) >= 2)
  dgList <- dgList[keep,]
  dgList <- calcNormFactors(dgList, method="TMM")
  design <- model.matrix(~condition, data = colData(fluidigm))
  dgList <- estimateDisp(dgList, design)
  fit <- glmFit(dgList, design)
  lrt <- glmLRT(fit, coef=2)
  edgeR <- as.data.frame(topTags(lrt,n=nrow(lrt$coefficients)))
  edgeR <- edgeR[rownames(counts),]
  
  
  
  #### edgeR zinbwave ####
  library(edgeR)
  dge <- DGEList(assay(fluidigm_zinb))
  dge <- calcNormFactors(dge)

  design <- model.matrix(~condition, data = colData(fluidigm))
  dge$weights <- weights
  dge <- estimateDisp(dge, design)
  fit <- glmFit(dge, design)
  lrt <- glmWeightedF(fit, coef = 2)
  edgeR.weighted <- as.data.frame(topTags(lrt,n=nrow(lrt$coefficients)))
  edgeR.weighted <- edgeR.weighted[rownames(counts),]
  # lrt <- lrt[which(lrt$PValue<0.05),]
  # pred <- rownames(counts)
  # pred[which(pred %in% rownames(lrt))] <- "DE"
  # pred[which(pred != "DE")] <- "common"
  #
  # pred <- as.factor(pred)
  # 
  # edgeR <- pred
  
  # edgeR <- table(truth,pred)
  # edgeR <- data.frame(edgeR0=edgeR[,1],edgeR1=edgeR[,2])


  ####DEsingle####
  library(DEsingle)
  library(stats)
  result <- DEsingle(counts,group,parallel = TRUE)
  DEsingle <- result
  DEsingle <- DEsingle[rownames(counts),]
  # results.classified <- DEtype(results = result, threshold = 0.05)
  # results.sig <- results.classified[results.classified$pvalue.adj.FDR < 0.05, ]
  # results.sig
  # pred <- rownames(counts)
  # pred[which(pred %in% rownames(results.sig))] <- "DE"
  # pred[which(pred != "DE")] <- "common"
  # 
  # pred <- as.factor(pred)
  # DEsingle <- pred
  # truth <- rownames(counts)
  # truth[DEgene] <- "DE"
  # truth[-DEgene] <- "common"
  # truth <- as.factor(truth)
  # 
  # DEsingle <- table(truth,pred)
  # DEsingle <- data.frame(DEsingle0=DEsingle[,1],DEsingle1=DEsingle[,2])



  #### limma voom
  library(edgeR)
  d0 <- DGEList(counts)
  d0 <- calcNormFactors(d0)
  
  mm <- model.matrix(~0 + group)
  y <- voom(d0, mm, plot = F)

  fit <- lmFit(y, mm)

  contr <- makeContrasts(group1 - group2, levels = colnames(coef(fit)))
  tmp <- contrasts.fit(fit, contr)
  tmp <- eBayes(tmp)
  top.table <- topTable(tmp, sort.by = "none", n = Inf)
  limma <- top.table
  limma <- limma[rownames(counts),]
  # top.table <- top.table[which(top.table$adj.P.Val < 0.05),]
  # top.table
  # pred <- rownames(counts)
  # pred[which(pred %in% rownames(top.table))] <- "DE"
  # pred[which(pred != "DE")] <- "common"
  # pred <- as.factor(pred)
  # limma <- pred
  # truth <- rownames(counts)
  # truth[DEgene] <- "DE"
  # truth[-DEgene] <- "common"
  # truth <- as.factor(truth)
  # 
  # limma <- table(truth,pred)
  # limma <- data.frame(limma0=limma[,1],limma1=limma[,2])
  
  #### limma voom ZINB WAVE
  library(edgeR)
  d0.weights <- DGEList(counts)
  d0.weights <- calcNormFactors(d0.weights)
  
  mm <- model.matrix(~0 + group)
  y.weights <- voom(d0.weights, mm, plot = F)
  y.weights$weights <- weights*y.weights$weights
  
  fit.weights <- lmFit(y.weights, mm)
  
  contr.weights <- makeContrasts(group1 - group2, levels = colnames(coef(fit.weights)))
  tmp.weights <- contrasts.fit(fit.weights, contr.weights)
  tmp.weights <- eBayes(tmp.weights)
  top.table.weights <- topTable(tmp.weights, sort.by = "none", n = Inf)
  limma.weights <- top.table.weights
  limma.weighted <- limma.weights[rownames(counts),]
  
  


  #### MAST
  group <- factor(c(rep("Group0",ncol(exps)),rep("Group1",ncol(expsDE))))
  group <- relevel(group,"Group0")
  coldata <- data.frame(condition=group)
  rownames(coldata) <- colnames(counts)
  fluidigm <- DESeqDataSetFromMatrix(countData = counts,
                                   colData = coldata,
                                   design = ~ condition)
  library(MAST)
  library(data.table)
  sca <- FromMatrix(counts,data.frame(condition=group),
                  data.frame(Geneid=rownames(counts)),check_sanity = FALSE)
  zlmCond <- zlm(~condition, sca)

  summaryCond <- summary(zlmCond, doLRT='conditionGroup1') 
  print(summaryCond, n=4)
  summaryDt <- summaryCond$datatable
  fcHurdle <- merge(summaryDt[contrast=='conditionGroup1' & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                    summaryDt[contrast=='conditionGroup1' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients
  MAST=fcHurdle
  MAST <- MAST[rownames(counts),]

  res1 <- list(edgeR=edgeR, DESeq2=DESeq2, edgeR.weighted=edgeR.weighted, DESeq2.weighted=DESeq2.weighted,
               limma=limma,  limma.weighted=limma.weighted, DEsingle=DEsingle,  MAST=MAST, truth=truth)
  return(res1)
}




# 
# 
# limma_wave <- function(sim, simDE, genenameDE){
# 
# library("Biobase")
# library(splatter)
# library(gtools)
# library(lcmix)  
# exps <- counts(sim)
# expsDE <- counts(simDE)
# counts <- cbind(exps,expsDE)
# sum=apply(counts,1,sum)
# counts=counts[which(sum>0),]
# 
# colnames(counts) <- paste0("cell",1:ncol(counts))
# rownames(counts) <- paste0("gene",1:nrow(counts))
# 
# truth <- rownames(counts)
# truth[which(rownames(counts) %in% genenameDE)] <- "DE"
# truth[-which(rownames(counts) %in% genenameDE)] <- "common"
# truth <- as.factor(truth)
# 
# ##################### DE methods ##############################
# #################### ZINB WAVE##################################
# library(zinbwave)
# library(DESeq2)
# group <- factor(c(rep(1,ncol(exps)),rep(2,ncol(expsDE))))
# coldata <- data.frame(condition=group)
# rownames(coldata) <- colnames(counts)
# fluidigm <- DESeqDataSetFromMatrix(countData = counts,
#                                    colData = coldata,
#                                    design = ~ condition)
# zinb <- zinbFit(fluidigm, K=2, epsilon=1000)
# fluidigm_zinb <- zinbwave(fluidigm, fitted_model = zinb, K = 2, epsilon=1000,
#                           observationalWeights = TRUE)
# weights <- assay(fluidigm_zinb, "weights")
# 
# 
# #### limma voom
# library(edgeR)
# d0 <- DGEList(counts)
# d0 <- calcNormFactors(d0)
# 
# mm <- model.matrix(~0 + group)
# y <- voom(d0, mm, plot = F)
# y$weights <- weights*y$weights
# 
# fit <- lmFit(y, mm)
# 
# contr <- makeContrasts(group1 - group2, levels = colnames(coef(fit)))
# tmp <- contrasts.fit(fit, contr)
# tmp <- eBayes(tmp)
# top.table <- topTable(tmp, sort.by = "none", n = Inf)
# limma <- top.table
# limma.weighted <- limma[rownames(counts),]
# 
# return(limma.weighted)
# }
# 



### Fold Change
setwd(paste0("E:\\DB\\Dropbox\\Qinfei\\Simulation of SC based on splatter\\DE\\DEdata\\Tung"))
res1.2 = get(load(file="simu data2000 fold1.2.RData"))
res1.2 <- fun_DE(sim=res1.2[[1]], simDE=res1.2[[2]], genenameDE=res1.2[[3]])

res1.4 = get(load(file="simu data2000 fold1.4.RData"))
res1.4 <- fun_DE(sim=res1.4[[1]], simDE=res1.4[[2]], genenameDE=res1.4[[3]])

res1.6 = get(load(file="simu data2000 fold1.6.RData"))
res1.6 <- fun_DE(sim=res1.6[[1]], simDE=res1.6[[2]], genenameDE=res1.6[[3]])

res1.8 = get(load(file="simu data2000 fold1.8.RData"))
res1.8 <- fun_DE(sim=res1.8[[1]], simDE=res1.8[[2]], genenameDE=res1.8[[3]])

load(paste0("DE2000 results fold1.4 ",name[1],".RData"))




# save(res1.2,file="DE2000 results fold1.2.RData")
# save(res1.4,file="DE2000 results fold1.4.RData")
# save(res1.6,file="DE2000 results fold1.6.RData")
# save(res1.8,file="DE2000 results fold1.8.RData")


setwd(paste0("E:\\DB\\Dropbox\\Qinfei\\Simulation of SC based on splatter\\DE\\DEdata\\Tung"))
res1.2 = get(load(file="simu data2000 fold1.2.RData"))
res1.4 = get(load(file="simu data2000 fold1.4.RData"))
res1.6 = get(load(file="simu data2000 fold1.6.RData"))
res1.8 = get(load(file="simu data2000 fold1.8.RData"))


FC1.2 <- limma_wave(sim=res1.2[[1]], simDE=res1.2[[2]], genenameDE=res1.2[[3]])
FC1.4 <- limma_wave(sim=res1.4[[1]], simDE=res1.4[[2]], genenameDE=res1.4[[3]])
FC1.6 <- limma_wave(sim=res1.6[[1]], simDE=res1.6[[2]], genenameDE=res1.6[[3]])
FC1.8 <- limma_wave(sim=res1.8[[1]], simDE=res1.8[[2]], genenameDE=res1.8[[3]])


setwd("E:\\DB\\Dropbox\\Qinfei\\Simulation of SC based on splatter\\DE\\DEdata\\Tung")
result1 <- get(load("DE2000 results fold1.2.RData"))
result2 <- get(load("DE2000 results fold1.4.RData"))
result3 <- get(load("DE2000 results fold1.6.RData"))
result4 <- get(load("DE2000 results fold1.8.RData"))


result1[["limma.weighted"]] <- FC1.2
save(result1,file=paste0("DE2000 results fold1.2 new.RData"))

result2[["limma.weighted"]] <- FC1.4
save(result2,file=paste0("DE2000 results fold1.4 new.RData"))

result3[["limma.weighted"]] <- FC1.6
save(result3,file=paste0("DE2000 results fold1.6 new.RData"))

result4[["limma.weighted"]] <- FC1.8
save(result4,file=paste0("DE2000 results fold1.8 new.RData"))









### Dropout;  
name=c("Drop0.1","Drop0.2","Drop0.3","Drop0.4")

setwd(paste0("E:\\DB\\Dropbox\\Qinfei\\Simulation of SC based on splatter\\DE\\DEdata\\Tung"))
DEres1 = get(load(file=paste0("simu data2000 fold1.4 ",name[1],".RData")))
DEres2 = get(load(file=paste0("simu data2000 fold1.4 ",name[2],".RData")))
DEres3 = get(load(file=paste0("simu data2000 fold1.4 ",name[3],".RData")))
DEres4 = get(load(file=paste0("simu data2000 fold1.4 ",name[4],".RData")))

result1 <- fun_DE(sim=DEres1[[1]], simDE=DEres1[[2]], genenameDE=DEres1[[3]])
setwd("E:\\DB\\Dropbox\\Qinfei\\Simulation of SC based on splatter\\DE\\DEdata\\Tung")
save(result1,file=paste0("DE2000 results fold1.4 ",name[1],".RData"))

result2 <- fun_DE(sim=DEres2[[1]], simDE=DEres2[[2]], genenameDE=DEres2[[3]])
setwd("E:\\DB\\Dropbox\\Qinfei\\Simulation of SC based on splatter\\DE\\DEdata\\Tung")
save(result2,file=paste0("DE2000 results fold1.4 ",name[2],".RData"))

result3 <- fun_DE(sim=DEres3[[1]], simDE=DEres3[[2]], genenameDE=DEres3[[3]])
setwd("E:\\DB\\Dropbox\\Qinfei\\Simulation of SC based on splatter\\DE\\DEdata\\Tung")
save(result3,file=paste0("DE2000 results fold1.4 ",name[3],".RData"))

result4 <- fun_DE(sim=DEres4[[1]], simDE=DEres4[[2]], genenameDE=DEres4[[3]])
setwd("E:\\DB\\Dropbox\\Qinfei\\Simulation of SC based on splatter\\DE\\DEdata\\Tung")
save(result4,file=paste0("DE2000 results fold1.4 ",name[4],".RData"))


name=c("Drop0.1","Drop0.2","Drop0.3","Drop0.4")
dropout1 <- limma_wave(sim=DEres1[[1]], simDE=DEres1[[2]], genenameDE=DEres1[[3]])
dropout2 <- limma_wave(sim=DEres2[[1]], simDE=DEres2[[2]], genenameDE=DEres2[[3]])
dropout3 <- limma_wave(sim=DEres3[[1]], simDE=DEres3[[2]], genenameDE=DEres3[[3]])
dropout4 <- limma_wave(sim=DEres4[[1]], simDE=DEres4[[2]], genenameDE=DEres4[[3]])
load(paste0("DE2000 results fold1.4 ",name[1],".RData"))
result1 <- append(result1,list(limma.weighted=dropout1))
save(result1,file=paste0("DE2000 results fold1.4 ",name[1],"new.RData"))

load(paste0("DE2000 results fold1.4 ",name[2],".RData"))
result2 <- append(result2,list(limma.weighted=dropout2))
save(result2,file=paste0("DE2000 results fold1.4 ",name[2],"new.RData"))

load(paste0("DE2000 results fold1.4 ",name[3],".RData"))
result3 <- append(result3,list(limma.weighted=dropout3))
save(result3,file=paste0("DE2000 results fold1.4 ",name[3],"new.RData"))

load(paste0("DE2000 results fold1.4 ",name[4],".RData"))
result4 <- append(result4,list(limma.weighted=dropout4))
save(result4,file=paste0("DE2000 results fold1.4 ",name[4],"new.RData"))



#### library size;
# name=c("lib5k","lib10k","lib20k","lib40k")
# 
# setwd(paste0("E:\\DB\\Dropbox\\Qinfei\\Simulation of SC based on splatter\\DE\\DEdata\\Tung"))
# DEres1 = get(load(file=paste0("simu data2000 fold1.4 ",name[1],".RData")))
# DEres2 = get(load(file=paste0("simu data2000 fold1.4 ",name[2],".RData")))
# DEres3 = get(load(file=paste0("simu data2000 fold1.4 ",name[3],".RData")))
# DEres4 = get(load(file=paste0("simu data2000 fold1.4 ",name[4],".RData")))
# 
# result1 <- fun_DE(sim=DEres1[[1]], simDE=DEres1[[2]], genenameDE=DEres1[[3]])
# setwd("E:\\DB\\Dropbox\\Qinfei\\Simulation of SC based on splatter\\DE\\DEdata\\Tung")
# save(result1,file=paste0("DE2000 results fold1.4 ",name[1],".RData"))
# 
# result2 <- fun_DE(sim=DEres2[[1]], simDE=DEres2[[2]], genenameDE=DEres2[[3]])
# setwd("E:\\DB\\Dropbox\\Qinfei\\Simulation of SC based on splatter\\DE\\DEdata\\Tung")
# save(result2,file=paste0("DE2000 results fold1.4 ",name[2],".RData"))
# 
# result3 <- fun_DE(sim=DEres3[[1]], simDE=DEres3[[2]], genenameDE=DEres3[[3]])
# setwd("E:\\DB\\Dropbox\\Qinfei\\Simulation of SC based on splatter\\DE\\DEdata\\Tung")
# save(result3,file=paste0("DE2000 results fold1.4 ",name[3],".RData"))
# 
# result4 <- fun_DE(sim=DEres4[[1]], simDE=DEres4[[2]], genenameDE=DEres4[[3]])
# setwd("E:\\DB\\Dropbox\\Qinfei\\Simulation of SC based on splatter\\DE\\DEdata\\Tung")
# save(result4,file=paste0("DE2000 results fold1.4 ",name[4],".RData"))



name=c("lib5k","lib10k","lib20k","lib40k")
libsiz1 <- limma_wave(sim=DEres1[[1]], simDE=DEres1[[2]], genenameDE=DEres1[[3]])
libsiz2 <- limma_wave(sim=DEres2[[1]], simDE=DEres2[[2]], genenameDE=DEres2[[3]])
libsiz3 <- limma_wave(sim=DEres3[[1]], simDE=DEres3[[2]], genenameDE=DEres3[[3]])
libsiz4 <- limma_wave(sim=DEres4[[1]], simDE=DEres4[[2]], genenameDE=DEres4[[3]])
load(paste0("DE2000 results fold1.4 ",name[1],".RData"))
result1 <- append(result1,list(limma.weighted=libsiz1))
save(result1,file=paste0("DE2000 results fold1.4 ",name[1],"new.RData"))

load(paste0("DE2000 results fold1.4 ",name[2],".RData"))
result2 <- append(result2, list(limma.weighted=libsiz2))
save(result2,file=paste0("DE2000 results fold1.4 ",name[2],"new.RData"))

load(paste0("DE2000 results fold1.4 ",name[3],".RData"))
result3 <- append(result3, list(limma.weighted=libsiz3))
save(result3,file=paste0("DE2000 results fold1.4 ",name[3],"new.RData"))

load(paste0("DE2000 results fold1.4 ",name[4],".RData"))
result4 <- append(result4, list(limma.weighted=libsiz4))
save(result4,file=paste0("DE2000 results fold1.4 ",name[4],"new.RData"))




#### BCV.df 
name=c("BCV.df2","BCV.df5","BCV.df10","BCV.df20")

setwd(paste0("E:\\DB\\Dropbox\\Qinfei\\Simulation of SC based on splatter\\DE\\DEdata\\Tung"))
DEres1 = get(load(file=paste0("simu data2000 fold1.4 ",name[1],".RData")))
DEres2 = get(load(file=paste0("simu data2000 fold1.4 ",name[2],".RData")))
DEres3 = get(load(file=paste0("simu data2000 fold1.4 ",name[3],".RData")))
DEres4 = get(load(file=paste0("simu data2000 fold1.4 ",name[4],".RData")))


result1 <- fun_DE(sim=DEres1[[1]], simDE=DEres1[[2]], genenameDE=DEres1[[3]])
setwd("E:\\DB\\Dropbox\\Qinfei\\Simulation of SC based on splatter\\DE\\DEdata\\Tung")
save(result1,file=paste0("DE2000 results fold1.4 ",name[1],"new.RData"))

result2 <- fun_DE(sim=DEres2[[1]], simDE=DEres2[[2]], genenameDE=DEres2[[3]])
setwd("E:\\DB\\Dropbox\\Qinfei\\Simulation of SC based on splatter\\DE\\DEdata\\Tung")
save(result2,file=paste0("DE2000 results fold1.4 ",name[2],"new.RData"))

result3 <- fun_DE(sim=DEres3[[1]], simDE=DEres3[[2]], genenameDE=DEres3[[3]])
setwd("E:\\DB\\Dropbox\\Qinfei\\Simulation of SC based on splatter\\DE\\DEdata\\Tung")
save(result3,file=paste0("DE2000 results fold1.4 ",name[3],"new.RData"))

result4 <- fun_DE(sim=DEres4[[1]], simDE=DEres4[[2]], genenameDE=DEres4[[3]])
setwd("E:\\DB\\Dropbox\\Qinfei\\Simulation of SC based on splatter\\DE\\DEdata\\Tung")
save(result4,file=paste0("DE2000 results fold1.4 ",name[4],"new.RData"))


# name=c("BCV.df2","BCV.df5","BCV.df10","BCV.df20")
# BCV.df1 <- limma_wave(sim=DEres1[[1]], simDE=DEres1[[2]], genenameDE=DEres1[[3]])
# BCV.df2 <- limma_wave(sim=DEres2[[1]], simDE=DEres2[[2]], genenameDE=DEres2[[3]])
# BCV.df3 <- limma_wave(sim=DEres3[[1]], simDE=DEres3[[2]], genenameDE=DEres3[[3]])
# BCV.df4 <- limma_wave(sim=DEres4[[1]], simDE=DEres4[[2]], genenameDE=DEres4[[3]])
# 
# result1 <- get(load(paste0("DE2000 results fold1.4 ",name[1],"new.RData")))
# result1[["limma.weighted"]] <- BCV.df1
# save(result1,file=paste0("DE2000 results fold1.4 ",name[1],"new.RData"))
# 
# result2 <- get(load(paste0("DE2000 results fold1.4 ",name[2],"new.RData")))
# result2[["limma.weighted"]] <- BCV.df2
# save(result2,file=paste0("DE2000 results fold1.4 ",name[2],"new.RData"))
# 
# result3 <- get(load(paste0("DE2000 results fold1.4 ",name[3],"new.RData")))
# result3[["limma.weighted"]] <- BCV.df3
# save(result3,file=paste0("DE2000 results fold1.4 ",name[3],"new.RData"))
# 
# result4 <- get(load(paste0("DE2000 results fold1.4 ",name[4],"new.RData")))
# result4[["limma.weighted"]] <- BCV.df4
# save(result4,file=paste0("DE2000 results fold1.4 ",name[4],"new.RData"))




#### BCV.shink 
name=c("bcv.shrink0.5","bcv.shrink1.0","bcv.shrink1.5","bcv.shrink2.0")

setwd(paste0("E:\\DB\\Dropbox\\Qinfei\\Simulation of SC based on splatter\\DE\\DEdata\\Tung"))
DEres1 = get(load(file=paste0("simu data2000 fold1.4 ",name[1],".RData")))
DEres2 = get(load(file=paste0("simu data2000 fold1.4 ",name[2],".RData")))
DEres3 = get(load(file=paste0("simu data2000 fold1.4 ",name[3],".RData")))
DEres4 = get(load(file=paste0("simu data2000 fold1.4 ",name[4],".RData")))


result1 <- fun_DE(sim=DEres1[[1]], simDE=DEres1[[2]], genenameDE=DEres1[[3]])
setwd("E:\\DB\\Dropbox\\Qinfei\\Simulation of SC based on splatter\\DE\\DEdata\\Tung")
save(result1,file=paste0("DE2000 results fold1.4 ",name[1],"new.RData"))

result2 <- fun_DE(sim=DEres2[[1]], simDE=DEres2[[2]], genenameDE=DEres2[[3]])
setwd("E:\\DB\\Dropbox\\Qinfei\\Simulation of SC based on splatter\\DE\\DEdata\\Tung")
save(result2,file=paste0("DE2000 results fold1.4 ",name[2],"new.RData"))

result3 <- fun_DE(sim=DEres3[[1]], simDE=DEres3[[2]], genenameDE=DEres3[[3]])
setwd("E:\\DB\\Dropbox\\Qinfei\\Simulation of SC based on splatter\\DE\\DEdata\\Tung")
save(result3,file=paste0("DE2000 results fold1.4 ",name[3],"new.RData"))

result4 <- fun_DE(sim=DEres4[[1]], simDE=DEres4[[2]], genenameDE=DEres4[[3]])
setwd("E:\\DB\\Dropbox\\Qinfei\\Simulation of SC based on splatter\\DE\\DEdata\\Tung")
save(result4,file=paste0("DE2000 results fold1.4 ",name[4],"new.RData"))



BCV.shink1 <- limma_wave(sim=DEres1[[1]], simDE=DEres1[[2]], genenameDE=DEres1[[3]])
BCV.shink2 <- limma_wave(sim=DEres2[[1]], simDE=DEres2[[2]], genenameDE=DEres2[[3]])
BCV.shink3 <- limma_wave(sim=DEres3[[1]], simDE=DEres3[[2]], genenameDE=DEres3[[3]])
BCV.shink4 <- limma_wave(sim=DEres4[[1]], simDE=DEres4[[2]], genenameDE=DEres4[[3]])

result1 <- get(load(paste0("DE2000 results fold1.4 ",name[1],"new.RData")))
result1[["limma.weighted"]] <- BCV.shink1
save(result1,file=paste0("DE2000 results fold1.4 ",name[1],"new.RData"))

result2 <- get(load(paste0("DE2000 results fold1.4 ",name[2],"new.RData")))
result2[["limma.weighted"]] <- BCV.shink2
save(result2,file=paste0("DE2000 results fold1.4 ",name[2],"new.RData"))

result3 <- get(load(paste0("DE2000 results fold1.4 ",name[3],"new.RData")))
result3[["limma.weighted"]] <- BCV.shink3
save(result3,file=paste0("DE2000 results fold1.4 ",name[3],"new.RData"))

result4 <- get(load(paste0("DE2000 results fold1.4 ",name[4],"new.RData")))
result4[["limma.weighted"]] <- BCV.shink4
save(result4,file=paste0("DE2000 results fold1.4 ",name[4],"new.RData"))











library(ROCit)
library(ROCR)
tpfp.gen <- function(res,nDE){
  
  truth <- res$truth

  TPFP <- NULL
  for (alpha in (seq(0,1000000,5000))/1000000){
  # for (alpha in 0.05){
  
    edgeR <- res$edgeR
    lrt <- edgeR[which(edgeR$PValue <= alpha),]
    pred <- rownames(edgeR)
    pred[which(pred %in% rownames(lrt))] <- "DE"
    pred[which(pred != "DE")] <- "common"
    tab <- table(pred,truth)
    if (dim(tab)[1]==2) {
    edgeR.tp <- tab[2,2]/nDE
    edgeR.fp <- tab[2,1]/(length(pred)-nDE)
    }
    if (dim(tab)[1]==1) {
      if (rownames(tab)=="common"){
        edgeR.tp <- 0
        edgeR.fp <- 0
      } 
      if (rownames(tab)=="DE"){
        edgeR.tp <- 1
        edgeR.fp <- 1
      }
    } 
    
    edgeR.zinb <- res$edgeR.weighted
    lrt <- edgeR.zinb[which(edgeR.zinb$padjFilter <= alpha),]
    pred <- rownames(edgeR)
    pred[which(pred %in% rownames(lrt))] <- "DE"
    pred[which(pred != "DE")] <- "common"
    tab <- table(pred,truth)
    if (dim(tab)[1]==2) {
      edgeR.zinb.tp <- tab[2,2]/nDE
      edgeR.zinb.fp <- tab[2,1]/(length(pred)-nDE)
    }
    if (dim(tab)[1]==1) {
      if (rownames(tab)=="common"){
        edgeR.zinb.tp <- 0
        edgeR.zinb.fp <- 0
      } 
      if (rownames(tab)=="DE"){
        edgeR.zinb.tp <- 1
        edgeR.zinb.fp <- 1
      }
    } 
    

    DESeq2 <- res$DESeq2
    lrt <- DESeq2[which(DESeq2$padj <= alpha),]
    pred <- rownames(DESeq2)
    pred[which(pred %in% rownames(lrt))] <- "DE"
    pred[which(pred != "DE")] <- "common"
    tab <- table(pred,truth)
    if (dim(tab)[1]==2) {
      DESeq2.tp <- tab[2,2]/nDE
      DESeq2.fp <- tab[2,1]/(length(pred)-nDE)
    }
    if (dim(tab)[1]==1) {
      if (rownames(tab)=="common"){
        DESeq2.tp <- 0
        DESeq2.fp <- 0
      } 
      if (rownames(tab)=="DE"){
        DESeq2.tp <- 1
        DESeq2.fp <- 1
      }
    } 
    
    
    DESeq2.zinb <- res$DESeq2.weighted
    lrt <- DESeq2.zinb[which(DESeq2.zinb$padj <= alpha),]
    pred <- rownames(DESeq2)
    pred[which(pred %in% rownames(lrt))] <- "DE"
    pred[which(pred != "DE")] <- "common"
    tab <- table(pred,truth)
    if (dim(tab)[1]==2) {
      DESeq2.zinb.tp <- tab[2,2]/nDE
      DESeq2.zinb.fp <- tab[2,1]/(length(pred)-nDE)
    }
    if (dim(tab)[1]==1) {
      if (rownames(tab)=="common"){
        DESeq2.zinb.tp <- 0
        DESeq2.zinb.fp <- 0
      } 
      if (rownames(tab)=="DE"){
        DESeq2.zinb.tp <- 1
        DESeq2.zinb.fp <- 1
      }
    } 
    
    DEsingle <- res$DEsingle
    DEsingle <- DEsingle[rownames(edgeR),]
    lrt <- DEsingle[DEsingle$pvalue.adj.FDR <= alpha, ]
    pred <- rownames(DEsingle)
    pred[which(pred %in% rownames(lrt))] <- "DE"
    pred[which(pred != "DE")] <- "common"
    tab <- table(pred,truth)
    if (dim(tab)[1]==2) {
      DEsingle.tp <- tab[2,2]/nDE
      DEsingle.fp <- tab[2,1]/(length(pred)-nDE)
    }
    if (dim(tab)[1]==1) {
      if (rownames(tab)=="common"){
        DEsingle.tp <- 0
        DEsingle.fp <- 0
      } 
      if (rownames(tab)=="DE"){
        DEsingle.tp <- 1
        DEsingle.fp <- 1
      }
    } 
    
    

    limma <- res$limma
    limma <- limma[rownames(edgeR),]
    lrt <- limma[which(limma$adj.P.Val <= alpha),]
    pred <- rownames(limma)
    pred[which(pred %in% rownames(lrt))] <- "DE"
    pred[which(pred != "DE")] <- "common"
    tab <- table(pred,truth)
    if (dim(tab)[1]==2) {
      limma.tp <- tab[2,2]/nDE
      limma.fp <- tab[2,1]/(length(pred)-nDE)
    }
    if (dim(tab)[1]==1) {
      if (rownames(tab)=="common"){
        limma.tp <- 0
        limma.fp <- 0
      } 
      if (rownames(tab)=="DE"){
        limma.tp <- 1
        limma.fp <- 1
      }
    } 
    
    
    limma.weighted <- res$limma.weighted
    limma.weighted <- limma.weighted[rownames(edgeR),]
    lrt <- limma.weighted[which(limma.weighted$adj.P.Val <= alpha),]
    pred <- rownames(limma.weighted)
    pred[which(pred %in% rownames(lrt))] <- "DE"
    pred[which(pred != "DE")] <- "common"
    tab <- table(pred,truth)
    if (dim(tab)[1]==2) {
      limma.weighted.tp <- tab[2,2]/nDE
      limma.weighted.fp <- tab[2,1]/(length(pred)-nDE)
    }
    if (dim(tab)[1]==1) {
      if (rownames(tab)=="common"){
        limma.weighted.tp <- 0
        limma.weighted.fp <- 0
      } 
      if (rownames(tab)=="DE"){
        limma.weighted.tp <- 1
        limma.weighted.fp <- 1
      }
    } 
    
    
    MAST <- res$MAST
    rownames(MAST) <- MAST$primerid
    lrt <- MAST[which(MAST$`Pr(>Chisq)` <= alpha),]
    pred <- MAST$primerid
    pred[which(pred %in% lrt$primerid)] <- "DE"
    pred[which(pred != "DE")] <- "common"
    tab <- table(pred,truth)
    if (dim(tab)[1]==2) {
      MAST.tp <- tab[2,2]/nDE
      MAST.fp <- tab[2,1]/(length(pred)-nDE)
    }
    if (dim(tab)[1]==1) {
      if (rownames(tab)=="common"){
        MAST.tp <- 0
        MAST.fp <- 0
      } 
      if (rownames(tab)=="DE"){
        MAST.tp <- 1
        MAST.fp <- 1
      }
    } 
    
    
    tpfp <- data.frame(edgeR.tp, edgeR.fp, 
                       DESeq2.tp, DESeq2.fp,
                       edgeR.zinb.tp, edgeR.zinb.fp,
                       DESeq2.zinb.tp, DESeq2.zinb.fp,
                       DEsingle.tp, DEsingle.fp, 
                       limma.tp, limma.fp, limma.weighted.tp, limma.weighted.fp,
                       MAST.tp, MAST.fp, 
                       alpha)
    TPFP <- rbind(TPFP,tpfp)
  }
return(TPFP)
}



library(zoo)
AUCdot <- function(TPFP,type){
  names <- c("edgeR","DESeq2","ZINB-Wave.edgeR","ZINB-Wave.DESeq2",
             "DEsingle",
             "limma_voom","ZINB-Wave.limma_voom","MAST")
  AUCdata <- NULL
  for (i in 1:8) {
    x=TPFP[,2*i]
    y=TPFP[,2*i-1]
    id <- order(x)
    AUC <- sum(diff(x[id])*rollmean(y[id],2))
    AUCdata <- c(AUCdata,AUC)
  }
  AUCdata <- data.frame(AUC=AUCdata)
  AUCdata$name <- names
  AUCdata$type <- type
  return(AUCdata)
}


setwd(paste0("E:\\DB\\Dropbox\\Qinfei\\Simulation of SC based on splatter\\DE\\DEdata\\Tung"))
res <- get(load(file="DE2000 results fold1.2 new.RData"))
TPFP1.2 <- tpfp.gen(res=res,nDE=200)
res <- get(load(file="DE2000 results fold1.4 new.RData"))
TPFP1.4 <- tpfp.gen(res=res,nDE=200)
res <- get(load(file="DE2000 results fold1.6 new.RData"))
TPFP1.6 <- tpfp.gen(res=res,nDE=200)
res <- get(load(file="DE2000 results fold1.8 new.RData"))
TPFP1.8 <- tpfp.gen(res=res,nDE=200)

AUCdata1 <- AUCdot(TPFP=TPFP1.2, type="1.2")
AUCdata2 <- AUCdot(TPFP=TPFP1.4, type="1.4")
AUCdata3 <- AUCdot(TPFP=TPFP1.6, type="1.6")
AUCdata4 <- AUCdot(TPFP=TPFP1.8, type="1.8")
AUCdata.FC <- rbind(AUCdata1,AUCdata2,AUCdata3,AUCdata4)
AUCdata.FC$type <- as.factor(AUCdata.FC$type)
AUCdata.FC$name <- factor(AUCdata.FC$name,level=c("edgeR","DESeq2","limma_voom","MAST","ZINB-Wave.edgeR","ZINB-Wave.DESeq2","ZINB-Wave.limma_voom"))


name=c("Drop0.1","Drop0.2","Drop0.3","Drop0.4")
res <- get(load(file=paste0("DE2000 results fold1.4 ",name[1],"new.RData")))
TPFP1 <- tpfp.gen(res=res,nDE=200)
res <- get(load(file=paste0("DE2000 results fold1.4 ",name[2],"new.RData")))
TPFP2 <- tpfp.gen(res=res,nDE=200)
res <- get(load(file=paste0("DE2000 results fold1.4 ",name[3],"new.RData")))
TPFP3 <- tpfp.gen(res=res,nDE=200)
res <- get(load(file=paste0("DE2000 results fold1.4 ",name[4],"new.RData")))
TPFP4 <- tpfp.gen(res=res,nDE=200)

AUCdata1 <- AUCdot(TPFP=TPFP1, type="0.1")
AUCdata2 <- AUCdot(TPFP=TPFP2, type="0.2")
AUCdata3 <- AUCdot(TPFP=TPFP3, type="0.3")
AUCdata4 <- AUCdot(TPFP=TPFP4, type="0.4")
AUCdata.DP <- rbind(AUCdata1,AUCdata2,AUCdata3,AUCdata4)
AUCdata.DP$type <- as.factor(AUCdata.DP$type)
AUCdata.DP$name <- factor(AUCdata.DP$name,level=c("edgeR","DESeq2","limma_voom","MAST","ZINB-Wave.edgeR","ZINB-Wave.DESeq2","ZINB-Wave.limma_voom"))


name=c("lib5k","lib10k","lib20k","lib40k")
res <- get(load(file=paste0("DE2000 results fold1.4 ",name[1],"new.RData")))
TPFP1 <- tpfp.gen(res=res,nDE=200)
res <- get(load(file=paste0("DE2000 results fold1.4 ",name[2],"new.RData")))
TPFP2 <- tpfp.gen(res=res,nDE=200)
res <- get(load(file=paste0("DE2000 results fold1.4 ",name[3],"new.RData")))
TPFP3 <- tpfp.gen(res=res,nDE=200)
res <- get(load(file=paste0("DE2000 results fold1.4 ",name[4],"new.RData")))
TPFP4 <- tpfp.gen(res=res,nDE=200)

AUCdata1 <- AUCdot(TPFP=TPFP1, type="5k")
AUCdata2 <- AUCdot(TPFP=TPFP2, type="10k")
AUCdata3 <- AUCdot(TPFP=TPFP3, type="20k")
AUCdata4 <- AUCdot(TPFP=TPFP4, type="40k")
AUCdata.lib <- rbind(AUCdata1,AUCdata2,AUCdata3,AUCdata4)
AUCdata.lib$type <- factor(AUCdata.lib$type,level=c("5k","10k","20k","40k"))
AUCdata.lib$name <- factor(AUCdata.lib$name,level=c("edgeR","DESeq2","limma_voom","MAST","ZINB-Wave.edgeR","ZINB-Wave.DESeq2","ZINB-Wave.limma_voom"))



name=c("BCV.df2","BCV.df5","BCV.df10","BCV.df20")
res <- get(load(file=paste0("DE2000 results fold1.4 ",name[1],"new.RData")))
TPFP1 <- tpfp.gen(res=res,nDE=200)
res <- get(load(file=paste0("DE2000 results fold1.4 ",name[2],"new.RData")))
TPFP2 <- tpfp.gen(res=res,nDE=200)
res <- get(load(file=paste0("DE2000 results fold1.4 ",name[3],"new.RData")))
TPFP3 <- tpfp.gen(res=res,nDE=200)
res <- get(load(file=paste0("DE2000 results fold1.4 ",name[4],"new.RData")))
TPFP4 <- tpfp.gen(res=res,nDE=200)

AUCdata1 <- AUCdot(TPFP=TPFP1, type="2")
AUCdata2 <- AUCdot(TPFP=TPFP2, type="5")
AUCdata3 <- AUCdot(TPFP=TPFP3, type="10")
AUCdata4 <- AUCdot(TPFP=TPFP4, type="20")
AUCdata.BCV.df <- rbind(AUCdata1,AUCdata2,AUCdata3,AUCdata4)
AUCdata.BCV.df$type <- factor(AUCdata.BCV.df$type,level=c("2","5","10","20"))
AUCdata.BCV.df$name <- factor(AUCdata.BCV.df$name,level=c("edgeR","DESeq2","limma_voom","MAST","ZINB-Wave.edgeR","ZINB-Wave.DESeq2","ZINB-Wave.limma_voom"))



name=c("bcv.shrink0.5","bcv.shrink1.0","bcv.shrink1.5","bcv.shrink2.0")
res <- get(load(file=paste0("DE2000 results fold1.4 ",name[1],"new.RData")))
TPFP1 <- tpfp.gen(res=res,nDE=200)
res <- get(load(file=paste0("DE2000 results fold1.4 ",name[2],"new.RData")))
TPFP2 <- tpfp.gen(res=res,nDE=200)
res <- get(load(file=paste0("DE2000 results fold1.4 ",name[3],"new.RData")))
TPFP3 <- tpfp.gen(res=res,nDE=200)
res <- get(load(file=paste0("DE2000 results fold1.4 ",name[4],"new.RData")))
TPFP4 <- tpfp.gen(res=res,nDE=200)

AUCdata1 <- AUCdot(TPFP=TPFP1, type="0.5")
AUCdata2 <- AUCdot(TPFP=TPFP2, type="1.0")
AUCdata3 <- AUCdot(TPFP=TPFP3, type="1.5")
AUCdata4 <- AUCdot(TPFP=TPFP4, type="2.0")
AUCdata.BCV.shrink <- rbind(AUCdata1,AUCdata2,AUCdata3,AUCdata4)
AUCdata.BCV.shrink$type <- as.factor(AUCdata.BCV.shrink$type)
AUCdata.BCV.shrink$name <- factor(AUCdata.BCV.shrink$name,level=c("edgeR","DESeq2","limma_voom","MAST","ZINB-Wave.edgeR","ZINB-Wave.DESeq2","ZINB-Wave.limma_voom"))















library(ggplot2)
setwd(paste0("E:\\DB\\Dropbox\\Qinfei\\Simulation of SC based on splatter\\Submition\\Bioinformatics\\Simulation data"))

data=AUCdata.FC; scenario="Fold change"
data <- data[which(data$name != "DEsingle"),]
FC <- qplot(type, AUC, data = data, color=name, group = name)+    
  geom_point(size=2,alpha=0.8)+ geom_line() +  
  ggtitle(scenario)+
  xlab(scenario)+  ylim(0.7, 1) +
  guides(color=guide_legend("DE methods"))+
  theme(plot.title=element_text(hjust=0.5,size=15,face="bold"),axis.text=element_text(size=12),
        axis.title.x = element_text(size=13),axis.title.y =element_text(size=13),
        strip.text.x=element_text(size=12),legend.text=element_text(size=12),
        legend.title = element_text(size=12), panel.background = element_blank(),
        panel.grid.major = element_line(colour = "grey93"),
        axis.line = element_line(colour = "black"),
        plot.margin = unit(c(0,0,0,0,0), "lines"))


data=AUCdata.DP; scenario="Drop out"
data <- data[which(data$name != "DEsingle"),]
DP <- qplot(type, AUC, data = data, color=name, group = name)+    
  geom_point(size=2,alpha=0.8)+ geom_line() +  
  ggtitle(scenario)+
  xlab(scenario)+  ylim(0.7, 1) +
  guides(color=guide_legend("DE methods"))+
  theme(plot.title=element_text(hjust=0.5,size=15,face="bold"),axis.text=element_text(size=12),
        axis.title.x = element_text(size=13),axis.title.y =element_text(size=13),
        strip.text.x=element_text(size=12),legend.text=element_text(size=12),
        legend.title = element_text(size=12), panel.background = element_blank(),
        panel.grid.major = element_line(colour = "grey93"),
        axis.line = element_line(colour = "black"),
        plot.margin = unit(c(0,0,0,0), "lines"))


data=AUCdata.lib; scenario="Library size"
data <- data[which(data$name != "DEsingle"),]
lib <- qplot(type, AUC, data = data, color=name, group = name)+    
  geom_point(size=2,alpha=0.8)+ geom_line() +  
  ggtitle(scenario)+
  xlab(scenario)+  ylim(0.7, 1) +
  guides(color=guide_legend("DE methods"))+
  theme(plot.title=element_text(hjust=0.5,size=15,face="bold"),axis.text=element_text(size=12),
        axis.title.x = element_text(size=13),axis.title.y =element_text(size=13),
        strip.text.x=element_text(size=12),legend.text=element_text(size=12),
        legend.title = element_text(size=12), panel.background = element_blank(),
        panel.grid.major = element_line(colour = "grey93"),
        axis.line = element_line(colour = "black"),
        plot.margin = unit(c(0,0,0,0), "lines"))


data=AUCdata.BCV.df; scenario="BCV.df"
data <- data[which(data$name != "DEsingle"),]
BCV.df <- qplot(type, AUC, data = data, color=name, group = name)+    
  geom_point(size=2,alpha=0.8)+ geom_line() +  
  ggtitle(scenario)+
  xlab(scenario)+  ylim(0.7, 1) +
  guides(color=guide_legend("DE methods"))+
  theme(plot.title=element_text(hjust=0.5,size=15,face="bold"),axis.text=element_text(size=12),
        axis.title.x = element_text(size=13),axis.title.y =element_text(size=13),
        strip.text.x=element_text(size=12),legend.text=element_text(size=12),
        legend.title = element_text(size=12), panel.background = element_blank(),
        panel.grid.major = element_line(colour = "grey93"),
        axis.line = element_line(colour = "black"),
        plot.margin = unit(c(0,0,0,0), "lines"))


data=AUCdata.BCV.shrink; scenario="BCV.shrink"
data <- data[which(data$name != "DEsingle"),]
BCV.shrink <- qplot(type, AUC, data = data, color=name, group = name)+    
  geom_point(size=2,alpha=0.8)+ geom_line() +  
  ggtitle(scenario)+
  xlab(scenario)+  ylim(0.7, 1) +
  guides(color=guide_legend("DE methods"))+
  theme(plot.title=element_text(hjust=0.5,size=15,face="bold"),axis.text=element_text(size=12),
        axis.title.x = element_text(size=13),axis.title.y =element_text(size=13),
        strip.text.x=element_text(size=12),legend.text=element_text(size=12),
        legend.title = element_text(size=12), panel.background = element_blank(),
        panel.grid.major = element_line(colour = "grey93"),
        axis.line = element_line(colour = "black"),
        plot.margin = unit(c(0,0,0,0), "lines"))




library(ggpubr)
pdf("DE.evaluation.AUCdata.BCV.shrink.pdf",height=6,width=10)
BCV.shrink
dev.off()


library(ggpubr)
pdf("DE.evaluation.updated.pdf",height=6,width=10)
ggarrange(FC, DP, lib, BCV.df, BCV.shrink ,nrow=2, ncol=3, common.legend = T, legend = "bottom",
          font.label = list(size = 13))
dev.off()












### Fixed alpha
setwd(paste0("E:\\DB\\Dropbox\\Qinfei\\Simulation of SC based on splatter\\DE\\DEdata\\Tung"))
res <- get(load(file="DE2000 results fold1.2 new.RData"))
TPFP1 <- tpfp.gen(res=res,nDE=200)
res <- get(load(file="DE2000 results fold1.4 new.RData"))
TPFP2 <- tpfp.gen(res=res,nDE=200)
res <- get(load(file="DE2000 results fold1.6 new.RData"))
TPFP3 <- tpfp.gen(res=res,nDE=200)
res <- get(load(file="DE2000 results fold1.8 new.RData"))
TPFP4 <- tpfp.gen(res=res,nDE=200)

TPFP <- rbind(TPFP1,TPFP2,TPFP3,TPFP4)
TPFP.FC <- data.frame(TP=c(TPFP[,1],TPFP[,3],TPFP[,5],TPFP[,7],TPFP[,9],TPFP[,11],TPFP[,13],TPFP[,15]),
                      FP=c(TPFP[,2],TPFP[,4],TPFP[,6],TPFP[,8],TPFP[,10],TPFP[,12],TPFP[,14],TPFP[,16]),
                      alpha=0.05,
                      type="FC",
                      num=rep(c(1.2,1.4,1.6,1.8),8),
                      name=rep(c("edgeR","DESeq2","ZINB-Wave.edgeR","ZINB-Wave.DESeq2","DEsingle","limma_voom","ZINB-Wave.limma_voom","MAST"),each=4))
TPFP.FC$num <- as.factor(TPFP.FC$num)
TPFP.FC$name <- factor(TPFP.FC$name,level=c("edgeR","DESeq2","limma_voom","MAST","ZINB-Wave.edgeR","ZINB-Wave.DESeq2","ZINB-Wave.limma_voom"))






name=c("Drop0.1","Drop0.2","Drop0.3","Drop0.4")
res <- get(load(file=paste0("DE2000 results fold1.4 ",name[1],"new.RData")))
TPFP1 <- tpfp.gen(res=res,nDE=200)
res <- get(load(file=paste0("DE2000 results fold1.4 ",name[2],"new.RData")))
TPFP2 <- tpfp.gen(res=res,nDE=200)
res <- get(load(file=paste0("DE2000 results fold1.4 ",name[3],"new.RData")))
TPFP3 <- tpfp.gen(res=res,nDE=200)
res <- get(load(file=paste0("DE2000 results fold1.4 ",name[4],"new.RData")))
TPFP4 <- tpfp.gen(res=res,nDE=200)

TPFP <- rbind(TPFP1,TPFP2,TPFP3,TPFP4)
TPFP.dropout <- data.frame(TP=c(TPFP[,1],TPFP[,3],TPFP[,5],TPFP[,7],TPFP[,9],TPFP[,11],TPFP[,13],TPFP[,15]),
                      FP=c(TPFP[,2],TPFP[,4],TPFP[,6],TPFP[,8],TPFP[,10],TPFP[,12],TPFP[,14],TPFP[,16]),
                      alpha=0.05,
                      type="Drop out",
                      num=rep(c(0.1,0.2,0.3,0.4),8),
                      name=rep(c("edgeR","DESeq2","ZINB-Wave.edgeR","ZINB-Wave.DESeq2","DEsingle","limma_voom","ZINB-Wave.limma_voom","MAST"),each=4))
TPFP.dropout$num <- as.factor(TPFP.dropout$num)
TPFP.dropout$name <- factor(TPFP.dropout$name,level=c("edgeR","DESeq2","limma_voom","MAST","ZINB-Wave.edgeR","ZINB-Wave.DESeq2","ZINB-Wave.limma_voom"))






name=c("lib5k","lib10k","lib20k","lib40k")
res <- get(load(file=paste0("DE2000 results fold1.4 ",name[1],"new.RData")))
TPFP1 <- tpfp.gen(res=res,nDE=200)
res <- get(load(file=paste0("DE2000 results fold1.4 ",name[2],"new.RData")))
TPFP2 <- tpfp.gen(res=res,nDE=200)
res <- get(load(file=paste0("DE2000 results fold1.4 ",name[3],"new.RData")))
TPFP3 <- tpfp.gen(res=res,nDE=200)
res <- get(load(file=paste0("DE2000 results fold1.4 ",name[4],"new.RData")))
TPFP4 <- tpfp.gen(res=res,nDE=200)

TPFP <- rbind(TPFP1,TPFP2,TPFP3,TPFP4)
TPFP.Libsize <- data.frame(TP=c(TPFP[,1],TPFP[,3],TPFP[,5],TPFP[,7],TPFP[,9],TPFP[,11],TPFP[,13],TPFP[,15]),
                      FP=c(TPFP[,2],TPFP[,4],TPFP[,6],TPFP[,8],TPFP[,10],TPFP[,12],TPFP[,14],TPFP[,16]),
                      alpha=0.05,
                      type="Library size",
                      num=rep(c("5k","10k","20k","40k"),8),
                      name=rep(c("edgeR","DESeq2","ZINB-Wave.edgeR","ZINB-Wave.DESeq2","DEsingle","limma_voom","ZINB-Wave.limma_voom","MAST"),each=4))
TPFP.Libsize$num <- factor(TPFP.Libsize$num,level=c("5k","10k","20k","40k"))
TPFP.Libsize$name <- factor(TPFP.Libsize$name,level=c("edgeR","DESeq2","limma_voom","MAST","ZINB-Wave.edgeR","ZINB-Wave.DESeq2","ZINB-Wave.limma_voom"))





name=c("BCV.df2","BCV.df5","BCV.df10","BCV.df20")
res <- get(load(file=paste0("DE2000 results fold1.4 ",name[1],"new.RData")))
TPFP1 <- tpfp.gen(res=res,nDE=200)
res <- get(load(file=paste0("DE2000 results fold1.4 ",name[2],"new.RData")))
TPFP2 <- tpfp.gen(res=res,nDE=200)
res <- get(load(file=paste0("DE2000 results fold1.4 ",name[3],"new.RData")))
TPFP3 <- tpfp.gen(res=res,nDE=200)
res <- get(load(file=paste0("DE2000 results fold1.4 ",name[4],"new.RData")))
TPFP4 <- tpfp.gen(res=res,nDE=200)

TPFP <- rbind(TPFP1,TPFP2,TPFP3,TPFP4)
TPFP.BCVdf <- data.frame(TP=c(TPFP[,1],TPFP[,3],TPFP[,5],TPFP[,7],TPFP[,9],TPFP[,11],TPFP[,13],TPFP[,15]),
                      FP=c(TPFP[,2],TPFP[,4],TPFP[,6],TPFP[,8],TPFP[,10],TPFP[,12],TPFP[,14],TPFP[,16]),
                      alpha=0.05,
                      type="BCV df",
                      num=rep(c(2,5,10,20),8),
                      name=rep(c("edgeR","DESeq2","ZINB-Wave.edgeR","ZINB-Wave.DESeq2","DEsingle","limma_voom","ZINB-Wave.limma_voom","MAST"),each=4))
TPFP.BCVdf$num <- as.factor(TPFP.BCVdf$num)
TPFP.BCVdf$name <- factor(TPFP.BCVdf$name,level=c("edgeR","DESeq2","limma_voom","MAST","ZINB-Wave.edgeR","ZINB-Wave.DESeq2","ZINB-Wave.limma_voom"))






name=c("bcv.shrink0.5","bcv.shrink1.0","bcv.shrink1.5","bcv.shrink2.0")
res <- get(load(file=paste0("DE2000 results fold1.4 ",name[1],"new.RData")))
TPFP1 <- tpfp.gen(res=res,nDE=200)
res <- get(load(file=paste0("DE2000 results fold1.4 ",name[2],"new.RData")))
TPFP2 <- tpfp.gen(res=res,nDE=200)
res <- get(load(file=paste0("DE2000 results fold1.4 ",name[3],"new.RData")))
TPFP3 <- tpfp.gen(res=res,nDE=200)
res <- get(load(file=paste0("DE2000 results fold1.4 ",name[4],"new.RData")))
TPFP4 <- tpfp.gen(res=res,nDE=200)


TPFP <- rbind(TPFP1,TPFP2,TPFP3,TPFP4)
TPFP.BCVshrink <- data.frame(TP=c(TPFP[,1],TPFP[,3],TPFP[,5],TPFP[,7],TPFP[,9],TPFP[,11],TPFP[,13],TPFP[,15]),
                      FP=c(TPFP[,2],TPFP[,4],TPFP[,6],TPFP[,8],TPFP[,10],TPFP[,12],TPFP[,14],TPFP[,16]),
                      alpha=0.05,
                      type="BCV shrink",
                      num=rep(c(0.5, 1.0, 1.5, 2.0),8),
                      name=rep(c("edgeR","DESeq2","ZINB-Wave.edgeR","ZINB-Wave.DESeq2","DEsingle","limma_voom","ZINB-Wave.limma_voom","MAST"),each=4))
TPFP.BCVshrink$num <- as.factor(TPFP.BCVshrink$num)
TPFP.BCVshrink$name <- factor(TPFP.BCVshrink$name,level=c("edgeR","DESeq2","limma_voom","MAST","ZINB-Wave.edgeR","ZINB-Wave.DESeq2","ZINB-Wave.limma_voom"))


library(ggplot2)
setwd(paste0("E:\\DB\\Dropbox\\Qinfei\\Simulation of SC based on splatter\\Submition\\Bioinformatics\\Simulation data"))

data=TPFP.FC; scenario="Fold change"
data <- data[which(data$name != "DEsingle"),]
FC <- qplot(num, TP, data = data, color=name, group = name)+    
  geom_point(size=2,alpha=0.8)+ geom_line() +  
  ggtitle(scenario)+
  xlab(scenario)+ ylim(0, 1) +
  guides(color=guide_legend("DE methods"))+
  theme(plot.title=element_text(hjust=0.5,size=15,face="bold"),axis.text=element_text(size=12),
        axis.title.x = element_text(size=13),axis.title.y =element_text(size=13),
        strip.text.x=element_text(size=12),legend.text=element_text(size=12),
        legend.title = element_text(size=12), panel.background = element_blank(),
        panel.grid.major = element_line(colour = "grey93"),
        axis.line = element_line(colour = "black"),
        plot.margin = unit(c(0.1,0.1,0.1,0.1,0.1), "lines"))


data=TPFP.dropout; scenario="Drop out"
data <- data[which(data$name != "DEsingle"),]
DP <- qplot(num, TP, data = data, color=name, group = name)+  
  geom_point(size=2,alpha=0.8)+ geom_line() +  
  ggtitle(scenario)+
  xlab(scenario) +ylim(0, 1) +
  guides(color=guide_legend("DE methods"))+
  theme(plot.title=element_text(hjust=0.5,size=15,face="bold"),axis.text=element_text(size=12),
        axis.title.x = element_text(size=13),axis.title.y =element_text(size=13),
        strip.text.x=element_text(size=12),legend.text=element_text(size=12),
        legend.title = element_text(size=12), panel.background = element_blank(),
        panel.grid.major = element_line(colour = "grey93"),
        axis.line = element_line(colour = "black"),
        plot.margin = unit(c(0.1,0.1,0.1,0.1,0.1), "lines"))



data=TPFP.Libsize; scenario="Library size"
data <- data[which(data$name != "DEsingle"),]
lib <- qplot(num, TP, data = data, color=name, group = name)+  
  geom_point(size=2,alpha=0.8)+ geom_line() +  
  ggtitle(scenario)+
  xlab(scenario)+  ylim(0, 1) +
  guides(color=guide_legend("DE methods"))+
  theme(plot.title=element_text(hjust=0.5,size=15,face="bold"),axis.text=element_text(size=12),
        axis.title.x = element_text(size=13),axis.title.y =element_text(size=13),
        strip.text.x=element_text(size=12),legend.text=element_text(size=12),
        legend.title = element_text(size=12), panel.background = element_blank(),
        panel.grid.major = element_line(colour = "grey93"),
        axis.line = element_line(colour = "black"),
        plot.margin = unit(c(0.1,0.1,0.1,0.1,0.1), "lines"))


data=TPFP.BCVdf; scenario="BCV.df"
data <- data[which(data$name != "DEsingle"),]
BCV.df <- qplot(num, TP, data = data, color=name, group = name)+  
  geom_point(size=2,alpha=0.8)+ geom_line() +  
  ggtitle(scenario)+
  xlab(scenario)+  ylim(0, 1) +
  guides(color=guide_legend("DE methods"))+
  theme(plot.title=element_text(hjust=0.5,size=15,face="bold"),axis.text=element_text(size=12),
        axis.title.x = element_text(size=13),axis.title.y =element_text(size=13),
        strip.text.x=element_text(size=12),legend.text=element_text(size=12),
        legend.title = element_text(size=12), panel.background = element_blank(),
        panel.grid.major = element_line(colour = "grey93"),
        axis.line = element_line(colour = "black"),
        plot.margin = unit(c(0.1,0.1,0.1,0.1,0.1), "lines"))


data=TPFP.BCVshrink; scenario="BCV.shrink"
data <- data[which(data$name != "DEsingle"),]
BCV.shrink <- qplot(num, TP, data = data, color=name, group = name)+  
  geom_point(size=2,alpha=0.8)+ geom_line() +  
  ggtitle(scenario)+
  xlab(scenario)+ ylim(0, 1) +
  guides(color=guide_legend("DE methods"))+
  theme(plot.title=element_text(hjust=0.5,size=15,face="bold"),axis.text=element_text(size=12),
        axis.title.x = element_text(size=13),axis.title.y =element_text(size=13),
        strip.text.x=element_text(size=12),legend.text=element_text(size=12),
        legend.title = element_text(size=12), panel.background = element_blank(),
        panel.grid.major = element_line(colour = "grey93"),
        axis.line = element_line(colour = "black"),
        plot.margin = unit(c(0.1,0.1,0.1,0.1,0.1), "lines"))


library(ggpubr)
pdf("DE.FC.fixed.alpha.pdf",height=6,width=8)
   FC
dev.off()


library(ggpubr)
pdf("DE.evaluation.fixed.alpha.pdf",height=6,width=9)
ggarrange(FC, DP, lib, BCV.df, BCV.shrink ,nrow=2, ncol=3, common.legend = T, legend = "bottom",
          font.label = list(size = 13))
dev.off()


