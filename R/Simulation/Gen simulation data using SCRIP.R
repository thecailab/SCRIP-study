setwd("/home/fqin/SCRIP/Klein")

Klein <- read.csv("Klein.csv",header=T,row.names = 1)
bursting <- read.csv("Klein_bursting_output.csv")
colnames(bursting) <- c("Gene", "kon","koff","s")
expre_data <- as.matrix(Klein)

library("Biobase")
library(splatter)
library(SymSim)

params <- splatEstimate(expre_data)

library(SCRIP)

start_time <- Sys.time()  
single_cellsample <- SCRIPsimu(data=expre_data, params=params,
                               mode="GP-commonBCV")
single_cellsample <- single_cellsample@assays@data$counts
end_time <- Sys.time()
print("the time for running GP-commonBCV is")
print(end_time-start_time)
save(single_cellsample,file="GP-commonBCV.RData")


start_time <- Sys.time()  
single_cellsample <- SCRIPsimu(data=expre_data, params=params,
                               mode="GP-trendedBCV")
single_cellsample <- single_cellsample@assays@data$counts
end_time <- Sys.time()
print("the time for running GP-trendedBCV is")
print(end_time-start_time)
save(single_cellsample,file="GP-trendedBCV.RData")

library(SymSim)
start_time <- Sys.time()  
single_cellsample <- SCRIPsimu(data=expre_data, params=params, kinetics=bursting, 
                               mode="BGP-commonBCV")
single_cellsample <- single_cellsample@assays@data$counts
end_time <- Sys.time()
print("the time for running BGP-commonBCV is")
print(end_time-start_time)
save(single_cellsample,file="BGP-commonBCV.RData")


start_time <- Sys.time()  
single_cellsample <- SCRIPsimu(data=expre_data, params=params, kinetics=bursting, 
                               mode="BGP-trendedBCV")
single_cellsample <- single_cellsample@assays@data$counts
end_time <- Sys.time()
print("the time for running BGP-trendedBCV is")
print(end_time-start_time)
save(single_cellsample,file="BGP-trendedBCV.RData")

start_time <- Sys.time()  
single_cellsample <- SCRIPsimu(data=expre_data, params=params, kinetics=bursting, 
                               mode="BP")
single_cellsample <- single_cellsample@assays@data$counts
end_time <- Sys.time()
print("the time for running BP is")
print(end_time-start_time)
save(single_cellsample,file="BP.RData")

