#--------------------------------------------------------------------------------------------------
#                                  PREPROCESSING ANALYSIS
#--------------------------------------------------------------------------------------------------



###################################################################################################
#                                   SCnorm
###################################################################################################

#source("https://bioconductor.org/biocLite.R")
#biocLite("SCnorm")
library(SCnorm)

#### SCnorm

dge_scnorm = SCnorm(Data = raw.data, Conditions = rep(c(1), each = ncol(norm.data)), PrintProgressPlots = T,
                    reportSF = T, FilterCellNum = 10, FilterExpression = 0,
                    Thresh = 0.1, K = NULL, NCores = NULL, ditherCounts = FALSE,
                    PropToUse = 0.25, Tau = 0.5, withinSample = NULL, useSpikes = FALSE,
                    useZerosToScale = FALSE)
# Setting up parallel computation using 15 cores
# Gene filter is applied within each condition.
# 248 genes in condition 1 will not be included in the normalization due to 
# the specified filter criteria.
# A list of these genes can be accessed in output, 
# see vignette for example.
# Finding K for Condition 1
# Trying K = 1
# Trying K = 2
# Trying K = 3
# Trying K = 4
# Done!
scnorm.data = dge_scnorm@metadata$NormalizedData
write.table(scnorm.data,file = "dge_scnorm_with_raw.txt",sep = "\t",row.names = T,col.names = T)
#didn't understand ???????????????????????????????????????????????????????????????????
scnorm_filtered = dge_scnorm@metadata$GenesFilteredOut
write.table(scnorm_filtered,file = "dge_scnorm_filtered_with_raw.txt",sep = "\t",row.names = T,col.names = T)
# filtered out genes from driver or not
intersect(scnorm_filtered, insitu.driver)

#### DistMap with scnorm.data

#scnorm.data = as.matrix(scnorm.data)
#input data into DistMap class
dm_scnorm = new("DistMap", raw.data=raw.data,data=scnorm.data, insitu.matrix=insitu.matrix, geometry=geometry)
#get binasied dge input dataset
dm_scnorm <- binarizeSingleCellData(dm_scnorm, seq(0.15, 0.5, 0.01))
write.table(dm_scnorm@binarized.data,file = "dge_scnorm_binarized_with_raw.csv",sep = ",",row.names = T,col.names = T)
#MCC counts of mapping of cells to each bins
dm_scnorm <- mapCells(dm_scnorm)
write.table(dm_scnorm@mcc.scores,file = "dge_scnorm_mccscores_with_raw.csv",sep = ",",row.names = T,col.names = T)
#didn't understand ???????????????????????????????????????????????????????????????????
pha_scnorm = computeVISH(dm_scnorm, 'sna', threshold=0.75)
computeGeneGradient(dm_scnorm, 'sna') 
# retrieve best quantile score
bq_scnorm <- BestQuantile(dm_scnorm, seq(0.15, 0.5, 0.01))

###################################################################################################
#                                   Comparison 
###################################################################################################

#### by means
mean(dm@mcc.scores)
mean(dm_scnorm@mcc.scores)

#### ttest
t.test(as.vector(dm@mcc.scores),as.vector(dm_scnorm@mcc.scores))

#### corr between mccs
cor(as.vector(dm@mcc.scores),as.vector(dm_scnorm@mcc.scores))

#### rmse
#install.packages("Metrics")
library(Metrics)
rmse(as.vector(dm@mcc.scores),as.vector(dm_scnorm@mcc.scores))

#### by euclidean
dm_max <- c()
for(i in 1:ncol(dm@mcc.scores)){dm_max <- c(dm_max,which.max(dm@mcc.scores[,i]))}
dm_scnorm_max <- c()
for(i in 1:ncol(dm_scnorm@mcc.scores)){dm_scnorm_max <- c(dm_scnorm_max,which.max(dm_scnorm@mcc.scores[,i]))}
geo_dm = geometry[dm_max,]
geo_dm_scnorm = geometry[dm_scnorm_max,]
eucl_dm_x = c()
for(i in 1:nrow(geo_dm)){
  p = (geo_dm[i,1] - geo_dm_scnorm[i,1])**2
  eucl_dm_x <- c(eucl_dm_x,p)
}
eucl_dm_y = c()
for(i in 1:nrow(geo_dm)){
  p = (geo_dm[i,2] - geo_dm_scnorm[i,2])**2
  eucl_dm_y <- c(eucl_dm_y,p)
}
eucl_dm_z = c()
for(i in 1:nrow(geo_dm)){
  p = (geo_dm[i,3] - geo_dm_scnorm[i,3])**2
  eucl_dm_z <- c(eucl_dm_z,p)
}
eucl_dm = c()
for(i in 1:nrow(geo_dm)){
  sum = sqrt(eucl_dm_x[i] + eucl_dm_y[i] + eucl_dm_z[i])
  eucl_dm <- c(eucl_dm, sum)
}



###################################################################################################
#                                   TO DO LIST
###################################################################################################

#https://www.synapse.org/#!Synapse:syn15665609/wiki/583230
#https://www.biostars.org/p/335555/
# For normalisation:
#   scran  ------------------------notuseful
#   SCnorm ------------------------Done-----------------notuseful
#   DESeq2 (uses RLE)--------------notuseful
#   Seurat - CCA algorithm---------only lognormalisation
#   baynorm------------------------not installing
library(devtools)
devtools::install_github("WT215/bayNorm")
install.packages("BiocManager")
library(BiocManager)
BiocManager::install("bayNorm")
#   SAVER-------------------------- 
#python DCA
library(SAVER)


###################################################################################################
#                                   SAVER
###################################################################################################

#install.packages("SAVER")
library(SAVER)

dge_saver = saver(raw.data, ncores = 12)
# Done!
#   Finish time: 2018-10-29 22:01:19
# Total time: 12.98048 mins

saver.data = dge_saver$estimate
write.table(saver.data,file = "dge_saver_with_raw.txt",sep = "\t",row.names = T,col.names = T)
saver_se = dge_saver$se
write.table(saver_se,file = "dge_saver_se_with_raw.txt",sep = "\t",row.names = T,col.names = T)
save(dge_saver,file = "dge_saver_complete_with_raw.RData")

#### DistMap with saver.data

saver.data = as.matrix(saver.data)

#input data into DistMap class
dm_saver = new("DistMap", raw.data=raw.data, data=saver.data, insitu.matrix=insitu.matrix, geometry=geometry)
#get binasied dge input dataset
dm_saver <- binarizeSingleCellData(dm_saver, seq(0.15, 0.5, 0.01))
write.table(dm_saver@binarized.data,file = "dge_saver_binarized_with_raw.csv",sep = ",",row.names = T,col.names = T)
#MCC counts of mapping of cells to each bins
dm_saver <- mapCells(dm_saver)
write.table(dm_saver@mcc.scores,file = "dge_saver_mccscores_with_raw.csv",sep = ",",row.names = T,col.names = T)
# retrieve best quantile score
bq_saver <- BestQuantile(dm_saver, seq(0.15, 0.5, 0.01))



###################################################################################################
#                                   scImpute
###################################################################################################

library(devtools)
#install_github("Vivianstats/scImpute")
library(scImpute)

dge_scimpute = scimpute(count_path = "./dge_raw.txt", infile = "txt", outfile = "txt",type = "count", out_dir = "./", labeled = FALSE, Kcluster = 9, ncores = 12, drop_thre = 0.7)


#### DistMap with scimpute.data

scimpute.data = read.table("scimpute_count.txt", header= TRUE, sep=" ", row.names = 1)
scimpute.data = as.matrix(scimpute.data)
#input data into DistMap class
dm_scimpute = new("DistMap", raw.data=raw.data, data=scimpute.data, insitu.matrix=insitu.matrix, geometry=geometry)
#get binasied dge input dataset
dm_scimpute <- binarizeSingleCellData(dm_scimpute, seq(0.15, 0.5, 0.01))
write.table(dm_scimpute@binarized.data,file = "dge_scimpute_binarized_with_raw.csv",sep = ",",row.names = T,col.names = T)
#MCC counts of mapping of cells to each bins
dm_scimpute <- mapCells(dm_scimpute)
write.table(dm_scimpute@mcc.scores,file = "dge_scimpute_mccscores_with_raw.csv",sep = ",",row.names = T,col.names = T)
# retrieve best quantile score
bq_scimpute <- BestQuantile(dm_scimpute, seq(0.15, 0.5, 0.01))





###################################################################################################
#                                   MAGiC
###################################################################################################

# #install.packages("Rmagic")
# library(Rmagic)
# if (!require(phateR)) install.packages("phateR")
# if (!require(viridis)) install.packages("viridis")
# if (!require(ggplot2)) install.packages("ggplot2")
# if (!require(readr)) install.packages("readr")
# ## Loading required package: Matrix
# library(ggplot2)
# library(readr)
# library(viridisLite)
# library(viridis)
# library(phateR)
# 
# bmmsc <- t(raw.data)
# keep_cols <- colSums(bmmsc > 0) > 10
# bmmsc <- bmmsc[,keep_cols]
# 
# keep_rows <- rowSums(bmmsc) > 1000
# bmmsc <- bmmsc[keep_rows,]
# bmmsc <- library.size.normalize(bmmsc)
# bmmsc <- sqrt(bmmsc)
# bmmsc_MAGIC <- magic(bmmsc, genes=colnames(bmmsc))
#                       # Error in py_module_import(module, convert = convert) : 
#                       #   ImportError: No module named Tkinter
#                       # Error: Python module magic was not found.
#                       # 
#                       # Detected Python configuration:
#                       #   
#                       #   python:         /usr/bin/python
#                       # libpython:      /usr/lib64/python2.7/config/libpython2.7.so
#                       # pythonhome:     /usr:/usr
#                       # version:        2.7.5 (default, Jul 13 2018, 13:06:57)  [GCC 4.8.5 20150623 (Red Hat 4.8.5-28)]
#                       # numpy:          /usr/lib64/python2.7/site-packages/numpy
#                       # numpy_version:  1.14.5
#                       # magic:          /home/shruti/.local/lib/python2.7/site-packages/magic
# 

#### DistMap with scimpute.data

dge_magic = read.table("dge_raw_magic_imputed.txt", header= TRUE, sep="\t", row.names = 1)
magic.data = as.matrix(dge_magic)
#input data into DistMap class
dm_magic = new("DistMap", raw.data=raw.data, data=magic.data, insitu.matrix=insitu.matrix, geometry=geometry)
#get binasied dge input dataset
dm_magic <- binarizeSingleCellData(dm_magic, seq(0.15, 0.5, 0.01))
write.table(dm_magic@binarized.data,file = "dge_magic_binarized_with_raw.csv",sep = ",",row.names = T,col.names = T)
#MCC counts of mapping of cells to each bins
dm_magic <- mapCells(dm_magic)
write.table(dm_magic@mcc.scores,file = "dge_magic_mccscores_with_raw.csv",sep = ",",row.names = T,col.names = T)
# retrieve best quantile score
bq_magic <- BestQuantile(dm_magic, seq(0.15, 0.5, 0.01))


###################################################################################################
#                                   TO DO LIST
###################################################################################################

#https://www.synapse.org/#!Synapse:syn15665609/wiki/583230

#https://www.biostars.org/p/335555/
# For normalisation/imputation:
#   scran  ------------------------notuseful
#   SCnorm ------------------------Done----------------- USEFUL !!
#   DESeq2 (uses RLE)--------------notuseful
#   Seurat - CCA algorithm---------only lognormalisation
#   SAVER--------------------------Done-----------------notuseful
#   scImpute--------------------------Done-----------------notuseful
#   python DCA---------------------notuseful
#   baynorm------------------------not installing
# library(devtools)
# devtools::install_github("WT215/bayNorm")
# install.packages("BiocManager")
# library(BiocManager)
# BiocManager::install("bayNorm")


###################################################################################################
#                                   *** END ***
###################################################################################################



###################################################################################################
#                                   SCnorm
###################################################################################################

#source("https://bioconductor.org/biocLite.R")
#biocLite("SCnorm")
library(SCnorm)

#### SCnorm

dge_scnorm = SCnorm(Data = raw.data, Conditions = rep(c(1), each = ncol(norm.data)), PrintProgressPlots = T,
                    reportSF = T, FilterCellNum = 10, FilterExpression = 0,
                    Thresh = 0.1, K = NULL, NCores = NULL, ditherCounts = FALSE,
                    PropToUse = 0.25, Tau = 0.5, withinSample = NULL, useSpikes = FALSE,
                    useZerosToScale = FALSE)
# Setting up parallel computation using 15 cores
# Gene filter is applied within each condition.
# 248 genes in condition 1 will not be included in the normalization due to 
# the specified filter criteria.
# A list of these genes can be accessed in output, 
# see vignette for example.
# Finding K for Condition 1
# Trying K = 1
# Trying K = 2
# Trying K = 3
# Trying K = 4
# Done!
scnorm.data = dge_scnorm@metadata$NormalizedData
write.table(scnorm.data,file = "dge_scnorm_with_raw.txt",sep = "\t",row.names = T,col.names = T)
#didn't understand ???????????????????????????????????????????????????????????????????
scnorm_filtered = dge_scnorm@metadata$GenesFilteredOut
write.table(scnorm_filtered,file = "dge_scnorm_filtered_with_raw.txt",sep = "\t",row.names = T,col.names = T)
# filtered out genes from driver or not
intersect(scnorm_filtered, insitu.driver)

#### DistMap with scnorm.data

#scnorm.data = as.matrix(scnorm.data)
#input data into DistMap class
dm_scnorm = new("DistMap", raw.data=raw.data,data=scnorm.data, insitu.matrix=insitu.matrix, geometry=geometry)
#get binasied dge input dataset
dm_scnorm <- binarizeSingleCellData(dm_scnorm, seq(0.15, 0.5, 0.01))
write.table(dm_scnorm@binarized.data,file = "dge_scnorm_binarized_with_raw.csv",sep = ",",row.names = T,col.names = T)
#MCC counts of mapping of cells to each bins
dm_scnorm <- mapCells(dm_scnorm)
write.table(dm_scnorm@mcc.scores,file = "dge_scnorm_mccscores_with_raw.csv",sep = ",",row.names = T,col.names = T)
#didn't understand ???????????????????????????????????????????????????????????????????
pha_scnorm = computeVISH(dm_scnorm, 'sna', threshold=0.75)
computeGeneGradient(dm_scnorm, 'sna') 
# retrieve best quantile score
bq_scnorm <- BestQuantile(dm_scnorm, seq(0.15, 0.5, 0.01))

###################################################################################################
#                                   Comparison 
###################################################################################################

#### by means
mean(dm@mcc.scores)
mean(dm_scnorm@mcc.scores)
mean(dm_saver@mcc.scores)
mean(dm_scimpute@mcc.scores)
mean(dm_magic@mcc.scores)


#### ttest
t.test(as.vector(dm@mcc.scores),as.vector(dm_scnorm@mcc.scores))
t.test(as.vector(dm@mcc.scores),as.vector(dm_saver@mcc.scores))
t.test(as.vector(dm@mcc.scores),as.vector(dm_scimpute@mcc.scores))


#### corr between mccs
cor(as.vector(dm@mcc.scores),as.vector(dm_scnorm@mcc.scores))
cor(as.vector(dm@mcc.scores),as.vector(dm_saver@mcc.scores))
cor(as.vector(dm@mcc.scores),as.vector(dm_scimpute@mcc.scores))
cor(as.vector(dm@mcc.scores),as.vector(dm_magic@mcc.scores))

#### rmse
#install.packages("Metrics")
library(Metrics)
rmse(as.vector(dm@mcc.scores),as.vector(dm_scnorm@mcc.scores))
rmse(as.vector(dm@mcc.scores),as.vector(dm_saver@mcc.scores))
rmse(as.vector(dm@mcc.scores),as.vector(dm_scimpute@mcc.scores))


#### by euclidean
dm_max <- c()
for(i in 1:ncol(dm@mcc.scores)){dm_max <- c(dm_max,which.max(dm@mcc.scores[,i]))}
dm_scnorm_max <- c()
for(i in 1:ncol(dm_scnorm@mcc.scores)){dm_scnorm_max <- c(dm_scnorm_max,which.max(dm_scnorm@mcc.scores[,i]))}
geo_dm = geometry[dm_max,]
geo_dm_scnorm = geometry[dm_scnorm_max,]
eucl_dm_x = c()
for(i in 1:nrow(geo_dm)){
  p = (geo_dm[i,1] - geo_dm_scnorm[i,1])**2
  eucl_dm_x <- c(eucl_dm_x,p)
}
eucl_dm_y = c()
for(i in 1:nrow(geo_dm)){
  p = (geo_dm[i,2] - geo_dm_scnorm[i,2])**2
  eucl_dm_y <- c(eucl_dm_y,p)
}
eucl_dm_z = c()
for(i in 1:nrow(geo_dm)){
  p = (geo_dm[i,3] - geo_dm_scnorm[i,3])**2
  eucl_dm_z <- c(eucl_dm_z,p)
}
eucl_dm = c()
for(i in 1:nrow(geo_dm)){
  sum = sqrt(eucl_dm_x[i] + eucl_dm_y[i] + eucl_dm_z[i])
  eucl_dm <- c(eucl_dm, sum)
}

