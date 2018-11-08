#--------------------------------------------------------------------------------------------------
#                                   DOWNLOAD LINKS
#--------------------------------------------------------------------------------------------------

if(!all(file.exists(c("dge_raw.txt.gz","dge_normalized.txt.gz","binarized_bdtnp.csv.gz","bdtnp.txt.gz","geometry.txt.gz")))){
  download.file("http://bimsbstatic.mdc-berlin.de/rajewsky/DVEX/dge_raw.txt.gz",destfile = "dge_raw.txt.gz")
  download.file("http://bimsbstatic.mdc-berlin.de/rajewsky/DVEX/dge_normalized.txt.gz",destfile = "dge_normalized.txt.gz")
  download.file("http://bimsbstatic.mdc-berlin.de/rajewsky/DVEX/binarized_bdtnp.csv.gz",destfile = "binarized_bdtnp.csv.gz")
  download.file("http://bimsbstatic.mdc-berlin.de/rajewsky/DVEX/bdtnp.txt.gz",destfile = "bdtnp.txt.gz")
  download.file("http://bimsbstatic.mdc-berlin.de/rajewsky/DVEX/geometry.txt.gz",destfile = "geometry.txt.gz")
}

library(R.utils)

gunzip("bdtnp.txt.gz")
gunzip("binarized_bdtnp.csv.gz")
gunzip("dge_normalized.txt.gz")
gunzip("dge_raw.txt.gz")#,"geometry.txt.gz")
gunzip("geometry.txt.gz")


#--------------------------------------------------------------------------------------------------
#                                   RUNNING DISTMAP
#--------------------------------------------------------------------------------------------------

#installation 
install.packages("devtools","ggplot2")
install.packages("ggplot2")
library(devtools)
install_github("rajewsky-lab/DistMap")

#IMPORTANT
#edited two names in vi (dge_raw&norm) before running this at pos 387- betaCOPsg and 7348 - PP2A-B

library(ggplot2)
library(DistMap)

raw.data = read.table("dge_raw.txt", header= FALSE, sep="\t", row.names = 1)
raw.data = as.matrix(raw.data)
raw.data.genes = rownames(raw.data)

norm.data = read.table("dge_normalized.txt", header= TRUE, sep="\t",row.names = 1)
norm.data = as.matrix(norm.data)
norm.data.genes = rownames(norm.data)
norm.data.cell = colnames(norm.data)

#check whether gene names match in order too in raw.data and norm.data
stopifnot(all(norm.data.genes == raw.data.genes))

#load 84 driver gene expression matrix in 3039 bins
insitu.matrix = read.table("binarized_bdtnp.csv", sep = ",",header = T)
insitu.matrix = as.matrix(insitu.matrix)
#changes in colnames in insitu.matrix
colnames(insitu.matrix)[36] = "E(spl)m5-HLH"
colnames(insitu.matrix)[6] = "Blimp-1"
#driver genes extraction and check for missing genes in dge_raw&norm
insitu.driver = colnames(insitu.matrix)
missingGenes = insitu.driver[which(!insitu.driver %in% norm.data.genes)]
print(missingGenes)
#check whether gene names match
stopifnot(all(insitu.driver %in% raw.data.genes))

#loading geometry data
geometry = read.csv("geometry.txt", sep= " ", header = TRUE)
geometry = as.matrix(geometry)

#### Running DistMap
#input data into DIstMap class
dm = new("DistMap", raw.data=raw.data,data=norm.data,
         insitu.matrix=insitu.matrix, geometry=geometry)
#get binasied dge input dataset
dm <- binarizeSingleCellData(dm, seq(0.15, 0.5, 0.01))
#MCC counts of mapping of cells to each bins
dm <- mapCells(dm)

#results
dm_mcc <- dm@mcc.scores
dm_actual <- c()
for (i in 1:ncol(dm_mcc)){
  dm_cell <- dm_mcc[,i]
  dm_cell_ordered <- order(dm_cell, decreasing = T)[1]
  dm_actual <- c(dm_actual, dm_cell_ordered)
}
#saving results
saveRDS(dm_actual,"dm_actual.rds")
saveRDS(dm_mcc,"dm_mcc.rds")
