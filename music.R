setwd("/sc/arion/projects/MetaDope/Teesta/Proj_Placenta/singlecell/singlecelldata/MUSIC")
library(SingleCellExperiment)
library(data.table)
library(dplyr)
library(tidyverse)
bulk = read.csv("../rat.placenta.RNAseq/counts_STARalign_jan142023/meanGreaterThan10_rawcounts_placenta_rnaseq_january172023.csv")
rownames(bulk) = bulk[,1]
bulk = bulk[,2:ncol(bulk)]
bulk.mtx <- as.matrix(bulk)
bulk.meta = read.table("../rat.placenta.RNAseq/metadata.txt")
head(bulk.meta)
rownames(bulk.meta) = bulk.meta[,1]
bulk.meta = bulk.meta[,2:ncol(bulk.meta)]
colnames(bulk.meta) = c("Group","Sex")
bulk.meta = bulk.meta[2:nrow(bulk.meta),]
bulk.meta.mtx = as.matrix(bulk.meta)
single_cell_counts = read.table("../rat_placenta_marsh/GSE152248_AllStages_AllNuclei_datamatrix.txt")
sc.exprs <- as.matrix(single_cell_counts)
marsh.meta = read.table("../rat_placenta_marsh/GSE152248_AllStages_AllNuclei_clusters.txt)
marsh.meta.sc = as.matrix(marsh.meta)  






single_cell_counts = read.table("rat_placenta_marsh/marsh_placenta_trophoblast_counts.txt")
sc.exprs <- as.matrix(single_cell_counts)
marsh.meta = read.table("rat_placenta_marsh/marsh_placenta_trophoblast_metadata.txt")
marsh.meta.sc = as.matrix(marsh.meta)  
all(rownames(marsh.meta.sc) %in% colnames(sc.exprs))
pheno_data <- as.data.frame(marsh.meta.sc)
library(Biobase)                                        

# Convert single-cell count data to a matrix
sc.exprs <- as.matrix(single_cell_counts)  # Replace with actual variable

# Convert single-cell metadata to data frame
pheno_data <- as.data.frame(marsh.meta)

# Ensure row names in metadata match column names in count matrix
rownames(pheno_data) <- colnames(sc.exprs)

# Create an ExpressionSet
sc.eset <- ExpressionSet(assayData = sc.exprs, phenoData = AnnotatedDataFrame(pheno_data))
head(pData(sc.eset))  # Check available metadata columns
pData(sc.eset)$sampleID <- rownames(marsh.meta)
library(SingleCellExperiment)
sc.sce <- SingleCellExperiment(assays = list(counts = as.matrix(exprs(sc.eset))))
is(sc.sce, "SingleCellExperiment")
colnames(colData(sc.sce))
colData(sc.sce) <- DataFrame(pData(sc.eset))
colnames(colData(sc.sce))
head(sc.sce)
table(colData(sc.sce)$cluster)   # Ensure clusters exist
table(colData(sc.sce)$sampleID)  # Ensure sample IDs exist
Est.prop <- music_prop(
  bulk.mtx = bulk.mtx,
  sc.sce = sc.sce,
  clusters = "cluster",
  samples = "sampleID",
  select.ct = NULL,
  verbose = FALSE
)
names(Est.prop)
[1] "Est.prop.weighted" "Est.prop.allgene"  "Weight.gene"      
[4] "r.squared.full"    "Var.prop"   
