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
