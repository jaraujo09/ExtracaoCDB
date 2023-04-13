if(!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

#if (!requireNamespace("ComplexHeatmap", quietly = TRUE))
#  BiocManager::install("ComplexHeatmap") 

if (!requireNamespace("TCGAbiolinks", quietly = TRUE))
  BiocManager::install("TCGAbiolinks")

BiocManager::install("qplot")

# BiocManager::install("maftools")

install.packages("ggplot2")

library(TCGAbiolinks)
library(MultiAssayExperiment)
#library(maftools)
#library(dplyr)
library("ggplot2")

projects <- getGDCprojects()
projects$id

## example DESeq
## https://rpubs.com/tiagochst/TCGAbiolinks_to_DESEq2

proj <- "TCGA-BRCA"
query <- GDCquery(
  project = proj,
  data.category = "Transcriptome Profiling", 
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)
GDCdownload(query)
data_rna_brca  <- GDCprepare(query)

class(data_rna_brca)
dim(data_rna_brca)

data_rna_brca$paper_BRCA_Pathology
data_rna_brca$age_at_index
data_rna_brca$gender
data_rna_brca$icd_10_code

meta_brca = colData(data_rna_brca)
dim(meta_brca)
meta_brca$patient
meta_brca$`paper_Pan-Gyn Clusters`

BiocManager::install("DESeq2")
library(DESeq2)

data_de <- data_rna_brca[,!is.na(data_rna_brca$definition)]

ddsSE <- DESeqDataSet(data_de, design = ~ definition)

keep <- rowSums(counts(ddsSE)) >= 10
ddsSE <- ddsSE[keep,]
ddsSE =DESeq(ddsSE)

resultsNames(ddsSE)

res1 <- results(ddsSE, name = "definition_Primary.solid.Tumor_vs_Metastatic")
res1
res2<- results(ddsSE, name = "definition_Solid.Tissue.Normal_vs_Metastatic")

dea1 <- as.data.frame(res1)
dea2 <- as.data.frame(res2)

plot(dea1)
plot(dea2)

summary(res1)
summary(res2)
sum(res1$padj < 0.1, na.rm=TRUE)  #1523 genes diferencialmente expressos dos quais 897 sao sobrexpressos e 626 sao subexpressos 
head(assay(data_rna_brca))
colSums(assay(data_rna_brca))
colData(data_rna_brca)

data_rna_brca$sample_type
DESeq2::plotMA(res1, main="DESeq2", ylim=c(-10,25))
DESeq2::plotMA(res2, main="DESeq2", ylim=c(-12,25))
## dados clinicos

query_clin <- GDCquery(project = "TCGA-GBM", 
                       data.category = "Clinical",
                       data.type = "Clinical Supplement", 
                       data.format = "BCR Biotab")
GDCdownload(query_clin)
clinical.gbm <- GDCprepare(query_clin)
names(clinical.gbm)

head (clinical.gbm$clinical_drug_gbm )
df = as.data.frame(clinical.gbm$clinical_patient_gbm)
View(df)


