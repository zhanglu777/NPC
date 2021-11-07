#/mnt/data/user_data/zhaolei/program/database/TCGA

library(TCGAbiolinks)
library(SummarizedExperiment)
library(stringr)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(stats4)
library(BiocGenerics)
library(parallel)


#normal
TCGAbiolinks:::getProjectSummary("TCGA-BLCA")
query <- GDCquery(project = "TCGA-BLCA", 
                 legacy = FALSE, 
                 experimental.strategy = "RNA-Seq", 
                 data.category = "Transcriptome Profiling", 
                 data.type = "Gene Expression Quantification", 
                 workflow.type = "HTSeq - Counts",
                 sample.type = c("Solid Tissue Normal"))
setwd("/mnt/data/user_data/zlu/TCGA_database/TCGA_BLCA")
GDCdownload(query)
data <- GDCprepare(query)
matrix <- assay(data)
write.csv(matrix,"TCGA_BLCA_normal_raw_count.csv")




TCGAbiolinks:::getProjectSummary("TCGA-ESCA")
query <- GDCquery(project = "TCGA-ESCA", 
                 legacy = FALSE, 
                 experimental.strategy = "RNA-Seq", 
                 data.category = "Transcriptome Profiling", 
                 data.type = "Gene Expression Quantification", 
                 workflow.type = "HTSeq - Counts",
                 sample.type = c("Solid Tissue Normal"))
setwd("/mnt/data/user_data/zlu/TCGA_database/TCGA_ESCA")
GDCdownload(query)
data <- GDCprepare(query)
matrix <- assay(data)
write.csv(matrix,"TCGA_ESCA_normal_raw_count.csv")



TCGAbiolinks:::getProjectSummary("TCGA-LIHC")
query <- GDCquery(project = "TCGA-LIHC", 
                 legacy = FALSE, 
                 experimental.strategy = "RNA-Seq", 
                 data.category = "Transcriptome Profiling", 
                 data.type = "Gene Expression Quantification", 
                 workflow.type = "HTSeq - Counts",
                 sample.type = c("Solid Tissue Normal"))
setwd("/mnt/data/user_data/zlu/TCGA_database/TCGA_LIHC/")
GDCdownload(query)
data <- GDCprepare(query)
matrix <- assay(data)
head(matrix)
write.csv(matrix,"TCGA_LIHC_normal_raw_count.csv")


TCGAbiolinks:::getProjectSummary("TCGA-STAD")
query <- GDCquery(project = "TCGA-STAD", 
                 legacy = FALSE, 
                 experimental.strategy = "RNA-Seq", 
                 data.category = "Transcriptome Profiling", 
                 data.type = "Gene Expression Quantification", 
                 workflow.type = "HTSeq - Counts",
                 sample.type = c("Solid Tissue Normal"))
setwd("/mnt/data/user_data/zlu/TCGA_database/TCGA_STAD/")
GDCdownload(query)
data <- GDCprepare(query)
matrix <- assay(data)
head(matrix)
write.csv(matrix,"TCGA_STAD_normal_raw_count.csv")


TCGAbiolinks:::getProjectSummary("TCGA-UCEC")
query <- GDCquery(project = "TCGA-UCEC", 
                 legacy = FALSE, 
                 experimental.strategy = "RNA-Seq", 
                 data.category = "Transcriptome Profiling", 
                 data.type = "Gene Expression Quantification", 
                 workflow.type = "HTSeq - Counts",
                 sample.type = c("Solid Tissue Normal"))
setwd("/mnt/data/user_data/zlu/TCGA_database/TCGA_UCEC/")
GDCdownload(query)
data <- GDCprepare(query)
matrix <- assay(data)
head(matrix)
write.csv(matrix,"TCGA_UCEC_normal_raw_count.csv")



TCGAbiolinks:::getProjectSummary("TCGA-LUAD")
query <- GDCquery(project = "TCGA-LUAD", 
                 legacy = FALSE, 
                 experimental.strategy = "RNA-Seq", 
                 data.category = "Transcriptome Profiling", 
                 data.type = "Gene Expression Quantification", 
                 workflow.type = "HTSeq - Counts",
                 sample.type = c("Solid Tissue Normal"))
setwd("/mnt/data/user_data/zlu/TCGA_database/TCGA_LUAD/")
GDCdownload(query)
data <- GDCprepare(query)
matrix <- assay(data)
head(matrix)
write.csv(matrix,"TCGA_LUAD_normal_raw_count.csv")

#每次做完都要清楚变量，并重新载入数据
cd /mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/nasopha_vs_others/human_tissue



suppressPackageStartupMessages({
        library(Rsamtools)
        library(GenomicFeatures)
        library(GenomicAlignments)
        library(BiocParallel)
        library(DESeq2)
        library(pheatmap)
        library(RColorBrewer)
        library(PoiClaClu)
        library(org.Mm.eg.db)
        library(AnnotationDbi)
        library(DOSE)
        library(clusterProfiler)
        library(topGO)
        library(pathview)
        library(org.Hs.eg.db)
        library(AnnotationDbi)
        library(DOSE)
        library(clusterProfiler)
        library(topGO)
        library(ggplot2)
})

load(file="/mnt/data/user_data/xiangyu/workshop/WORKFLOW_RNAseq/ebg_GRCh38.RData")
load(file="/mnt/data/user_data/xiangyu/workshop/WORKFLOW_RNAseq/txdb_GRCh38.RData")

library(future)
library(future.apply)
options(future.globals.maxSize = 3000 * 1024^2)
plan("multiprocess", workers = 8)
plan()
library("ggplot2")
library("ggrepel")
library("dplyr")
library("biomaRt")
library("remotes")
library("nichenetr")

***********************************************************
***********************************************************
***********************************************************
nasopha <- read.csv(file = "nasopha_normal_raw_count.csv",row.names = 1)
blader <- read.csv(file = "/mnt/data/user_data/zlu/TCGA_database/TCGA_BLCA/TCGA_BLCA_normal_raw_count.csv",row.names = 1)

nasopha_vs_blader <- merge(blader,nasopha,by = "row.names", all=FALSE)
rownames(nasopha_vs_blader) <-  nasopha_vs_blader[,1]
nasopha_vs_blader <-  nasopha_vs_blader[,-1]
colnames(nasopha_vs_blader)
id <- colnames(nasopha_vs_blader)
deal <- c(rep('blader',19),rep('nasopha',3))
nasopha_vs_blader_df <- data.frame(id,deal)
dds <- DESeqDataSetFromMatrix(countData = nasopha_vs_blader, colData = nasopha_vs_blader_df, design = ~ deal)

countdata <- nasopha_vs_blader
coldata <- nasopha_vs_blader_df
DESeq_counts <- DESeq(dds)
normalized_DEseq <- counts(DESeq_counts, normalized=TRUE)
DEseq_res <- results(DESeq_counts,contrast=c("deal","nasopha","blader"))
res_1 <- cbind(normalized_DEseq,DEseq_res)

res_1$symbol <- mapIds(x = org.Hs.eg.db,
                        keys = rownames(res_1),
                        keytype ="ENSEMBL",
                        column ="SYMBOL",
                        multiVals="first")
res_1$entrez <- mapIds(x = org.Hs.eg.db,
                        keys = rownames(res_1),
                        keytype ="ENSEMBL",
                        column ="ENTREZID",
                        multiVals="first")
AA <- res_1$symbol
AA <- as.character(AA)
res_1$GENENAME <- mapIds(x = org.Hs.eg.db,
                        keys = AA,
                        keytype ="SYMBOL",
                        column ="GENENAME",
                        multiVals="first")
all_summry <- cbind(countdata,res_1)
head(all_summry)
write.csv(all_summry, "nasopha_vs_blader_normal_tissue_DESeq.csv")

all_TF <- all_summry
aa <- as.character(all_TF$symbol) %>% convert_human_to_mouse_symbols()
all_TF$mouse_symbol <- as.character(aa)
matrix <- all_TF[!duplicated(all_TF$mouse_symbol),]
head(matrix)
write.csv(matrix,"nasopha_vs_blader_normal_tissue_DESeq_mmsymbol.csv")
rm(list = ls())

***********************************************************
***********************************************************
***********************************************************
nasopha <- read.csv(file = "nasopha_normal_raw_count.csv",row.names = 1)
esophag <- read.csv(file = "/mnt/data/user_data/zlu/TCGA_database/TCGA_ESCA/TCGA_ESCA_normal_raw_count.csv",row.names = 1)

nasopha_vs_esophag <- merge(esophag,nasopha,by = "row.names", all=FALSE)
rownames(nasopha_vs_esophag) <-  nasopha_vs_esophag[,1]
nasopha_vs_esophag <-  nasopha_vs_esophag[,-1]
colnames(nasopha_vs_esophag)
id <- colnames(nasopha_vs_esophag)
deal <- c(rep('esophag',11),rep('nasopha',3))
nasopha_vs_esophag_df <- data.frame(id,deal)
dds <- DESeqDataSetFromMatrix(countData = nasopha_vs_esophag, colData = nasopha_vs_esophag_df, design = ~ deal)

countdata <- nasopha_vs_esophag
coldata <- nasopha_vs_esophag_df
DESeq_counts <- DESeq(dds)
normalized_DEseq <- counts(DESeq_counts, normalized=TRUE)

DEseq_res <- results(DESeq_counts,contrast=c("deal","nasopha","esophag"))
res_1 <- cbind(normalized_DEseq,DEseq_res)

res_1$symbol <- mapIds(x = org.Hs.eg.db,
                        keys = rownames(res_1),
                        keytype ="ENSEMBL",
                        column ="SYMBOL",
                        multiVals="first")
res_1$entrez <- mapIds(x = org.Hs.eg.db,
                        keys = rownames(res_1),
                        keytype ="ENSEMBL",
                        column ="ENTREZID",
                        multiVals="first")
AA <- res_1$symbol
AA <- as.character(AA)
res_1$GENENAME <- mapIds(x = org.Hs.eg.db,
                        keys = AA,
                        keytype ="SYMBOL",
                        column ="GENENAME",
                        multiVals="first")
all_summry <- cbind(countdata,res_1)
head(all_summry)
write.csv(all_summry, "nasopha_vs_esophag_normal_tissue_DESeq.csv")
all_TF <- all_summry
aa <- as.character(all_TF$symbol) %>% convert_human_to_mouse_symbols()
all_TF$mouse_symbol <- as.character(aa)
matrix <- all_TF[!duplicated(all_TF$mouse_symbol),]
head(matrix)
write.csv(matrix,"nasopha_vs_esophag_normal_tissue_DESeq_mmsymbol.csv")
rm(list = ls())

***********************************************************
***********************************************************
***********************************************************
nasopha <- read.csv(file = "nasopha_normal_raw_count.csv",row.names = 1)
liver <- read.csv(file = "/mnt/data/user_data/zlu/TCGA_database/TCGA_LIHC/TCGA_LIHC_normal_raw_count.csv",row.names = 1)

nasopha_vs_liver <- merge(liver,nasopha,by = "row.names", all=FALSE)
rownames(nasopha_vs_liver) <-  nasopha_vs_liver[,1]
nasopha_vs_liver <-  nasopha_vs_liver[,-1]
colnames(nasopha_vs_liver)
id <- colnames(nasopha_vs_liver)
deal <- c(rep('liver',50),rep('nasopha',3))
nasopha_vs_liver_df <- data.frame(id,deal)
dds <- DESeqDataSetFromMatrix(countData = nasopha_vs_liver, colData = nasopha_vs_liver_df, design = ~ deal)

countdata <- nasopha_vs_liver
coldata <- nasopha_vs_liver_df
DESeq_counts <- DESeq(dds)
normalized_DEseq <- counts(DESeq_counts, normalized=TRUE)

DEseq_res <- results(DESeq_counts,contrast=c("deal","nasopha","liver"))
res_1 <- cbind(normalized_DEseq,DEseq_res)

res_1$symbol <- mapIds(x = org.Hs.eg.db,
                        keys = rownames(res_1),
                        keytype ="ENSEMBL",
                        column ="SYMBOL",
                        multiVals="first")
res_1$entrez <- mapIds(x = org.Hs.eg.db,
                        keys = rownames(res_1),
                        keytype ="ENSEMBL",
                        column ="ENTREZID",
                        multiVals="first")
AA <- res_1$symbol
AA <- as.character(AA)
res_1$GENENAME <- mapIds(x = org.Hs.eg.db,
                        keys = AA,
                        keytype ="SYMBOL",
                        column ="GENENAME",
                        multiVals="first")
all_summry <- cbind(countdata,res_1)
head(all_summry)
write.csv(all_summry, "nasopha_vs_liver_normal_tissue_DESeq.csv")

all_TF <- all_summry
aa <- as.character(all_TF$symbol) %>% convert_human_to_mouse_symbols()
all_TF$mouse_symbol <- as.character(aa)
matrix <- all_TF[!duplicated(all_TF$mouse_symbol),]
head(matrix)
write.csv(matrix,"nasopha_vs_liver_normal_tissue_DESeq_mmsymbol.csv")
rm(list = ls())

***********************************************************
***********************************************************
***********************************************************
nasopha <- read.csv(file = "nasopha_normal_raw_count.csv",row.names = 1)
stomach <- read.csv(file = "/mnt/data/user_data/zlu/TCGA_database/TCGA_STAD/TCGA_STAD_normal_raw_count.csv",row.names = 1)

nasopha_vs_stomach <- merge(stomach,nasopha,by = "row.names", all=FALSE)
rownames(nasopha_vs_stomach) <-  nasopha_vs_stomach[,1]
nasopha_vs_stomach <-  nasopha_vs_stomach[,-1]
colnames(nasopha_vs_stomach)
id <- colnames(nasopha_vs_stomach)
deal <- c(rep('stomach',32),rep('nasopha',3))
nasopha_vs_stomach_df <- data.frame(id,deal)
dds <- DESeqDataSetFromMatrix(countData = nasopha_vs_stomach, colData = nasopha_vs_stomach_df, design = ~ deal)

countdata <- nasopha_vs_stomach
coldata <- nasopha_vs_stomach_df
DESeq_counts <- DESeq(dds)
normalized_DEseq <- counts(DESeq_counts, normalized=TRUE)

DEseq_res <- results(DESeq_counts,contrast=c("deal","nasopha","stomach"))
res_1 <- cbind(normalized_DEseq,DEseq_res)

res_1$symbol <- mapIds(x = org.Hs.eg.db,
                        keys = rownames(res_1),
                        keytype ="ENSEMBL",
                        column ="SYMBOL",
                        multiVals="first")
res_1$entrez <- mapIds(x = org.Hs.eg.db,
                        keys = rownames(res_1),
                        keytype ="ENSEMBL",
                        column ="ENTREZID",
                        multiVals="first")
AA <- res_1$symbol
AA <- as.character(AA)
res_1$GENENAME <- mapIds(x = org.Hs.eg.db,
                        keys = AA,
                        keytype ="SYMBOL",
                        column ="GENENAME",
                        multiVals="first")
all_summry <- cbind(countdata,res_1)
head(all_summry)
write.csv(all_summry, "nasopha_vs_stomach_normal_tissue_DESeq.csv")
all_TF <- all_summry
aa <- as.character(all_TF$symbol) %>% convert_human_to_mouse_symbols()
all_TF$mouse_symbol <- as.character(aa)
matrix <- all_TF[!duplicated(all_TF$mouse_symbol),]
head(matrix)
write.csv(matrix,"nasopha_vs_stomach_normal_tissue_DESeq_mmsymbol.csv")
rm(list = ls())
***********************************************************
***********************************************************
***********************************************************
nasopha <- read.csv(file = "nasopha_normal_raw_count.csv",row.names = 1)
uteri <- read.csv(file = "/mnt/data/user_data/zlu/TCGA_database/TCGA_UCEC/TCGA_UCEC_normal_raw_count.csv",row.names = 1)

nasopha_vs_uteri <- merge(uteri,nasopha,by = "row.names", all=FALSE)
rownames(nasopha_vs_uteri) <-  nasopha_vs_uteri[,1]
nasopha_vs_uteri <-  nasopha_vs_uteri[,-1]
colnames(nasopha_vs_uteri)
id <- colnames(nasopha_vs_uteri)
deal <- c(rep('uteri',35),rep('nasopha',3))
nasopha_vs_uteri_df <- data.frame(id,deal)
dds <- DESeqDataSetFromMatrix(countData = nasopha_vs_uteri, colData = nasopha_vs_uteri_df, design = ~ deal)

countdata <- nasopha_vs_uteri
coldata <- nasopha_vs_uteri_df
DESeq_counts <- DESeq(dds)
normalized_DEseq <- counts(DESeq_counts, normalized=TRUE)

DEseq_res <- results(DESeq_counts,contrast=c("deal","nasopha","uteri"))
res_1 <- cbind(normalized_DEseq,DEseq_res)

res_1$symbol <- mapIds(x = org.Hs.eg.db,
                        keys = rownames(res_1),
                        keytype ="ENSEMBL",
                        column ="SYMBOL",
                        multiVals="first")
res_1$entrez <- mapIds(x = org.Hs.eg.db,
                        keys = rownames(res_1),
                        keytype ="ENSEMBL",
                        column ="ENTREZID",
                        multiVals="first")
AA <- res_1$symbol
AA <- as.character(AA)
res_1$GENENAME <- mapIds(x = org.Hs.eg.db,
                        keys = AA,
                        keytype ="SYMBOL",
                        column ="GENENAME",
                        multiVals="first")
all_summry <- cbind(countdata,res_1)
head(all_summry)
write.csv(all_summry, "nasopha_vs_uteri_normal_tissue_DESeq.csv")
all_TF <- all_summry
aa <- as.character(all_TF$symbol) %>% convert_human_to_mouse_symbols()
all_TF$mouse_symbol <- as.character(aa)
matrix <- all_TF[!duplicated(all_TF$mouse_symbol),]
head(matrix)
write.csv(matrix,"nasopha_vs_uteri_normal_tissue_DESeq_mmsymbol.csv")
rm(list = ls())

***********************************************************
***********************************************************
***********************************************************
nasopha <- read.csv(file = "nasopha_normal_raw_count.csv",row.names = 1)
lung <- read.csv(file = "/mnt/data/user_data/zlu/TCGA_database/TCGA_LUAD/TCGA_LUAD_normal_raw_count.csv",row.names = 1)
nasopha_vs_lung <- merge(lung,nasopha,by = "row.names", all=FALSE)
rownames(nasopha_vs_lung) <-  nasopha_vs_lung[,1]
nasopha_vs_lung <-  nasopha_vs_lung[,-1]
colnames(nasopha_vs_lung)
id <- colnames(nasopha_vs_lung)
deal <- c(rep('lung',59),rep('nasopha',3))
nasopha_vs_lung_df <- data.frame(id,deal)
dds <- DESeqDataSetFromMatrix(countData = nasopha_vs_lung, colData = nasopha_vs_lung_df, design = ~ deal)

countdata <- nasopha_vs_lung
coldata <- nasopha_vs_lung_df
DESeq_counts <- DESeq(dds)
normalized_DEseq <- counts(DESeq_counts, normalized=TRUE)

DEseq_res <- results(DESeq_counts,contrast=c("deal","nasopha","lung"))
res_1 <- cbind(normalized_DEseq,DEseq_res)

res_1$symbol <- mapIds(x = org.Hs.eg.db,
                        keys = rownames(res_1),
                        keytype ="ENSEMBL",
                        column ="SYMBOL",
                        multiVals="first")
res_1$entrez <- mapIds(x = org.Hs.eg.db,
                        keys = rownames(res_1),
                        keytype ="ENSEMBL",
                        column ="ENTREZID",
                        multiVals="first")
AA <- res_1$symbol
AA <- as.character(AA)
res_1$GENENAME <- mapIds(x = org.Hs.eg.db,
                        keys = AA,
                        keytype ="SYMBOL",
                        column ="GENENAME",
                        multiVals="first")
all_summry <- cbind(countdata,res_1)
head(all_summry)
write.csv(all_summry, "nasopha_vs_lung_normal_tissue_DESeq.csv")
all_TF <- all_summry
aa <- as.character(all_TF$symbol) %>% convert_human_to_mouse_symbols()
all_TF$mouse_symbol <- as.character(aa)
matrix <- all_TF[!duplicated(all_TF$mouse_symbol),]
head(matrix)
write.csv(matrix,"nasopha_vs_lung_normal_tissue_DESeq_mmsymbol.csv")
rm(list = ls())

******************************************************************others
nasopha <- read.csv(file = "nasopha_normal_raw_count.csv",row.names = 1)
blader <- read.csv(file = "/mnt/data/user_data/zlu/TCGA_database/TCGA_BLCA/TCGA_BLCA_normal_raw_count.csv",row.names = 1)
esophag <- read.csv(file = "/mnt/data/user_data/zlu/TCGA_database/TCGA_ESCA/TCGA_ESCA_normal_raw_count.csv",row.names = 1)
liver <- read.csv(file = "/mnt/data/user_data/zlu/TCGA_database/TCGA_LIHC/TCGA_LIHC_normal_raw_count.csv",row.names = 1)
stomach <- read.csv(file = "/mnt/data/user_data/zlu/TCGA_database/TCGA_STAD/TCGA_STAD_normal_raw_count.csv",row.names = 1)
uteri <- read.csv(file = "/mnt/data/user_data/zlu/TCGA_database/TCGA_UCEC/TCGA_UCEC_normal_raw_count.csv",row.names = 1)
lung <- read.csv(file = "/mnt/data/user_data/zlu/TCGA_database/TCGA_LUAD/TCGA_LUAD_normal_raw_count.csv",row.names = 1)

others <- cbind(blader,esophag,liver,stomach,uteri,lung)


nasopha_vs_others <- merge(others,nasopha,by = "row.names", all=FALSE)
rownames(nasopha_vs_others) <-  nasopha_vs_others[,1]
nasopha_vs_others <-  nasopha_vs_others[,-1]
colnames(nasopha_vs_others)

id <- colnames(nasopha_vs_others)
class <- c(rep('blader',19),rep('esophag',11),rep('liver',50),rep('stomach',32),rep('uteri',35),rep('lung',59),rep('nasopha',3))
deal <- c(rep('others',206),rep('nasopha',3))
nasopha_vs_others_df <- data.frame(id,class,deal)

dds <- DESeqDataSetFromMatrix(countData = nasopha_vs_others, colData = nasopha_vs_others_df, design = ~ deal)

countdata <- nasopha_vs_others
coldata <- nasopha_vs_others_df
DESeq_counts <- DESeq(dds)
normalized_DEseq <- counts(DESeq_counts, normalized=TRUE)

DEseq_res <- results(DESeq_counts,contrast=c("deal","nasopha","others"))
res_1 <- cbind(normalized_DEseq,DEseq_res)

res_1$symbol <- mapIds(x = org.Hs.eg.db,
                        keys = rownames(res_1),
                        keytype ="ENSEMBL",
                        column ="SYMBOL",
                        multiVals="first")
res_1$entrez <- mapIds(x = org.Hs.eg.db,
                        keys = rownames(res_1),
                        keytype ="ENSEMBL",
                        column ="ENTREZID",
                        multiVals="first")
AA <- res_1$symbol
AA <- as.character(AA)
res_1$GENENAME <- mapIds(x = org.Hs.eg.db,
                        keys = AA,
                        keytype ="SYMBOL",
                        column ="GENENAME",
                        multiVals="first")
all_summry <- cbind(countdata,res_1)
head(all_summry)
write.csv(all_summry, "nasopha_vs_others_normal_tissue_DESeq.csv")
all_TF <- all_summry
aa <- as.character(all_TF$symbol) %>% convert_human_to_mouse_symbols()
all_TF$mouse_symbol <- as.character(aa)
matrix <- all_TF[!duplicated(all_TF$mouse_symbol),]
head(matrix)
write.csv(matrix,"nasopha_vs_others_normal_tissue_DESeq_mmsymbol.csv")
rm(list = ls())

******************************************************************3 tissue
nasopha <- read.csv(file = "nasopha_normal_raw_count.csv",row.names = 1)
liver <- read.csv(file = "/mnt/data/user_data/zlu/TCGA_database/TCGA_LIHC/TCGA_LIHC_normal_raw_count.csv",row.names = 1)
uteri <- read.csv(file = "/mnt/data/user_data/zlu/TCGA_database/TCGA_UCEC/TCGA_UCEC_normal_raw_count.csv",row.names = 1)
lung <- read.csv(file = "/mnt/data/user_data/zlu/TCGA_database/TCGA_LUAD/TCGA_LUAD_normal_raw_count.csv",row.names = 1)

three_tissue <- cbind(liver,uteri,lung)


nasopha_vs_three_tissue <- merge(three_tissue,nasopha,by = "row.names", all=FALSE)
rownames(nasopha_vs_three_tissue) <-  nasopha_vs_three_tissue[,1]
nasopha_vs_three_tissue <-  nasopha_vs_three_tissue[,-1]
colnames(nasopha_vs_three_tissue)

id <- colnames(nasopha_vs_three_tissue)
class <- c(rep('liver',50),rep('uteri',35),rep('lung',59),rep('nasopha',3))
deal <- c(rep('others',144),rep('nasopha',3))
nasopha_vs_three_tissue_df <- data.frame(id,class,deal)

dds <- DESeqDataSetFromMatrix(countData = nasopha_vs_three_tissue, colData = nasopha_vs_three_tissue_df, design = ~ deal)

countdata <- nasopha_vs_three_tissue
coldata <- nasopha_vs_three_tissue_df
DESeq_counts <- DESeq(dds)
normalized_DEseq <- counts(DESeq_counts, normalized=TRUE)

DEseq_res <- results(DESeq_counts,contrast=c("deal","nasopha","others"))      #最后一个是CTRL
res_1 <- cbind(normalized_DEseq,DEseq_res)

res_1$symbol <- mapIds(x = org.Hs.eg.db,
                        keys = rownames(res_1),
                        keytype ="ENSEMBL",
                        column ="SYMBOL",
                        multiVals="first")
res_1$entrez <- mapIds(x = org.Hs.eg.db,
                        keys = rownames(res_1),
                        keytype ="ENSEMBL",
                        column ="ENTREZID",
                        multiVals="first")

all_summry <- cbind(countdata,res_1)
head(all_summry)
write.csv(all_summry, "nasopha_vs_three_tissue_normal_tissue_DESeq.csv")
all_TF <- all_summry
aa <- as.character(all_TF$symbol) %>% convert_human_to_mouse_symbols()
all_TF$mouse_symbol <- as.character(aa)
matrix <- all_TF[!duplicated(all_TF$mouse_symbol),]
head(matrix)
write.csv(matrix,"nasopha_vs_three_tissue_normal_tissue_DESeq_mmsymbol.csv")
rm(list = ls())
 
******************************************************************3 tissue 加上esophag
nasopha <- read.csv(file = "nasopha_normal_raw_count.csv",row.names = 1)
liver <- read.csv(file = "/mnt/data/user_data/zlu/TCGA_database/TCGA_LIHC/TCGA_LIHC_normal_raw_count.csv",row.names = 1)
uteri <- read.csv(file = "/mnt/data/user_data/zlu/TCGA_database/TCGA_UCEC/TCGA_UCEC_normal_raw_count.csv",row.names = 1)
lung <- read.csv(file = "/mnt/data/user_data/zlu/TCGA_database/TCGA_LUAD/TCGA_LUAD_normal_raw_count.csv",row.names = 1)
esophag <- read.csv(file = "/mnt/data/user_data/zlu/TCGA_database/TCGA_ESCA/TCGA_ESCA_normal_raw_count.csv",row.names = 1)

four_tissue <- cbind(liver,uteri,lung)


nasopha_vs_four_tissue <- merge(four_tissue,nasopha,by = "row.names", all=FALSE)
rownames(nasopha_vs_four_tissue) <-  nasopha_vs_four_tissue[,1]
nasopha_vs_four_tissue <-  nasopha_vs_four_tissue[,-1]
colnames(nasopha_vs_four_tissue)

id <- colnames(nasopha_vs_four_tissue)
class <- c(rep('liver',50),rep('uteri',35),rep('lung',59),rep('nasopha',3))
deal <- c(rep('others',144),rep('nasopha',3))
nasopha_vs_four_tissue_df <- data.frame(id,class,deal)

dds <- DESeqDataSetFromMatrix(countData = nasopha_vs_four_tissue, colData = nasopha_vs_four_tissue_df, design = ~ deal)

countdata <- nasopha_vs_four_tissue
coldata <- nasopha_vs_four_tissue_df
DESeq_counts <- DESeq(dds)
normalized_DEseq <- counts(DESeq_counts, normalized=TRUE)

DEseq_res <- results(DESeq_counts,contrast=c("deal","nasopha","others"))      #最后一个是CTRL
res_1 <- cbind(normalized_DEseq,DEseq_res)

res_1$symbol <- mapIds(x = org.Hs.eg.db,
                        keys = rownames(res_1),
                        keytype ="ENSEMBL",
                        column ="SYMBOL",
                        multiVals="first")
res_1$entrez <- mapIds(x = org.Hs.eg.db,
                        keys = rownames(res_1),
                        keytype ="ENSEMBL",
                        column ="ENTREZID",
                        multiVals="first")

all_summry <- cbind(countdata,res_1)
head(all_summry)
write.csv(all_summry, "nasopha_vs_four_tissue_normal_tissue_DESeq.csv")
all_TF <- all_summry
aa <- as.character(all_TF$symbol) %>% convert_human_to_mouse_symbols()
all_TF$mouse_symbol <- as.character(aa)
matrix <- all_TF[!duplicated(all_TF$mouse_symbol),]
head(matrix)
write.csv(matrix,"nasopha_vs_four_tissue_normal_tissue_DESeq_mmsymbol.csv")
rm(list = ls())

$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$mouse tissue$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE138103
suppressPackageStartupMessages({
        library(Rsamtools)
        library(GenomicFeatures)
        library(GenomicAlignments)
        library(BiocParallel)
        library(DESeq2)
        library(pheatmap)
        library(RColorBrewer)
        library(PoiClaClu)
        library(org.Mm.eg.db)
        library(AnnotationDbi)
        library(DOSE)
        library(clusterProfiler)
        library(topGO)
        library(pathview)
        library(org.Hs.eg.db)
        library(AnnotationDbi)
        library(DOSE)
        library(clusterProfiler)
        library(topGO)
        library(ggplot2)
})



nasopha <- read.csv(file = "/mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/nasopha_vs_others/mouse_tissue/nasophag_normal_tissue.csv")
three_tissue <- read.csv(file = "/mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/nasopha_vs_others/mouse_tissue/GSE138103_uteri_lung_liver_normal_tissue.csv")

nasopha_vs_three_tissue <- merge(three_tissue,nasopha,by = "X", all=FALSE)
nasopha_vs_three_tissue <- nasopha_vs_three_tissue[!duplicated(nasopha_vs_three_tissue$X),]
rownames(nasopha_vs_three_tissue) <-  nasopha_vs_three_tissue[,1]
nasopha_vs_three_tissue <-  nasopha_vs_three_tissue[,-1]
colnames(nasopha_vs_three_tissue)
id <- colnames(nasopha_vs_three_tissue)
deal <- c(rep('others',12),rep('nasopha',2))
nasopha_vs_three_tissue_df <- data.frame(id,deal)
dds <- DESeqDataSetFromMatrix(countData = nasopha_vs_three_tissue, colData = nasopha_vs_three_tissue_df, design = ~ deal)

countdata <- nasopha_vs_three_tissue
coldata <- nasopha_vs_three_tissue_df
DESeq_counts <- DESeq(dds)
normalized_DEseq <- counts(DESeq_counts, normalized=TRUE)

DEseq_res <- results(DESeq_counts,contrast=c("deal","nasopha","others"))
res_1 <- cbind(normalized_DEseq,DEseq_res)

all_summry <- cbind(countdata,res_1)
head(all_summry)
write.csv(all_summry, "nasopha_vs_three_tissue_normal_tissue_DESeq.csv")

*************************************organoid  GRCm38
nasopha_vs_three_organoid <- read.csv(file = "/mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/nasopha_vs_others/mouse_tissue/organoid_nasopha_vs_3tissue_raw_count.csv",row.names = 1)


colnames(nasopha_vs_three_organoid)
id <- colnames(nasopha_vs_three_organoid)
deal <- c(rep('others',9),rep('nasopha',3))
nasopha_vs_three_organoid_df <- data.frame(id,deal)
dds <- DESeqDataSetFromMatrix(countData = nasopha_vs_three_organoid, colData = nasopha_vs_three_organoid_df, design = ~ deal)

countdata <- nasopha_vs_three_organoid
coldata <- nasopha_vs_three_organoid_df
DESeq_counts <- DESeq(dds)
normalized_DEseq <- counts(DESeq_counts, normalized=TRUE)

DEseq_res <- results(DESeq_counts,contrast=c("deal","nasopha","others"))
res_1 <- cbind(normalized_DEseq,DEseq_res)


load(file="/mnt/data/user_data/xiangyu/workshop/WORKFLOW_RNAseq/ebg_GRCm38.RData")
load(file="/mnt/data/user_data/xiangyu/workshop/WORKFLOW_RNAseq/txdb_GRCm38.RData")
res_1$symbol <- mapIds(x = org.Mm.eg.db,
                        keys = rownames(res_1),
                        keytype ="ENSEMBL",
                        column ="SYMBOL",
                        multiVals="first")


all_summry <- cbind(countdata,res_1)
head(all_summry)
write.csv(all_summry, "nasopha_vs_three_organoid_normal_organoid_DESeq.csv")

*************************************organoid  GRCm38
nasopha_vs_four_organoid <- read.csv(file = "/mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/nasopha_vs_others/organoids/nasopha_vs_others4_count_symbol.csv",row.names = 1)


colnames(nasopha_vs_four_organoid)
id <- colnames(nasopha_vs_four_organoid)
deal <- c(rep('others',12),rep('nasopha',3))
nasopha_vs_four_organoid_df <- data.frame(id,deal)
dds <- DESeqDataSetFromMatrix(countData = nasopha_vs_four_organoid, colData = nasopha_vs_four_organoid_df, design = ~ deal)

countdata <- nasopha_vs_four_organoid
coldata <- nasopha_vs_four_organoid_df
DESeq_counts <- DESeq(dds)
normalized_DEseq <- counts(DESeq_counts, normalized=TRUE)

DEseq_res <- results(DESeq_counts,contrast=c("deal","nasopha","others"))
res_1 <- cbind(normalized_DEseq,DEseq_res)
write.csv(res_1, "nasopha_vs_four_organoid_normal_organoid_DESeq.csv")


res_1 <- na.omit(res_1)
y <- res_1[c(1:nrow(res_1)),20]
res_1 <- res_1[with(res_1,y<0.05),]
y <- res_1[c(1:nrow(res_1)),17]
upres_1 <- res_1[with(res_1,y>=1),]
upres_1 <- upres_1[order(upres_1$log2FoldChange,decreasing = TRUE),]



zscore <- t(apply(all_gene, 1, function(x) (x-mean(x))/sd(x)))
zscore <- na.omit(zscore)
head(zscore)

library(Seurat)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
require(RColorBrewer)
source("/mnt/data/user_data/xiangyu/programme/R_PACKAGES/my_code/MyBestFunction_scRNA.R")
source("/mnt/data/user_data/xiangyu/programme/R_PACKAGES/my_code/Pseudo_CNV_series.R")

SeuratObject <- CreateSeuratObject(counts = zscore, project = "WXD")
SeuratObject$orig.ident <-factor(SeuratObject$orig.ident,levels = c("lung","liver","uteri","esophag","nasopha"))
levels(SeuratObject$orig.ident) 
gene <- rownames(zscore)
pdf("nasopha_marker_heatmap_add_esophag.pdf")
XY_heatmap(seurat_obj=SeuratObject,group="orig.ident", genes=gene,all_num=FALSE,new_names=NULL,labels_rot=90,assay_sel="RNA",color=colorRampPalette(brewer.pal(10, "RdBu"))(101),min_and_max_cut=1.2,show_row_names=FALSE,scale=FALSE,label_size=0,mark_gene=c("Krt5","Krt13","Krt14","Krt15","Krt16","Krt17","Tuba1c","Cfap45","Nes","Muc15","Muc16","Trp63"))
dev.off()



&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&换不同器官





$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$organoid  GRCm38$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

organoids_nasopha_vs_liver_stomach_esophagus <- read.csv(file = "/mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/nasopha_vs_others/nasopha_vs_liver_stomach_esophagus/mouse_organoids_liver_esophag_stomach_counts_mmsymbol.csv",row.names = 1)

organoids_nasopha_vs_liver_stomach_esophagus <- na.omit(organoids_nasopha_vs_liver_stomach_esophagus)
organoids_nasopha_vs_liver_stomach_esophagus <- organoids_nasopha_vs_liver_stomach_esophagus[!duplicated(organoids_nasopha_vs_liver_stomach_esophagus$symbol),]
rownames(organoids_nasopha_vs_liver_stomach_esophagus) <- organoids_nasopha_vs_liver_stomach_esophagus[,13]
organoids_nasopha_vs_liver_stomach_esophagus <- organoids_nasopha_vs_liver_stomach_esophagus[,-13]


colnames(organoids_nasopha_vs_liver_stomach_esophagus)
id <- colnames(organoids_nasopha_vs_liver_stomach_esophagus)
deal <- c(rep('others',9),rep('nasopha',3))
organoids_nasopha_vs_liver_stomach_esophagus_df <- data.frame(id,deal)
dds <- DESeqDataSetFromMatrix(countData = organoids_nasopha_vs_liver_stomach_esophagus, colData = organoids_nasopha_vs_liver_stomach_esophagus_df, design = ~ deal)

countdata <- organoids_nasopha_vs_liver_stomach_esophagus

DESeq_counts <- DESeq(dds)
normalized_DEseq <- counts(DESeq_counts, normalized=TRUE)

DEseq_res <- results(DESeq_counts,contrast=c("deal","nasopha","others"))
res_1 <- cbind(normalized_DEseq,DEseq_res)

all_summry <- cbind(countdata,res_1)
head(all_summry)
write.csv(all_summry, "organoids_nasopha_vs_liver_stomach_esophagus_count_normalized_DESeq.csv")


$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$mouse tissue$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#https://www.nature.com/articles/sdata2017185#Sec8
suppressPackageStartupMessages({
        library(Rsamtools)
        library(GenomicFeatures)
        library(GenomicAlignments)
        library(BiocParallel)
        library(DESeq2)
        library(pheatmap)
        library(RColorBrewer)
        library(PoiClaClu)
        library(org.Mm.eg.db)
        library(AnnotationDbi)
        library(DOSE)
        library(clusterProfiler)
        library(topGO)
        library(pathview)
        library(org.Hs.eg.db)
        library(AnnotationDbi)
        library(DOSE)
        library(clusterProfiler)
        library(topGO)
        library(ggplot2)
})



nasopha <- read.csv(file = "/mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/nasopha_vs_others/mouse_tissue/nasophag_normal_tissue.csv")
liver_stomach_esophagus <- read.csv(file = "/mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/nasopha_vs_others/nasopha_vs_liver_stomach_esophagus/mouse_tissue_liver_esophag_stomach_counts_GRCm38.csv",row.names=1)

load(file="/mnt/data/user_data/xiangyu/workshop/WORKFLOW_RNAseq/ebg_GRCm38.RData")
load(file="/mnt/data/user_data/xiangyu/workshop/WORKFLOW_RNAseq/txdb_GRCm38.RData")
liver_stomach_esophagus$symbol <- mapIds(x = org.Mm.eg.db,
                        keys = rownames(liver_stomach_esophagus),
                        keytype ="ENSEMBL",
                        column ="SYMBOL",
                        multiVals="first")

liver_stomach_esophagus <- na.omit(liver_stomach_esophagus)
liver_stomach_esophagus <- liver_stomach_esophagus[!duplicated(liver_stomach_esophagus$symbol),]
rownames(liver_stomach_esophagus) <- liver_stomach_esophagus[,10]
liver_stomach_esophagus <- liver_stomach_esophagus[,-10]


mus_tissue_nasopha_vs_liver_stomach_esophagus <- merge(liver_stomach_esophagus,nasopha,by.x = "row.names",by.y="X" )

rownames(mus_tissue_nasopha_vs_liver_stomach_esophagus) <-  mus_tissue_nasopha_vs_liver_stomach_esophagus[,1]
mus_tissue_nasopha_vs_liver_stomach_esophagus <-  mus_tissue_nasopha_vs_liver_stomach_esophagus[,-1]
colnames(mus_tissue_nasopha_vs_liver_stomach_esophagus)
id <- colnames(mus_tissue_nasopha_vs_liver_stomach_esophagus)
deal <- c(rep('others',9),rep('nasopha',2))
mus_tissue_nasopha_vs_liver_stomach_esophagus_df <- data.frame(id,deal)
dds <- DESeqDataSetFromMatrix(countData = mus_tissue_nasopha_vs_liver_stomach_esophagus, colData = mus_tissue_nasopha_vs_liver_stomach_esophagus_df, design = ~ deal)

countdata <- mus_tissue_nasopha_vs_liver_stomach_esophagus

DESeq_counts <- DESeq(dds)
normalized_DEseq <- counts(DESeq_counts, normalized=TRUE)

DEseq_res <- results(DESeq_counts,contrast=c("deal","nasopha","others"))
res_1 <- cbind(normalized_DEseq,DEseq_res)

all_summry <- cbind(countdata,res_1)
head(all_summry)
write.csv(all_summry, "mus_tissue_nasopha_vs_liver_stomach_esophagus_count_normalized_DESeq.csv")

$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$human tissue$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
nasopha <- read.csv(file = "/mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/nasopha_vs_others/human_tissue/human_nasopha_normal_raw_count.csv",row.names = 1)
esophag <- read.csv(file = "/mnt/data/user_data/zlu/TCGA_database/TCGA_ESCA/TCGA_ESCA_normal_raw_count.csv",row.names = 1)
liver <- read.csv(file = "/mnt/data/user_data/zlu/TCGA_database/TCGA_LIHC/TCGA_LIHC_normal_raw_count.csv",row.names = 1)
stomach <- read.csv(file = "/mnt/data/user_data/zlu/TCGA_database/TCGA_STAD/TCGA_STAD_normal_raw_count.csv",row.names = 1)


liver_stomach_esophagus <- cbind(esophag,liver,stomach)


human_TCGA_nasopha_vs_liver_stomach_esophagus <- merge(liver_stomach_esophagus,nasopha,by = "row.names", all=FALSE)

load(file="/mnt/data/user_data/xiangyu/workshop/WORKFLOW_RNAseq/ebg_GRCh38.RData")
load(file="/mnt/data/user_data/xiangyu/workshop/WORKFLOW_RNAseq/txdb_GRCh38.RData")

human_TCGA_nasopha_vs_liver_stomach_esophagus$symbol <- mapIds(x = org.Hs.eg.db,
                        keys = rownames(human_TCGA_nasopha_vs_liver_stomach_esophagus),
                        keytype ="ENSEMBL",
                        column ="SYMBOL",
                        multiVals="first")

library(future)
library(future.apply)
options(future.globals.maxSize = 3000 * 1024^2)
plan("multiprocess", workers = 8)
plan()
library("ggplot2")
library("ggrepel")
library("dplyr")
library("biomaRt")
library("remotes")
library("nichenetr")

all_TF <- human_TCGA_nasopha_vs_liver_stomach_esophagus
aa <- as.character(all_TF$symbol) %>% convert_human_to_mouse_symbols()
all_TF$mouse_symbol <- as.character(aa)

human_TCGA_nasopha_vs_liver_stomach_esophagus <- all_TF[!duplicated(all_TF$mouse_symbol),]

human_TCGA_nasopha_vs_liver_stomach_esophagus <-  human_TCGA_nasopha_vs_liver_stomach_esophagus[,c(1:96,98)]
human_TCGA_nasopha_vs_liver_stomach_esophagus <- na.omit(human_TCGA_nasopha_vs_liver_stomach_esophagus)

rownames(human_TCGA_nasopha_vs_liver_stomach_esophagus) <-  human_TCGA_nasopha_vs_liver_stomach_esophagus[,97]
human_TCGA_nasopha_vs_liver_stomach_esophagus_mmID <-  human_TCGA_nasopha_vs_liver_stomach_esophagus[,-97]
colnames(human_TCGA_nasopha_vs_liver_stomach_esophagus_mmID)

id <- colnames(human_TCGA_nasopha_vs_liver_stomach_esophagus_mmID)
deal <- c(rep('others',93),rep('nasopha',3))
human_TCGA_nasopha_vs_liver_stomach_esophagus_mmID_df <- data.frame(id,deal)

dds <- DESeqDataSetFromMatrix(countData = human_TCGA_nasopha_vs_liver_stomach_esophagus_mmID, colData = human_TCGA_nasopha_vs_liver_stomach_esophagus_mmID_df, design = ~ deal)

countdata <- human_TCGA_nasopha_vs_liver_stomach_esophagus_mmID

DESeq_counts <- DESeq(dds)
normalized_DEseq <- counts(DESeq_counts, normalized=TRUE)

DEseq_res <- results(DESeq_counts,contrast=c("deal","nasopha","others"))
res_1 <- cbind(normalized_DEseq,DEseq_res)

all_summry <- cbind(countdata,res_1)
head(all_summry)

write.csv(all_summry, "mmID_human_TCGA_nasopha_vs_liver_stomach_esophagus_count_normalized_DESeq.csv")


***********************做GSEA的grp文件
cd /mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/nasopha_vs_others/nasopha_vs_liver_stomach_esophagus


library("dplyr")
mouse_tissue <- read.csv(file = "/mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/nasopha_vs_others/nasopha_vs_liver_stomach_esophagus/mus_tissue_nasopha_vs_liver_stomach_esophagus_count_normalized_DESeq.csv")
mouse_tissue <- mouse_tissue[,c(1,25,28,29)]
mouse_tissue_padj0.05_up <- mouse_tissue %>% filter(padj <= 0.05) %>% filter(log2FoldChange >= 1)
mouse_tissue_padj0.05_up <- mouse_tissue_padj0.05_up[order(mouse_tissue_padj0.05_up$log2FoldChange,decreasing = TRUE),]

mouse_tissue_padj0.05_up_top100 <- data.frame(mouse_tissue_padj0.05_up[1:100,1])
colnames(mouse_tissue_padj0.05_up_top100) <- "#mus_tissue_nasopha_vs_liver_stomach_esophagus_up_top100"
write.table(mouse_tissue_padj0.05_up_top100, "mus_tissue_nasopha_vs_liver_stomach_esophagus_up_top100.grp", quote = FALSE, row.names=FALSE)

mouse_tissue_padj0.05_up_top200 <- data.frame(mouse_tissue_padj0.05_up[1:200,1])
colnames(mouse_tissue_padj0.05_up_top200) <- "#mus_tissue_nasopha_vs_liver_stomach_esophagus_up_top200"
write.table(mouse_tissue_padj0.05_up_top200, "mus_tissue_nasopha_vs_liver_stomach_esophagus_up_top200.grp", quote = FALSE, row.names=FALSE)

mouse_tissue_padj0.05_up_top300 <- data.frame(mouse_tissue_padj0.05_up[1:300,1])
colnames(mouse_tissue_padj0.05_up_top300) <- "#mus_tissue_nasopha_vs_liver_stomach_esophagus_up_top300"
write.table(mouse_tissue_padj0.05_up_top300, "mus_tissue_nasopha_vs_liver_stomach_esophagus_up_top300.grp", quote = FALSE, row.names=FALSE)

mouse_tissue_padj0.05_up_top400 <- data.frame(mouse_tissue_padj0.05_up[1:400,1])
colnames(mouse_tissue_padj0.05_up_top400) <- "#mus_tissue_nasopha_vs_liver_stomach_esophagus_up_top400"
write.table(mouse_tissue_padj0.05_up_top400, "mus_tissue_nasopha_vs_liver_stomach_esophagus_up_top400.grp", quote = FALSE, row.names=FALSE)

mouse_tissue_padj0.05_up_top500 <- data.frame(mouse_tissue_padj0.05_up[1:500,1])
colnames(mouse_tissue_padj0.05_up_top500) <- "#mus_tissue_nasopha_vs_liver_stomach_esophagus_up_top500"
write.table(mouse_tissue_padj0.05_up_top500, "mus_tissue_nasopha_vs_liver_stomach_esophagus_up_top500.grp", quote = FALSE, row.names=FALSE)

mouse_tissue_padj0.05_up_top600 <- data.frame(mouse_tissue_padj0.05_up[1:600,1])
colnames(mouse_tissue_padj0.05_up_top600) <- "#mus_tissue_nasopha_vs_liver_stomach_esophagus_up_top600"
write.table(mouse_tissue_padj0.05_up_top600, "mus_tissue_nasopha_vs_liver_stomach_esophagus_up_top600.grp", quote = FALSE, row.names=FALSE)

mouse_tissue_padj0.05_up_top700 <- data.frame(mouse_tissue_padj0.05_up[1:700,1])
colnames(mouse_tissue_padj0.05_up_top700) <- "#mus_tissue_nasopha_vs_liver_stomach_esophagus_up_top700"
write.table(mouse_tissue_padj0.05_up_top700, "mus_tissue_nasopha_vs_liver_stomach_esophagus_up_top700.grp", quote = FALSE, row.names=FALSE)

mouse_tissue_padj0.05_up_top800 <- data.frame(mouse_tissue_padj0.05_up[1:800,1])
colnames(mouse_tissue_padj0.05_up_top800) <- "#mus_tissue_nasopha_vs_liver_stomach_esophagus_up_top800"
write.table(mouse_tissue_padj0.05_up_top800, "mus_tissue_nasopha_vs_liver_stomach_esophagus_up_top800.grp", quote = FALSE, row.names=FALSE)

mouse_tissue_padj0.05_up_top900 <- data.frame(mouse_tissue_padj0.05_up[1:900,1])
colnames(mouse_tissue_padj0.05_up_top900) <- "#mus_tissue_nasopha_vs_liver_stomach_esophagus_up_top900"
write.table(mouse_tissue_padj0.05_up_top900, "mus_tissue_nasopha_vs_liver_stomach_esophagus_up_top900.grp", quote = FALSE, row.names=FALSE)

mouse_tissue_padj0.05_up_top1000 <- data.frame(mouse_tissue_padj0.05_up[1:1000,1])
colnames(mouse_tissue_padj0.05_up_top1000) <- "#mus_tissue_nasopha_vs_liver_stomach_esophagus_up_top1000"
write.table(mouse_tissue_padj0.05_up_top1000, "mus_tissue_nasopha_vs_liver_stomach_esophagus_up_top1000.grp", quote = FALSE, row.names=FALSE)


mouse_tissue_padj0.05_dn <- mouse_tissue %>% filter(padj <= 0.05) %>% filter(log2FoldChange <= -1)
mouse_tissue_padj0.05_dn <- mouse_tissue_padj0.05_dn[order(mouse_tissue_padj0.05_dn$log2FoldChange,decreasing = FALSE),]


mouse_tissue_padj0.05_dn_top100 <- data.frame(mouse_tissue_padj0.05_dn[1:100,1])
colnames(mouse_tissue_padj0.05_dn_top100) <- "#mus_tissue_nasopha_vs_liver_stomach_esophagus_dn_top100"
write.table(mouse_tissue_padj0.05_dn_top100, "mus_tissue_nasopha_vs_liver_stomach_esophagus_dn_top100.grp", quote = FALSE, row.names=FALSE)

mouse_tissue_padj0.05_dn_top200 <- data.frame(mouse_tissue_padj0.05_dn[1:200,1])
colnames(mouse_tissue_padj0.05_dn_top200) <- "#mus_tissue_nasopha_vs_liver_stomach_esophagus_dn_top200"
write.table(mouse_tissue_padj0.05_dn_top200, "mus_tissue_nasopha_vs_liver_stomach_esophagus_dn_top200.grp", quote = FALSE, row.names=FALSE)

mouse_tissue_padj0.05_dn_top300 <- data.frame(mouse_tissue_padj0.05_dn[1:300,1])
colnames(mouse_tissue_padj0.05_dn_top300) <- "#mus_tissue_nasopha_vs_liver_stomach_esophagus_dn_top300"
write.table(mouse_tissue_padj0.05_dn_top300, "mus_tissue_nasopha_vs_liver_stomach_esophagus_dn_top300.grp", quote = FALSE, row.names=FALSE)

mouse_tissue_padj0.05_dn_top400 <- data.frame(mouse_tissue_padj0.05_dn[1:400,1])
colnames(mouse_tissue_padj0.05_dn_top400) <- "#mus_tissue_nasopha_vs_liver_stomach_esophagus_dn_top400"
write.table(mouse_tissue_padj0.05_dn_top400, "mus_tissue_nasopha_vs_liver_stomach_esophagus_dn_top400.grp", quote = FALSE, row.names=FALSE)

mouse_tissue_padj0.05_dn_top500 <- data.frame(mouse_tissue_padj0.05_dn[1:500,1])
colnames(mouse_tissue_padj0.05_dn_top500) <- "#mus_tissue_nasopha_vs_liver_stomach_esophagus_dn_top500"
write.table(mouse_tissue_padj0.05_dn_top500, "mus_tissue_nasopha_vs_liver_stomach_esophagus_dn_top500.grp", quote = FALSE, row.names=FALSE)

mouse_tissue_padj0.05_dn_top600 <- data.frame(mouse_tissue_padj0.05_dn[1:600,1])
colnames(mouse_tissue_padj0.05_dn_top600) <- "#mus_tissue_nasopha_vs_liver_stomach_esophagus_dn_top600"
write.table(mouse_tissue_padj0.05_dn_top600, "mus_tissue_nasopha_vs_liver_stomach_esophagus_dn_top600.grp", quote = FALSE, row.names=FALSE)

mouse_tissue_padj0.05_dn_top700 <- data.frame(mouse_tissue_padj0.05_dn[1:700,1])
colnames(mouse_tissue_padj0.05_dn_top700) <- "#mus_tissue_nasopha_vs_liver_stomach_esophagus_dn_top700"
write.table(mouse_tissue_padj0.05_dn_top700, "mus_tissue_nasopha_vs_liver_stomach_esophagus_dn_top700.grp", quote = FALSE, row.names=FALSE)

mouse_tissue_padj0.05_dn_top800 <- data.frame(mouse_tissue_padj0.05_dn[1:800,1])
colnames(mouse_tissue_padj0.05_dn_top800) <- "#mus_tissue_nasopha_vs_liver_stomach_esophagus_dn_top800"
write.table(mouse_tissue_padj0.05_dn_top800, "mus_tissue_nasopha_vs_liver_stomach_esophagus_dn_top800.grp", quote = FALSE, row.names=FALSE)

mouse_tissue_padj0.05_dn_top900 <- data.frame(mouse_tissue_padj0.05_dn[1:900,1])
colnames(mouse_tissue_padj0.05_dn_top900) <- "#mus_tissue_nasopha_vs_liver_stomach_esophagus_dn_top900"
write.table(mouse_tissue_padj0.05_dn_top900, "mus_tissue_nasopha_vs_liver_stomach_esophagus_dn_top900.grp", quote = FALSE, row.names=FALSE)

mouse_tissue_padj0.05_dn_top1000 <- data.frame(mouse_tissue_padj0.05_dn[1:1000,1])
colnames(mouse_tissue_padj0.05_dn_top1000) <- "#mus_tissue_nasopha_vs_liver_stomach_esophagus_dn_top1000"
write.table(mouse_tissue_padj0.05_dn_top1000, "mus_tissue_nasopha_vs_liver_stomach_esophagus_dn_top1000.grp", quote = FALSE, row.names=FALSE)


mouse_tissue_padj0.05_up_top50 <- data.frame(mouse_tissue_padj0.05_up[1:50,1])
colnames(mouse_tissue_padj0.05_up_top50) <- "mouse_tissue_UP_50    "
def <- data.frame(c("> Genes in the mouse_tissue up."))
colnames(def) <- "mouse_tissue_UP_50    "
mouse_tissue_padj0.05_up_top50 <- rbind(def,mouse_tissue_padj0.05_up_top50)


mouse_tissue_padj0.05_dn_top50 <- data.frame(mouse_tissue_padj0.05_dn[1:50,1])
colnames(mouse_tissue_padj0.05_dn_top50) <- "mouse_tissue_DN_50    "
def <- data.frame(c("> Genes in the mouse_tissue dn."))
colnames(def) <- "mouse_tissue_DN_50    "
mouse_tissue_padj0.05_dn_top50 <- rbind(def,mouse_tissue_padj0.05_dn_top50)


mouse_tissue_padj0.05_top50 <- t(cbind(mouse_tissue_padj0.05_up_top50,mouse_tissue_padj0.05_dn_top50))
write.table(mouse_tissue_padj0.05_top50, "mouse_tissue_nasopha_vs_liver_stomach_esophagus_top50.gmt", quote = FALSE, col.names=FALSE)



mmID_human_TCGA <- read.csv(file = "/mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/nasopha_vs_others/nasopha_vs_liver_stomach_esophagus/mmID_human_TCGA_nasopha_vs_liver_stomach_esophagus_count_normalized_DESeq.csv")

mmID_human_TCGA <- mmID_human_TCGA[,c(1,195,198,199)]
mmID_human_TCGA_padj0.05_up <- mmID_human_TCGA %>% filter(padj <= 0.05) %>% filter(log2FoldChange >= 1)
mmID_human_TCGA_padj0.05_up <- mmID_human_TCGA_padj0.05_up[order(mmID_human_TCGA_padj0.05_up$log2FoldChange,decreasing = TRUE),]

mmID_human_TCGA_padj0.05_up_top100 <- data.frame(mmID_human_TCGA_padj0.05_up[1:100,1])
colnames(mmID_human_TCGA_padj0.05_up_top100) <- "#mmID_human_TCGA_nasopha_vs_liver_stomach_esophagus_up_top100"
write.table(mmID_human_TCGA_padj0.05_up_top100, "mmID_human_TCGA_nasopha_vs_liver_stomach_esophagus_up_top100.grp", quote = FALSE, row.names=FALSE)

mmID_human_TCGA_padj0.05_up_top200 <- data.frame(mmID_human_TCGA_padj0.05_up[1:200,1])
colnames(mmID_human_TCGA_padj0.05_up_top200) <- "#mmID_human_TCGA_nasopha_vs_liver_stomach_esophagus_up_top200"
write.table(mmID_human_TCGA_padj0.05_up_top200, "mmID_human_TCGA_nasopha_vs_liver_stomach_esophagus_up_top200.grp", quote = FALSE, row.names=FALSE)

mmID_human_TCGA_padj0.05_up_top300 <- data.frame(mmID_human_TCGA_padj0.05_up[1:300,1])
colnames(mmID_human_TCGA_padj0.05_up_top300) <- "#mmID_human_TCGA_nasopha_vs_liver_stomach_esophagus_up_top300"
write.table(mmID_human_TCGA_padj0.05_up_top300, "mmID_human_TCGA_nasopha_vs_liver_stomach_esophagus_up_top300.grp", quote = FALSE, row.names=FALSE)

mmID_human_TCGA_padj0.05_up_top400 <- data.frame(mmID_human_TCGA_padj0.05_up[1:400,1])
colnames(mmID_human_TCGA_padj0.05_up_top400) <- "#mmID_human_TCGA_nasopha_vs_liver_stomach_esophagus_up_top400"
write.table(mmID_human_TCGA_padj0.05_up_top400, "mmID_human_TCGA_nasopha_vs_liver_stomach_esophagus_up_top400.grp", quote = FALSE, row.names=FALSE)

mmID_human_TCGA_padj0.05_up_top500 <- data.frame(mmID_human_TCGA_padj0.05_up[1:500,1])
colnames(mmID_human_TCGA_padj0.05_up_top500) <- "#mmID_human_TCGA_nasopha_vs_liver_stomach_esophagus_up_top500"
write.table(mmID_human_TCGA_padj0.05_up_top500, "mmID_human_TCGA_nasopha_vs_liver_stomach_esophagus_up_top500.grp", quote = FALSE, row.names=FALSE)

mmID_human_TCGA_padj0.05_up_top600 <- data.frame(mmID_human_TCGA_padj0.05_up[1:600,1])
colnames(mmID_human_TCGA_padj0.05_up_top600) <- "#mmID_human_TCGA_nasopha_vs_liver_stomach_esophagus_up_top600"
write.table(mmID_human_TCGA_padj0.05_up_top600, "mmID_human_TCGA_nasopha_vs_liver_stomach_esophagus_up_top600.grp", quote = FALSE, row.names=FALSE)

mmID_human_TCGA_padj0.05_up_top700 <- data.frame(mmID_human_TCGA_padj0.05_up[1:700,1])
colnames(mmID_human_TCGA_padj0.05_up_top700) <- "#mmID_human_TCGA_nasopha_vs_liver_stomach_esophagus_up_top700"
write.table(mmID_human_TCGA_padj0.05_up_top700, "mmID_human_TCGA_nasopha_vs_liver_stomach_esophagus_up_top700.grp", quote = FALSE, row.names=FALSE)

mmID_human_TCGA_padj0.05_up_top800 <- data.frame(mmID_human_TCGA_padj0.05_up[1:800,1])
colnames(mmID_human_TCGA_padj0.05_up_top800) <- "#mmID_human_TCGA_nasopha_vs_liver_stomach_esophagus_up_top800"
write.table(mmID_human_TCGA_padj0.05_up_top800, "mmID_human_TCGA_nasopha_vs_liver_stomach_esophagus_up_top800.grp", quote = FALSE, row.names=FALSE)

mmID_human_TCGA_padj0.05_up_top900 <- data.frame(mmID_human_TCGA_padj0.05_up[1:900,1])
colnames(mmID_human_TCGA_padj0.05_up_top900) <- "#mmID_human_TCGA_nasopha_vs_liver_stomach_esophagus_up_top900"
write.table(mmID_human_TCGA_padj0.05_up_top900, "mmID_human_TCGA_nasopha_vs_liver_stomach_esophagus_up_top900.grp", quote = FALSE, row.names=FALSE)

mmID_human_TCGA_padj0.05_up_top1000 <- data.frame(mmID_human_TCGA_padj0.05_up[1:1000,1])
colnames(mmID_human_TCGA_padj0.05_up_top1000) <- "#mmID_human_TCGA_nasopha_vs_liver_stomach_esophagus_up_top1000"
write.table(mmID_human_TCGA_padj0.05_up_top1000, "mmID_human_TCGA_nasopha_vs_liver_stomach_esophagus_up_top1000.grp", quote = FALSE, row.names=FALSE)



mmID_human_TCGA_padj0.05_dn <- mmID_human_TCGA %>% filter(padj <= 0.05) %>% filter(log2FoldChange <= -1)
mmID_human_TCGA_padj0.05_dn <- mmID_human_TCGA_padj0.05_dn[order(mmID_human_TCGA_padj0.05_dn$log2FoldChange,decreasing = FALSE),]


mmID_human_TCGA_padj0.05_dn_top100 <- data.frame(mmID_human_TCGA_padj0.05_dn[1:100,1])
colnames(mmID_human_TCGA_padj0.05_dn_top100) <- "#mmID_human_TCGA_nasopha_vs_liver_stomach_esophagus_dn_top100"
write.table(mmID_human_TCGA_padj0.05_dn_top100, "mmID_human_TCGA_nasopha_vs_liver_stomach_esophagus_dn_top100.grp", quote = FALSE, row.names=FALSE)

mmID_human_TCGA_padj0.05_dn_top200 <- data.frame(mmID_human_TCGA_padj0.05_dn[1:200,1])
colnames(mmID_human_TCGA_padj0.05_dn_top200) <- "#mmID_human_TCGA_nasopha_vs_liver_stomach_esophagus_dn_top200"
write.table(mmID_human_TCGA_padj0.05_dn_top200, "mmID_human_TCGA_nasopha_vs_liver_stomach_esophagus_dn_top200.grp", quote = FALSE, row.names=FALSE)

mmID_human_TCGA_padj0.05_dn_top300 <- data.frame(mmID_human_TCGA_padj0.05_dn[1:300,1])
colnames(mmID_human_TCGA_padj0.05_dn_top300) <- "#mmID_human_TCGA_nasopha_vs_liver_stomach_esophagus_dn_top300"
write.table(mmID_human_TCGA_padj0.05_dn_top300, "mmID_human_TCGA_nasopha_vs_liver_stomach_esophagus_dn_top300.grp", quote = FALSE, row.names=FALSE)

mmID_human_TCGA_padj0.05_dn_top400 <- data.frame(mmID_human_TCGA_padj0.05_dn[1:400,1])
colnames(mmID_human_TCGA_padj0.05_dn_top400) <- "#mmID_human_TCGA_nasopha_vs_liver_stomach_esophagus_dn_top400"
write.table(mmID_human_TCGA_padj0.05_dn_top400, "mmID_human_TCGA_nasopha_vs_liver_stomach_esophagus_dn_top400.grp", quote = FALSE, row.names=FALSE)

mmID_human_TCGA_padj0.05_dn_top500 <- data.frame(mmID_human_TCGA_padj0.05_dn[1:500,1])
colnames(mmID_human_TCGA_padj0.05_dn_top500) <- "#mmID_human_TCGA_nasopha_vs_liver_stomach_esophagus_dn_top500"
write.table(mmID_human_TCGA_padj0.05_dn_top500, "mmID_human_TCGA_nasopha_vs_liver_stomach_esophagus_dn_top500.grp", quote = FALSE, row.names=FALSE)

mmID_human_TCGA_padj0.05_dn_top600 <- data.frame(mmID_human_TCGA_padj0.05_dn[1:600,1])
colnames(mmID_human_TCGA_padj0.05_dn_top600) <- "#mmID_human_TCGA_nasopha_vs_liver_stomach_esophagus_dn_top600"
write.table(mmID_human_TCGA_padj0.05_dn_top600, "mmID_human_TCGA_nasopha_vs_liver_stomach_esophagus_dn_top600.grp", quote = FALSE, row.names=FALSE)

mmID_human_TCGA_padj0.05_dn_top700 <- data.frame(mmID_human_TCGA_padj0.05_dn[1:700,1])
colnames(mmID_human_TCGA_padj0.05_dn_top700) <- "#mmID_human_TCGA_nasopha_vs_liver_stomach_esophagus_dn_top700"
write.table(mmID_human_TCGA_padj0.05_dn_top700, "mmID_human_TCGA_nasopha_vs_liver_stomach_esophagus_dn_top700.grp", quote = FALSE, row.names=FALSE)

mmID_human_TCGA_padj0.05_dn_top800 <- data.frame(mmID_human_TCGA_padj0.05_dn[1:800,1])
colnames(mmID_human_TCGA_padj0.05_dn_top800) <- "#mmID_human_TCGA_nasopha_vs_liver_stomach_esophagus_dn_top800"
write.table(mmID_human_TCGA_padj0.05_dn_top800, "mmID_human_TCGA_nasopha_vs_liver_stomach_esophagus_dn_top800.grp", quote = FALSE, row.names=FALSE)

mmID_human_TCGA_padj0.05_dn_top900 <- data.frame(mmID_human_TCGA_padj0.05_dn[1:900,1])
colnames(mmID_human_TCGA_padj0.05_dn_top900) <- "#mmID_human_TCGA_nasopha_vs_liver_stomach_esophagus_dn_top900"
write.table(mmID_human_TCGA_padj0.05_dn_top900, "mmID_human_TCGA_nasopha_vs_liver_stomach_esophagus_dn_top900.grp", quote = FALSE, row.names=FALSE)

mmID_human_TCGA_padj0.05_dn_top1000 <- data.frame(mmID_human_TCGA_padj0.05_dn[1:1000,1])
colnames(mmID_human_TCGA_padj0.05_dn_top1000) <- "#mmID_human_TCGA_nasopha_vs_liver_stomach_esophagus_dn_top1000"
write.table(mmID_human_TCGA_padj0.05_dn_top1000, "mmID_human_TCGA_nasopha_vs_liver_stomach_esophagus_dn_top1000.grp", quote = FALSE, row.names=FALSE)



mmID_human_TCGA_padj0.05_up_top50 <- data.frame(mmID_human_TCGA_padj0.05_up[1:50,1])
colnames(mmID_human_TCGA_padj0.05_up_top50) <- "human_TCGA_UP_50    "
def <- data.frame(c("> Genes in the TCGA up."))
colnames(def) <- "human_TCGA_UP_50    "
mmID_human_TCGA_padj0.05_up_top50 <- rbind(def,mmID_human_TCGA_padj0.05_up_top50)


mmID_human_TCGA_padj0.05_dn_top50 <- data.frame(mmID_human_TCGA_padj0.05_dn[1:50,1])
colnames(mmID_human_TCGA_padj0.05_dn_top50) <- "human_TCGA_DN_50    "
def <- data.frame(c("> Genes in the TCGA dn."))
colnames(def) <- "human_TCGA_DN_50    "
mmID_human_TCGA_padj0.05_dn_top50 <- rbind(def,mmID_human_TCGA_padj0.05_dn_top50)


mmID_human_TCGA_padj0.05_top50 <- t(cbind(mmID_human_TCGA_padj0.05_up_top50,mmID_human_TCGA_padj0.05_dn_top50))
write.table(mmID_human_TCGA_padj0.05_top50, "mmID_human_TCGA_nasopha_vs_liver_stomach_esophagus_top50.gmt", quote = FALSE, col.names=FALSE)




************************************************************************************鼻咽组织marker

organoids_nasopha_vs_liver_stomach_esophagus <- read.csv("/mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/nasopha_vs_others/nasopha_vs_liver_stomach_esophagus/organoids_nasopha_vs_liver_stomach_esophagus_count_normalized_DESeq.csv",row.names=1)
organoids_nasopha_vs_liver_stomach_esophagus <- organoids_nasopha_vs_liver_stomach_esophagus[,c(13:24,26,29:30)]


organoids_nasopha_vs_liver_stomach_esophagus_pvalue0.05_log2FC1_up <- organoids_nasopha_vs_liver_stomach_esophagus %>% filter(padj <= 0.05) %>% filter(log2FoldChange >= 1)
nrow(organoids_nasopha_vs_liver_stomach_esophagus_pvalue0.05_log2FC1_up)

organoids_nasopha_vs_liver_stomach_esophagus_pvalue0.05_log2FC1_up <- organoids_nasopha_vs_liver_stomach_esophagus_pvalue0.05_log2FC1_up[order(organoids_nasopha_vs_liver_stomach_esophagus_pvalue0.05_log2FC1_up$log2FoldChange,decreasing = TRUE),]
head(organoids_nasopha_vs_liver_stomach_esophagus_pvalue0.05_log2FC1_up)
tail(organoids_nasopha_vs_liver_stomach_esophagus_pvalue0.05_log2FC1_up)
organoids_nasopha_vs_liver_stomach_esophagus_pvalue0.05_log2FC1_up_heatmap <- organoids_nasopha_vs_liver_stomach_esophagus_pvalue0.05_log2FC1_up[,1:12]
head(organoids_nasopha_vs_liver_stomach_esophagus_pvalue0.05_log2FC1_up_heatmap)

write.csv(organoids_nasopha_vs_liver_stomach_esophagus_pvalue0.05_log2FC1_up_heatmap, "/mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/nasopha_vs_others/nasopha_vs_liver_stomach_esophagus/organoids_nasopha_vs_liver_stomach_esophagus_pvalue0.05_log2FC1_up_heatmap.csv")


all_gene <- organoids_nasopha_vs_liver_stomach_esophagus_pvalue0.05_log2FC1_up_heatmap
head(all_gene)

all_gene <- na.omit(all_gene)
#all_gene_log2 <- log(all_gene+1,2)

zscore <- t(apply(all_gene, 1, function(x) (x-mean(x))/sd(x)))
zscore <- na.omit(zscore)
head(zscore)

library(Seurat)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
require(RColorBrewer)
source("/mnt/data/user_data/xiangyu/programme/R_PACKAGES/my_code/MyBestFunction_scRNA.R")
source("/mnt/data/user_data/xiangyu/programme/R_PACKAGES/my_code/Pseudo_CNV_series.R")

#range(zscore)
#zscore[zscore > 2] <- 2
#zscore[zscore < -2] <- -2

pheatmap(zscore,color = colorRampPalette(c("navy", "white", "firebrick3"))(50),show_rownames=FALSE,cluster_row = FALSE,cluster_col= FALSE,border=FALSE)


SeuratObject <- CreateSeuratObject(counts = zscore, project = "WXD")
SeuratObject$orig.ident <-factor(SeuratObject$orig.ident,levels = c("liver","stomach","esophag","nasopha"))   
levels(SeuratObject$orig.ident) 
gene <- rownames(zscore)
XY_heatmap(seurat_obj=SeuratObject,group="orig.ident", genes=gene,all_num=FALSE,new_names=NULL,labels_rot=90,assay_sel="RNA",color=colorRampPalette(brewer.pal(10, "RdBu"))(101),min_and_max_cut=1.2,show_row_names=FALSE,scale=FALSE,label_size=0,mark_gene=c("Krt5","Krt13","Krt14","Krt15","Krt16","Krt17","Tuba1c","Cfap45","Nes","Muc15","Muc16","Trp63"))

pdf("/mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/nasopha_vs_others/nasopha_vs_liver_stomach_esophagus/nasopha_marker_heatmap_esophag.pdf")
XY_heatmap(seurat_obj=SeuratObject,group="orig.ident", genes=gene,all_num=FALSE,new_names=NULL,labels_rot=90,assay_sel="RNA",color=colorRampPalette(brewer.pal(10, "RdBu"))(101),min_and_max_cut=1.2,show_row_names=FALSE,scale=FALSE,label_size=0,mark_gene=c("Krt5","Krt13","Krt14","Krt15","Krt16","Krt17","Tuba1c","Cfap45","Nes","Muc15","Muc16","Trp63"))
dev.off()
