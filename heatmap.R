*********************************************************************************************************************************
**************************************************************heatmap************************************************************
*********************************************************************************************************************************
###############nasopharynx marker
all_gene <- read.csv("/mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/nasopha_vs_others/nasopha_vs_three_organoid_heatmap_data_up_log2fFC1.csv")

head(all_gene)
all_gene <- all_gene[!duplicated(all_gene$X),]
all_gene <- na.omit(all_gene)
rownames(all_gene) <- all_gene[,1]
all_gene <- all_gene[,-1]
colnames(all_gene)
all_gene <- all_gene[,c(1:3,10:12)]

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

SeuratObject <- CreateSeuratObject(counts = zscore, project = "NPC")
gene <- rownames(zscore)

pdf("nasopha_marker_compared_lung_up_log2fFC1_heatmap.pdf")
XY_heatmap(seurat_obj=SeuratObject,group="orig.ident", genes=gene,all_num=FALSE,new_names=NULL,labels_rot=90,assay_sel="RNA",color=colorRampPalette(brewer.pal(10, "RdBu"))(101),min_and_max_cut=1.2,show_row_names=FALSE,scale=FALSE,label_size=0,mark_gene=c("Krt5","Krt13","Krt14","Krt15","Krt16","Krt17","Tuba1c","Cfap45","Nes","Muc15","Muc16","Trp63"))
dev.off()

#####################NPC marker
all_gene <- read.csv("/mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/P_ctrl_vs_nomal_organoid/Primary_TMP_v_normal_organoid_heatmap_data.csv")
head(all_gene)
all_gene <- all_gene[!duplicated(all_gene$X),]
all_gene <- na.omit(all_gene)
rownames(all_gene) <- all_gene[,1]
all_gene <- all_gene[,-1]

colnames(all_gene)
all_gene <- all_gene[,-c(7:9)]
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
SeuratObject$orig.ident <-factor(SeuratObject$orig.ident,levels = c("normal","ctrl"))
gene <- rownames(zscore)

pdf("NPC_marker_heatmap.pdf")
XY_heatmap(seurat_obj=SeuratObject,group="orig.ident",genes=gene,all_num=FALSE,new_names=NULL,labels_rot=90,assay_sel="RNA",color=colorRampPalette(brewer.pal(10, "RdBu"))(101),min_and_max_cut=1.2,show_row_names=FALSE,scale=FALSE,label_size=0,mark_gene=c("Bmi1","Tubb2a","Tubb2b","Muc5ac","Krt5","Krt13","Krt14","Krt16","Notch1","Egr2","Cd44","Cav1","Cd109"))
dev.off()
