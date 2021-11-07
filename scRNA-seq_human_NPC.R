smart-seq2 样本信息 https://bigd.big.ac.cn/gsa-human/browse/HRA000087

ls *_r1.fq.gz |while read id;
do  
STAR  --outSAMtype BAM SortedByCoordinate \
--runThreadN 20 \
--readFilesCommand zcat \
--readFilesIn /mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/GSA_scNPC/fq_file/$id \
/mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/GSA_scNPC/fq_file/${id%_*}_r2.fq.gz \
--genomeDir /mnt/data/public_data/genome_index/genome_index/Human_GRCh38_readlength150_index/ \
--outFileNamePrefix /mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/GSA_scNPC/star_out/${id%_*}.
done


vi sample_sampletable.csv

ids,sample,deal,order
1,HRR026852.Aligned.sortedByCoord.out.bam,Pri_1,primary,1
2,HRR026854.Aligned.sortedByCoord.out.bam,Pri_2,primary,2
3,HRR026856.Aligned.sortedByCoord.out.bam,Pri_3,primary,3
4,HRR026851.Aligned.sortedByCoord.out.bam,Met_1,metastasis,4
5,HRR026853.Aligned.sortedByCoord.out.bam,Met_2,metastasis,5
6,HRR026855.Aligned.sortedByCoord.out.bam,Met_3,metastasis,6

vi sample_sampletable_without_P2.csv

ids,sample,deal,order
1,HRR026852.Aligned.sortedByCoord.out.bam,Pri_1,primary,1
2,HRR026856.Aligned.sortedByCoord.out.bam,Pri_3,primary,2
3,HRR026851.Aligned.sortedByCoord.out.bam,Met_1,metastasis,3
4,HRR026853.Aligned.sortedByCoord.out.bam,Met_2,metastasis,4
5,HRR026855.Aligned.sortedByCoord.out.bam,Met_3,metastasis,5

vi sample_sampletable_without_P3.csv

ids,sample,deal,order
1,HRR026852.Aligned.sortedByCoord.out.bam,Pri_1,primary,1
2,HRR026854.Aligned.sortedByCoord.out.bam,Pri_2,primary,2
3,HRR026851.Aligned.sortedByCoord.out.bam,Met_1,metastasis,3
4,HRR026853.Aligned.sortedByCoord.out.bam,Met_2,metastasis,4
5,HRR026855.Aligned.sortedByCoord.out.bam,Met_3,metastasis,5


project <- c("scNPC")
deal_design <- c("metastasis","primary")
significant_cutoff <- c(1)
organism <- "human"

file_path <- "/mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/GSA_scNPC/work_file/"
sample_sampletable.path <- "/mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/GSA_scNPC/star_out"

deal_case <- deal_design[1]
dela_control <- deal_design[2]
deal_names <- paste(deal_case,dela_control,sep="_v_")
project_pathway <- paste(".",project,sep="/")
se.save_RData <- paste(project_pathway,"00_RNAseq_se.RData",sep="_")
cluster_sample_heatmap_rld <- paste(project_pathway,"01_cluster_similarity_rld_pheatmap.pdf",sep="_")
cluster_sample_heatmap_vsd <- paste(project_pathway,"01_cluster_similarity_vsd_pheatmap.pdf",sep="_")
cluster_sample_pca <- paste(project_pathway,"01_cluster_similarity_pca.pdf",sep="_")
tpm_csv <- paste(project_pathway,"02_tpm.csv",sep="_")
res_1_file.csv <- paste(project_pathway,"03",deal_names,"result.csv",sep="_")
res_1_filter <- c("deal",deal_design)
colnames_1 <- deal_names
res_1_file_sym_name <- paste(project_pathway,"04",deal_names,"symbol_and_anno.csv",sep="_")
ALL.CSV_FILE <- paste(project_pathway,"05",deal_names,"count_tpm_symbol_and_anno.csv",sep="_")
KEGGupres_1_file.csv <- paste(project_pathway,"06",deal_names,"KEGG_UP.csv",sep="_")
KEGGdownres_1_file.csv <- paste(project_pathway,"07",deal_names,"KEGG_DOWN.csv",sep="_")
KEGGupres_1_pdf <- paste(project_pathway,"08",deal_names,"KEGG_UP.pdf",sep="_")
KEGGdownres_1_pdf <- paste(project_pathway,"09",deal_names,"KEGG_DOWN.pdf",sep="_")
KEGGres_1_all_pdf <- paste(project_pathway,"10",deal_names,"KEGG_all.pdf",sep="_")
KEGGres_1_all_title <- paste(deal_names,"KEGG Pathway Enrichment",sep=" ")
GOres_1_all_UP_csv <- paste(project_pathway,"11",deal_names,"GO_UP.csv",sep="_")
GOres_1_all_DOWN_csv <- paste(project_pathway,"12",deal_names,"GO_DOWN.csv",sep="_")
GOres_1_all_UP_pdf <- paste(project_pathway,"13",deal_names,"GO_UP.pdf",sep="_")
GOres_1_all_DOWN_pdf <- paste(project_pathway,"14",deal_names,"GO_DOWN.pdf",sep="_")

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


if (organism %in% c("mouse")) {
	anno_data=org.Mm.eg.db
	print("organism is mouse")
	load(file="/mnt/data/user_data/xiangyu/workshop/WORKFLOW_RNAseq/ebg_mm10.RData")
	load(file="/mnt/data/user_data/xiangyu/workshop/WORKFLOW_RNAseq/txdb_mm10.RData")
} else {
	anno_data=org.Hs.eg.db
	print("organism is human")
	load(file="/mnt/data/user_data/xiangyu/workshop/WORKFLOW_RNAseq/ebg_GRCh38.RData")
	load(file="/mnt/data/user_data/xiangyu/workshop/WORKFLOW_RNAseq/txdb_GRCh38.RData")
}


indir <- file.path(sample_sampletable.path)
csvfile <- file.path(indir, "sample_sampletable.csv")
csvfile <- file.path(indir, "sample_sampletable_without_P2.csv")
csvfile <- file.path(indir, "sample_sampletable_without_P3.csv")

sampleTable <- read.csv(csvfile, row.names = 1)
sampleTable
filenames <- file.path(indir, sampleTable$ids)
file.exists(filenames)

bamfiles <- BamFileList(filenames, yieldSize=2000000)
seqinfo(bamfiles[1])
#需要调用BiocParallel，windows用户可能报错，可以先运行一句：register(SerialParam())
register(SerialParam())
#常用的计数软件包括HTseq-count, featureCounts,以及summarizeOverlaps等。后两者是bioconductor中的软件包.
#summarizeOverlaps其实是HTseq-count中R包装而已。在计算速度上可以和HTseq-count相较，
#它的优点在于统计结果可以直接封装在RangedSummarizedExperiment对象中，下游的DESeq2对它的支持非常好。
se <- summarizeOverlaps(features=ebg, reads=bamfiles,
mode="Union",   # describes what kind of read overlaps will be counted.读段覆盖的模式
singleEnd=FALSE,  ##表示是pair-end
ignore.strand=TRUE, ##是否是链特异性,True 表示忽略±链的限制
fragments=FALSE)  ##只应用于双末端测序,true表示非成对的对端应该被计数
setwd(file_path)

se <- summarizeOverlaps(features=ebg, reads=bamfiles,
mode="Union",   
singleEnd=FALSE,  
ignore.strand=TRUE,
fragments=FALSE)  
setwd(file_path)
save(se,file=se.save_RData)
save(se,file="./scNPC_00_RNAseq_se_kickP2.RData")
save(se,file="./scNPC_00_RNAseq_se_kickP3.RData")

colData(se) <- DataFrame(sampleTable)

dds <- DESeqDataSet(se, design = ~deal)
countdata <- assay(se)
coldata <- colData(se)

rld <- rlog(dds, blind = FALSE)
sampleDists <- dist(t(assay(rld)))
samplePoisDistMatrix <- as.matrix(sampleDists)

rownames(samplePoisDistMatrix) <- paste(rld$sample)
colnames(samplePoisDistMatrix) <- NULL
pdf(file=cluster_sample_heatmap_rld)
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(samplePoisDistMatrix,
clustering_distance_rows = sampleDists,
clustering_distance_cols = sampleDists,
col = colors)
dev.off()

vsd <- vst(dds, blind = FALSE)
sampleDists <- dist(t(assay(vsd)))
samplePoisDistMatrix <- as.matrix(sampleDists)
rownames(samplePoisDistMatrix) <- paste(vsd$sample)
colnames(samplePoisDistMatrix) <- NULL
pdf(file=cluster_sample_heatmap_vsd)
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(samplePoisDistMatrix,
clustering_distance_rows = sampleDists,
clustering_distance_cols = sampleDists,
col = colors)
dev.off()

pdf(file=cluster_sample_pca)
pdf(file="./scNPC_01_cluster_similarity_pca_kickP2.pdf")
pdf(file="./scNPC_01_cluster_similarity_pca_kickP3.pdf")

library("ggplot2")
pcaData <- plotPCA(rld, intgroup = c("sample","sample"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(x = PC1, y = PC2, color = sample)) +
geom_point(size =1) +
geom_text(aes(label=sample))+
xlab(paste0("PC1: ", percentVar[1], "% variance")) +
ylab(paste0("PC2: ", percentVar[2], "% variance")) +
coord_fixed()
dev.off()


exons<-as.data.frame(ebg)
exon_length<-data.frame(gene_name=c(0),length=c(0))
for(i in 1:nrow(assay(se))){
  exonsfromgene<-exons[which(exons$group==i),]
  length<-sum(exonsfromgene$width)
  gene_name<-exonsfromgene[1,2]
  exon_length<-rbind(exon_length,c(gene_name,length))
}
exon_length<-exon_length[-1,]
whole_length<-as.data.frame(round(colSums(assay(dds))/1e6,1))
countM<-as.data.frame(assay(se))
FPKM<-countM
for(i in 1:nrow(FPKM)){
  gene<-rownames(FPKM[i,])
  length<-as.numeric(exon_length[which(exon_length$gene_name==gene),]$length)
  for(j in 1:ncol(FPKM)){
    FPKM[i,j]<-FPKM[i,j]/(whole_length[j,]*length*0.001) 
  }
}
fpkmToTpm <- function(FPKM){
    exp(log(FPKM) - log(sum(FPKM)) + log(1e6))
}
tpm <- fpkmToTpm(FPKM)
count <- assay(dds)
colnames(tpm) <- sampleTable$sample
colnames(count) <- sampleTable$sample
count_and_tpm <- cbind(count,tpm)
#cbind： 根据列进行合并，即叠加所有列，m列的矩阵与n列的矩阵cbind()最后变成m+n列，合并前提：cbind(a, b)中矩阵a、b的行数必需相符
#rbind： 根据行进行合并，就是行的叠加，m行的矩阵与n行的矩阵rbind()最后变成m+n行，合并前提：rbind(a, b)中矩阵a、b的列数必需相符
write.csv(tpm,file=tpm_csv)

dds <- DESeq(dds)
res_1 <- results(dds, contrast=res_1_filter)
colnames(res_1) = paste(colnames_1,colnames(res_1),sep="_")
write.csv(res_1, file = res_1_file.csv)

res_1$symbol <- mapIds(x = anno_data,
                        keys = rownames(res_1),
                        keytype ="ENSEMBL",
                        column ="SYMBOL",
                        multiVals="first")
res_1$entrez <- mapIds(x = anno_data,
                        keys = rownames(res_1),
                        keytype ="ENSEMBL",
                        column ="ENTREZID",
                        multiVals="first")
AA <- res_1$symbol
AA <- as.character(AA)
res_1$GENENAME <- mapIds(x = anno_data,
                        keys = AA,
                        keytype ="SYMBOL",
                        column ="GENENAME",
                        multiVals="first")


write.csv(res_1, file = res_1_file_sym_name)
all_summry <- cbind(count_and_tpm,res_1)
write.csv(all_summry, file = ALL.CSV_FILE)


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

all_TF <- all_summry
aa <- as.character(all_TF$symbol) %>% convert_human_to_mouse_symbols()
all_TF$mouse_symbol <- as.character(aa)
matrix <- all_TF[!duplicated(all_TF$mouse_symbol),]
head(matrix)
write.csv(matrix,"scNPC_05_Met_v_Pri_count_tpm_symbol_and_anno_mmsymbol.csv")


*******************************************************************消除LN和nasopharyngeal的器官差异
正常LN @@@@@@此数据存疑
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6148066/
/mnt/data/user_data/xuelan/project/11_WML_bladder/2_metastasis/1_download_data/1_bam

LN_vs_nasopha <- read.csv(file = "/mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/GSA_scNPC/LN_vs_nasopha/LN_vs_nasopha_raw_data_symbol.csv")
LN_vs_nasopha <- LN_vs_nasopha[!duplicated(LN_vs_nasopha$X),]
LN_vs_nasopha <- na.omit(LN_vs_nasopha)
rownames(LN_vs_nasopha) <- LN_vs_nasopha[,1]
LN_vs_nasopha <-  LN_vs_nasopha[,-1]
colnames(LN_vs_nasopha)

id <- colnames(LN_vs_nasopha)
deal <- c(rep('nasopha',3),rep('LN',5))
LN_vs_nasopha_df <- data.frame(id,deal)

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


dds <- DESeqDataSetFromMatrix(countData = LN_vs_nasopha, colData = LN_vs_nasopha_df, design = ~ deal)

countdata <- LN_vs_nasopha
coldata <- LN_vs_nasopha_df
DESeq_counts <- DESeq(dds)
normalized_DEseq <- counts(DESeq_counts, normalized=TRUE)

DEseq_res <- results(DESeq_counts,contrast=c("deal","LN","nasopha")) #最后一个是CTRL
res_1 <- cbind(normalized_DEseq,DEseq_res)


all_summry <- cbind(countdata,res_1)
head(all_summry)
write.csv(all_summry, "LN_vs_nasopha_normal_human_tissue_DESeq.csv")

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


all_TF <- all_summry
aa <- as.character(rownames(all_TF)) %>% convert_human_to_mouse_symbols()
all_TF$mouse_symbol <- as.character(aa)
matrix <- all_TF[!duplicated(all_TF$mouse_symbol),]
head(matrix)
write.csv(matrix,"LN_vs_nasopha_normal_human_tissue_DESeq_mmsymbol.csv")

@@@@@@@@@@@@@@@@@换数据   二者参考转录本不同 合并后就只剩19000多gene

https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4277406/#pone.0115911.s003
   GRCh37

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

nasopha <- read.csv(file = "/mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/ncbi_human_NPC/tumor_vs_normal/humanNPC_tumor_VS_normal_02_tpm.csv")
nasopha <- nasopha[,c(1,5:7)]
colnames(nasopha) <- c("X","nasopha_1","nasopha_2","nasopha_3")

load(file="/mnt/data/user_data/xiangyu/workshop/WORKFLOW_RNAseq/ebg_GRCh38.RData")
load(file="/mnt/data/user_data/xiangyu/workshop/WORKFLOW_RNAseq/txdb_GRCh38.RData")
nasopha$symbol <- mapIds(x = org.Hs.eg.db,
                        keys = as.character(nasopha$X),
						keytype ="ENSEMBL",
						column ="SYMBOL",
						multiVals="first")

lymph_node <- read.csv(file = "/mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/GSA_scNPC/LN_vs_nasopha/E_MTAB_1733_FPKM_allsamples_gooddata.csv",row.names=1)
lymph_node <- lymph_node[,34:38]
colnames(lymph_node) <- c("lymph_node_1","lymph_node_2","lymph_node_3","lymph_node_4","lymph_node_5")
fpkmToTpm <- function(FPKM){
    exp(log(lymph_node) - log(sum(lymph_node)) + log(1e6))
}
lymph_node <- fpkmToTpm(lymph_node)

lymph_node$X<- rownames(lymph_node)
lymph_node <- lymph_node[,c(6,1:5)]
rownames(lymph_node)=NULL

load(file="/mnt/data/user_data/xiangyu/workshop/WORKFLOW_RNAseq/ebg_hg19.RData")
load(file="/mnt/data/user_data/xiangyu/workshop/WORKFLOW_RNAseq/txdb_hg19.RData")
lymph_node$symbol <- mapIds(x = org.Hs.eg.db,
                        keys = as.character(lymph_node$X),
						keytype ="ENSEMBL",
						column ="SYMBOL",
						multiVals="first")


common <- as.data.frame(intersect(x = lymph_node$symbol, y = nasopha$symbol))

lymph_node <- na.omit(lymph_node)
nrow(lymph_node)
lymph_node <- lymph_node[!duplicated(lymph_node$symbol),]
nasopha <- na.omit(nasopha)
nrow(nasopha)
nasopha <- nasopha[!duplicated(nasopha$symbol),]

lymphnode_vs_nasopha_bysymbol <- merge(lymph_node,nasopha,by = "symbol", all=FALSE)
18845
write.csv(lymphnode_vs_nasopha_bysymbol,"/mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/GSA_scNPC/LN_vs_nasopha/lymphnode_vs_nasopha_merge_GRCh37_GRCh38_merge_bysymbol.csv")

lymphnode_vs_nasopha_byENSEMBL <- merge(lymph_node,nasopha,by = "X", all=FALSE)
18828
write.csv(lymphnode_vs_nasopha_byENSEMBL,"/mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/GSA_scNPC/LN_vs_nasopha/lymphnode_vs_nasopha_merge_GRCh37_GRCh38_merge_byENSEMBL.csv")
#两种merge方法，比较发现，会有两个ensemble id 对应同一个symbol id的情况，所以将二者都转换为symbol id再合并较好
#这些ARL17A、BPY2、DNAJC9-AS1、FBXO10、GDF1、HLA-DQA1、IFNAR2、IRS4、NPFF、PALM2AKAP2、POLR2J2、PRH1、PRORP、PRY、SMN1、TRIM74、VCY，在GRCh37、GRCh38中对应不同的ensembl


lymphnode_vs_nasopha_bysymbol <- na.omit(lymphnode_vs_nasopha_bysymbol)
colnames(lymphnode_vs_nasopha_bysymbol) 
lymphnode_vs_nasopha_bysymbol <- lymphnode_vs_nasopha_bysymbol[,c(1,3:7,9:11)]
rownames(lymphnode_vs_nasopha_bysymbol) <- lymphnode_vs_nasopha_bysymbol[,1]
lymphnode_vs_nasopha_bysymbol <- lymphnode_vs_nasopha_bysymbol[,-1]

nasopha_tpm <- c("nasopha_1","nasopha_2","nasopha_3")
lymphnode_tpm <- c("lymph_node_1","lymph_node_2","lymph_node_3","lymph_node_4","lymph_node_5")

lymphnode_tpm <- lymphnode_vs_nasopha_bysymbol[,colnames(lymphnode_vs_nasopha_bysymbol) %in% lymphnode_tpm]
nasopha_tpm <- lymphnode_vs_nasopha_bysymbol[,colnames(lymphnode_vs_nasopha_bysymbol) %in% nasopha_tpm]

logFC <- log2((rowMeans(lymphnode_tpm)+1) / (rowMeans(nasopha_tpm)+1))
library("future.apply")
p_values <- future_lapply(seq(1,nrow(lymphnode_tpm)), function(x){
  res <- wilcox.test(x = t(lymphnode_tpm[x,])[,1], y = t(nasopha_tpm[x,])[,1] )
  res$p.value
})
p <- unlist(p_values)
p.adj <- p.adjust(p, method = "fdr")
genelist <- as.data.frame(logFC)
genelist$p_values <- p
genelist$p_adj <- p.adj
head(genelist)
all_summry <- cbind(lymphnode_vs_nasopha_bysymbol,genelist)


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

all_TF <- all_summry
aa <- as.character(rownames(all_TF)) %>% convert_human_to_mouse_symbols()
all_TF$mouse_symbol <- as.character(aa)
matrix <- all_TF[!duplicated(all_TF$mouse_symbol),]
head(matrix)
write.csv(matrix,"/mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/GSA_scNPC/LN_vs_nasopha/lymphnode_vs_nasopha_merge_GRCh37_GRCh38_bycalculate_mmsymbol.csv")


@@@@@@@@@@@@@@@将之前的normal nasopharyngeal数据换hg19重新跑一次  合并后就只剩16000多gene

nasopha <- read.csv(file = "/mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/GSA_scNPC/LN_vs_nasopha/hg19_humanNPC_nasopha_tpm_mmsymbol.csv")

#更换指定列名    colnames(nasopha)[1] <- NULL

lymph_node <- read.csv(file = "/mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/GSA_scNPC/LN_vs_nasopha/E_MTAB_1733_FPKM_allsamples_gooddata.csv",row.names=1)
lymph_node <- lymph_node[,34:38]
colnames(lymph_node) <- c("lymph_node_1","lymph_node_2","lymph_node_3","lymph_node_4","lymph_node_5")

fpkmToTpm <- function(FPKM){
    exp(log(lymph_node) - log(sum(lymph_node)) + log(1e6))
}
lymph_node <- fpkmToTpm(lymph_node)
lymph_node$X<- rownames(lymph_node)
lymph_node <- lymph_node[,c(6,1:5)]
rownames(lymph_node)=NULL



lymphnode_vs_nasopha_byENSEMBL <- merge(lymph_node,nasopha,by = "X", all=FALSE)


lymphnode_vs_nasopha_byENSEMBL <- na.omit(lymphnode_vs_nasopha_byENSEMBL)

nasopha_tpm <- c("nasopha_1","nasopha_2","nasopha_3")
lymphnode_tpm <- c("lymph_node_1","lymph_node_2","lymph_node_3","lymph_node_4","lymph_node_5")

lymphnode_tpm <- lymphnode_vs_nasopha_byENSEMBL[,colnames(lymphnode_vs_nasopha_byENSEMBL) %in% lymphnode_tpm]
nasopha_tpm <- lymphnode_vs_nasopha_byENSEMBL[,colnames(lymphnode_vs_nasopha_byENSEMBL) %in% nasopha_tpm]

logFC <- log2((rowMeans(lymphnode_tpm)+1) / (rowMeans(nasopha_tpm)+1))
library("future.apply")
p_values <- future_lapply(seq(1,nrow(lymphnode_tpm)), function(x){
  res <- wilcox.test(x = t(lymphnode_tpm[x,])[,1], y = t(nasopha_tpm[x,])[,1] )
  res$p.value
})
p <- unlist(p_values)
p.adj <- p.adjust(p, method = "fdr")
genelist <- as.data.frame(logFC)
genelist$p_values <- p
genelist$p_adj <- p.adj
head(genelist)
all_summry <- cbind(lymphnode_vs_nasopha_byENSEMBL,genelist)

write.csv(all_summry,"/mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/GSA_scNPC/LN_vs_nasopha/lymphnode_vs_nasopha_merge_GRCh37_bycalculate_mmsymbol.csv")


*******************************************************************消除lung和nasopharyngeal organoids的器官差异

lung_vs_nasopha <- read.csv(file = "/mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/GSA_scNPC/lung_vs_nasopha_organoids/lung_vs_nasopha__mouse_organoids_symbol.csv")
lung_vs_nasopha <- lung_vs_nasopha[!duplicated(lung_vs_nasopha$X),]
lung_vs_nasopha <- na.omit(lung_vs_nasopha)
rownames(lung_vs_nasopha) <- lung_vs_nasopha[,1]
lung_vs_nasopha <-  lung_vs_nasopha[,-1]
colnames(lung_vs_nasopha)

id <- colnames(lung_vs_nasopha)
deal <- c(rep('nasopha',3),rep('lung',3))
lung_vs_nasopha_df <- data.frame(id,deal)

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


dds <- DESeqDataSetFromMatrix(countData = lung_vs_nasopha, colData = lung_vs_nasopha_df, design = ~ deal)

countdata <- lung_vs_nasopha
coldata <- lung_vs_nasopha_df
DESeq_counts <- DESeq(dds)
normalized_DEseq <- counts(DESeq_counts, normalized=TRUE)

DEseq_res <- results(DESeq_counts,contrast=c("deal","lung","nasopha")) #最后一个是CTRL
res_1 <- cbind(normalized_DEseq,DEseq_res)


all_summry <- cbind(countdata,res_1)
head(all_summry)
write.csv(all_summry, "lung_vs_nasopha_normal_mouse_organoids_DESeq.csv")



*******************************************GSEA 做相关性*******************mimical***************
UP_DOWN <- read.csv(file = "/mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/GSA_scNPC/LN_vs_nasopha/up_down_gene_GRCh37_GRCh38_bycalculate_mmsymbol.csv")
colnames(UP_DOWN) <- c("up","down")

all <- read.csv(file = "/mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/GSA_scNPC/work_file/scNPC_Met_v_Pri_count_tpm_symbol_and_anno_mmsymbol.csv")

up_gene <- merge(all,UP_DOWN, by.x = "mouse_symbol", by.y = "up", all=FALSE)
write.csv(up_gene, "humanNPC_up_without_organeffect.csv")

down_gene <- merge(all,UP_DOWN, by.x = "mouse_symbol", by.y = "down", all=FALSE)
write.csv(down_gene, "humanNPC_down_without_organeffect.csv")


********************************************************一对一单独分析，差异基因取交集确定common差异基因*********************************
all <- read.csv(file = "/mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/GSA_scNPC/1vs1_analyze/scNPC_Met_v_Pri_tpm_mmsymbol.csv")

aa <- all[,c(1,2,5,8)]
head(aa)
Met <- aa[,3]
Pri <- aa[,2]
logFC <- log2((Met + 1) / (Pri + 1))
aa$logFC <- logFC

nrow(aa)
aa <- na.omit(aa)
y <- aa[c(1:nrow(aa)),5]
aa_up <- aa[with(aa, y >= 0.5),]
nrow(aa_up)
aa_down <- aa[with(aa, y <= -0.5),]
nrow(aa_down)

write.csv(aa, "/mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/GSA_scNPC/1vs1_analyze/scNPC_Met1_v_Pri1_tpm_log2FC_mmsymbol.csv")
﻿
bb <- all[,c(1,3,6,8)]
head(bb)
Met <- bb[,3]
Pri <- bb[,2]
logFC <- log2((Met + 1) / (Pri + 1))
bb$logFC <- logFC
nrow(bb)
bb <- na.omit(bb)
y <- bb[c(1:nrow(bb)),5]
bb_up <- bb[with(bb, y >= 0.5),]
nrow(bb_up)
bb_down <- bb[with(bb, y <= -0.5),]
nrow(bb_down)

write.csv(bb, "/mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/GSA_scNPC/1vs1_analyze/scNPC_Met2_v_Pri2_tpm_log2FC_mmsymbol.csv")


cc <- all[,c(1,4,7,8)]
head(cc)
Met <- cc[,3]
Pri <- cc[,2]
logFC <- log2((Met + 1) / (Pri + 1))
cc$logFC <- logFC
nrow(cc)
cc <- na.omit(cc)
y <- cc[c(1:nrow(cc)),5]
cc_up <- cc[with(cc, y >= 0.5),]
nrow(cc_up)
cc_down <- cc[with(cc, y <= -0.5),]
nrow(cc_down)

write.csv(cc, "/mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/GSA_scNPC/1vs1_analyze/scNPC_Met3_v_Pri3_tpm_log2FC_mmsymbol.csv")


ab_up <- merge(aa_up,bb_up, by.x = "X",by.y = "X", all=FALSE)
nrow(ab_up)
abc_up <- merge(ab_up,cc_up, by.x = "X",by.y = "X", all=FALSE)
nrow(abc_up)
colnames(abc_up)



ab_down <- merge(aa_down,bb_down, by.x = "X",by.y = "X", all=FALSE)
nrow(ab_down)
abc_down <- merge(ab_down,cc_down, by.x = "X",by.y = "X", all=FALSE)
nrow(abc_down)

write.csv(abc_up, "/mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/GSA_scNPC/1vs1_analyze/scNPC_Met_v_Pri_up_mmsymbol.csv")
write.csv(abc_down, "/mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/GSA_scNPC/1vs1_analyze/scNPC_Met_v_Pri_down_mmsymbol.csv")
