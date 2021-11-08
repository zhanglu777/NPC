#smart-seq2 样本信息 https://bigd.big.ac.cn/gsa-human/browse/HRA000087
*************************************************star 比对
*************************************************star 比对
*************************************************star 比对
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

*************************************************DESeq2 差异分析
*************************************************DESeq2 差异分析
*************************************************DESeq2 差异分析

vi sample_sampletable.csv

ids,sample,deal,order
1,HRR026852.Aligned.sortedByCoord.out.bam,Pri_1,primary,1
2,HRR026854.Aligned.sortedByCoord.out.bam,Pri_2,primary,2
3,HRR026856.Aligned.sortedByCoord.out.bam,Pri_3,primary,3
4,HRR026851.Aligned.sortedByCoord.out.bam,Met_1,metastasis,4
5,HRR026853.Aligned.sortedByCoord.out.bam,Met_2,metastasis,5
6,HRR026855.Aligned.sortedByCoord.out.bam,Met_3,metastasis,6

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

sampleTable <- read.csv(csvfile, row.names = 1)
sampleTable
filenames <- file.path(indir, sampleTable$ids)
file.exists(filenames)

bamfiles <- BamFileList(filenames, yieldSize=2000000)
seqinfo(bamfiles[1])
register(SerialParam())

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

#转换mouse ID
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
