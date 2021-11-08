*********************************************************************************************************************************
*********************************************************1.检验数据完整性*********************************************************
*********************************************************************************************************************************
cd /mnt/data/sequencedata/RNAseq/RNAseq_76_ZJP_bladder_20210105_8samples/CP2018091207789/H101SC19062343/RSCR0102/X101SC19062343-Z01/X101SC19062343-Z01-J075
md5sum -c md5.txt
md5sum -c MD5_Ctr-1_FRAS210005133-1r.txt
*********************************************************************************************************************************
************************************************************2.star比对***********************************************************
*********************************************************************************************************************************
vi config 

ctrl_1 /mnt/data/sequencedata/RNAseq/RNAseq_77_WXD_NPC_20210108_11samples/CP2018091214925/H101SC20030543/RSCS7000/X101SC20030543-Z01/X101SC20030543-Z01-J090/2.cleandata/Ctr-1_FRAS210005133-1r/Ctr-1_FRAS210005133-1r_1.clean.fq.gz /mnt/data/sequencedata/RNAseq/RNAseq_77_WXD_NPC_20210108_11samples/CP2018091214925/H101SC20030543/RSCS7000/X101SC20030543-Z01/X101SC20030543-Z01-J090/2.cleandata/Ctr-1_FRAS210005133-1r/Ctr-1_FRAS210005133-1r_2.clean.fq.gz
ctrl_2 /mnt/data/sequencedata/RNAseq/RNAseq_77_WXD_NPC_20210108_11samples/CP2018091214925/H101SC20030543/RSCS7000/X101SC20030543-Z01/X101SC20030543-Z01-J090/2.cleandata/Ctr-2_FRAS210005134-1r/Ctr-2_FRAS210005134-1r_1.clean.fq.gz /mnt/data/sequencedata/RNAseq/RNAseq_77_WXD_NPC_20210108_11samples/CP2018091214925/H101SC20030543/RSCS7000/X101SC20030543-Z01/X101SC20030543-Z01-J090/2.cleandata/Ctr-2_FRAS210005134-1r/Ctr-2_FRAS210005134-1r_2.clean.fq.gz
ctrl_3 /mnt/data/sequencedata/RNAseq/RNAseq_77_WXD_NPC_20210108_11samples/CP2018091214925/H101SC20030543/RSCS7000/X101SC20030543-Z01/X101SC20030543-Z01-J090/2.cleandata/Ctr-3_FRAS210005135-1r/Ctr-3_FRAS210005135-1r_1.clean.fq.gz /mnt/data/sequencedata/RNAseq/RNAseq_77_WXD_NPC_20210108_11samples/CP2018091214925/H101SC20030543/RSCS7000/X101SC20030543-Z01/X101SC20030543-Z01-J090/2.cleandata/Ctr-3_FRAS210005135-1r/Ctr-3_FRAS210005135-1r_2.clean.fq.gz
LMP1_1 /mnt/data/sequencedata/RNAseq/RNAseq_77_WXD_NPC_20210108_11samples/CP2018091214925/H101SC20030543/RSCS7000/X101SC20030543-Z01/X101SC20030543-Z01-J090/2.cleandata/LMP1-1_FRAS210005127-2r/LMP1-1_FRAS210005127-2r_1.clean.fq.gz /mnt/data/sequencedata/RNAseq/RNAseq_77_WXD_NPC_20210108_11samples/CP2018091214925/H101SC20030543/RSCS7000/X101SC20030543-Z01/X101SC20030543-Z01-J090/2.cleandata/LMP1-1_FRAS210005127-2r/LMP1-1_FRAS210005127-2r_2.clean.fq.gz
LMP1_2 /mnt/data/sequencedata/RNAseq/RNAseq_77_WXD_NPC_20210108_11samples/CP2018091214925/H101SC20030543/RSCS7000/X101SC20030543-Z01/X101SC20030543-Z01-J090/2.cleandata/LMP1-2_FRAS210005128-2r/LMP1-2_FRAS210005128-2r_1.clean.fq.gz /mnt/data/sequencedata/RNAseq/RNAseq_77_WXD_NPC_20210108_11samples/CP2018091214925/H101SC20030543/RSCS7000/X101SC20030543-Z01/X101SC20030543-Z01-J090/2.cleandata/LMP1-2_FRAS210005128-2r/LMP1-2_FRAS210005128-2r_2.clean.fq.gz
LMP1_3 /mnt/data/sequencedata/RNAseq/RNAseq_77_WXD_NPC_20210108_11samples/CP2018091214925/H101SC20030543/RSCS7000/X101SC20030543-Z01/X101SC20030543-Z01-J090/2.cleandata/LMP1-3_FRAS210005129-2r/LMP1-3_FRAS210005129-2r_1.clean.fq.gz /mnt/data/sequencedata/RNAseq/RNAseq_77_WXD_NPC_20210108_11samples/CP2018091214925/H101SC20030543/RSCS7000/X101SC20030543-Z01/X101SC20030543-Z01-J090/2.cleandata/LMP1-3_FRAS210005129-2r/LMP1-3_FRAS210005129-2r_2.clean.fq.gz
Lung_ctrl_1 /mnt/data/sequencedata/RNAseq/RNAseq_77_WXD_NPC_20210108_11samples/CP2018091214925/H101SC20030543/RSCS7000/X101SC20030543-Z01/X101SC20030543-Z01-J090/2.cleandata/Lung-Ctr-1_FRAS210005136-2r/Lung-Ctr-1_FRAS210005136-2r_1.clean.fq.gz /mnt/data/sequencedata/RNAseq/RNAseq_77_WXD_NPC_20210108_11samples/CP2018091214925/H101SC20030543/RSCS7000/X101SC20030543-Z01/X101SC20030543-Z01-J090/2.cleandata/Lung-Ctr-1_FRAS210005136-2r/Lung-Ctr-1_FRAS210005136-2r_2.clean.fq.gz
Lung_ctrl_2 /mnt/data/sequencedata/RNAseq/RNAseq_77_WXD_NPC_20210108_11samples/CP2018091214925/H101SC20030543/RSCS7000/X101SC20030543-Z01/X101SC20030543-Z01-J090/2.cleandata/Lung-Ctr-2_FRAS210005137-2r/Lung-Ctr-2_FRAS210005137-2r_1.clean.fq.gz /mnt/data/sequencedata/RNAseq/RNAseq_77_WXD_NPC_20210108_11samples/CP2018091214925/H101SC20030543/RSCS7000/X101SC20030543-Z01/X101SC20030543-Z01-J090/2.cleandata/Lung-Ctr-2_FRAS210005137-2r/Lung-Ctr-2_FRAS210005137-2r_2.clean.fq.gz
Lung_LMP1_1 /mnt/data/sequencedata/RNAseq/RNAseq_77_WXD_NPC_20210108_11samples/CP2018091214925/H101SC20030543/RSCS7000/X101SC20030543-Z01/X101SC20030543-Z01-J090/2.cleandata/Lung_-LMP1-1_FRAS210005130-1r/Lung_-LMP1-1_FRAS210005130-1r_1.clean.fq.gz /mnt/data/sequencedata/RNAseq/RNAseq_77_WXD_NPC_20210108_11samples/CP2018091214925/H101SC20030543/RSCS7000/X101SC20030543-Z01/X101SC20030543-Z01-J090/2.cleandata/Lung_-LMP1-1_FRAS210005130-1r/Lung_-LMP1-1_FRAS210005130-1r_2.clean.fq.gz
Lung_LMP1_2 /mnt/data/sequencedata/RNAseq/RNAseq_77_WXD_NPC_20210108_11samples/CP2018091214925/H101SC20030543/RSCS7000/X101SC20030543-Z01/X101SC20030543-Z01-J090/2.cleandata/Lung-LMP1-2_FRAS210005131-1r/Lung-LMP1-2_FRAS210005131-1r_1.clean.fq.gz /mnt/data/sequencedata/RNAseq/RNAseq_77_WXD_NPC_20210108_11samples/CP2018091214925/H101SC20030543/RSCS7000/X101SC20030543-Z01/X101SC20030543-Z01-J090/2.cleandata/Lung-LMP1-2_FRAS210005131-1r/Lung-LMP1-2_FRAS210005131-1r_2.clean.fq.gz
Lung_LMP1_3 /mnt/data/sequencedata/RNAseq/RNAseq_77_WXD_NPC_20210108_11samples/CP2018091214925/H101SC20030543/RSCS7000/X101SC20030543-Z01/X101SC20030543-Z01-J090/2.cleandata/Lung-LMP1-3_FRAS210005132-1r/Lung-LMP1-3_FRAS210005132-1r_1.clean.fq.gz /mnt/data/sequencedata/RNAseq/RNAseq_77_WXD_NPC_20210108_11samples/CP2018091214925/H101SC20030543/RSCS7000/X101SC20030543-Z01/X101SC20030543-Z01-J090/2.cleandata/Lung-LMP1-3_FRAS210005132-1r/Lung-LMP1-3_FRAS210005132-1r_2.clean.fq.gz
normal_organoid_1 /mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/organoid_raw_data/nasopha_mouse1_organoid_1.fq.gz /mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/organoid_raw_data/nasopha_mouse1_organoid_2.fq.gz
normal_organoid_2 /mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/organoid_raw_data/nasopha_mouse2_organoid_1.fq.gz /mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/organoid_raw_data/nasopha_mouse2_organoid_2.fq.gz
normal_organoid_3 /mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/organoid_raw_data/nasopha_mouse3_organoid_1.fq.gz /mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/organoid_raw_data/nasopha_mouse3_organoid_2.fq.gz

cat config |while read id;
do echo $id
arr=($id)
fq1=${arr[1]}
fq2=${arr[2]}
sample=${arr[0]}
echo $fq1
echo $fq2
echo $sample
STAR  --readFilesCommand zcat \
--outSAMtype BAM SortedByCoordinate \
--runThreadN 20 \
--readFilesIn $fq1 $fq2 \
--genomeDir /mnt/data/public_data/genome_index/genome_index/mm10_genome_index_150/ \
--outFileNamePrefix /mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/star_mm10_out/$sample. ;done

#构建index
ls *.sortedByCoord.out.bam |while read id
do
samtools index ./$id
done

#igv
bash /mnt/data/program/IGV_2.4.10/igv.sh



*********************************************************************************************************************************
**********************************************************3.差异分析和KEGG、GO、GSEA************************************************
*********************************************************************************************************************************

****************************************************primary_TMPL vs primary_TMP*********************************************************
#move to bamfiles path

cd /mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/star_mm10_out

vi sample_sampletable_primary.csv

ids,sample,deal,order
1,ctrl_1.Aligned.sortedByCoord.out.bam,ctrl_1,ctrl,1
2,ctrl_2.Aligned.sortedByCoord.out.bam,ctrl_2,ctrl,2
3,ctrl_3.Aligned.sortedByCoord.out.bam,ctrl_3,ctrl,3
4,LMP1_2.Aligned.sortedByCoord.out.bam,LMP1_2,LMP1,4
5,LMP1_3.Aligned.sortedByCoord.out.bam,LMP1_3,LMP1,5


############# R environment  #####
project <- c("primary_LMP1")
deal_design <- c("LMP1","ctrl")
significant_cutoff <- c(0.5)
organism <- "mouse"
##物种为人就是hsa，是老鼠就是mouse
file_path <- "/mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/primary_LMP1_vs_ctrl/"
sample_sampletable.path <- "/mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/star_mm10_out"
load(file="/mnt/data/user_data/xiangyu/workshop/WORKFLOW_RNAseq/ebg_mm10.RData")
load(file="/mnt/data/user_data/xiangyu/workshop/WORKFLOW_RNAseq/txdb_mm10.RData")  
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
csvfile <- file.path(indir, "sample_sampletable_primary.csv")
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
library("ggplot2")
pcaData <- plotPCA(rld, intgroup = c("sample","order"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(x = PC1, y = PC2, color = sample)) +
geom_point(size =1) +
geom_text(aes(label=order))+
xlab(paste0("PC1: ", percentVar[1], "% variance")) +
ylab(paste0("PC2: ", percentVar[2], "% variance")) +
coord_fixed()
dev.off()

#normalize by TPM
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
#DESeq2 差异分析
dds <- DESeq(dds)
res_1 <- results(dds, contrast=res_1_filter)
colnames(res_1) = paste(colnames_1,colnames(res_1),sep="_")
write.csv(res_1, file = res_1_file.csv)

res_1$ensembl <- mapIds(x = anno_data,
                        keys = rownames(res_1),
						keytype ="SYMBOL",
						column ="ENSEMBL",
						multiVals="first")
res_1$entrez <- mapIds(x = anno_data,
                        keys = rownames(res_1),
						keytype ="SYMBOL",
						column ="ENTREZID",
						multiVals="first")
AA <- rownames(res_1)
AA <- as.character(AA)
res_1$GENENAME <- mapIds(x = anno_data,
                        keys = AA,
						keytype ="SYMBOL",
						column ="GENENAME",
						multiVals="first")
write.csv(res_1, file = res_1_file_sym_name)
all_summry <- cbind(count_and_tpm,res_1)
write.csv(all_summry, file = ALL.CSV_FILE)

res_1 <- na.omit(res_1)
y <- res_1[c(1:nrow(res_1)),6]
res_1 <- res_1[with(res_1,y<0.05),]
y <- res_1[c(1:nrow(res_1)),2]
#####KEGG 富集
upres_1 <- res_1[with(res_1,y>=significant_cutoff),]
ee	<-as.matrix(upres_1$entrez)
	dd <- as.vector(ee)
KEGGupres_1 <- enrichKEGG(gene =dd, 
					organism = organism, 
					keyType = "ncbi-geneid",
					 pvalueCutoff = 0.05,
				       pAdjustMethod = "BH", 
				       minGSSize = 10, 
				       maxGSSize = 500,
				       qvalueCutoff = 0.2, 
				       use_internal_data = FALSE)
KEGG_pro_enhan_up <- setReadable(KEGGupres_1, "org.Mm.eg.db", keyType="ENTREZID")
write.csv(KEGG_pro_enhan_up, file = "primary_LMP1_07_LMP1_v_ctrl_KEGG_UP_genes.csv")
pdf(file = KEGGupres_1_pdf)
dotplot(KEGGupres_1, showCategory=20)
dev.off()

downres_1 <- res_1[with(res_1,y<=-significant_cutoff),]
hh	<-as.matrix(downres_1$entrez)
	ii <- as.vector(hh)
KEGGdownres_1 <- enrichKEGG(gene =ii, 
					organism = organism, 
					keyType = "ncbi-geneid",
					 pvalueCutoff = 0.05,
				       pAdjustMethod = "BH", 
				       minGSSize = 10, 
				       maxGSSize = 500,
				       qvalueCutoff = 0.2, 
				       use_internal_data = FALSE)

KEGG_pro_enhan_down <- setReadable(KEGGdownres_1, "org.Mm.eg.db", keyType="ENTREZID")
write.csv(KEGG_pro_enhan_down, file = "primary_LMP1_07_LMP1_v_ctrl_KEGG_DOWN_genes.csv")
pdf(file = "primary_LMP1_v_ctrl_KEGG_UP_top15.pdf")
dotplot(KEGGupres_1, showCategory=15)
dev.off()

pdf(file = KEGGres_1_all_pdf)
KEGGupres_1uporder <- KEGGupres_1[order((KEGGupres_1$Count),decreasing = TRUE),]
KEGGdownres_1downorder <- KEGGdownres_1[order((KEGGdownres_1$Count),decreasing = TRUE),]
down_kegg<-KEGGdownres_1downorder[KEGGdownres_1downorder$p.adjust<0.05,];
down_kegg$group=-1
up_kegg<-KEGGupres_1uporder[KEGGupres_1uporder$p.adjust<0.05,];
up_kegg$group=1
down_kegg1 <- down_kegg[-(21:nrow(down_kegg)),]
up_kegg1 <- up_kegg[-(21:nrow(up_kegg)),]
  dat=rbind(up_kegg1,down_kegg1)
  colnames(dat)
    dat$p.adjust = -log10(dat$p.adjust)
  dat$p.adjust=dat$p.adjust*dat$group 
  dat=dat[order(dat$p.adjust,decreasing = F),]
  library("ggplot2")
  g_kegg<- ggplot(dat, aes(x=reorder(Description,order(p.adjust, decreasing = F)), y=p.adjust, fill=group)) + 
    geom_bar(stat="identity") + 
    scale_fill_gradient(low="blue",high="red",guide = FALSE) + 
    scale_x_discrete(name ="Pathway names") +
    scale_y_continuous(name ="-log10p_adjust") +
    coord_flip() +
    ggtitle(KEGGres_1_all_title)
  print(g_kegg)
dev.off()

#####GO 富集
ee	<-as.matrix(upres_1$entrez)
	dd <- as.vector(ee)
GOupres_1_all <- enrichGO(gene = dd, 
	           OrgDb = anno_data,
				ont = "all", 
		             pvalueCutoff = 0.01, 
                     pAdjustMethod = "BH", 
                     qvalueCutoff = 0.05,
                     minGSSize = 10, 
                     maxGSSize = 500, 
                     readable = T, 
                     pool = FALSE)
write.csv(GOupres_1_all,file=GOres_1_all_UP_csv)

ff <- barplot(GOupres_1_all, split="ONTOLOGY",showCategory=15) + facet_grid(ONTOLOGY~., scale="free")
ff1 <- ff + theme(axis.text.y = element_text(size = 8)) +labs(title = NULL)
pdf(file = GOres_1_all_UP_pdf)
ff1
dev.off()

hh	<-as.matrix(downres_1$entrez)
	ii <- as.vector(hh)
GOdownres_1_all <- enrichGO(gene = ii, 
	           OrgDb = anno_data,
				ont = "all", 
		             pvalueCutoff = 0.01, 
                     pAdjustMethod = "BH", 
                     qvalueCutoff = 0.05,
                     minGSSize = 10, 
                     maxGSSize = 500, 
                     readable = T, 
                     pool = FALSE)
write.csv(GOdownres_1_all,file=GOres_1_all_DOWN_csv)

ff <- barplot(GOdownres_1_all, split="ONTOLOGY",showCategory=15) + facet_grid(ONTOLOGY~., scale="free")
ff1 <- ff + theme(axis.text.y = element_text(size = 8)) +labs(title = NULL)
pdf(file = GOres_1_all_DOWN_pdf)
ff1
dev.off()

ff <- barplot(GOdownres_1_all, split="ONTOLOGY",showCategory=5) + facet_grid(ONTOLOGY~., scale="free")
ff1 <- ff + theme(axis.text.y = element_text(size = 8)) +labs(title = NULL)
pdf(file = "primary_LMP1_v_ctrl_GO_DOWN_top5.pdf")
ff1
dev.off()

#####GSEA 富集
gsea_data <- read.csv("/mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/primary_LMP1_vs_ctrl/primary_LMP1_05_LMP1_v_ctrl_count_tpm_symbol_and_anno.csv")
gsea_data$Description <- c("NULL")
gsea_data <- gsea_data[,c(18,21,10,11,7,8,9)]
nrow(gsea_data)
gsea_data <- na.omit(gsea_data)
nrow(gsea_data)
colnames(gsea_data)
colnames(gsea_data) <- c("NAME","Description","LMP1_2","LMP1_3","ctrl_1","ctrl_2","ctrl_3")
write.csv(gsea_data,file="/mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/primary_LMP1_vs_ctrl/gsea_data.csv")
#在excel中改好格式(就以Ensembl id为行名)，然后直接保存为txt格式，再导入R修改文件名为.gct即可
#加上  
       #1.2
       #总基因行数 总样本数

vi gsea.cls
5 2 1
# LMP1 ctrl
0 0 1 1 1

#样本多的时候，-Xmx1024m改为-Xmx10240 \       9022      注意\后面不能有空格
vi gsea_shell.sh

java -cp /mnt/data/user_data/xiangyu/programme/gsea/gsea-3.0.jar -Xmx1024m xtools.gsea.Gsea \
-res gsea_data.gct \
-cls gsea.cls \
-gmx /mnt/data/user_data/xiangyu/programme/gsea/msigdb_v6.1_GMTs/h.all.v6.1.symbols.gmt \
-chip /mnt/data/user_data/xiangyu/programme/gsea/chip/ENSEMBL_mouse_gene.chip \
-collapse true -mode Max_probe -norm meandiv -nperm 1000 -permute gene_set \
-rnd_type no_balance -scoring_scheme weighted -rpt_label my_analysis \
-metric Signal2Noise -sort real -order descending -include_only_symbols true \
-make_sets true -median false -num 100 -plot_top_x 50 -rnd_seed timestamp \
-save_rnd_lists false -set_max 500 -set_min 15 -zip_report false \
-out ./LMP1_vs_ctrl_h -gui false

java -cp /mnt/data/user_data/xiangyu/programme/gsea/gsea-3.0.jar -Xmx1024m xtools.gsea.Gsea \
-res gsea_data.gct \
-cls gsea.cls \
-gmx /mnt/data/user_data/xiangyu/programme/gsea/msigdb_v6.1_GMTs/c1.all.v6.1.symbols.gmt \
-chip /mnt/data/user_data/xiangyu/programme/gsea/chip/ENSEMBL_mouse_gene.chip \
-collapse true -mode Max_probe -norm meandiv -nperm 1000 -permute gene_set \
-rnd_type no_balance -scoring_scheme weighted -rpt_label my_analysis \
-metric Signal2Noise -sort real -order descending -include_only_symbols true \
-make_sets true -median false -num 100 -plot_top_x 50 -rnd_seed timestamp \
-save_rnd_lists false -set_max 500 -set_min 15 -zip_report false \
-out ./LMP1_vs_ctrl_c1 -gui false

java -cp /mnt/data/user_data/xiangyu/programme/gsea/gsea-3.0.jar -Xmx1024m xtools.gsea.Gsea \
-res gsea_data.gct \
-cls gsea.cls \
-gmx /mnt/data/user_data/xiangyu/programme/gsea/msigdb_v6.1_GMTs/c2.all.v6.1.symbols.gmt \
-chip /mnt/data/user_data/xiangyu/programme/gsea/chip/ENSEMBL_mouse_gene.chip \
-collapse true -mode Max_probe -norm meandiv -nperm 1000 -permute gene_set \
-rnd_type no_balance -scoring_scheme weighted -rpt_label my_analysis \
-metric Signal2Noise -sort real -order descending -include_only_symbols true \
-make_sets true -median false -num 100 -plot_top_x 50 -rnd_seed timestamp \
-save_rnd_lists false -set_max 500 -set_min 15 -zip_report false \
-out ./LMP1_vs_ctrl_c2 -gui false

java -cp /mnt/data/user_data/xiangyu/programme/gsea/gsea-3.0.jar -Xmx1024m xtools.gsea.Gsea \
-res gsea_data.gct \
-cls gsea.cls \
-gmx /mnt/data/user_data/xiangyu/programme/gsea/msigdb_v6.1_GMTs/c3.all.v6.1.symbols.gmt \
-chip /mnt/data/user_data/xiangyu/programme/gsea/chip/ENSEMBL_mouse_gene.chip \
-collapse true -mode Max_probe -norm meandiv -nperm 1000 -permute gene_set \
-rnd_type no_balance -scoring_scheme weighted -rpt_label my_analysis \
-metric Signal2Noise -sort real -order descending -include_only_symbols true \
-make_sets true -median false -num 100 -plot_top_x 50 -rnd_seed timestamp \
-save_rnd_lists false -set_max 500 -set_min 15 -zip_report false \
-out ./LMP1_vs_ctrl_c3 -gui false

java -cp /mnt/data/user_data/xiangyu/programme/gsea/gsea-3.0.jar -Xmx1024m xtools.gsea.Gsea \
-res gsea_data.gct \
-cls gsea.cls \
-gmx /mnt/data/user_data/xiangyu/programme/gsea/msigdb_v6.1_GMTs/c4.all.v6.1.symbols.gmt \
-chip /mnt/data/user_data/xiangyu/programme/gsea/chip/ENSEMBL_mouse_gene.chip \
-collapse true -mode Max_probe -norm meandiv -nperm 1000 -permute gene_set \
-rnd_type no_balance -scoring_scheme weighted -rpt_label my_analysis \
-metric Signal2Noise -sort real -order descending -include_only_symbols true \
-make_sets true -median false -num 100 -plot_top_x 50 -rnd_seed timestamp \
-save_rnd_lists false -set_max 500 -set_min 15 -zip_report false \
-out ./LMP1_vs_ctrl_c4 -gui false

java -cp /mnt/data/user_data/xiangyu/programme/gsea/gsea-3.0.jar -Xmx1024m xtools.gsea.Gsea \
-res gsea_data.gct \
-cls gsea.cls \
-gmx /mnt/data/user_data/xiangyu/programme/gsea/msigdb_v6.1_GMTs/c5.all.v6.1.symbols.gmt \
-chip /mnt/data/user_data/xiangyu/programme/gsea/chip/ENSEMBL_mouse_gene.chip \
-collapse true -mode Max_probe -norm meandiv -nperm 1000 -permute gene_set \
-rnd_type no_balance -scoring_scheme weighted -rpt_label my_analysis \
-metric Signal2Noise -sort real -order descending -include_only_symbols true \
-make_sets true -median false -num 100 -plot_top_x 50 -rnd_seed timestamp \
-save_rnd_lists false -set_max 500 -set_min 15 -zip_report false \
-out ./LMP1_vs_ctrl_c5 -gui false

java -cp /mnt/data/user_data/xiangyu/programme/gsea/gsea-3.0.jar -Xmx1024m xtools.gsea.Gsea \
-res gsea_data.gct \
-cls gsea.cls \
-gmx /mnt/data/user_data/xiangyu/programme/gsea/msigdb_v6.1_GMTs/c6.all.v6.1.symbols.gmt \
-chip /mnt/data/user_data/xiangyu/programme/gsea/chip/ENSEMBL_mouse_gene.chip \
-collapse true -mode Max_probe -norm meandiv -nperm 1000 -permute gene_set \
-rnd_type no_balance -scoring_scheme weighted -rpt_label my_analysis \
-metric Signal2Noise -sort real -order descending -include_only_symbols true \
-make_sets true -median false -num 100 -plot_top_x 50 -rnd_seed timestamp \
-save_rnd_lists false -set_max 500 -set_min 15 -zip_report false \
-out ./LMP1_vs_ctrl_c6 -gui false

java -cp /mnt/data/user_data/xiangyu/programme/gsea/gsea-3.0.jar -Xmx1024m xtools.gsea.Gsea \
-res gsea_data.gct \
-cls gsea.cls \
-gmx /mnt/data/user_data/xiangyu/programme/gsea/msigdb_v6.1_GMTs/c7.all.v6.1.symbols.gmt \
-chip /mnt/data/user_data/xiangyu/programme/gsea/chip/ENSEMBL_mouse_gene.chip \
-collapse true -mode Max_probe -norm meandiv -nperm 1000 -permute gene_set \
-rnd_type no_balance -scoring_scheme weighted -rpt_label my_analysis \
-metric Signal2Noise -sort real -order descending -include_only_symbols true \
-make_sets true -median false -num 100 -plot_top_x 50 -rnd_seed timestamp \
-save_rnd_lists false -set_max 500 -set_min 15 -zip_report false \
-out ./LMP1_vs_ctrl_c7 -gui false


bash gsea_shell.sh


****************************************************metastasis_TMPL vs metastasis_TPM*******************************************

############################
#move to bamfiles path  #####
#move to bamfiles path  #####
#############################
cd /mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/star_mm10_out

vi sample_sampletable_metastasis.csv

ids,sample,deal,order
1,Lung_ctrl_1.Aligned.sortedByCoord.out.bam,ctrl_1,ctrl,1
2,Lung_ctrl_2.Aligned.sortedByCoord.out.bam,ctrl_2,ctrl,2
3,Lung_LMP1_1.Aligned.sortedByCoord.out.bam,LMP1_1,LMP1,3
4,Lung_LMP1_2.Aligned.sortedByCoord.out.bam,LMP1_2,LMP1,4
5,Lung_LMP1_3.Aligned.sortedByCoord.out.bam,LMP1_3,LMP1,5

############################
#work in R environment  #####
#work in R environment  #####
#############################



project <- c("metastasis_LMP1")
deal_design <- c("LMP1","ctrl")
significant_cutoff <- c(1)
organism <- "mouse"
##物种为人就是hsa，是老鼠就是mouse

file_path <- "/mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/metastasis_LMP1_vs_ctrl/"
sample_sampletable.path <- "/mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/star_mm10_out"
load(file="/mnt/data/user_data/xiangyu/workshop/WORKFLOW_RNAseq/ebg_mm10.RData")
load(file="/mnt/data/user_data/xiangyu/workshop/WORKFLOW_RNAseq/txdb_mm10.RData")  


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
csvfile <- file.path(indir, "sample_sampletable_metastasis.csv")
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
library("ggplot2")
pcaData <- plotPCA(rld, intgroup = c("sample","order"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(x = PC1, y = PC2, color = sample)) +
geom_point(size =1) +
geom_text(aes(label=order))+
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

res_1$ensembl <- mapIds(x = anno_data,
                        keys = rownames(res_1),
						keytype ="SYMBOL",
						column ="ENSEMBL",
						multiVals="first")
res_1$entrez <- mapIds(x = anno_data,
                        keys = rownames(res_1),
						keytype ="SYMBOL",
						column ="ENTREZID",
						multiVals="first")
AA <- rownames(res_1)
AA <- as.character(AA)
res_1$GENENAME <- mapIds(x = anno_data,
                        keys = AA,
						keytype ="SYMBOL",
						column ="GENENAME",
						multiVals="first")
write.csv(res_1, file = res_1_file_sym_name)
all_summry <- cbind(count_and_tpm,res_1)
write.csv(all_summry, file = ALL.CSV_FILE)

res_1 <- na.omit(res_1)
y <- res_1[c(1:nrow(res_1)),6]
res_1 <- res_1[with(res_1,y<0.05),]
y <- res_1[c(1:nrow(res_1)),2]

upres_1 <- res_1[with(res_1,y>=significant_cutoff),]
ee	<-as.matrix(upres_1$entrez)
	dd <- as.vector(ee)
KEGGupres_1 <- enrichKEGG(gene =dd, 
					organism = organism, 
					keyType = "ncbi-geneid",
					 pvalueCutoff = 0.05,
				       pAdjustMethod = "BH", 
				       minGSSize = 10, 
				       maxGSSize = 500,
				       qvalueCutoff = 0.2, 
				       use_internal_data = FALSE)

downres_1  <- res_1[with(res_1,y<=-significant_cutoff),]
downres_1  <- res_1[with(res_1,y<=-0.5),]
hh	<-as.matrix(downres_1$entrez)
	ii <- as.vector(hh)
KEGGdownres_1 <- enrichKEGG(gene =ii, 
					organism = organism, 
					keyType = "ncbi-geneid",
					 pvalueCutoff = 0.05,
				       pAdjustMethod = "BH", 
				       minGSSize = 10, 
				       maxGSSize = 500,
				       qvalueCutoff = 0.2, 
				       use_internal_data = FALSE)

write.csv(KEGGupres_1, file = KEGGupres_1_file.csv)
write.csv(KEGGdownres_1, file = KEGGdownres_1_file.csv)

pdf(file = KEGGupres_1_pdf)
dotplot(KEGGupres_1, showCategory=20)
dev.off()
pdf(file = KEGGdownres_1_pdf)
dotplot(KEGGdownres_1, showCategory=20)
dev.off()

pdf(file = KEGGres_1_all_pdf)
KEGGupres_1uporder <- KEGGupres_1[order((KEGGupres_1$Count),decreasing = TRUE),]
KEGGdownres_1downorder <- KEGGdownres_1[order((KEGGdownres_1$Count),decreasing = TRUE),]
down_kegg<-KEGGdownres_1downorder[KEGGdownres_1downorder$p.adjust<0.05,];
down_kegg$group=-1
up_kegg<-KEGGupres_1uporder[KEGGupres_1uporder$p.adjust<0.05,];
up_kegg$group=1
down_kegg1 <- down_kegg[-(21:nrow(down_kegg)),]
up_kegg1 <- up_kegg[-(21:nrow(up_kegg)),]
  dat=rbind(up_kegg1,down_kegg1)
  colnames(dat)
    dat$p.adjust = -log10(dat$p.adjust)
  dat$p.adjust=dat$p.adjust*dat$group 
  dat=dat[order(dat$p.adjust,decreasing = F),]
  library("ggplot2")
  g_kegg<- ggplot(dat, aes(x=reorder(Description,order(p.adjust, decreasing = F)), y=p.adjust, fill=group)) + 
    geom_bar(stat="identity") + 
    scale_fill_gradient(low="blue",high="red",guide = FALSE) + 
    scale_x_discrete(name ="Pathway names") +
    scale_y_continuous(name ="-log10p_adjust") +
    coord_flip() +
    ggtitle(KEGGres_1_all_title)
  print(g_kegg)
dev.off()

GOupres_1_all <- enrichGO(gene = dd, 
	           OrgDb = anno_data,
				ont = "all", 
		             pvalueCutoff = 0.01, 
                     pAdjustMethod = "BH", 
                     qvalueCutoff = 0.05,
                     minGSSize = 10, 
                     maxGSSize = 500, 
                     readable = T, 
                     pool = FALSE)
write.csv(GOupres_1_all,file=GOres_1_all_UP_csv)



ff <- barplot(GOupres_1_all, split="ONTOLOGY",showCategory=15) + facet_grid(ONTOLOGY~., scale="free")
ff1 <- ff + theme(axis.text.y = element_text(size = 8)) +labs(title = NULL)
pdf(file = GOres_1_all_UP_pdf)
ff1
dev.off()

hh	<-as.matrix(downres_1$entrez)
	ii <- as.vector(hh)
GOdownres_1_all <- enrichGO(gene = ii, 
	           OrgDb = anno_data,
				ont = "all", 
		             pvalueCutoff = 0.01, 
                     pAdjustMethod = "BH", 
                     qvalueCutoff = 0.05,
                     minGSSize = 10, 
                     maxGSSize = 500, 
                     readable = T, 
                     pool = FALSE)
write.csv(GOdownres_1_all,file=GOres_1_all_DOWN_csv)

ff <- barplot(GOdownres_1_all, split="ONTOLOGY",showCategory=15) + facet_grid(ONTOLOGY~., scale="free")
ff1 <- ff + theme(axis.text.y = element_text(size = 8)) +labs(title = NULL)
pdf(file = GOres_1_all_DOWN_pdf)
ff1
dev.off()


*********gsea输入文件准备***********
gsea_data <- read.csv("/mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/metastasis_LMP1_vs_ctrl/metastasis_LMP1_05_LMP1_v_ctrl_count_tpm_symbol_and_anno.csv")
gsea_data$Description <- c("NULL")
gsea_data <- gsea_data[,c(18,21,9,10,11,7,8)]
nrow(gsea_data)
gsea_data <- na.omit(gsea_data)
nrow(gsea_data)
colnames(gsea_data)
colnames(gsea_data) <- c("NAME","Description","Lung_LMP1_1","Lung_LMP1_2","Lung_LMP1_3","Lung_ctrl_1","Lung_ctrl_2")
write.csv(gsea_data,file="/mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/metastasis_LMP1_vs_ctrl/gsea_data.csv")
#在excel中改好格式(就以Ensembl id为行名)，然后直接保存为txt格式，再导入R修改文件名为.gct即可
#加上  
       #1.2
       #总基因行数 总样本数

vi gsea.cls
5 2 1
# Lung_LMP1 Lung_ctrl
0 0 0 1 1



#样本多的时候，-Xmx1024m改为-Xmx10240 \       9022      注意\后面不能有空格
vi gsea_shell.sh

java -cp /mnt/data/user_data/xiangyu/programme/gsea/gsea-3.0.jar -Xmx1024m xtools.gsea.Gsea \
-res gsea_data.gct \
-cls gsea.cls \
-gmx /mnt/data/user_data/xiangyu/programme/gsea/msigdb_v6.1_GMTs/h.all.v6.1.symbols.gmt \
-chip /mnt/data/user_data/xiangyu/programme/gsea/chip/ENSEMBL_mouse_gene.chip \
-collapse true -mode Max_probe -norm meandiv -nperm 1000 -permute gene_set \
-rnd_type no_balance -scoring_scheme weighted -rpt_label my_analysis \
-metric Signal2Noise -sort real -order descending -include_only_symbols true \
-make_sets true -median false -num 100 -plot_top_x 50 -rnd_seed timestamp \
-save_rnd_lists false -set_max 500 -set_min 15 -zip_report false \
-out ./Lung_LMP1_vs_ctrl_h -gui false

java -cp /mnt/data/user_data/xiangyu/programme/gsea/gsea-3.0.jar -Xmx1024m xtools.gsea.Gsea \
-res gsea_data.gct \
-cls gsea.cls \
-gmx /mnt/data/user_data/xiangyu/programme/gsea/msigdb_v6.1_GMTs/c1.all.v6.1.symbols.gmt \
-chip /mnt/data/user_data/xiangyu/programme/gsea/chip/ENSEMBL_mouse_gene.chip \
-collapse true -mode Max_probe -norm meandiv -nperm 1000 -permute gene_set \
-rnd_type no_balance -scoring_scheme weighted -rpt_label my_analysis \
-metric Signal2Noise -sort real -order descending -include_only_symbols true \
-make_sets true -median false -num 100 -plot_top_x 50 -rnd_seed timestamp \
-save_rnd_lists false -set_max 500 -set_min 15 -zip_report false \
-out ./Lung_LMP1_vs_ctrl_c1 -gui false

java -cp /mnt/data/user_data/xiangyu/programme/gsea/gsea-3.0.jar -Xmx1024m xtools.gsea.Gsea \
-res gsea_data.gct \
-cls gsea.cls \
-gmx /mnt/data/user_data/xiangyu/programme/gsea/msigdb_v6.1_GMTs/c2.all.v6.1.symbols.gmt \
-chip /mnt/data/user_data/xiangyu/programme/gsea/chip/ENSEMBL_mouse_gene.chip \
-collapse true -mode Max_probe -norm meandiv -nperm 1000 -permute gene_set \
-rnd_type no_balance -scoring_scheme weighted -rpt_label my_analysis \
-metric Signal2Noise -sort real -order descending -include_only_symbols true \
-make_sets true -median false -num 100 -plot_top_x 50 -rnd_seed timestamp \
-save_rnd_lists false -set_max 500 -set_min 15 -zip_report false \
-out ./Lung_LMP1_vs_ctrl_c2 -gui false

java -cp /mnt/data/user_data/xiangyu/programme/gsea/gsea-3.0.jar -Xmx1024m xtools.gsea.Gsea \
-res gsea_data.gct \
-cls gsea.cls \
-gmx /mnt/data/user_data/xiangyu/programme/gsea/msigdb_v6.1_GMTs/c3.all.v6.1.symbols.gmt \
-chip /mnt/data/user_data/xiangyu/programme/gsea/chip/ENSEMBL_mouse_gene.chip \
-collapse true -mode Max_probe -norm meandiv -nperm 1000 -permute gene_set \
-rnd_type no_balance -scoring_scheme weighted -rpt_label my_analysis \
-metric Signal2Noise -sort real -order descending -include_only_symbols true \
-make_sets true -median false -num 100 -plot_top_x 50 -rnd_seed timestamp \
-save_rnd_lists false -set_max 500 -set_min 15 -zip_report false \
-out ./Lung_LMP1_vs_ctrl_c3 -gui false

java -cp /mnt/data/user_data/xiangyu/programme/gsea/gsea-3.0.jar -Xmx1024m xtools.gsea.Gsea \
-res gsea_data.gct \
-cls gsea.cls \
-gmx /mnt/data/user_data/xiangyu/programme/gsea/msigdb_v6.1_GMTs/c4.all.v6.1.symbols.gmt \
-chip /mnt/data/user_data/xiangyu/programme/gsea/chip/ENSEMBL_mouse_gene.chip \
-collapse true -mode Max_probe -norm meandiv -nperm 1000 -permute gene_set \
-rnd_type no_balance -scoring_scheme weighted -rpt_label my_analysis \
-metric Signal2Noise -sort real -order descending -include_only_symbols true \
-make_sets true -median false -num 100 -plot_top_x 50 -rnd_seed timestamp \
-save_rnd_lists false -set_max 500 -set_min 15 -zip_report false \
-out ./Lung_LMP1_vs_ctrl_c4 -gui false

java -cp /mnt/data/user_data/xiangyu/programme/gsea/gsea-3.0.jar -Xmx1024m xtools.gsea.Gsea \
-res gsea_data.gct \
-cls gsea.cls \
-gmx /mnt/data/user_data/xiangyu/programme/gsea/msigdb_v6.1_GMTs/c5.all.v6.1.symbols.gmt \
-chip /mnt/data/user_data/xiangyu/programme/gsea/chip/ENSEMBL_mouse_gene.chip \
-collapse true -mode Max_probe -norm meandiv -nperm 1000 -permute gene_set \
-rnd_type no_balance -scoring_scheme weighted -rpt_label my_analysis \
-metric Signal2Noise -sort real -order descending -include_only_symbols true \
-make_sets true -median false -num 100 -plot_top_x 50 -rnd_seed timestamp \
-save_rnd_lists false -set_max 500 -set_min 15 -zip_report false \
-out ./Lung_LMP1_vs_ctrl_c5 -gui false

java -cp /mnt/data/user_data/xiangyu/programme/gsea/gsea-3.0.jar -Xmx1024m xtools.gsea.Gsea \
-res gsea_data.gct \
-cls gsea.cls \
-gmx /mnt/data/user_data/xiangyu/programme/gsea/msigdb_v6.1_GMTs/c6.all.v6.1.symbols.gmt \
-chip /mnt/data/user_data/xiangyu/programme/gsea/chip/ENSEMBL_mouse_gene.chip \
-collapse true -mode Max_probe -norm meandiv -nperm 1000 -permute gene_set \
-rnd_type no_balance -scoring_scheme weighted -rpt_label my_analysis \
-metric Signal2Noise -sort real -order descending -include_only_symbols true \
-make_sets true -median false -num 100 -plot_top_x 50 -rnd_seed timestamp \
-save_rnd_lists false -set_max 500 -set_min 15 -zip_report false \
-out ./Lung_LMP1_vs_ctrl_c6 -gui false

java -cp /mnt/data/user_data/xiangyu/programme/gsea/gsea-3.0.jar -Xmx1024m xtools.gsea.Gsea \
-res gsea_data.gct \
-cls gsea.cls \
-gmx /mnt/data/user_data/xiangyu/programme/gsea/msigdb_v6.1_GMTs/c7.all.v6.1.symbols.gmt \
-chip /mnt/data/user_data/xiangyu/programme/gsea/chip/ENSEMBL_mouse_gene.chip \
-collapse true -mode Max_probe -norm meandiv -nperm 1000 -permute gene_set \
-rnd_type no_balance -scoring_scheme weighted -rpt_label my_analysis \
-metric Signal2Noise -sort real -order descending -include_only_symbols true \
-make_sets true -median false -num 100 -plot_top_x 50 -rnd_seed timestamp \
-save_rnd_lists false -set_max 500 -set_min 15 -zip_report false \
-out ./Lung_LMP1_vs_ctrl_c7 -gui false


bash gsea_shell.sh

****************************************************TMPL_metastasis vs TMPL_primary********************************************

############################
#move to bamfiles path  #####
#move to bamfiles path  #####
#############################
cd /mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/star_mm10_out


vi sample_sampletable_LMP1.csv

ids,sample,source,order
1,LMP1_2.Aligned.sortedByCoord.out.bam,LMP1_2,primary,1
2,LMP1_3.Aligned.sortedByCoord.out.bam,LMP1_3,primary,2
3,Lung_LMP1_1.Aligned.sortedByCoord.out.bam,Lung_LMP1_1,metastasis,3
4,Lung_LMP1_2.Aligned.sortedByCoord.out.bam,Lung_LMP1_2,metastasis,4
5,Lung_LMP1_3.Aligned.sortedByCoord.out.bam,Lung_LMP1_3,metastasis,5

############################
#work in R environment  #####
#work in R environment  #####
#############################


project <- c("LMP1_metastasis_vs_primary")
deal_design <- c("metastasis","primary")
significant_cutoff <- c(1)
organism <- "mouse"
##物种为人就是hsa，是老鼠就是mouse

file_path <- "/mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/LMP1_metastasis_vs_primary/"
sample_sampletable.path <- "/mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/star_mm10_out"
load(file="/mnt/data/user_data/xiangyu/workshop/WORKFLOW_RNAseq/ebg_mm10.RData")
load(file="/mnt/data/user_data/xiangyu/workshop/WORKFLOW_RNAseq/txdb_mm10.RData")  


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
res_1_filter <- c("source",deal_design)
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
csvfile <- file.path(indir, "sample_sampletable_LMP1.csv")
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
getwd()

colData(se) <- DataFrame(sampleTable)
dds <- DESeqDataSet(se, design = ~source)
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
library("ggplot2")
pcaData <- plotPCA(rld, intgroup = c("sample","order"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(x = PC1, y = PC2, color = sample)) +
geom_point(size =1) +
geom_text(aes(label=order))+
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

res_1$ensembl <- mapIds(x = anno_data,
                        keys = rownames(res_1),
						keytype ="SYMBOL",
						column ="ENSEMBL",
						multiVals="first")
res_1$entrez <- mapIds(x = anno_data,
                        keys = rownames(res_1),
						keytype ="SYMBOL",
						column ="ENTREZID",
						multiVals="first")
AA <- rownames(res_1)
AA <- as.character(AA)
res_1$GENENAME <- mapIds(x = anno_data,
                        keys = AA,
						keytype ="SYMBOL",
						column ="GENENAME",
						multiVals="first")
write.csv(res_1, file = res_1_file_sym_name)
all_summry <- cbind(count_and_tpm,res_1)
write.csv(all_summry, file = ALL.CSV_FILE)

res_1 <- na.omit(res_1)
y <- res_1[c(1:nrow(res_1)),6]
res_1 <- res_1[with(res_1,y<0.05),]
y <- res_1[c(1:nrow(res_1)),2]

upres_1 <- res_1[with(res_1,y>=significant_cutoff),]
ee	<-as.matrix(upres_1$entrez)
	dd <- as.vector(ee)
KEGGupres_1 <- enrichKEGG(gene =dd, 
					organism = organism, 
					keyType = "ncbi-geneid",
					 pvalueCutoff = 0.05,
				       pAdjustMethod = "BH", 
				       minGSSize = 10, 
				       maxGSSize = 500,
				       qvalueCutoff = 0.2, 
				       use_internal_data = FALSE)

downres_1 <- res_1[with(res_1,y<=-significant_cutoff),]
hh	<-as.matrix(downres_1$entrez)
	ii <- as.vector(hh)
KEGGdownres_1 <- enrichKEGG(gene =ii, 
					organism = organism, 
					keyType = "ncbi-geneid",
					 pvalueCutoff = 0.05,
				       pAdjustMethod = "BH", 
				       minGSSize = 10, 
				       maxGSSize = 500,
				       qvalueCutoff = 0.2, 
				       use_internal_data = FALSE)

write.csv(KEGGupres_1, file = KEGGupres_1_file.csv)
write.csv(KEGGdownres_1, file = KEGGdownres_1_file.csv)

pdf(file = KEGGupres_1_pdf)
dotplot(KEGGupres_1, showCategory=20)
dev.off()
pdf(file = KEGGdownres_1_pdf)
dotplot(KEGGdownres_1, showCategory=20)
dev.off()

pdf(file = KEGGres_1_all_pdf)
KEGGupres_1uporder <- KEGGupres_1[order((KEGGupres_1$Count),decreasing = TRUE),]
KEGGdownres_1downorder <- KEGGdownres_1[order((KEGGdownres_1$Count),decreasing = TRUE),]
down_kegg<-KEGGdownres_1downorder[KEGGdownres_1downorder$p.adjust<0.05,];
down_kegg$group=-1
up_kegg<-KEGGupres_1uporder[KEGGupres_1uporder$p.adjust<0.05,];
up_kegg$group=1
down_kegg1 <- down_kegg[-(21:nrow(down_kegg)),]
up_kegg1 <- up_kegg[-(21:nrow(up_kegg)),]
  dat=rbind(up_kegg1,down_kegg1)
  colnames(dat)
    dat$p.adjust = -log10(dat$p.adjust)
  dat$p.adjust=dat$p.adjust*dat$group 
  dat=dat[order(dat$p.adjust,decreasing = F),]
  library("ggplot2")
  g_kegg<- ggplot(dat, aes(x=reorder(Description,order(p.adjust, decreasing = F)), y=p.adjust, fill=group)) + 
    geom_bar(stat="identity") + 
    scale_fill_gradient(low="blue",high="red",guide = FALSE) + 
    scale_x_discrete(name ="Pathway names") +
    scale_y_continuous(name ="-log10p_adjust") +
    coord_flip() +
    ggtitle(KEGGres_1_all_title)
  print(g_kegg)
dev.off()

ee	<-as.matrix(upres_1$entrez)
	dd <- as.vector(ee)
GOupres_1_all <- enrichGO(gene = dd, 
	           OrgDb = anno_data,
				ont = "all", 
		             pvalueCutoff = 0.01, 
                     pAdjustMethod = "BH", 
                     qvalueCutoff = 0.05,
                     minGSSize = 10, 
                     maxGSSize = 500, 
                     readable = T, 
                     pool = FALSE)
write.csv(GOupres_1_all,file=GOres_1_all_UP_csv)

ff <- barplot(GOupres_1_all, split="ONTOLOGY",showCategory=15) + facet_grid(ONTOLOGY~., scale="free")
ff1 <- ff + theme(axis.text.y = element_text(size = 8)) +labs(title = NULL)
pdf(file = GOres_1_all_UP_pdf)
ff1
dev.off()

hh	<-as.matrix(downres_1$entrez)
	ii <- as.vector(hh)
GOdownres_1_all <- enrichGO(gene = ii, 
	           OrgDb = anno_data,
				ont = "all", 
		             pvalueCutoff = 0.01, 
                     pAdjustMethod = "BH", 
                     qvalueCutoff = 0.05,
                     minGSSize = 10, 
                     maxGSSize = 500, 
                     readable = T, 
                     pool = FALSE)
write.csv(GOdownres_1_all,file=GOres_1_all_DOWN_csv)


ff <- barplot(GOdownres_1_all, split="ONTOLOGY",showCategory=15) + facet_grid(ONTOLOGY~., scale="free")
ff1 <- ff + theme(axis.text.y = element_text(size = 8)) +labs(title = NULL)
pdf(file = GOres_1_all_DOWN_pdf)
ff1
dev.off()


*********gsea输入文件准备***********
gsea_data <- read.csv("/mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/LMP1_metastasis_vs_primary/LMP1_metastasis_vs_primary_05_metastasis_v_primary_count_tpm_symbol_and_anno.csv")
gsea_data$Description <- c("NULL")
head(gsea_data)
gsea_data <- gsea_data[,c(18,21,9,10,11,7,8)]
nrow(gsea_data)
gsea_data <- na.omit(gsea_data)
nrow(gsea_data)
colnames(gsea_data)
colnames(gsea_data) <- c("NAME","Description","Lung_LMP1_1","Lung_LMP1_2","Lung_LMP1_3","LMP1_2","LMP1_3")
write.csv(gsea_data,file="/mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/LMP1_metastasis_vs_primary/gsea_data.csv")
#在excel中改好格式(就以Ensembl id为行名)，然后直接保存为txt格式，再导入R修改文件名为.gct即可
#加上  
       #1.2
       #总基因行数 总样本数

vi gsea.cls
5 2 1
# Lung_LMP1 LMP1
0 0 0 1 1



#样本多的时候，-Xmx1024m改为-Xmx10240 \       9022      注意\后面不能有空格
vi gsea_shell.sh

java -cp /mnt/data/user_data/xiangyu/programme/gsea/gsea-3.0.jar -Xmx1024m xtools.gsea.Gsea \
-res gsea_data.gct \
-cls gsea.cls \
-gmx /mnt/data/user_data/xiangyu/programme/gsea/msigdb_v6.1_GMTs/h.all.v6.1.symbols.gmt \
-chip /mnt/data/user_data/xiangyu/programme/gsea/chip/ENSEMBL_mouse_gene.chip \
-collapse true -mode Max_probe -norm meandiv -nperm 1000 -permute gene_set \
-rnd_type no_balance -scoring_scheme weighted -rpt_label my_analysis \
-metric Signal2Noise -sort real -order descending -include_only_symbols true \
-make_sets true -median false -num 100 -plot_top_x 50 -rnd_seed timestamp \
-save_rnd_lists false -set_max 500 -set_min 15 -zip_report false \
-out ./GSEA_out/LMP1_metastasis_vs_primary_h -gui false

java -cp /mnt/data/user_data/xiangyu/programme/gsea/gsea-3.0.jar -Xmx1024m xtools.gsea.Gsea \
-res gsea_data.gct \
-cls gsea.cls \
-gmx /mnt/data/user_data/xiangyu/programme/gsea/msigdb_v6.1_GMTs/c1.all.v6.1.symbols.gmt \
-chip /mnt/data/user_data/xiangyu/programme/gsea/chip/ENSEMBL_mouse_gene.chip \
-collapse true -mode Max_probe -norm meandiv -nperm 1000 -permute gene_set \
-rnd_type no_balance -scoring_scheme weighted -rpt_label my_analysis \
-metric Signal2Noise -sort real -order descending -include_only_symbols true \
-make_sets true -median false -num 100 -plot_top_x 50 -rnd_seed timestamp \
-save_rnd_lists false -set_max 500 -set_min 15 -zip_report false \
-out ./GSEA_out/LMP1_metastasis_vs_primary_c1 -gui false

java -cp /mnt/data/user_data/xiangyu/programme/gsea/gsea-3.0.jar -Xmx1024m xtools.gsea.Gsea \
-res gsea_data.gct \
-cls gsea.cls \
-gmx /mnt/data/user_data/xiangyu/programme/gsea/msigdb_v6.1_GMTs/c2.all.v6.1.symbols.gmt \
-chip /mnt/data/user_data/xiangyu/programme/gsea/chip/ENSEMBL_mouse_gene.chip \
-collapse true -mode Max_probe -norm meandiv -nperm 1000 -permute gene_set \
-rnd_type no_balance -scoring_scheme weighted -rpt_label my_analysis \
-metric Signal2Noise -sort real -order descending -include_only_symbols true \
-make_sets true -median false -num 100 -plot_top_x 50 -rnd_seed timestamp \
-save_rnd_lists false -set_max 500 -set_min 15 -zip_report false \
-out ./GSEA_out/LMP1_metastasis_vs_primary_c2 -gui false

java -cp /mnt/data/user_data/xiangyu/programme/gsea/gsea-3.0.jar -Xmx1024m xtools.gsea.Gsea \
-res gsea_data.gct \
-cls gsea.cls \
-gmx /mnt/data/user_data/xiangyu/programme/gsea/msigdb_v6.1_GMTs/c3.all.v6.1.symbols.gmt \
-chip /mnt/data/user_data/xiangyu/programme/gsea/chip/ENSEMBL_mouse_gene.chip \
-collapse true -mode Max_probe -norm meandiv -nperm 1000 -permute gene_set \
-rnd_type no_balance -scoring_scheme weighted -rpt_label my_analysis \
-metric Signal2Noise -sort real -order descending -include_only_symbols true \
-make_sets true -median false -num 100 -plot_top_x 50 -rnd_seed timestamp \
-save_rnd_lists false -set_max 500 -set_min 15 -zip_report false \
-out ./GSEA_out/LMP1_metastasis_vs_primary_c3 -gui false

java -cp /mnt/data/user_data/xiangyu/programme/gsea/gsea-3.0.jar -Xmx1024m xtools.gsea.Gsea \
-res gsea_data.gct \
-cls gsea.cls \
-gmx /mnt/data/user_data/xiangyu/programme/gsea/msigdb_v6.1_GMTs/c4.all.v6.1.symbols.gmt \
-chip /mnt/data/user_data/xiangyu/programme/gsea/chip/ENSEMBL_mouse_gene.chip \
-collapse true -mode Max_probe -norm meandiv -nperm 1000 -permute gene_set \
-rnd_type no_balance -scoring_scheme weighted -rpt_label my_analysis \
-metric Signal2Noise -sort real -order descending -include_only_symbols true \
-make_sets true -median false -num 100 -plot_top_x 50 -rnd_seed timestamp \
-save_rnd_lists false -set_max 500 -set_min 15 -zip_report false \
-out ./GSEA_out/LMP1_metastasis_vs_primary_c4 -gui false

java -cp /mnt/data/user_data/xiangyu/programme/gsea/gsea-3.0.jar -Xmx1024m xtools.gsea.Gsea \
-res gsea_data.gct \
-cls gsea.cls \
-gmx /mnt/data/user_data/xiangyu/programme/gsea/msigdb_v6.1_GMTs/c5.all.v6.1.symbols.gmt \
-chip /mnt/data/user_data/xiangyu/programme/gsea/chip/ENSEMBL_mouse_gene.chip \
-collapse true -mode Max_probe -norm meandiv -nperm 1000 -permute gene_set \
-rnd_type no_balance -scoring_scheme weighted -rpt_label my_analysis \
-metric Signal2Noise -sort real -order descending -include_only_symbols true \
-make_sets true -median false -num 100 -plot_top_x 50 -rnd_seed timestamp \
-save_rnd_lists false -set_max 500 -set_min 15 -zip_report false \
-out ./GSEA_out/LMP1_metastasis_vs_primary_c5 -gui false

java -cp /mnt/data/user_data/xiangyu/programme/gsea/gsea-3.0.jar -Xmx1024m xtools.gsea.Gsea \
-res gsea_data.gct \
-cls gsea.cls \
-gmx /mnt/data/user_data/xiangyu/programme/gsea/msigdb_v6.1_GMTs/c6.all.v6.1.symbols.gmt \
-chip /mnt/data/user_data/xiangyu/programme/gsea/chip/ENSEMBL_mouse_gene.chip \
-collapse true -mode Max_probe -norm meandiv -nperm 1000 -permute gene_set \
-rnd_type no_balance -scoring_scheme weighted -rpt_label my_analysis \
-metric Signal2Noise -sort real -order descending -include_only_symbols true \
-make_sets true -median false -num 100 -plot_top_x 50 -rnd_seed timestamp \
-save_rnd_lists false -set_max 500 -set_min 15 -zip_report false \
-out ./GSEA_out/LMP1_metastasis_vs_primary_c6 -gui false

java -cp /mnt/data/user_data/xiangyu/programme/gsea/gsea-3.0.jar -Xmx1024m xtools.gsea.Gsea \
-res gsea_data.gct \
-cls gsea.cls \
-gmx /mnt/data/user_data/xiangyu/programme/gsea/msigdb_v6.1_GMTs/c7.all.v6.1.symbols.gmt \
-chip /mnt/data/user_data/xiangyu/programme/gsea/chip/ENSEMBL_mouse_gene.chip \
-collapse true -mode Max_probe -norm meandiv -nperm 1000 -permute gene_set \
-rnd_type no_balance -scoring_scheme weighted -rpt_label my_analysis \
-metric Signal2Noise -sort real -order descending -include_only_symbols true \
-make_sets true -median false -num 100 -plot_top_x 50 -rnd_seed timestamp \
-save_rnd_lists false -set_max 500 -set_min 15 -zip_report false \
-out ./GSEA_out/LMP1_metastasis_vs_primary_c7 -gui false


bash gsea_shell.sh


****************************************************TMP_metastasis vs TMP_primary****************************************


############################
#move to bamfiles path  #####
#move to bamfiles path  #####
#############################
cd /mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/star_mm10_out

vi sample_sampletable_ctrl.csv

ids,sample,source,order
1,ctrl_1.Aligned.sortedByCoord.out.bam,ctrl_1,primary,1
2,ctrl_2.Aligned.sortedByCoord.out.bam,ctrl_2,primary,2
3,ctrl_3.Aligned.sortedByCoord.out.bam,ctrl_3,primary,3
4,Lung_ctrl_1.Aligned.sortedByCoord.out.bam,Lung_ctrl_1,metastasis,4
5,Lung_ctrl_2.Aligned.sortedByCoord.out.bam,Lung_ctrl_2,metastasis,5



############################
#work in R environment  #####
#work in R environment  #####
#############################


project <- c("ctrl_metastasis_vs_primary")
deal_design <- c("metastasis","primary")
significant_cutoff <- c(0.5)
organism <- "mouse"
##物种为人就是hsa，是老鼠就是mouse

file_path <- "/mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/ctrl_metastasis_vs_primary/"
sample_sampletable.path <- "/mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/star_mm10_out"
load(file="/mnt/data/user_data/xiangyu/workshop/WORKFLOW_RNAseq/ebg_mm10.RData")
load(file="/mnt/data/user_data/xiangyu/workshop/WORKFLOW_RNAseq/txdb_mm10.RData")  


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
res_1_filter <- c("source",deal_design)
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
csvfile <- file.path(indir, "sample_sampletable_ctrl.csv")
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
getwd()

colData(se) <- DataFrame(sampleTable)
dds <- DESeqDataSet(se, design = ~source)
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
library("ggplot2")
pcaData <- plotPCA(rld, intgroup = c("sample","order"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(x = PC1, y = PC2, color = sample)) +
geom_point(size =1) +
geom_text(aes(label=order))+
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

res_1$ensembl <- mapIds(x = anno_data,
                        keys = rownames(res_1),
						keytype ="SYMBOL",
						column ="ENSEMBL",
						multiVals="first")
res_1$entrez <- mapIds(x = anno_data,
                        keys = rownames(res_1),
						keytype ="SYMBOL",
						column ="ENTREZID",
						multiVals="first")
AA <- rownames(res_1)
AA <- as.character(AA)
res_1$GENENAME <- mapIds(x = anno_data,
                        keys = AA,
						keytype ="SYMBOL",
						column ="GENENAME",
						multiVals="first")
write.csv(res_1, file = res_1_file_sym_name)
all_summry <- cbind(count_and_tpm,res_1)
write.csv(all_summry, file = ALL.CSV_FILE)

res_1 <- na.omit(res_1)
y <- res_1[c(1:nrow(res_1)),5]
res_1 <- res_1[with(res_1,y<0.05),]
y <- res_1[c(1:nrow(res_1)),2]

upres_1 <- res_1[with(res_1,y>=significant_cutoff),]
ee	<-as.matrix(upres_1$entrez)
	dd <- as.vector(ee)
KEGGupres_1 <- enrichKEGG(gene =dd, 
					organism = organism, 
					keyType = "ncbi-geneid",
					 pvalueCutoff = 0.05,
				       pAdjustMethod = "BH", 
				       minGSSize = 10, 
				       maxGSSize = 500,
				       qvalueCutoff = 0.2, 
				       use_internal_data = FALSE)

downres_1 <- res_1[with(res_1,y<=-significant_cutoff),]
hh	<-as.matrix(downres_1$entrez)
	ii <- as.vector(hh)
KEGGdownres_1 <- enrichKEGG(gene =ii, 
					organism = organism, 
					keyType = "ncbi-geneid",
					 pvalueCutoff = 0.05,
				       pAdjustMethod = "BH", 
				       minGSSize = 10, 
				       maxGSSize = 500,
				       qvalueCutoff = 0.2, 
				       use_internal_data = FALSE)

write.csv(KEGGupres_1, file = KEGGupres_1_file.csv)
write.csv(KEGGdownres_1, file = KEGGdownres_1_file.csv)


pdf(file = KEGGupres_1_pdf)
dotplot(KEGGupres_1, showCategory=20)
dev.off()
pdf(file = KEGGdownres_1_pdf)
dotplot(KEGGdownres_1, showCategory=20)
dev.off()

pdf(file = KEGGres_1_all_pdf)q
KEGGupres_1uporder <- KEGGupres_1[order((KEGGupres_1$Count),decreasing = TRUE),]
KEGGdownres_1downorder <- KEGGdownres_1[order((KEGGdownres_1$Count),decreasing = TRUE),]
down_kegg<-KEGGdownres_1downorder[KEGGdownres_1downorder$pvalue<0.05,];
down_kegg$group=-1
up_kegg<-KEGGupres_1uporder[KEGGupres_1uporder$pvalue<0.05,];
up_kegg$group=1
down_kegg1 <- down_kegg[-(21:nrow(down_kegg)),]
up_kegg1 <- up_kegg[-(21:nrow(up_kegg)),]
  dat=rbind(up_kegg1,down_kegg1)
  colnames(dat)
  dat$pvalue = -log10(dat$pvalue)
  dat$pvalue=dat$pvalue*dat$group 
  dat=dat[order(dat$pvalue,decreasing = F),]
  library("ggplot2")
  g_kegg<- ggplot(dat, aes(x=reorder(Description,order(pvalue, decreasing = F)), y=pvalue, fill=group)) + 
    geom_bar(stat="identity") + 
    scale_fill_gradient(low="blue",high="red",guide = FALSE) + 
    scale_x_discrete(name ="Pathway names") +
    scale_y_continuous(name ="-log10pvalue") +
    coord_flip() +
    ggtitle(KEGGres_1_all_title)
  print(g_kegg)
dev.off()


ee	<-as.matrix(upres_1$entrez)
	dd <- as.vector(ee)
GOupres_1_all <- enrichGO(gene = dd, 
	           OrgDb = anno_data,
				ont = "all", 
		             pvalueCutoff = 0.01, 
                     pAdjustMethod = "BH", 
                     qvalueCutoff = 0.05,
                     minGSSize = 10, 
                     maxGSSize = 500, 
                     readable = T, 
                     pool = FALSE)
write.csv(GOupres_1_all,file=GOres_1_all_UP_csv)

ff <- barplot(GOupres_1_all, split="ONTOLOGY",showCategory=15) + facet_grid(ONTOLOGY~., scale="free")
ff1 <- ff + theme(axis.text.y = element_text(size = 8)) +labs(title = NULL)
pdf(file = GOres_1_all_UP_pdf)
ff1
dev.off()

hh	<-as.matrix(downres_1$entrez)
	ii <- as.vector(hh)
GOdownres_1_all <- enrichGO(gene = ii, 
	           OrgDb = anno_data,
				ont = "all", 
		             pvalueCutoff = 0.01, 
                     pAdjustMethod = "BH", 
                     qvalueCutoff = 0.05,
                     minGSSize = 10, 
                     maxGSSize = 500, 
                     readable = T, 
                     pool = FALSE)
write.csv(GOdownres_1_all,file=GOres_1_all_DOWN_csv)

ff <- barplot(GOdownres_1_all, split="ONTOLOGY",showCategory=15) + facet_grid(ONTOLOGY~., scale="free")
ff1 <- ff + theme(axis.text.y = element_text(size = 8)) +labs(title = NULL)
pdf(file = GOres_1_all_DOWN_pdf)
ff1
dev.off()




*********gsea输入文件准备***********
gsea_data <- read.csv("/mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/ctrl_metastasis_vs_primary/ctrl_metastasis_vs_primary_05_metastasis_v_primary_count_tpm_symbol_and_anno.csv")
gsea_data$Description <- c("NULL")
gsea_data <- gsea_data[,c(18,21,10,11,7,8,9)]
nrow(gsea_data)
gsea_data <- na.omit(gsea_data)
nrow(gsea_data)
colnames(gsea_data)
colnames(gsea_data) <- c("NAME","Description","Lung_ctrl_1","Lung_ctrl_2","ctrl_1","ctrl_2","ctrl_3")
write.csv(gsea_data,file="/mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/ctrl_metastasis_vs_primary/gsea_data.csv")
#在excel中改好格式(就以Ensembl id为行名)，然后直接保存为txt格式，再导入R修改文件名为.gct即可
#加上  
       #1.2
       #总基因行数 总样本数

vi gsea.cls
5 2 1
# Lung_ctrl ctrl
0 0 1 1 1



#样本多的时候，-Xmx1024m改为-Xmx10240 \       9022      注意\后面不能有空格
vi gsea_shell.sh

java -cp /mnt/data/user_data/xiangyu/programme/gsea/gsea-3.0.jar -Xmx1024m xtools.gsea.Gsea \
-res gsea_data.gct \
-cls gsea.cls \
-gmx /mnt/data/user_data/xiangyu/programme/gsea/msigdb_v6.1_GMTs/h.all.v6.1.symbols.gmt \
-chip /mnt/data/user_data/xiangyu/programme/gsea/chip/ENSEMBL_mouse_gene.chip \
-collapse true -mode Max_probe -norm meandiv -nperm 1000 -permute gene_set \
-rnd_type no_balance -scoring_scheme weighted -rpt_label my_analysis \
-metric Signal2Noise -sort real -order descending -include_only_symbols true \
-make_sets true -median false -num 100 -plot_top_x 50 -rnd_seed timestamp \
-save_rnd_lists false -set_max 500 -set_min 15 -zip_report false \
-out ./ctrl_metastasis_vs_primary_h -gui false

java -cp /mnt/data/user_data/xiangyu/programme/gsea/gsea-3.0.jar -Xmx1024m xtools.gsea.Gsea \
-res gsea_data.gct \
-cls gsea.cls \
-gmx /mnt/data/user_data/xiangyu/programme/gsea/msigdb_v6.1_GMTs/c1.all.v6.1.symbols.gmt \
-chip /mnt/data/user_data/xiangyu/programme/gsea/chip/ENSEMBL_mouse_gene.chip \
-collapse true -mode Max_probe -norm meandiv -nperm 1000 -permute gene_set \
-rnd_type no_balance -scoring_scheme weighted -rpt_label my_analysis \
-metric Signal2Noise -sort real -order descending -include_only_symbols true \
-make_sets true -median false -num 100 -plot_top_x 50 -rnd_seed timestamp \
-save_rnd_lists false -set_max 500 -set_min 15 -zip_report false \
-out ./ctrl_metastasis_vs_primary_c1 -gui false

java -cp /mnt/data/user_data/xiangyu/programme/gsea/gsea-3.0.jar -Xmx1024m xtools.gsea.Gsea \
-res gsea_data.gct \
-cls gsea.cls \
-gmx /mnt/data/user_data/xiangyu/programme/gsea/msigdb_v6.1_GMTs/c2.all.v6.1.symbols.gmt \
-chip /mnt/data/user_data/xiangyu/programme/gsea/chip/ENSEMBL_mouse_gene.chip \
-collapse true -mode Max_probe -norm meandiv -nperm 1000 -permute gene_set \
-rnd_type no_balance -scoring_scheme weighted -rpt_label my_analysis \
-metric Signal2Noise -sort real -order descending -include_only_symbols true \
-make_sets true -median false -num 100 -plot_top_x 50 -rnd_seed timestamp \
-save_rnd_lists false -set_max 500 -set_min 15 -zip_report false \
-out ./ctrl_metastasis_vs_primary_c2 -gui false

java -cp /mnt/data/user_data/xiangyu/programme/gsea/gsea-3.0.jar -Xmx1024m xtools.gsea.Gsea \
-res gsea_data.gct \
-cls gsea.cls \
-gmx /mnt/data/user_data/xiangyu/programme/gsea/msigdb_v6.1_GMTs/c3.all.v6.1.symbols.gmt \
-chip /mnt/data/user_data/xiangyu/programme/gsea/chip/ENSEMBL_mouse_gene.chip \
-collapse true -mode Max_probe -norm meandiv -nperm 1000 -permute gene_set \
-rnd_type no_balance -scoring_scheme weighted -rpt_label my_analysis \
-metric Signal2Noise -sort real -order descending -include_only_symbols true \
-make_sets true -median false -num 100 -plot_top_x 50 -rnd_seed timestamp \
-save_rnd_lists false -set_max 500 -set_min 15 -zip_report false \
-out ./ctrl_metastasis_vs_primary_c3 -gui false

java -cp /mnt/data/user_data/xiangyu/programme/gsea/gsea-3.0.jar -Xmx1024m xtools.gsea.Gsea \
-res gsea_data.gct \
-cls gsea.cls \
-gmx /mnt/data/user_data/xiangyu/programme/gsea/msigdb_v6.1_GMTs/c4.all.v6.1.symbols.gmt \
-chip /mnt/data/user_data/xiangyu/programme/gsea/chip/ENSEMBL_mouse_gene.chip \
-collapse true -mode Max_probe -norm meandiv -nperm 1000 -permute gene_set \
-rnd_type no_balance -scoring_scheme weighted -rpt_label my_analysis \
-metric Signal2Noise -sort real -order descending -include_only_symbols true \
-make_sets true -median false -num 100 -plot_top_x 50 -rnd_seed timestamp \
-save_rnd_lists false -set_max 500 -set_min 15 -zip_report false \
-out ./ctrl_metastasis_vs_primary_c4 -gui false

java -cp /mnt/data/user_data/xiangyu/programme/gsea/gsea-3.0.jar -Xmx1024m xtools.gsea.Gsea \
-res gsea_data.gct \
-cls gsea.cls \
-gmx /mnt/data/user_data/xiangyu/programme/gsea/msigdb_v6.1_GMTs/c5.all.v6.1.symbols.gmt \
-chip /mnt/data/user_data/xiangyu/programme/gsea/chip/ENSEMBL_mouse_gene.chip \
-collapse true -mode Max_probe -norm meandiv -nperm 1000 -permute gene_set \
-rnd_type no_balance -scoring_scheme weighted -rpt_label my_analysis \
-metric Signal2Noise -sort real -order descending -include_only_symbols true \
-make_sets true -median false -num 100 -plot_top_x 50 -rnd_seed timestamp \
-save_rnd_lists false -set_max 500 -set_min 15 -zip_report false \
-out ./ctrl_metastasis_vs_primary_c5 -gui false

java -cp /mnt/data/user_data/xiangyu/programme/gsea/gsea-3.0.jar -Xmx1024m xtools.gsea.Gsea \
-res gsea_data.gct \
-cls gsea.cls \
-gmx /mnt/data/user_data/xiangyu/programme/gsea/msigdb_v6.1_GMTs/c6.all.v6.1.symbols.gmt \
-chip /mnt/data/user_data/xiangyu/programme/gsea/chip/ENSEMBL_mouse_gene.chip \
-collapse true -mode Max_probe -norm meandiv -nperm 1000 -permute gene_set \
-rnd_type no_balance -scoring_scheme weighted -rpt_label my_analysis \
-metric Signal2Noise -sort real -order descending -include_only_symbols true \
-make_sets true -median false -num 100 -plot_top_x 50 -rnd_seed timestamp \
-save_rnd_lists false -set_max 500 -set_min 15 -zip_report false \
-out ./ctrl_metastasis_vs_primary_c6 -gui false

java -cp /mnt/data/user_data/xiangyu/programme/gsea/gsea-3.0.jar -Xmx1024m xtools.gsea.Gsea \
-res gsea_data.gct \
-cls gsea.cls \
-gmx /mnt/data/user_data/xiangyu/programme/gsea/msigdb_v6.1_GMTs/c7.all.v6.1.symbols.gmt \
-chip /mnt/data/user_data/xiangyu/programme/gsea/chip/ENSEMBL_mouse_gene.chip \
-collapse true -mode Max_probe -norm meandiv -nperm 1000 -permute gene_set \
-rnd_type no_balance -scoring_scheme weighted -rpt_label my_analysis \
-metric Signal2Noise -sort real -order descending -include_only_symbols true \
-make_sets true -median false -num 100 -plot_top_x 50 -rnd_seed timestamp \
-save_rnd_lists false -set_max 500 -set_min 15 -zip_report false \
-out ./ctrl_metastasis_vs_primary_c7 -gui false

bash gsea_shell.sh

java -cp /mnt/data/user_data/xiangyu/programme/gsea/gsea-3.0.jar -Xmx1024m xtools.gsea.Gsea \
-res gsea_data.gct \
-cls gsea.cls \
-gmx /mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/metastasis/GSEA_gmt/metastasis.gmt \
-chip /mnt/data/user_data/xiangyu/programme/gsea/chip/ENSEMBL_mouse_gene.chip \
-collapse true -mode Max_probe -norm meandiv -nperm 1000 -permute gene_set \
-rnd_type no_balance -scoring_scheme weighted -rpt_label my_analysis \
-metric Signal2Noise -sort real -order descending -include_only_symbols true \
-make_sets true -median false -num 100 -plot_top_x 50 -rnd_seed timestamp \
-save_rnd_lists false -set_max 500 -set_min 15 -zip_report false \
-out ./ctrl_metastasis_vs_primary_metastasis -gui false


********************************************************primary_TMP vs normal******************************
vi sample_sampletable_ctrl_normal_organoid.csv

ids,sample,deal,order
1,ctrl_1.Aligned.sortedByCoord.out.bam,ctrl_1,P_ctrl,1
2,ctrl_2.Aligned.sortedByCoord.out.bam,ctrl_2,P_ctrl,2
3,ctrl_3.Aligned.sortedByCoord.out.bam,ctrl_3,P_ctrl,3
4,normal_organoid_1.Aligned.sortedByCoord.out.bam,normal_organoid_1,normal_organoid,4
5,normal_organoid_2.Aligned.sortedByCoord.out.bam,normal_organoid_2,normal_organoid,5
6,normal_organoid_3.Aligned.sortedByCoord.out.bam,normal_organoid_3,normal_organoid,6

############################
#work in R environment  #####
#work in R environment  #####
#############################


project <- c("P_ctrl_VS_normal_organoid")
deal_design <- c("P_ctrl","normal_organoid")
significant_cutoff <- c(1)
organism <- "mouse"
##物种为人就是hsa，是老鼠就是mouse

file_path <- "/mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/P_ctrl_vs_nomal_organoid/"
sample_sampletable.path <- "/mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/star_mm10_out"
load(file="/mnt/data/user_data/xiangyu/workshop/WORKFLOW_RNAseq/ebg_mm10.RData")
load(file="/mnt/data/user_data/xiangyu/workshop/WORKFLOW_RNAseq/txdb_mm10.RData")  


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
csvfile <- file.path(indir, "sample_sampletable_ctrl_normal_organoid.csv")
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

setwd(file_path)
getwd()

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
library("ggplot2")
pcaData <- plotPCA(rld, intgroup = c("sample","order"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(x = PC1, y = PC2, color = sample)) +
geom_point(size =1) +
geom_text(aes(label=order))+
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

#DESeq2 normalize
colData(se) <- DataFrame(sampleTable)
library("DESeq2")
dds <- DESeqDataSet(se, design = ~ deal)
countdata <- assay(se)
coldata <- colData(se)
DESeq_counts <- DESeq(dds)
normalized_DEseq <- counts(DESeq_counts, normalized=TRUE)
DEseq_res <- results(DESeq_counts,contrast=c("deal","P_ctrl","normal_organoid"))    #最后一个是CTRL
res_1 <- cbind(normalized_DEseq,DEseq_res)
write.csv(tpm,file=tpm_csv)


dds <- DESeq(dds)
res_1 <- results(dds, contrast=res_1_filter)
colnames(res_1) = paste(colnames_1,colnames(res_1),sep="_")
write.csv(res_1, file = res_1_file.csv)

res_1$ensembl <- mapIds(x = anno_data,
                        keys = rownames(res_1),
						keytype ="SYMBOL",
						column ="ENSEMBL",
						multiVals="first")
res_1$entrez <- mapIds(x = anno_data,
                        keys = rownames(res_1),
						keytype ="SYMBOL",
						column ="ENTREZID",
						multiVals="first")
AA <- rownames(res_1)
AA <- as.character(AA)
res_1$GENENAME <- mapIds(x = anno_data,
                        keys = AA,
						keytype ="SYMBOL",
						column ="GENENAME",
						multiVals="first")
write.csv(res_1, file = res_1_file_sym_name)
all_summry <- cbind(count_and_tpm,res_1)
write.csv(all_summry, file = ALL.CSV_FILE)


res_1 <- na.omit(res_1)
y <- res_1[c(1:nrow(res_1)),6]
res_1 <- res_1[with(res_1,y<0.05),]
y <- res_1[c(1:nrow(res_1)),2]

upres_1 <- res_1[with(res_1,y>=significant_cutoff),]
ee	<-as.matrix(upres_1$entrez)
	dd <- as.vector(ee)
KEGGupres_1 <- enrichKEGG(gene =dd, 
					organism = organism, 
					keyType = "ncbi-geneid",
					 pvalueCutoff = 0.05,
				       pAdjustMethod = "BH", 
				       minGSSize = 10, 
				       maxGSSize = 500,
				       qvalueCutoff = 0.2, 
				       use_internal_data = FALSE)

downres_1  <- res_1[with(res_1,y<=-significant_cutoff),]
hh	<-as.matrix(downres_1$entrez)
	ii <- as.vector(hh)
KEGGdownres_1 <- enrichKEGG(gene =ii, 
					organism = organism, 
					keyType = "ncbi-geneid",
					 pvalueCutoff = 0.05,
				       pAdjustMethod = "BH", 
				       minGSSize = 10, 
				       maxGSSize = 500,
				       qvalueCutoff = 0.2, 
				       use_internal_data = FALSE)

KEGG_pro_enhan_up <- setReadable(KEGGupres_1, "org.Mm.eg.db", keyType="ENTREZID")
KEGG_pro_enhan_down <- setReadable(KEGGdownres_1, "org.Mm.eg.db", keyType="ENTREZID")

write.csv(KEGG_pro_enhan_up, file = KEGGupres_1_file.csv)
write.csv(KEGG_pro_enhan_down, file = KEGGdownres_1_file.csv)

pdf(file = KEGGupres_1_pdf)
dotplot(KEGGupres_1, showCategory=20)
dev.off()
pdf(file = KEGGdownres_1_pdf)
dotplot(KEGGdownres_1, showCategory=20)
dev.off()


pdf(file = KEGGres_1_all_pdf)
KEGGupres_1uporder <- KEGGupres_1[order((KEGGupres_1$Count),decreasing = TRUE),]
KEGGdownres_1downorder <- KEGGdownres_1[order((KEGGdownres_1$Count),decreasing = TRUE),]
down_kegg<-KEGGdownres_1downorder[KEGGdownres_1downorder$p.adjust<0.05,];
down_kegg$group=-1
up_kegg<-KEGGupres_1uporder[KEGGupres_1uporder$p.adjust<0.05,];
up_kegg$group=1
down_kegg1 <- down_kegg[-(21:nrow(down_kegg)),]
up_kegg1 <- up_kegg[-(21:nrow(up_kegg)),]
  dat=rbind(up_kegg1,down_kegg1)
  colnames(dat)
    dat$p.adjust = -log10(dat$p.adjust)
  dat$p.adjust=dat$p.adjust*dat$group 
  dat=dat[order(dat$p.adjust,decreasing = F),]
  library("ggplot2")
  g_kegg<- ggplot(dat, aes(x=reorder(Description,order(p.adjust, decreasing = F)), y=p.adjust, fill=group)) + 
    geom_bar(stat="identity") + 
    scale_fill_gradient(low="blue",high="red",guide = FALSE) + 
    scale_x_discrete(name ="Pathway names") +
    scale_y_continuous(name ="-log10p_adjust") +
    coord_flip() +
    ggtitle(KEGGres_1_all_title)
  print(g_kegg)
dev.off()

GOupres_1_all <- enrichGO(gene = dd, 
	           OrgDb = anno_data,
				ont = "all", 
		             pvalueCutoff = 0.01, 
                     pAdjustMethod = "BH", 
                     qvalueCutoff = 0.05,
                     minGSSize = 10, 
                     maxGSSize = 500, 
                     readable = T, 
                     pool = FALSE)
write.csv(GOupres_1_all,file=GOres_1_all_UP_csv)

ff <- barplot(GOupres_1_all, split="ONTOLOGY",showCategory=15) + facet_grid(ONTOLOGY~., scale="free")
ff1 <- ff + theme(axis.text.y = element_text(size = 8)) +labs(title = NULL)
pdf(file = GOres_1_all_UP_pdf)
ff1
dev.off()

hh	<-as.matrix(downres_1$entrez)
	ii <- as.vector(hh)
GOdownres_1_all <- enrichGO(gene = ii, 
	           OrgDb = anno_data,
				ont = "all", 
		             pvalueCutoff = 0.01, 
                     pAdjustMethod = "BH", 
                     qvalueCutoff = 0.05,
                     minGSSize = 10, 
                     maxGSSize = 500, 
                     readable = T, 
                     pool = FALSE)
write.csv(GOdownres_1_all,file=GOres_1_all_DOWN_csv)

ff <- barplot(GOdownres_1_all, split="ONTOLOGY",showCategory=15) + facet_grid(ONTOLOGY~., scale="free")
ff1 <- ff + theme(axis.text.y = element_text(size = 8)) +labs(title = NULL)
pdf(file = GOres_1_all_DOWN_pdf)
ff1
dev.off()


*********gsea输入文件准备***********
gsea_data <- read.csv("/mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/P_ctrl_vs_nomal_organoid/P_ctrl_VS_normal_organoid_05_P_ctrl_v_normal_organoid_count_tpm_symbol_and_anno.csv")
gsea_data$Description <- c("NULL")
gsea_data <- gsea_data[,c(20,23,8,9,10,11,12,13)]
nrow(gsea_data)
gsea_data <- na.omit(gsea_data)
nrow(gsea_data)
colnames(gsea_data)
colnames(gsea_data) <- c("NAME","Description","P_ctrl_1","P_ctrl_2","P_ctrl_3","normal_organoid_1","normal_organoid_2","normal_organoid_3")
write.csv(gsea_data,file="/mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/P_ctrl_vs_nomal_organoid/gsea_data.csv")
#在excel中改好格式(就以Ensembl id为行名)，然后直接保存为txt格式，再导入R修改文件名为.gct即可
#加上  
       #1.2
       #总基因行数 总样本数

vi gsea.cls
6 2 1
# P_ctrl normal_organoid
0 0 0 1 1 1



#样本多的时候，-Xmx1024m改为-Xmx10240 \       9022      注意\后面不能有空格
vi gsea_shell.sh

java -cp /mnt/data/user_data/xiangyu/programme/gsea/gsea-3.0.jar -Xmx1024m xtools.gsea.Gsea \
-res gsea_data.gct \
-cls gsea.cls \
-gmx /mnt/data/user_data/xiangyu/programme/gsea/msigdb_v6.1_GMTs/h.all.v6.1.symbols.gmt \
-chip /mnt/data/user_data/xiangyu/programme/gsea/chip/ENSEMBL_mouse_gene.chip \
-collapse true -mode Max_probe -norm meandiv -nperm 1000 -permute gene_set \
-rnd_type no_balance -scoring_scheme weighted -rpt_label my_analysis \
-metric Signal2Noise -sort real -order descending -include_only_symbols true \
-make_sets true -median false -num 100 -plot_top_x 50 -rnd_seed timestamp \
-save_rnd_lists false -set_max 500 -set_min 15 -zip_report false \
-out ./P_ctrl_VS_normal_organoid_h -gui false

java -cp /mnt/data/user_data/xiangyu/programme/gsea/gsea-3.0.jar -Xmx1024m xtools.gsea.Gsea \
-res gsea_data.gct \
-cls gsea.cls \
-gmx /mnt/data/user_data/xiangyu/programme/gsea/msigdb_v6.1_GMTs/c1.all.v6.1.symbols.gmt \
-chip /mnt/data/user_data/xiangyu/programme/gsea/chip/ENSEMBL_mouse_gene.chip \
-collapse true -mode Max_probe -norm meandiv -nperm 1000 -permute gene_set \
-rnd_type no_balance -scoring_scheme weighted -rpt_label my_analysis \
-metric Signal2Noise -sort real -order descending -include_only_symbols true \
-make_sets true -median false -num 100 -plot_top_x 50 -rnd_seed timestamp \
-save_rnd_lists false -set_max 500 -set_min 15 -zip_report false \
-out ./P_ctrl_VS_normal_organoid_c1 -gui false

java -cp /mnt/data/user_data/xiangyu/programme/gsea/gsea-3.0.jar -Xmx1024m xtools.gsea.Gsea \
-res gsea_data.gct \
-cls gsea.cls \
-gmx /mnt/data/user_data/xiangyu/programme/gsea/msigdb_v6.1_GMTs/c2.all.v6.1.symbols.gmt \
-chip /mnt/data/user_data/xiangyu/programme/gsea/chip/ENSEMBL_mouse_gene.chip \
-collapse true -mode Max_probe -norm meandiv -nperm 1000 -permute gene_set \
-rnd_type no_balance -scoring_scheme weighted -rpt_label my_analysis \
-metric Signal2Noise -sort real -order descending -include_only_symbols true \
-make_sets true -median false -num 100 -plot_top_x 50 -rnd_seed timestamp \
-save_rnd_lists false -set_max 500 -set_min 15 -zip_report false \
-out ./P_ctrl_VS_normal_organoid_c2 -gui false

java -cp /mnt/data/user_data/xiangyu/programme/gsea/gsea-3.0.jar -Xmx1024m xtools.gsea.Gsea \
-res gsea_data.gct \
-cls gsea.cls \
-gmx /mnt/data/user_data/xiangyu/programme/gsea/msigdb_v6.1_GMTs/c3.all.v6.1.symbols.gmt \
-chip /mnt/data/user_data/xiangyu/programme/gsea/chip/ENSEMBL_mouse_gene.chip \
-collapse true -mode Max_probe -norm meandiv -nperm 1000 -permute gene_set \
-rnd_type no_balance -scoring_scheme weighted -rpt_label my_analysis \
-metric Signal2Noise -sort real -order descending -include_only_symbols true \
-make_sets true -median false -num 100 -plot_top_x 50 -rnd_seed timestamp \
-save_rnd_lists false -set_max 500 -set_min 15 -zip_report false \
-out ./P_ctrl_VS_normal_organoid_c3 -gui false

java -cp /mnt/data/user_data/xiangyu/programme/gsea/gsea-3.0.jar -Xmx1024m xtools.gsea.Gsea \
-res gsea_data.gct \
-cls gsea.cls \
-gmx /mnt/data/user_data/xiangyu/programme/gsea/msigdb_v6.1_GMTs/c4.all.v6.1.symbols.gmt \
-chip /mnt/data/user_data/xiangyu/programme/gsea/chip/ENSEMBL_mouse_gene.chip \
-collapse true -mode Max_probe -norm meandiv -nperm 1000 -permute gene_set \
-rnd_type no_balance -scoring_scheme weighted -rpt_label my_analysis \
-metric Signal2Noise -sort real -order descending -include_only_symbols true \
-make_sets true -median false -num 100 -plot_top_x 50 -rnd_seed timestamp \
-save_rnd_lists false -set_max 500 -set_min 15 -zip_report false \
-out ./P_ctrl_VS_normal_organoid_c4 -gui false

java -cp /mnt/data/user_data/xiangyu/programme/gsea/gsea-3.0.jar -Xmx1024m xtools.gsea.Gsea \
-res gsea_data.gct \
-cls gsea.cls \
-gmx /mnt/data/user_data/xiangyu/programme/gsea/msigdb_v6.1_GMTs/c5.all.v6.1.symbols.gmt \
-chip /mnt/data/user_data/xiangyu/programme/gsea/chip/ENSEMBL_mouse_gene.chip \
-collapse true -mode Max_probe -norm meandiv -nperm 1000 -permute gene_set \
-rnd_type no_balance -scoring_scheme weighted -rpt_label my_analysis \
-metric Signal2Noise -sort real -order descending -include_only_symbols true \
-make_sets true -median false -num 100 -plot_top_x 50 -rnd_seed timestamp \
-save_rnd_lists false -set_max 500 -set_min 15 -zip_report false \
-out ./P_ctrl_VS_normal_organoid_c5 -gui false

java -cp /mnt/data/user_data/xiangyu/programme/gsea/gsea-3.0.jar -Xmx1024m xtools.gsea.Gsea \
-res gsea_data.gct \
-cls gsea.cls \
-gmx /mnt/data/user_data/xiangyu/programme/gsea/msigdb_v6.1_GMTs/c6.all.v6.1.symbols.gmt \
-chip /mnt/data/user_data/xiangyu/programme/gsea/chip/ENSEMBL_mouse_gene.chip \
-collapse true -mode Max_probe -norm meandiv -nperm 1000 -permute gene_set \
-rnd_type no_balance -scoring_scheme weighted -rpt_label my_analysis \
-metric Signal2Noise -sort real -order descending -include_only_symbols true \
-make_sets true -median false -num 100 -plot_top_x 50 -rnd_seed timestamp \
-save_rnd_lists false -set_max 500 -set_min 15 -zip_report false \
-out ./P_ctrl_VS_normal_organoid_c6 -gui false

java -cp /mnt/data/user_data/xiangyu/programme/gsea/gsea-3.0.jar -Xmx1024m xtools.gsea.Gsea \
-res gsea_data.gct \
-cls gsea.cls \
-gmx /mnt/data/user_data/xiangyu/programme/gsea/msigdb_v6.1_GMTs/c7.all.v6.1.symbols.gmt \
-chip /mnt/data/user_data/xiangyu/programme/gsea/chip/ENSEMBL_mouse_gene.chip \
-collapse true -mode Max_probe -norm meandiv -nperm 1000 -permute gene_set \
-rnd_type no_balance -scoring_scheme weighted -rpt_label my_analysis \
-metric Signal2Noise -sort real -order descending -include_only_symbols true \
-make_sets true -median false -num 100 -plot_top_x 50 -rnd_seed timestamp \
-save_rnd_lists false -set_max 500 -set_min 15 -zip_report false \
-out ./P_ctrl_VS_normal_organoid_c7 -gui false


bash gsea_shell.sh


****************************************normal_nasopha vs normal_lung***************************************
vi sampletable_nasopha_vs_lung.csv

ids,sample,deal,order
1,lung111_mouse1_NFF.bam,lung_1,lung,1
2,lung111_mouse2_NFF.bam,lung_2,lung,2
3,lung111_mouse3_NFF.bam,lung_3,lung,3
4,nasopha_mouse1_LLH.bam,nasopha_1,nasopha,4
5,nasopha_mouse2_LLH.bam,nasopha_2,nasopha,5
6,nasopha_mouse3_LLH.bam,nasopha_3,nasopha,6

project <- c("nasopha_vs_lung")
deal_design <- c("nasopha","lung")
organism <- "mouse"
##物种为人就是hsa，是老鼠就是mouse

file_path <- "/mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/nasopha_vs_others/"
sample_sampletable.path <- "/mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/GRCm38_bam"
load(file="/mnt/data/user_data/xiangyu/workshop/WORKFLOW_RNAseq/ebg_GRCm38.RData")
load(file="/mnt/data/user_data/xiangyu/workshop/WORKFLOW_RNAseq/txdb_GRCm38.RData")  

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


indir <- file.path(sample_sampletable.path)
csvfile <- file.path(indir, "sampletable_nasopha_vs_lung.csv")
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
getwd()

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

