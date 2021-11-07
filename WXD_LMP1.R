#数据位置
/mnt/data/sequencedata/RNAseq/RNAseq_77_WXD_NPC_20210108_11samples/CP2018091214925/H101SC20030543/RSCS7000/X101SC20030543-Z01/X101SC20030543-Z01-J090/2.cleandata

#检验数据完整性
cd /mnt/data/sequencedata/RNAseq/RNAseq_76_ZJP_bladder_20210105_8samples/CP2018091207789/H101SC19062343/RSCR0102/X101SC19062343-Z01/X101SC19062343-Z01-J075
md5sum -c md5.txt
md5sum -c MD5_Ctr-1_FRAS210005133-1r.txt

#star比对
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

#test
cat config |while read id;
do echo $id
arr=($id)
fq1=${arr[1]}
fq2=${arr[2]}
sample=${arr[0]}
echo $fq1
echo $fq2
echo $sample
done

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

cd /mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/star_mm10_out/primary
cd /mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/star_mm10_out/metastasis
#构建index
ls *.sortedByCoord.out.bam |while read id
do
samtools index ./$id
done

#igv
bash /mnt/data/program/IGV_2.4.10/igv.sh



**************************************************外源LMP1可视化************************************************************************
cd /mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/bwa_LMP1_out

vim WXD_LMP1.fasta
>LMP1
atggaacgcgaccttgagaggggcccaccgggcccgccacggccccctctaggaccccccctctcctcttccataggccttgctctccttctcctgctcttggcgctactgttctggctgtatatcgttatgagtaactggactggaggagcgctccttgtcctctattcctttgctctcatgcttattattatcattctcatcatctttatcttcagaagagaccttctctgtccacttggaggccttggtctactcctactgatgatcaccctcctactcatcgctctctggaatttgcacggacaggcattgtaccttggaattgtgctgttcatctttggctgcttacttgtcttaggtctctggatctacttcttggagattctctggcggcttggtgccaccatctggcagcttttggccttcatcctagccttcttcctagccatcatcctgcttattattgctctctatctacaacaaaactggtggactctattggttgatctcctttggctcctcctgtttatggccattttaatctggatgtattatcatggaccacgacacactgatgaacaccaccacgatgactccctcccgcaccctcaacaagctaccgacgattctagccatgaatctgactctaactccaacgagggcagacaccacctgctcgtgagtggggccggcgacggacccccactctgctctcaaaacctaggcgcacctggaggtggtcctgacaatggcccacaggaccctgacaacactgatgacaatggcccacaggaccctgacaacactgatgacaatggcccacatgacccgctgcctcaggaccctgacaacactgatgacaatggcccacaggaccctgacaacactgatgacaatggcccacatgacccgctgcctcataaccctagcgactctgctggaaatgatggaggccctccaaatttgacggaagaggttgaaaacaaaggaggtgaccgggacccgccttcgatgacagacggtggcggcggtgatccacaccttcctacgctgcttttgggtacttctggttccggtggagatgatgacgacccccacggcccagttcagctaagctactatgactaa

$$$config
$$$config
$$$config

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


/mnt/data/user_data/xiangyu/programme/bwa-0.7.17/bwa index WXD_LMP1.fasta

###bwa mem

cat config  |while read id;
do echo $id
arr=($id)
fq2=${arr[2]}
fq1=${arr[1]}
sample=${arr[0]}
echo $fq1
echo $fq2
echo $sample
/mnt/data/user_data/xiangyu/programme/bwa-0.7.17/bwa mem -t 20 -o ${sample}'_virus'.sam \
WXD_LMP1.fasta $fq1 $fq2 ;done

GATK=/mnt/data/user_data/xiangyu/programme/gatk-4.1.3.0/gatk
cat config  |while read id;
do echo $id
arr=($id)
sample=${arr[0]}
echo $sample
$GATK  --java-options "-Xmx20G -Djava.io.tmpdir=./"  SortSam -SO coordinate  -I ${sample}'_virus'.sam -O ${sample}'_virus'.bam > ${sample}'_virus'.log.sort
samtools index ${sample}'_virus'.bam ;done

*******************************************************外源BFP可视化**************************************************************************
cd /mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/bwa_BFP_out

vim WXD_BFP.fasta
>BFP
atgagcgagctgattaaggagaacatgcacatgaagctgtacatggagggcaccgtggacaaccatcacttcaagtgcacatccgagggcgaaggcaagccctacgagggcacccagaccatgagaatcaaggtggtcgagggcggccctctccccttcgccttcgacatcctggctactagcttcctctacggcagcaagaccttcatcaaccacacccagggcatccccgacttcttcaagcagtccttccctgagggcttcacatgggagagagtcaccacatacgaagacgggggcgtgctgaccgctacccaggacaccagcctccaggacggctgcctcatctacaacgtcaagatcagaggggtgaacttcacatccaacggccctgtgatgcagaagaaaacactcggctgggaggccttcaccgagacgctgtaccccgctgacggcggcctggaaggcagaaacgacatggccctgaagctcgtgggcgggagccatctgatcgcaaacatcaagaccacatatagatccaagaaacccgctaagaacctcaagatgcctggcgtctactatgtggactacagactggaaagaatcaaggaggccaacaacgagacctacgtcgagcagcacgaggtggcagtggccagatactgcgacctccctagcaaactggggcacaagcttaattga

/mnt/data/user_data/xiangyu/programme/bwa-0.7.17/bwa index WXD_BFP.fasta
samtools faidx WXD_BFP.fasta

###bwa mem

cat /mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/bwa_LMP1_out/config  |while read id;
do echo $id
arr=($id)
fq2=${arr[2]}
fq1=${arr[1]}
sample=${arr[0]}
echo $fq1
echo $fq2
echo $sample
/mnt/data/user_data/xiangyu/programme/bwa-0.7.17/bwa mem -t 20 -o ${sample}'_virus'.sam \
WXD_BFP.fasta $fq1 $fq2 ;done

GATK=/mnt/data/user_data/xiangyu/programme/gatk-4.1.3.0/gatk
cat /mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/bwa_LMP1_out/config |while read id;
do echo $id
arr=($id)
sample=${arr[0]}
echo $sample
$GATK  --java-options "-Xmx20G -Djava.io.tmpdir=./"  SortSam -SO coordinate  -I ${sample}'_virus'.sam -O ${sample}'_virus'.bam > ${sample}'_virus'.log.sort
samtools index ${sample}'_virus'.bam ;done

******************************************************外源LMP1-EFS-BFP可视化**************************************************************************
cd /mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/bwa_LMP1-EFS-BFP_out

vim WXD_LMP1-EFS-BFP.fasta
>LMP1_EFS_BFP
atggaacgcgaccttgagaggggcccaccgggcccgccacggccccctctaggaccccccctctcctcttccataggccttgctctccttctcctgctcttggcgctactgttctggctgtatatcgttatgagtaactggactggaggagcgctccttgtcctctattcctttgctctcatgcttattattatcattctcatcatctttatcttcagaagagaccttctctgtccacttggaggccttggtctactcctactgatgatcaccctcctactcatcgctctctggaatttgcacggacaggcattgtaccttggaattgtgctgttcatctttggctgcttacttgtcttaggtctctggatctacttcttggagattctctggcggcttggtgccaccatctggcagcttttggccttcatcctagccttcttcctagccatcatcctgcttattattgctctctatctacaacaaaactggtggactctattggttgatctcctttggctcctcctgtttatggccattttaatctggatgtattatcatggaccacgacacactgatgaacaccaccacgatgactccctcccgcaccctcaacaagctaccgacgattctagccatgaatctgactctaactccaacgagggcagacaccacctgctcgtgagtggggccggcgacggacccccactctgctctcaaaacctaggcgcacctggaggtggtcctgacaatggcccacaggaccctgacaacactgatgacaatggcccacaggaccctgacaacactgatgacaatggcccacatgacccgctgcctcaggaccctgacaacactgatgacaatggcccacaggaccctgacaacactgatgacaatggcccacatgacccgctgcctcataaccctagcgactctgctggaaatgatggaggccctccaaatttgacggaagaggttgaaaacaaaggaggtgaccgggacccgccttcgatgacagacggtggcggcggtgatccacaccttcctacgctgcttttgggtacttctggttccggtggagatgatgacgacccccacggcccagttcagctaagctactatgactaataggtcttgaaaggagtgggaattggctccggtgcccgtcagtgggcagagcgcacatcgcccacagtccccgagaagttggggggaggggtcggcaattgatccggtgcctagagaaggtggcgcggggtaaactgggaaagtgatgtcgtgtactggctccgcctttttcccgagggtgggggagaaccgtatataagtgcagtagtcgccgtgaacgttctttttcgcaacgggtttgccgccagaacacaggatgagcgagctgattaaggagaacatgcacatgaagctgtacatggagggcaccgtggacaaccatcacttcaagtgcacatccgagggcgaaggcaagccctacgagggcacccagaccatgagaatcaaggtggtcgagggcggccctctccccttcgccttcgacatcctggctactagcttcctctacggcagcaagaccttcatcaaccacacccagggcatccccgacttcttcaagcagtccttccctgagggcttcacatgggagagagtcaccacatacgaagacgggggcgtgctgaccgctacccaggacaccagcctccaggacggctgcctcatctacaacgtcaagatcagaggggtgaacttcacatccaacggccctgtgatgcagaagaaaacactcggctgggaggccttcaccgagacgctgtaccccgctgacggcggcctggaaggcagaaacgacatggccctgaagctcgtgggcgggagccatctgatcgcaaacatcaagaccacatatagatccaagaaacccgctaagaacctcaagatgcctggcgtctactatgtggactacagactggaaagaatcaaggaggccaacaacgagacctacgtcgagcagcacgaggtggcagtggccagatactgcgacctccctagcaaactggggcacaagcttaattga

/mnt/data/user_data/xiangyu/programme/bwa-0.7.17/bwa index WXD_LMP1-EFS-BFP.fasta
samtools faidx WXD_LMP1-EFS-BFP.fasta

###bwa mem

cat /mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/bwa_LMP1_out/config  |while read id;
do echo $id
arr=($id)
fq2=${arr[2]}
fq1=${arr[1]}
sample=${arr[0]}
echo $fq1
echo $fq2
echo $sample
/mnt/data/user_data/xiangyu/programme/bwa-0.7.17/bwa mem -t 20 -o ${sample}'_virus'.sam \
WXD_LMP1-EFS-BFP.fasta $fq1 $fq2 ;done

GATK=/mnt/data/user_data/xiangyu/programme/gatk-4.1.3.0/gatk
cat /mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/bwa_LMP1_out/config  |while read id;
do echo $id
arr=($id)
sample=${arr[0]}
echo $sample
$GATK  --java-options "-Xmx20G -Djava.io.tmpdir=./"  SortSam -SO coordinate  -I ${sample}'_virus'.sam -O ${sample}'_virus'.bam > ${sample}'_virus'.log.sort
samtools index ${sample}'_virus'.bam ;done



*********************************************************************************************************************************


*********************************************************************************************************************************
****************************************************primary_LMP1_vs_ctrl*****************************************************************************
*********************************************************************************************************************************

############################
#move to bamfiles path  #####
#move to bamfiles path  #####
#############################
cd /mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/star_mm10_out

vi sample_sampletable_primary.csv

ids,sample,deal,order
1,ctrl_1.Aligned.sortedByCoord.out.bam,ctrl_1,ctrl,1
2,ctrl_2.Aligned.sortedByCoord.out.bam,ctrl_2,ctrl,2
3,ctrl_3.Aligned.sortedByCoord.out.bam,ctrl_3,ctrl,3
4,LMP1_2.Aligned.sortedByCoord.out.bam,LMP1_2,LMP1,4
5,LMP1_3.Aligned.sortedByCoord.out.bam,LMP1_3,LMP1,5

############################
#work in R environment  #####
#work in R environment  #####
#############################


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

#发现收得较早的样LMP1_1和ctrl聚到一起，故决定踢掉LMP1_1这个样。

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

#res_1 <- read.csv("/mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/ctrl_metastasis_vs_primary/ctrl_metastasis_vs_primary_04_metastasis_v_primary_symbol_and_anno.csv",row.names=1)
#padj<0.05 gene，3791
#significant_cutoff <- c(0.5)
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
#up的gene为0，将logFC阈值设置为0.5的时候则能够做出图

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
************************************************KEGG down GENEid
res_1 <- read.csv("/mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/primary_LMP1_vs_ctrl/primary_LMP1_04_LMP1_v_ctrl_symbol_and_anno.csv")
significant_cutoff <- c(0.5)

rownames(res_1) <- res_1[,1]
res_1 <- res_1[,-1]
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


KEGG_pro_enhan_up <- setReadable(KEGGupres_1, "org.Mm.eg.db", keyType="ENTREZID")
#up的gene为0，将logFC阈值设置为0.5的时候则能够做出图

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

write.csv(KEGG_pro_enhan_up, file = "primary_LMP1_07_LMP1_v_ctrl_KEGG_UP_genes.csv")
write.csv(KEGG_pro_enhan_down, file = "primary_LMP1_07_LMP1_v_ctrl_KEGG_DOWN_genes.csv")
****************************************************************

pdf(file = KEGGupres_1_pdf)
dotplot(KEGGupres_1, showCategory=20)
dev.off()

pdf(file = "primary_LMP1_v_ctrl_KEGG_UP_top15.pdf")
dotplot(KEGGupres_1, showCategory=15)
dev.off()

KEGGupres_1_top20 <-  KEGGupres_1[1:15,]
x <- KEGGupres_1_top20[order(KEGGupres_1_top20$GeneRatio,decreasing = TRUE),]
x <- x[c(2:15,1),]
y=new("enrichResult",
      result=x)
y
dotplot(y,orderBy ='GeneRatio',showCategory=10)

f <- dotplot(y, orderBy ='GeneRatio',showCategory=11)
ggsave("/mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/primary_LMP1_vs_ctrl/primary_LMP1_v_ctrl_KEGG_UP_top11.svg", plot=f,width = 6,height = 5,dpi=1080)


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

ff <- barplot(GOdownres_1_all, split="ONTOLOGY",showCategory=5) + facet_grid(ONTOLOGY~., scale="free")
ff1 <- ff + theme(axis.text.y = element_text(size = 8)) +labs(title = NULL)
pdf(file = "primary_LMP1_v_ctrl_GO_DOWN_top5.pdf")
ff1
dev.off()



*********gsea输入文件准备***********
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

$$$$ motif 分析
$$$$ motif 分析

ESCA_up <- ESCA_squm %>% filter(p.adj <= 0.01) %>% filter(logFC >= 1)
ESCA_up <- na.omit(ESCA_up)

library(Homo.sapiens)
all_genes <- genes(Homo.sapiens)
all_gene_TSS <- resize(all_genes,1)
TSS_2k <-promoters(all_gene_TSS, 2000, 2000)
TSS_2k <- as.data.frame(TSS_2k)
TSS_2k$entrez <- row.names(TSS_2k)
all_gene_TSS <- as.data.frame(all_gene_TSS)

ESCA_up$entrez <- mapIds(x = anno_data,
            keys = rownames(ESCA_up),
            keytype ="ENSEMBL",
            column ="ENTREZID",
            multiVals="first")

ESCA_up_TSS <- merge(ESCA_up,TSS_2k,by = "entrez", all = FALSE)

setwd("/mnt/data/user_data/zhaolei/project/zms_EZH2/human_data/Ezh2_ESCA")
write.csv(ESCA_up_TSS,"ESCA_up_TSS.csv")
#看一下/mnt/data/user_data/zhaolei/project/zms_EZH2/human_data/路径下的bed格式，把文件改成对应的格式
findMotifsGenome.pl ESCA_up_TSS.bed hg19 ESCA_suqa_vs_other_up_motif -len 8,10,12


*********************************************************************************************************************************
****************************************************metastasis_LMP1_vs_ctrl*****************************************************************************
*********************************************************************************************************************************
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

#发现收得较早的样LMP1_1和ctrl聚到一起，故决定踢掉LMP1_1这个样。

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

#res_1 <- read.csv("/mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/primary_LMP1_vs_ctrl/primary_LMP1_04_LMP1_v_ctrl_symbol_and_anno.csv")
#padj<0.05 gene，436     pvalue<0.05 gene，2146
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



*******************************************************************human_data***********************************************************************************
#数据来源 https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi    GSE102349
#https://www.ncbi.nlm.nih.gov/Traces/study/?page=3&acc=PRJNA397538&o=acc_s%3Aa
#Genomic Analysis of Nasopharyngeal Carcinoma Reveals TME-Based Subtypes
#数据信息 https://www.ebi.ac.uk/ena/browser/view/PRJNA397538

vi ncbi_NPC.txt           #pair

/vol1/fastq/SRR590/007/SRR5908747/SRR5908747_1.fastq.gz
/vol1/fastq/SRR590/008/SRR5908748/SRR5908748_1.fastq.gz
/vol1/fastq/SRR590/009/SRR5908749/SRR5908749_1.fastq.gz
/vol1/fastq/SRR590/000/SRR5908750/SRR5908750_1.fastq.gz
/vol1/fastq/SRR590/001/SRR5908751/SRR5908751_1.fastq.gz
/vol1/fastq/SRR590/002/SRR5908752/SRR5908752_1.fastq.gz
/vol1/fastq/SRR590/003/SRR5908753/SRR5908753_1.fastq.gz
/vol1/fastq/SRR590/004/SRR5908754/SRR5908754_1.fastq.gz
/vol1/fastq/SRR590/005/SRR5908755/SRR5908755_1.fastq.gz
/vol1/fastq/SRR590/006/SRR5908756/SRR5908756_1.fastq.gz
/vol1/fastq/SRR590/007/SRR5908757/SRR5908757_1.fastq.gz
/vol1/fastq/SRR590/008/SRR5908758/SRR5908758_1.fastq.gz
/vol1/fastq/SRR590/009/SRR5908759/SRR5908759_1.fastq.gz
/vol1/fastq/SRR590/000/SRR5908760/SRR5908760_1.fastq.gz
/vol1/fastq/SRR590/001/SRR5908761/SRR5908761_1.fastq.gz
/vol1/fastq/SRR590/002/SRR5908762/SRR5908762_1.fastq.gz
/vol1/fastq/SRR590/003/SRR5908763/SRR5908763_1.fastq.gz
/vol1/fastq/SRR590/004/SRR5908764/SRR5908764_1.fastq.gz
/vol1/fastq/SRR590/005/SRR5908765/SRR5908765_1.fastq.gz
/vol1/fastq/SRR590/006/SRR5908766/SRR5908766_1.fastq.gz
/vol1/fastq/SRR590/007/SRR5908767/SRR5908767_1.fastq.gz
/vol1/fastq/SRR590/008/SRR5908768/SRR5908768_1.fastq.gz
/vol1/fastq/SRR590/009/SRR5908769/SRR5908769_1.fastq.gz
/vol1/fastq/SRR590/000/SRR5908770/SRR5908770_1.fastq.gz
/vol1/fastq/SRR590/001/SRR5908771/SRR5908771_1.fastq.gz
/vol1/fastq/SRR590/002/SRR5908772/SRR5908772_1.fastq.gz
/vol1/fastq/SRR590/003/SRR5908773/SRR5908773_1.fastq.gz
/vol1/fastq/SRR590/004/SRR5908774/SRR5908774_1.fastq.gz
/vol1/fastq/SRR590/005/SRR5908775/SRR5908775_1.fastq.gz
/vol1/fastq/SRR590/006/SRR5908776/SRR5908776_1.fastq.gz
/vol1/fastq/SRR590/007/SRR5908777/SRR5908777_1.fastq.gz
/vol1/fastq/SRR590/008/SRR5908778/SRR5908778_1.fastq.gz
/vol1/fastq/SRR590/009/SRR5908779/SRR5908779_1.fastq.gz
/vol1/fastq/SRR590/000/SRR5908780/SRR5908780_1.fastq.gz
/vol1/fastq/SRR590/001/SRR5908781/SRR5908781_1.fastq.gz
/vol1/fastq/SRR590/002/SRR5908782/SRR5908782_1.fastq.gz
/vol1/fastq/SRR590/003/SRR5908783/SRR5908783_1.fastq.gz
/vol1/fastq/SRR590/004/SRR5908784/SRR5908784_1.fastq.gz
/vol1/fastq/SRR590/005/SRR5908785/SRR5908785_1.fastq.gz
/vol1/fastq/SRR590/006/SRR5908786/SRR5908786_1.fastq.gz
/vol1/fastq/SRR590/007/SRR5908787/SRR5908787_1.fastq.gz
/vol1/fastq/SRR590/008/SRR5908788/SRR5908788_1.fastq.gz
/vol1/fastq/SRR590/009/SRR5908789/SRR5908789_1.fastq.gz
/vol1/fastq/SRR590/000/SRR5908790/SRR5908790_1.fastq.gz
/vol1/fastq/SRR590/001/SRR5908791/SRR5908791_1.fastq.gz
/vol1/fastq/SRR590/002/SRR5908792/SRR5908792_1.fastq.gz
/vol1/fastq/SRR590/003/SRR5908793/SRR5908793_1.fastq.gz
/vol1/fastq/SRR590/004/SRR5908794/SRR5908794_1.fastq.gz
/vol1/fastq/SRR590/005/SRR5908795/SRR5908795_1.fastq.gz
/vol1/fastq/SRR590/006/SRR5908796/SRR5908796_1.fastq.gz
/vol1/fastq/SRR590/007/SRR5908797/SRR5908797_1.fastq.gz
/vol1/fastq/SRR590/008/SRR5908798/SRR5908798_1.fastq.gz
/vol1/fastq/SRR590/009/SRR5908799/SRR5908799_1.fastq.gz
/vol1/fastq/SRR590/000/SRR5908800/SRR5908800_1.fastq.gz
/vol1/fastq/SRR590/001/SRR5908801/SRR5908801_1.fastq.gz
/vol1/fastq/SRR590/002/SRR5908802/SRR5908802_1.fastq.gz
/vol1/fastq/SRR590/003/SRR5908803/SRR5908803_1.fastq.gz
/vol1/fastq/SRR590/004/SRR5908804/SRR5908804_1.fastq.gz
/vol1/fastq/SRR590/005/SRR5908805/SRR5908805_1.fastq.gz
/vol1/fastq/SRR590/006/SRR5908806/SRR5908806_1.fastq.gz
/vol1/fastq/SRR590/007/SRR5908807/SRR5908807_1.fastq.gz
/vol1/fastq/SRR590/008/SRR5908808/SRR5908808_1.fastq.gz
/vol1/fastq/SRR590/009/SRR5908809/SRR5908809_1.fastq.gz
/vol1/fastq/SRR590/000/SRR5908810/SRR5908810_1.fastq.gz
/vol1/fastq/SRR590/001/SRR5908811/SRR5908811_1.fastq.gz
/vol1/fastq/SRR590/002/SRR5908812/SRR5908812_1.fastq.gz
/vol1/fastq/SRR590/003/SRR5908813/SRR5908813_1.fastq.gz
/vol1/fastq/SRR590/004/SRR5908814/SRR5908814_1.fastq.gz
/vol1/fastq/SRR590/005/SRR5908815/SRR5908815_1.fastq.gz
/vol1/fastq/SRR590/006/SRR5908816/SRR5908816_1.fastq.gz
/vol1/fastq/SRR590/007/SRR5908817/SRR5908817_1.fastq.gz
/vol1/fastq/SRR590/008/SRR5908818/SRR5908818_1.fastq.gz
/vol1/fastq/SRR590/009/SRR5908819/SRR5908819_1.fastq.gz
/vol1/fastq/SRR590/000/SRR5908820/SRR5908820_1.fastq.gz
/vol1/fastq/SRR590/001/SRR5908821/SRR5908821_1.fastq.gz
/vol1/fastq/SRR590/002/SRR5908822/SRR5908822_1.fastq.gz
/vol1/fastq/SRR590/003/SRR5908823/SRR5908823_1.fastq.gz
/vol1/fastq/SRR590/004/SRR5908824/SRR5908824_1.fastq.gz
/vol1/fastq/SRR590/005/SRR5908825/SRR5908825_1.fastq.gz
/vol1/fastq/SRR590/006/SRR5908826/SRR5908826_1.fastq.gz
/vol1/fastq/SRR590/007/SRR5908827/SRR5908827_1.fastq.gz
/vol1/fastq/SRR590/008/SRR5908828/SRR5908828_1.fastq.gz
/vol1/fastq/SRR590/009/SRR5908829/SRR5908829_1.fastq.gz
/vol1/fastq/SRR590/000/SRR5908830/SRR5908830_1.fastq.gz
/vol1/fastq/SRR590/001/SRR5908831/SRR5908831_1.fastq.gz
/vol1/fastq/SRR590/002/SRR5908832/SRR5908832_1.fastq.gz
/vol1/fastq/SRR590/003/SRR5908833/SRR5908833_1.fastq.gz
/vol1/fastq/SRR590/004/SRR5908834/SRR5908834_1.fastq.gz
/vol1/fastq/SRR590/005/SRR5908835/SRR5908835_1.fastq.gz
/vol1/fastq/SRR590/006/SRR5908836/SRR5908836_1.fastq.gz
/vol1/fastq/SRR590/007/SRR5908837/SRR5908837_1.fastq.gz
/vol1/fastq/SRR590/008/SRR5908838/SRR5908838_1.fastq.gz
/vol1/fastq/SRR590/009/SRR5908839/SRR5908839_1.fastq.gz
/vol1/fastq/SRR590/000/SRR5908840/SRR5908840_1.fastq.gz
/vol1/fastq/SRR590/001/SRR5908841/SRR5908841_1.fastq.gz
/vol1/fastq/SRR590/002/SRR5908842/SRR5908842_1.fastq.gz
/vol1/fastq/SRR590/003/SRR5908843/SRR5908843_1.fastq.gz
/vol1/fastq/SRR590/004/SRR5908844/SRR5908844_1.fastq.gz
/vol1/fastq/SRR590/005/SRR5908845/SRR5908845_1.fastq.gz
/vol1/fastq/SRR590/006/SRR5908846/SRR5908846_1.fastq.gz
/vol1/fastq/SRR590/007/SRR5908847/SRR5908847_1.fastq.gz
/vol1/fastq/SRR590/008/SRR5908848/SRR5908848_1.fastq.gz
/vol1/fastq/SRR590/009/SRR5908849/SRR5908849_1.fastq.gz
/vol1/fastq/SRR590/000/SRR5908850/SRR5908850_1.fastq.gz
/vol1/fastq/SRR590/001/SRR5908851/SRR5908851_1.fastq.gz
/vol1/fastq/SRR590/002/SRR5908852/SRR5908852_1.fastq.gz
/vol1/fastq/SRR590/003/SRR5908853/SRR5908853_1.fastq.gz
/vol1/fastq/SRR590/004/SRR5908854/SRR5908854_1.fastq.gz
/vol1/fastq/SRR590/005/SRR5908855/SRR5908855_1.fastq.gz
/vol1/fastq/SRR590/006/SRR5908856/SRR5908856_1.fastq.gz
/vol1/fastq/SRR590/007/SRR5908857/SRR5908857_1.fastq.gz
/vol1/fastq/SRR590/008/SRR5908858/SRR5908858_1.fastq.gz
/vol1/fastq/SRR590/009/SRR5908859/SRR5908859_1.fastq.gz
/vol1/fastq/SRR590/007/SRR5908747/SRR5908747_2.fastq.gz
/vol1/fastq/SRR590/008/SRR5908748/SRR5908748_2.fastq.gz
/vol1/fastq/SRR590/009/SRR5908749/SRR5908749_2.fastq.gz
/vol1/fastq/SRR590/000/SRR5908750/SRR5908750_2.fastq.gz
/vol1/fastq/SRR590/001/SRR5908751/SRR5908751_2.fastq.gz
/vol1/fastq/SRR590/002/SRR5908752/SRR5908752_2.fastq.gz
/vol1/fastq/SRR590/003/SRR5908753/SRR5908753_2.fastq.gz
/vol1/fastq/SRR590/004/SRR5908754/SRR5908754_2.fastq.gz
/vol1/fastq/SRR590/005/SRR5908755/SRR5908755_2.fastq.gz
/vol1/fastq/SRR590/006/SRR5908756/SRR5908756_2.fastq.gz
/vol1/fastq/SRR590/007/SRR5908757/SRR5908757_2.fastq.gz
/vol1/fastq/SRR590/008/SRR5908758/SRR5908758_2.fastq.gz
/vol1/fastq/SRR590/009/SRR5908759/SRR5908759_2.fastq.gz
/vol1/fastq/SRR590/000/SRR5908760/SRR5908760_2.fastq.gz
/vol1/fastq/SRR590/001/SRR5908761/SRR5908761_2.fastq.gz
/vol1/fastq/SRR590/002/SRR5908762/SRR5908762_2.fastq.gz
/vol1/fastq/SRR590/003/SRR5908763/SRR5908763_2.fastq.gz
/vol1/fastq/SRR590/004/SRR5908764/SRR5908764_2.fastq.gz
/vol1/fastq/SRR590/005/SRR5908765/SRR5908765_2.fastq.gz
/vol1/fastq/SRR590/006/SRR5908766/SRR5908766_2.fastq.gz
/vol1/fastq/SRR590/007/SRR5908767/SRR5908767_2.fastq.gz
/vol1/fastq/SRR590/008/SRR5908768/SRR5908768_2.fastq.gz
/vol1/fastq/SRR590/009/SRR5908769/SRR5908769_2.fastq.gz
/vol1/fastq/SRR590/000/SRR5908770/SRR5908770_2.fastq.gz
/vol1/fastq/SRR590/001/SRR5908771/SRR5908771_2.fastq.gz
/vol1/fastq/SRR590/002/SRR5908772/SRR5908772_2.fastq.gz
/vol1/fastq/SRR590/003/SRR5908773/SRR5908773_2.fastq.gz
/vol1/fastq/SRR590/004/SRR5908774/SRR5908774_2.fastq.gz
/vol1/fastq/SRR590/005/SRR5908775/SRR5908775_2.fastq.gz
/vol1/fastq/SRR590/006/SRR5908776/SRR5908776_2.fastq.gz
/vol1/fastq/SRR590/007/SRR5908777/SRR5908777_2.fastq.gz
/vol1/fastq/SRR590/008/SRR5908778/SRR5908778_2.fastq.gz
/vol1/fastq/SRR590/009/SRR5908779/SRR5908779_2.fastq.gz
/vol1/fastq/SRR590/000/SRR5908780/SRR5908780_2.fastq.gz
/vol1/fastq/SRR590/001/SRR5908781/SRR5908781_2.fastq.gz
/vol1/fastq/SRR590/002/SRR5908782/SRR5908782_2.fastq.gz
/vol1/fastq/SRR590/003/SRR5908783/SRR5908783_2.fastq.gz
/vol1/fastq/SRR590/004/SRR5908784/SRR5908784_2.fastq.gz
/vol1/fastq/SRR590/005/SRR5908785/SRR5908785_2.fastq.gz
/vol1/fastq/SRR590/006/SRR5908786/SRR5908786_2.fastq.gz
/vol1/fastq/SRR590/007/SRR5908787/SRR5908787_2.fastq.gz
/vol1/fastq/SRR590/008/SRR5908788/SRR5908788_2.fastq.gz
/vol1/fastq/SRR590/009/SRR5908789/SRR5908789_2.fastq.gz
/vol1/fastq/SRR590/000/SRR5908790/SRR5908790_2.fastq.gz
/vol1/fastq/SRR590/001/SRR5908791/SRR5908791_2.fastq.gz
/vol1/fastq/SRR590/002/SRR5908792/SRR5908792_2.fastq.gz
/vol1/fastq/SRR590/003/SRR5908793/SRR5908793_2.fastq.gz
/vol1/fastq/SRR590/004/SRR5908794/SRR5908794_2.fastq.gz
/vol1/fastq/SRR590/005/SRR5908795/SRR5908795_2.fastq.gz
/vol1/fastq/SRR590/006/SRR5908796/SRR5908796_2.fastq.gz
/vol1/fastq/SRR590/007/SRR5908797/SRR5908797_2.fastq.gz
/vol1/fastq/SRR590/008/SRR5908798/SRR5908798_2.fastq.gz
/vol1/fastq/SRR590/009/SRR5908799/SRR5908799_2.fastq.gz
/vol1/fastq/SRR590/000/SRR5908800/SRR5908800_2.fastq.gz
/vol1/fastq/SRR590/001/SRR5908801/SRR5908801_2.fastq.gz
/vol1/fastq/SRR590/002/SRR5908802/SRR5908802_2.fastq.gz
/vol1/fastq/SRR590/003/SRR5908803/SRR5908803_2.fastq.gz
/vol1/fastq/SRR590/004/SRR5908804/SRR5908804_2.fastq.gz
/vol1/fastq/SRR590/005/SRR5908805/SRR5908805_2.fastq.gz
/vol1/fastq/SRR590/006/SRR5908806/SRR5908806_2.fastq.gz
/vol1/fastq/SRR590/007/SRR5908807/SRR5908807_2.fastq.gz
/vol1/fastq/SRR590/008/SRR5908808/SRR5908808_2.fastq.gz
/vol1/fastq/SRR590/009/SRR5908809/SRR5908809_2.fastq.gz
/vol1/fastq/SRR590/000/SRR5908810/SRR5908810_2.fastq.gz
/vol1/fastq/SRR590/001/SRR5908811/SRR5908811_2.fastq.gz
/vol1/fastq/SRR590/002/SRR5908812/SRR5908812_2.fastq.gz
/vol1/fastq/SRR590/003/SRR5908813/SRR5908813_2.fastq.gz
/vol1/fastq/SRR590/004/SRR5908814/SRR5908814_2.fastq.gz
/vol1/fastq/SRR590/005/SRR5908815/SRR5908815_2.fastq.gz
/vol1/fastq/SRR590/006/SRR5908816/SRR5908816_2.fastq.gz
/vol1/fastq/SRR590/007/SRR5908817/SRR5908817_2.fastq.gz
/vol1/fastq/SRR590/008/SRR5908818/SRR5908818_2.fastq.gz
/vol1/fastq/SRR590/009/SRR5908819/SRR5908819_2.fastq.gz
/vol1/fastq/SRR590/000/SRR5908820/SRR5908820_2.fastq.gz
/vol1/fastq/SRR590/001/SRR5908821/SRR5908821_2.fastq.gz
/vol1/fastq/SRR590/002/SRR5908822/SRR5908822_2.fastq.gz
/vol1/fastq/SRR590/003/SRR5908823/SRR5908823_2.fastq.gz
/vol1/fastq/SRR590/004/SRR5908824/SRR5908824_2.fastq.gz
/vol1/fastq/SRR590/005/SRR5908825/SRR5908825_2.fastq.gz
/vol1/fastq/SRR590/006/SRR5908826/SRR5908826_2.fastq.gz
/vol1/fastq/SRR590/007/SRR5908827/SRR5908827_2.fastq.gz
/vol1/fastq/SRR590/008/SRR5908828/SRR5908828_2.fastq.gz
/vol1/fastq/SRR590/009/SRR5908829/SRR5908829_2.fastq.gz
/vol1/fastq/SRR590/000/SRR5908830/SRR5908830_2.fastq.gz
/vol1/fastq/SRR590/001/SRR5908831/SRR5908831_2.fastq.gz
/vol1/fastq/SRR590/002/SRR5908832/SRR5908832_2.fastq.gz
/vol1/fastq/SRR590/003/SRR5908833/SRR5908833_2.fastq.gz
/vol1/fastq/SRR590/004/SRR5908834/SRR5908834_2.fastq.gz
/vol1/fastq/SRR590/005/SRR5908835/SRR5908835_2.fastq.gz
/vol1/fastq/SRR590/006/SRR5908836/SRR5908836_2.fastq.gz
/vol1/fastq/SRR590/007/SRR5908837/SRR5908837_2.fastq.gz
/vol1/fastq/SRR590/008/SRR5908838/SRR5908838_2.fastq.gz
/vol1/fastq/SRR590/009/SRR5908839/SRR5908839_2.fastq.gz
/vol1/fastq/SRR590/000/SRR5908840/SRR5908840_2.fastq.gz
/vol1/fastq/SRR590/001/SRR5908841/SRR5908841_2.fastq.gz
/vol1/fastq/SRR590/002/SRR5908842/SRR5908842_2.fastq.gz
/vol1/fastq/SRR590/003/SRR5908843/SRR5908843_2.fastq.gz
/vol1/fastq/SRR590/004/SRR5908844/SRR5908844_2.fastq.gz
/vol1/fastq/SRR590/005/SRR5908845/SRR5908845_2.fastq.gz
/vol1/fastq/SRR590/006/SRR5908846/SRR5908846_2.fastq.gz
/vol1/fastq/SRR590/007/SRR5908847/SRR5908847_2.fastq.gz
/vol1/fastq/SRR590/008/SRR5908848/SRR5908848_2.fastq.gz
/vol1/fastq/SRR590/009/SRR5908849/SRR5908849_2.fastq.gz
/vol1/fastq/SRR590/000/SRR5908850/SRR5908850_2.fastq.gz
/vol1/fastq/SRR590/001/SRR5908851/SRR5908851_2.fastq.gz
/vol1/fastq/SRR590/002/SRR5908852/SRR5908852_2.fastq.gz
/vol1/fastq/SRR590/003/SRR5908853/SRR5908853_2.fastq.gz
/vol1/fastq/SRR590/004/SRR5908854/SRR5908854_2.fastq.gz
/vol1/fastq/SRR590/005/SRR5908855/SRR5908855_2.fastq.gz
/vol1/fastq/SRR590/006/SRR5908856/SRR5908856_2.fastq.gz
/vol1/fastq/SRR590/007/SRR5908857/SRR5908857_2.fastq.gz
/vol1/fastq/SRR590/008/SRR5908858/SRR5908858_2.fastq.gz
/vol1/fastq/SRR590/009/SRR5908859/SRR5908859_2.fastq.gz


ascp -v -QT -l 400m -P33001 -k1 -i /mnt/data/user_data/zlu/mytools/miniconda3/envs/sra-tools/etc/asperaweb_id_dsa.openssh --mode recv --host fasp.sra.ebi.ac.uk --user era-fasp --file-list  ncbi_NPC.txt  ./raw_data


*********************************************************************************************************************************
****************************************************ALL_PCA*****************************************************************************
*********************************************************************************************************************************

############################
#move to bamfiles path  #####
#move to bamfiles path  #####
#############################
cd /mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/star_mm10_out


vi sample_sampletable_all.csv

ids,sample,deal,source,order
1,ctrl_1.Aligned.sortedByCoord.out.bam,ctrl_1,ctrl,primary,1
2,ctrl_2.Aligned.sortedByCoord.out.bam,ctrl_2,ctrl,primary,2
3,ctrl_3.Aligned.sortedByCoord.out.bam,ctrl_3,ctrl,primary,3
4,LMP1_1.Aligned.sortedByCoord.out.bam,LMP1_1,LMP1,primary,4
5,LMP1_2.Aligned.sortedByCoord.out.bam,LMP1_2,LMP1,primary,5
6,LMP1_3.Aligned.sortedByCoord.out.bam,LMP1_3,LMP1,primary,6
7,Lung_ctrl_1.Aligned.sortedByCoord.out.bam,Lung_ctrl_1,ctrl,metastasis,7
8,Lung_ctrl_2.Aligned.sortedByCoord.out.bam,Lung_ctrl_2,ctrl,metastasis,8
9,Lung_LMP1_1.Aligned.sortedByCoord.out.bam,Lung_LMP1_1,LMP1,metastasis,9
10,Lung_LMP1_2.Aligned.sortedByCoord.out.bam,Lung_LMP1_2,LMP1,metastasis,10
11,Lung_LMP1_3.Aligned.sortedByCoord.out.bam,Lung_LMP1_3,LMP1,metastasis,11

############################
#work in R environment  #####
#work in R environment  #####
#############################


project <- c("all_PCA")
deal_design <- c("LMP1","ctrl")
significant_cutoff <- c(1)
organism <- "mouse"
##物种为人就是hsa，是老鼠就是mouse

file_path <- "/mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/all_PCA/"
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
csvfile <- file.path(indir, "sample_sampletable_all.csv")
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
getwd()

colData(se) <- DataFrame(sampleTable)
dds <- DESeqDataSet(se, design = ~deal + source)
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

#发现收得较早的样LMP1_1和ctrl聚到一起，故决定踢掉LMP1_1这个样。

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

*********************************************************************************************************************************
****************************************************LMP1_metastasis_vs_primary*****************************************************************************
*********************************************************************************************************************************

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

#发现收得较早的样LMP1_1和ctrl聚到一起，故决定踢掉LMP1_1这个样。

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

#padj<0.05 3733
res_1 <- na.omit(res_1)
y <- res_1[c(1:nrow(res_1)),6]
res_1 <- res_1[with(res_1,y<0.05),]
y <- res_1[c(1:nrow(res_1)),2]

#significant_cutoff <- c(0.5)

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
#up的gene为0，将logFC阈值设置为0.5的时候则能够做出图

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


*********************************************************************************************************************************
****************************************************ctrl_metastasis_vs_primary*****************************************************************************
*********************************************************************************************************************************

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


ids,sample,source,order
1,ctrl_2.Aligned.sortedByCoord.out.bam,ctrl_2,primary,1
2,ctrl_3.Aligned.sortedByCoord.out.bam,ctrl_3,primary,2
3,Lung_ctrl_1.Aligned.sortedByCoord.out.bam,Lung_ctrl_1,metastasis,3
4,Lung_ctrl_2.Aligned.sortedByCoord.out.bam,Lung_ctrl_2,metastasis,4


ids,sample,source,order
1,ctrl_1.Aligned.sortedByCoord.out.bam,ctrl_1,primary,1
2,ctrl_3.Aligned.sortedByCoord.out.bam,ctrl_3,primary,2
3,Lung_ctrl_1.Aligned.sortedByCoord.out.bam,Lung_ctrl_1,metastasis,3
4,Lung_ctrl_2.Aligned.sortedByCoord.out.bam,Lung_ctrl_2,metastasis,4


############################
#work in R environment  #####
#work in R environment  #####
#############################


project <- c("ctrl_metastasis_vs_primary")
deal_design <- c("metastasis","primary")
significant_cutoff <- c(1)
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
#cbind： 根据列进行合并，即叠加所有列，m列的矩阵与n列的矩阵cbind()最后变成m+n列，合并前提：cbind(a, b)中矩阵a、b的行数必需相符
#rbind： 根据行进行合并，就是行的叠加，m行的矩阵与n行的矩阵rbind()最后变成m+n行，合并前提：rbind(a, b)中矩阵a、b的列数必需相符
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

#padj<0.05的只有28 genes, pvalue<0.05 有579
res_1 <- na.omit(res_1)
y <- res_1[c(1:nrow(res_1)),5]
res_1 <- res_1[with(res_1,y<0.05),]
y <- res_1[c(1:nrow(res_1)),2]
#FC卡0.5
significant_cutoff <- c(0.5)

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
#up的gene为0，将logFC阈值设置为0.5的时候则能够做出图

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

@@@@@@@@@@@@@@@@@@@KEGG up top10

res_1 <- read.csv("/mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/ctrl_metastasis_vs_primary/ctrl_metastasis_vs_primary_04_metastasis_v_primary_symbol_and_anno.csv",row.names=1)
significant_cutoff <- c(0.5)

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

x <- KEGGupres_1[order(KEGGupres_1$GeneRatio,decreasing = TRUE),]
x <- x[c(10:15,1:9),]
y=new("enrichResult",
      result=x)
y
dotplot(y,orderBy ='GeneRatio',showCategory=10)

f <- dotplot(y, orderBy ='GeneRatio',showCategory=10)
ggsave("/mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/ctrl_metastasis_vs_primary/ctrl_metastasis_vs_primaryl_KEGG_UP_top10.svg", plot=f,width = 6,height = 5,dpi=1080)

y
dotplot(y,orderBy ='GeneRatio',showCategory=5)











@@@@@@@@@@@@@@@@@@@





#卡pvalue
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




*******************卡pvalue<0.05的基因list   579*******************
res_1 <- read.csv("/mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/ctrl_metastasis_vs_primary/ctrl_metastasis_vs_primary_pvalue0.05_578.csv")
res_1 <- na.omit(res_1)
y <- res_1[c(1:nrow(res_1)),2]
significant_cutoff <- c(0)
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

write.csv(KEGGupres_1, file = "ctrl_metastasis_vs_primary_pvalue0.05_578_KEGGup.csv")
write.csv(KEGGdownres_1, file = "ctrl_metastasis_vs_primary_pvalue0.05_578_KEGGdown.csv")


pdf(file = "ctrl_metastasis_vs_primary_pvalue0.05_578_KEGGup.pdf")
dotplot(KEGGupres_1, showCategory=20)
dev.off()
pdf(file = "ctrl_metastasis_vs_primary_pvalue0.05_578_KEGGdwon.pdf")
dotplot(KEGGdownres_1, showCategory=20)
dev.off()


#卡pvalue
pdf(file = KEGGres_1_all_pdf)
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
#无

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
write.csv(GOupres_1_all,file="ctrl_metastasis_vs_primary_pvalue0.05_578_GOup.csv")



ff <- barplot(GOupres_1_all, split="ONTOLOGY",showCategory=15) + facet_grid(ONTOLOGY~., scale="free")
ff1 <- ff + theme(axis.text.y = element_text(size = 8)) +labs(title = NULL)
pdf(file = "ctrl_metastasis_vs_primary_pvalue0.05_578_GOup.pdf")
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
write.csv(GOdownres_1_all,file="ctrl_metastasis_vs_primary_pvalue0.05_578_GOdwon.csv")


ff <- barplot(GOdownres_1_all, split="ONTOLOGY",showCategory=15) + facet_grid(ONTOLOGY~., scale="free")
ff1 <- ff + theme(axis.text.y = element_text(size = 8)) +labs(title = NULL)
pdf(file = "ctrl_metastasis_vs_primary_pvalue0.05_578_GOdown.pdf")
ff1
dev.off()
#无



*********gsea输入文件准备***********
setwd("/mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/ctrl_metastasis_vs_primary/ctrl_metastasis_vs_primary_pvalue0.05_578")
gsea_data <- read.csv("/mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/ctrl_metastasis_vs_primary/ctrl_metastasis_vs_primary_pvalue0.05_578.csv")
gsea_data$Description <- c("NULL")
gsea_data <- gsea_data[,c(18,21,10,11,7,8,9)]
nrow(gsea_data)
gsea_data <- na.omit(gsea_data)
nrow(gsea_data)
colnames(gsea_data)
colnames(gsea_data) <- c("NAME","Description","Lung_ctrl_1","Lung_ctrl_2","ctrl_1","ctrl_2","ctrl_3")
write.csv(gsea_data,file="/mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/ctrl_metastasis_vs_primary/ctrl_metastasis_vs_primary_pvalue0.05_578/gsea_data.csv")
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
-out ./598_ctrl_metastasis_vs_primary_h -gui false

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
-out ./598_ctrl_metastasis_vs_primary_c1 -gui false

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
-out ./598_ctrl_metastasis_vs_primary_c2 -gui false

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
-out ./598_ctrl_metastasis_vs_primary_c3 -gui false

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
-out ./598_ctrl_metastasis_vs_primary_c4 -gui false

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
-out ./598_ctrl_metastasis_vs_primary_c5 -gui false

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
-out ./598_ctrl_metastasis_vs_primary_c6 -gui false

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
-out ./598_ctrl_metastasis_vs_primary_c7 -gui false


bash gsea_shell.sh

**************************************************heatmap*****************************************************************************
$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$common_49genes
commonup_49genes <- read.csv("/mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/ctrl_metastasis_vs_primary/49genes_heatmap.csv",row.names = 1)
head(commonup_49genes)
#删掉log2FC这一列
commonup_49genes <- commonup_49genes[,-6]
colnames(commonup_49genes) <- c("ctrl_1","ctrl_2","ctrl_3","LMP1_2","LMP1_3")
zscore <- t(apply(commonup_49genes, 1, function(x) (x-mean(x))/sd(x)))
zscore <- na.omit(zscore)
head(zscore)

library(Seurat)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
require(RColorBrewer)
source("/mnt/data/user_data/xiangyu/programme/R_PACKAGES/my_code/MyBestFunction_scRNA.R")
source("/mnt/data/user_data/xiangyu/programme/R_PACKAGES/my_code/Pseudo_CNV_series.R")

pdf("commonup_49genes_heatmap.pdf",height = 14)
pheatmap(zscore,color = colorRampPalette(c("navy", "white", "firebrick3"))(50),show_rownames=TRUE,cluster_row = FALSE,cluster_col= FALSE,border=FALSE)
dev.off()

$$$$$$$$$$$$$$$$$$$$$$$$$$$all_gene
all_gene <- read.csv("/mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/ctrl_metastasis_vs_primary/all_gene_orderbylog2FC_heatmap.csv")
head(all_gene)
all_gene <- all_gene[!duplicated(all_gene$X),]
rownames(all_gene) <- all_gene[,1]
all_gene <- all_gene[,-1]
#删掉log2FC这一列
all_gene <- all_gene[,-6]


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

pdf("all_gene.pdf",height = 14)
pheatmap(zscore,color = colorRampPalette(c("navy", "white", "firebrick3"))(50),show_rownames=FALSE,cluster_row = FALSE,cluster_col= FALSE,border=FALSE)
dev.off()

SeuratObject <- CreateSeuratObject(counts = all_gene, project = "WXD")
gene <- rownames(all_gene)
pdf("all_genes_heatmap.pdf")
XY_heatmap(seurat_obj=SeuratObject,group="orig.ident",genes=gene,all_num=FALSE,new_names=NULL,labels_rot=90,assay_sel="RNA",color=colorRampPalette(brewer.pal(10, "RdBu"))(101),min_and_max_cut=1.2,show_row_names=FALSE,scale=FALSE,label_size=0,mark_gene=c("Cldn2","Palm2","Sulf2","Prex1","Spata13","Rbp1","Ctgf","Ffar4","Bmp4","Aqp5","Dapk1","Frmd4a","Amotl2","Cdh13","Lbh","Uap1l1","Grina","Lamtor2","Igsf8","Traf3ip2"))
dev.off()
 





********************************************************************tissue+organoids*******************************
#8922  /mnt/data/sequencedata/RNAseq/RNAseq_6_chenxuelan_organoids_25yangben/bam  #GRcm38
#scp -r /mnt/data/sequencedata/RNAseq/RNAseq_6_chenxuelan_organoids_25yangben/bam/nasopha_mouse1*.bam zhanglu@172.16.1.21:/mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/star_mm10_out/
vi config 

normal_organoid_1 /mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/organoid_raw_data/nasopha_mouse1_organoid_1.fq.gz /mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/organoid_raw_data/nasopha_mouse1_organoid_2.fq.gz
normal_organoid_2 /mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/organoid_raw_data/nasopha_mouse2_organoid_1.fq.gz /mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/organoid_raw_data/nasopha_mouse2_organoid_2.fq.gz
normal_organoid_3 /mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/organoid_raw_data/nasopha_mouse3_organoid_1.fq.gz /mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/organoid_raw_data/nasopha_mouse3_organoid_2.fq.gz
normal_tissue_1 /mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/organoid_raw_data/nasopha_mouse1_tissue_1.fq.gz /mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/organoid_raw_data/nasopha_mouse1_tissue_2.fq.gz
normal_tissue_2 /mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/organoid_raw_data/nasopha_mouse2_tissue_1.fq.gz /mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/organoid_raw_data/nasopha_mouse2_tissue_2.fq.gz
 

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

cd /mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/star_mm10_out
#构建index
ls normal*.sortedByCoord.out.bam |while read id
do
samtools index ./$id
done

$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$tissue

cd /mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/star_mm10_out

vi sample_sampletable_ctrl_normal_tissue.csv

ids,sample,deal,order
1,ctrl_1.Aligned.sortedByCoord.out.bam,ctrl_1,P_ctrl,1
2,ctrl_2.Aligned.sortedByCoord.out.bam,ctrl_2,P_ctrl,2
3,ctrl_3.Aligned.sortedByCoord.out.bam,ctrl_3,P_ctrl,3
4,normal_tissue_1.Aligned.sortedByCoord.out.bam,normal_tissue_1,normal_tissue,4
5,normal_tissue_2.Aligned.sortedByCoord.out.bam,normal_tissue_2,normal_tissue,5

############################
#work in R environment  #####
#work in R environment  #####
#############################


project <- c("P_ctrl_VS_normal_tissue")
deal_design <- c("P_ctrl","normal_tissue")
significant_cutoff <- c(1)
organism <- "mouse"
##物种为人就是hsa，是老鼠就是mouse

file_path <- "/mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/P_ctrl_vs_nomal_tissue/"
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
csvfile <- file.path(indir, "sample_sampletable_ctrl_normal_tissue.csv")
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

all_TF <- read.csv(file = "P_ctrl_VS_normal_organoid_05_P_ctrl_v_normal_organoid_count_tpm_symbol_and_anno.csv",row.names = 1)
aa <- as.character(rownames(all_TF)) %>% convert_mouse_to_human_symbols()
all_TF$human_symbol <- as.character(aa)
matrix <- all_TF[!duplicated(all_TF$human_symbol),]
write.csv(matrix,"P_ctrl_VS_normal_organoid_05_P_ctrl_v_normal_organoid_count_tpm_symbol_and_anno_human.csv")





res_1 <- read.csv("/mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/P_ctrl_vs_nomal_tissue/P_ctrl_VS_normal_tissue_04_P_ctrl_v_normal_tissue_symbol_and_anno.csv")

rownames(res_1) <- res_1[,1]
res_1 <- res_1[,-1]
significant_cutoff <- c(1)
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
gsea_data <- read.csv("/mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/P_ctrl_vs_nomal_tissue/P_ctrl_VS_normal_tissue_05_P_ctrl_v_normal_tissue_count_tpm_symbol_and_anno.csv")
gsea_data$Description <- c("NULL")
gsea_data <- gsea_data[,c(18,21,7,8,9,10,11)]
nrow(gsea_data)
gsea_data <- na.omit(gsea_data)
nrow(gsea_data)
colnames(gsea_data)
colnames(gsea_data) <- c("NAME","Description","P_ctrl_1","P_ctrl_2","P_ctrl_3","normal_tissue_1","normal_tissue_2")
write.csv(gsea_data,file="/mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/P_ctrl_vs_nomal_tissue/gsea_data.csv")
#在excel中改好格式(就以Ensembl id为行名)，然后直接保存为txt格式，再导入R修改文件名为.gct即可
#加上  
       #1.2
       #总基因行数 总样本数

vi gsea.cls
5 2 1
# P_ctrl normal_tissue
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
-out ./P_ctrl_VS_normal_tissue_h -gui false

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
-out ./P_ctrl_VS_normal_tissue_c1 -gui false

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
-out ./P_ctrl_VS_normal_tissue_c2 -gui false

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
-out ./P_ctrl_VS_normal_tissue_c3 -gui false

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
-out ./P_ctrl_VS_normal_tissue_c4 -gui false

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
-out ./P_ctrl_VS_normal_tissue_c5 -gui false

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
-out ./P_ctrl_VS_normal_tissue_c6 -gui false

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
-out ./P_ctrl_VS_normal_tissue_c7 -gui false


bash gsea_shell.sh


$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$organoids  vs P_ctrl

cd /mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/star_mm10_out

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
*********************DESeq2 normalize***********
colData(se) <- DataFrame(sampleTable)
library("DESeq2")
dds <- DESeqDataSet(se, design = ~ deal)
countdata <- assay(se)
coldata <- colData(se)
DESeq_counts <- DESeq(dds)
normalized_DEseq <- counts(DESeq_counts, normalized=TRUE)
DEseq_res <- results(DESeq_counts,contrast=c("deal","P_ctrl","normal_organoid"))    #最后一个是CTRL
res_1 <- cbind(normalized_DEseq,DEseq_res)

	load(file="/mnt/data/user_data/xiangyu/workshop/WORKFLOW_RNAseq/ebg_mm10.RData")
	load(file="/mnt/data/user_data/xiangyu/workshop/WORKFLOW_RNAseq/txdb_mm10.RData")
res_1$ENSEMBL <- mapIds(x = anno_data,
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
all_summry <- cbind(countdata,res_1)
names(all_summry) <- c("ctrl_1","ctrl_2","ctrl_3","normal_1","normal_2","normal_3","DESeq2_ctrl_1","DESeq2_ctrl_2","DESeq2_ctrl_3","DESeq2_normal_1","DESeq2_normal_2","DESeq2_normal_3","baseMean", "log2FoldChange", "lfcSE", "stat","pvalue","padj", "ENSEMBL","entrez","GENENAME")
write.csv(all_summry, "renew_DEseq2normalized_P_ctrl_VS_normal_organoid_allsummry.csv")

gsea_data <- read.csv("/mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/P_ctrl_vs_nomal_organoid/DESeq2_normalize_gsea/renew_DEseq2normalized_P_ctrl_VS_normal_organoid_allsummry.csv")
gsea_data$Description <- c("NULL")
gsea_data <- gsea_data[,c(20,23,8,9,10,11,12,13)]
nrow(gsea_data)
gsea_data <- na.omit(gsea_data)
nrow(gsea_data)
colnames(gsea_data)
colnames(gsea_data) <- c("NAME","Description","DESeq2_ctrl_1","DESeq2_ctrl_2","DESeq2_ctrl_3","DESeq2_normal_1","DESeq2_normal_2","DESeq2_normal_3")
write.csv(gsea_data,file="/mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/P_ctrl_vs_nomal_organoid/DESeq2_normalize_gsea/gsea_data.csv")
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
-out ./DESeq2_P_ctrl_VS_normal_organoid_h -gui false

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
-out ./DESeq2_P_ctrl_VS_normal_organoid_c2 -gui false

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
-out ./DESeq2_P_ctrl_VS_normal_organoid_c5 -gui false


bash gsea_shell.sh

***********************************************

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



significant_cutoff <- c(0.5)
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


***********************************************************************************************************KEGG top5 画图
res_1 <- read.csv("/mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/P_ctrl_vs_nomal_organoid/P_ctrl_VS_normal_organoid_04_P_ctrl_v_normal_organoid_symbol_and_anno.csv")
significant_cutoff <- c(0.5)

rownames(res_1) <- res_1[,1]
res_1 <- res_1[,-1]
res_1 <- na.omit(res_1)
y <- res_1[c(1:nrow(res_1)),6]
res_1 <- res_1[with(res_1,y<0.05),]
y <- res_1[c(1:nrow(res_1)),2]


x <- KEGGupres_1[order(KEGGupres_1$p.adjust,decreasing = FALSE),]
z <- x[1:5,]

dotplot(y,orderBy ='GeneRatio',showCategory=5)

z$p.adjust = -log10(z$p.adjust)
  z$p.adjust=z$p.adjust*1 
  z=z[order(z$p.adjust,decreasing = F),]
  z$group=1
  library("ggplot2")
g_kegg<- ggplot(z, aes(x=reorder(Description,order(p.adjust, decreasing = F)), y=p.adjust,fill=group)) + 
    geom_bar(stat="identity") + 
    scale_fill_gradient(low="blue",high="red",guide = "none") + 
    scale_x_discrete(name ="Pathway names") +
    scale_y_continuous(name ="-log10p_adjust") +
    coord_flip() +
    ggtitle("TMP vs. Normal up KEGG")
  print(g_kegg)
ggsave("/mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/P_ctrl_vs_nomal_organoid/P_TMP_v_normal_KEGG_UP_top5.svg", plot=g_kegg,width = 6,height = 5,dpi=1080)
***************************************************************************************************************

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





$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$organoids  vs P_LMP1

cd /mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/star_mm10_out

vi sample_sampletable_LMP1_normal_organoid.csv

ids,sample,deal,order
1,LMP1_1.Aligned.sortedByCoord.out.bam,LMP1_1,P_LMP1,1
2,LMP1_2.Aligned.sortedByCoord.out.bam,LMP1_2,P_LMP1,2
3,LMP1_3.Aligned.sortedByCoord.out.bam,LMP1_3,P_LMP1,3
4,normal_organoid_1.Aligned.sortedByCoord.out.bam,normal_organoid_1,normal_organoid,4
5,normal_organoid_2.Aligned.sortedByCoord.out.bam,normal_organoid_2,normal_organoid,5
6,normal_organoid_3.Aligned.sortedByCoord.out.bam,normal_organoid_3,normal_organoid,6

############################
#work in R environment  #####
#work in R environment  #####
#############################


project <- c("P_LMP1_VS_normal_organoid")
deal_design <- c("P_LMP1","normal_organoid")
organism <- "mouse"
##物种为人就是hsa，是老鼠就是mouse

file_path <- "/mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/P_LMP1_vs_nomal/"
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
csvfile <- file.path(indir, "sample_sampletable_LMP1_normal_organoid.csv")
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


***********************************************

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

*********************DESeq2 normalize***********
colData(se) <- DataFrame(sampleTable)
library("DESeq2")
dds <- DESeqDataSet(se, design = ~ deal)
countdata <- assay(se)
coldata <- colData(se)
DESeq_counts <- DESeq(dds)
normalized_DEseq <- counts(DESeq_counts, normalized=TRUE)
DEseq_res <- results(DESeq_counts,contrast=c("deal","P_LMP1","normal_organoid"))    #最后一个是CTRL
res_1 <- cbind(normalized_DEseq,DEseq_res)

res_1$ENSEMBL <- mapIds(x = anno_data,
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
all_summry <- cbind(countdata,res_1)
names(all_summry) <- c("LMP1_1","LMP1_2","LMP1_3","normal_1","normal_2","normal_3","DESeq2_LMP1_1","DESeq2_LMP1_2","DESeq2_LMP1_3","DESeq2_normal_1","DESeq2_normal_2","DESeq2_normal_3","baseMean", "log2FoldChange", "lfcSE", "stat","pvalue","padj", "ENSEMBL","entrez","GENENAME")
write.csv(all_summry, "renew_DEseq2normalized_P_LMP1_VS_normal_organoid_allsummry.csv")

res_1 <- read.csv(file = "renew_DEseq2normalized_P_LMP1_VS_normal_organoid_allsummry.csv",row.names = 1)
res_1 <- res_1[,13:21]
significant_cutoff <- c(0.5)
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
gsea_data <- read.csv("/mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/P_LMP1_vs_nomal/renew_DEseq2normalized_P_LMP1_VS_normal_organoid_allsummry.csv")
gsea_data$Description <- c("NULL")
gsea_data <- gsea_data[,c(20,23,8,9,10,11,12,13)]
nrow(gsea_data)
gsea_data <- na.omit(gsea_data)
nrow(gsea_data)
colnames(gsea_data)
colnames(gsea_data) <- c("NAME","Description","DESeq2_LMP1_1","DESeq2_LMP1_2","DESeq2_LMP1_3","DESeq2_normal_1","DESeq2_normal_2","DESeq2_normal_3")
write.csv(gsea_data,file="/mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/P_LMP1_vs_nomal/gsea_data.csv")
#在excel中改好格式(就以Ensembl id为行名)，然后直接保存为txt格式，再导入R修改文件名为.gct即可
#加上  
       #1.2
       #总基因行数 总样本数

vi gsea.cls

6 2 1
# P_LMP1 normal_organoid
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
-out ./P_LMP1_VS_normal_organoid_h -gui false

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
-out ./P_LMP1_VS_normal_organoid_c1 -gui false

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
-out ./P_LMP1_VS_normal_organoid_c2 -gui false

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
-out ./P_LMP1_VS_normal_organoid_c3 -gui false

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
-out ./P_LMP1_VS_normal_organoid_c4 -gui false

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
-out ./P_LMP1_VS_normal_organoid_c5 -gui false

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
-out ./P_LMP1_VS_normal_organoid_c6 -gui false

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
-out ./P_LMP1_VS_normal_organoid_c7 -gui false


bash gsea_shell.sh















**************************************************motif**************************************************************************
$$$$ motif 分析
$$$$ motif 分析
cd /mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/ctrl_metastasis_vs_primary/
ctrl_M_vs_P <- read.csv("/mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/ctrl_metastasis_vs_primary/ctrl_metastasis_vs_primary_05_metastasis_v_primary_count_tpm_symbol_and_anno.csv",row.names = 1)
ctrl_M_vs_P$HALLMARK_NOTCH_SIGNALING <- rownames(ctrl_M_vs_P)

library("dplyr")
ctrl_M_vs_P_up <- ctrl_M_vs_P %>% filter(metastasis_v_primary_pvalue <= 0.05) %>% filter(metastasis_v_primary_log2FoldChange >= 1)
ctrl_M_vs_P_down <- ctrl_M_vs_P %>% filter(metastasis_v_primary_pvalue <= 0.05) %>% filter(metastasis_v_primary_log2FoldChange <= -1)
ctrl_M_vs_P_up <- na.omit(ctrl_M_vs_P_up)
ctrl_M_vs_P_down <- na.omit(ctrl_M_vs_P_down)
library(Mus.musculus)
all_genes <- genes(Mus.musculus)
all_gene_TSS <- resize(all_genes,1)
TSS_2k <-promoters(all_gene_TSS, 2000, 2000)
TSS_2k <- as.data.frame(TSS_2k)
TSS_2k$entrez <- row.names(TSS_2k)
all_gene_TSS <- as.data.frame(all_gene_TSS)

#load(file="/mnt/data/user_data/xiangyu/workshop/WORKFLOW_RNAseq/ebg_mm10.RData")
#load(file="/mnt/data/user_data/xiangyu/workshop/WORKFLOW_RNAseq/txdb_mm10.RData")

#ctrl_M_vs_P_up$entrez <- mapIds(x = org.Mm.eg.db,
#            keys = rownames(ctrl_M_vs_P_up),
#            keytype ="ENSEMBL",
#            column ="ENTREZID",
#            multiVals="first")

ctrl_M_vs_P_up_TSS <- merge(ctrl_M_vs_P_up,TSS_2k,by = "entrez", all = FALSE)
head(ctrl_M_vs_P_up_TSS)
ctrl_M_vs_P_up_TSS <- ctrl_M_vs_P_up_TSS[,c(1,20:26)]
nrow(ctrl_M_vs_P_up_TSS)  #196
write.csv(ctrl_M_vs_P_up_TSS,"ctrl_M_vs_P_up_TSS.csv")

ctrl_M_vs_P_down_TSS <- merge(ctrl_M_vs_P_down,TSS_2k,by = "entrez", all = FALSE)
head(ctrl_M_vs_P_down_TSS)
ctrl_M_vs_P_down_TSS <- ctrl_M_vs_P_down_TSS[,c(1,20:26)]
head(ctrl_M_vs_P_down_TSS)
nrow(ctrl_M_vs_P_down_TSS) #144
write.csv(ctrl_M_vs_P_down_TSS,"ctrl_M_vs_P_down_TSS.csv")


mv ctrl_M_vs_P_up_TSS.txt ctrl_M_vs_P_up_TSS.bed
mv ctrl_M_vs_P_down_TSS.txt ctrl_M_vs_P_down_TSS.bed
#先看一下/mnt/data/user_data/zhaolei/project/zms_EZH2/human_data/路径下的bed格式，把文件改成对应的格式。在excel中删掉不需要的列，列名，存为txt格式，在重命名位.bed即可。
#下一步用homer 切换阿爆的账号使用
. /mnt/data/user_data/zhaolei/program/yes/bin/activate
findMotifsGenome.pl /mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/ctrl_metastasis_vs_primary/ctrl_M_vs_P_up_TSS.bed mm10 ctrl_M_vs_P_up_motif -len 8,10,12
findMotifsGenome.pl /mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/ctrl_metastasis_vs_primary/ctrl_M_vs_P_down_TSS.bed mm10 ctrl_M_vs_P_down_motif -len 8,10,12

****************************************
cd /mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/ctrl_metastasis_vs_primary/
primary_LMP1_vs_ctrl <- read.csv("/mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/primary_LMP1_vs_ctrl/primary_LMP1_05_LMP1_v_ctrl_count_tpm_symbol_and_anno.csv",row.names = 1)
primary_LMP1_vs_ctrl$HALLMARK_NOTCH_SIGNALING <- rownames(primary_LMP1_vs_ctrl)

library("dplyr")
primary_LMP1_vs_ctrl_up <- primary_LMP1_vs_ctrl %>% filter(LMP1_v_ctrl_padj <= 0.05) %>% filter(LMP1_v_ctrl_log2FoldChange >= 1)
primary_LMP1_vs_ctrl_down <- primary_LMP1_vs_ctrl %>% filter(LMP1_v_ctrl_padj <= 0.05) %>% filter(LMP1_v_ctrl_log2FoldChange <= -1)
primary_LMP1_vs_ctrl_up <- na.omit(primary_LMP1_vs_ctrl_up)
primary_LMP1_vs_ctrl_down <- na.omit(primary_LMP1_vs_ctrl_down)
library(Mus.musculus)
all_genes <- genes(Mus.musculus)
all_gene_TSS <- resize(all_genes,1)
TSS_2k <-promoters(all_gene_TSS, 2000, 2000)
TSS_2k <- as.data.frame(TSS_2k)
TSS_2k$entrez <- row.names(TSS_2k)
all_gene_TSS <- as.data.frame(all_gene_TSS)

primary_LMP1_vs_ctrl_up_TSS <- merge(primary_LMP1_vs_ctrl_up,TSS_2k,by = "entrez", all = FALSE)
head(primary_LMP1_vs_ctrl_up_TSS)
primary_LMP1_vs_ctrl_up_TSS <- primary_LMP1_vs_ctrl_up_TSS[,c(1,20:26)]
head(primary_LMP1_vs_ctrl_up_TSS)
nrow(primary_LMP1_vs_ctrl_up_TSS)  #834
write.csv(primary_LMP1_vs_ctrl_up_TSS,"primary_LMP1_vs_ctrl_up_TSS.csv")

primary_LMP1_vs_ctrl_down_TSS <- merge(primary_LMP1_vs_ctrl_down,TSS_2k,by = "entrez", all = FALSE)
head(primary_LMP1_vs_ctrl_down_TSS)
primary_LMP1_vs_ctrl_down_TSS <- primary_LMP1_vs_ctrl_down_TSS[,c(1,20:26)]
head(primary_LMP1_vs_ctrl_down_TSS)
nrow(primary_LMP1_vs_ctrl_down_TSS)  #1237
write.csv(primary_LMP1_vs_ctrl_down_TSS,"primary_LMP1_vs_ctrl_down_TSS.csv")


mv primary_LMP1_vs_ctrl_up_TSS.txt primary_LMP1_vs_ctrl_up_TSS.bed
mv primary_LMP1_vs_ctrl_down_TSS.txt primary_LMP1_vs_ctrl_down_TSS.bed

. /mnt/data/user_data/zhaolei/program/yes/bin/activate
findMotifsGenome.pl /mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/primary_LMP1_vs_ctrl/primary_LMP1_vs_ctrl_up_TSS.bed mm10 primary_LMP1_vs_ctrl_up_motif -len 8,10,12

findMotifsGenome.pl /mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/primary_LMP1_vs_ctrl/primary_LMP1_vs_ctrl_down_TSS.bed mm10 primary_LMP1_vs_ctrl_down_motif -len 8,10,12



***********************************************************************消除鼻咽和肺 器官之间的差异***********************************************************************
***********************************************************************消除鼻咽和肺 器官之间的差异***********************************************************************
***********************************************************************消除鼻咽和肺 器官之间的差异***********************************************************************
scp -r /mnt/data/sequencedata/RNAseq/RNAseq_6_chenxuelan_organoids_25yangben/bam/lung111_mouse* zhanglu@172.16.1.21:/mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/GRCm38_bam
scp -r /mnt/data/sequencedata/RNAseq/RNAseq_6_chenxuelan_organoids_25yangben/bam/nasopha_mouse*_LLH.bam zhanglu@172.16.1.21:/mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/GRCm38_bam

cd /mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/GRCm38_bam

vi sampletable_nasopha_vs_lung.csv

ids,sample,deal,order
1,lung111_mouse1_NFF.bam,lung_1,lung,1
2,lung111_mouse2_NFF.bam,lung_2,lung,2
3,lung111_mouse3_NFF.bam,lung_3,lung,3
4,nasopha_mouse1_LLH.bam,nasopha_1,nasopha,4
5,nasopha_mouse2_LLH.bam,nasopha_2,nasopha,5
6,nasopha_mouse3_LLH.bam,nasopha_3,nasopha,6

############################
#work in R environment  #####
#work in R environment  #####
#############################


project <- c("nasopha_vs_lung")
deal_design <- c("nasopha","lung")
organism <- "mouse"
##物种为人就是hsa，是老鼠就是mouse

file_path <- "/mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/nasopha_vs_lung/"
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
	load(file="/mnt/data/user_data/xiangyu/workshop/WORKFLOW_RNAseq/ebg_GRCm38.RData")
	load(file="/mnt/data/user_data/xiangyu/workshop/WORKFLOW_RNAseq/txdb_GRCm38.RData")
} else {
	anno_data=org.Hs.eg.db
	print("organism is human")
	load(file="/mnt/data/user_data/xiangyu/workshop/WORKFLOW_RNAseq/ebg_GRCh38.RData")
	load(file="/mnt/data/user_data/xiangyu/workshop/WORKFLOW_RNAseq/txdb_GRCh38.RData")
}


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

#发现收得较早的样LMP1_1和ctrl聚到一起，故决定踢掉LMP1_1这个样。

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
-out ./ctrl_M_vs_P_withoutlung&nasophary_h -gui false

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
-out ./ctrl_M_vs_P_withoutlung&nasophary_c1 -gui false

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
-out ./ctrl_M_vs_P_withoutlung&nasophary_c2 -gui false

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
-out ./ctrl_M_vs_P_withoutlung&nasophary_c3 -gui false

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
-out ./ctrl_M_vs_P_withoutlung&nasophary_c4 -gui false

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
-out ./ctrl_M_vs_P_withoutlung&nasophary_c5 -gui false

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
-out ./ctrl_M_vs_P_withoutlung&nasophary_c6 -gui false

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
-out ./ctrl_M_vs_P_withoutlung&nasophary_c7 -gui false

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
-out ./ctrl_M_vs_P_withoutlung&nasophary_metastasis -gui false

bash gsea_shell.sh


*******************************************************************************证明鼻咽organoid的可靠性
scp -r /mnt/data/sequencedata/RNAseq/RNAseq_6_chenxuelan_organoids_25yangben/bam/*.bam zhanglu@172.16.1.21:/mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/GRCm38_bam
******************************lung
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


*********************************************************************************************************************************
****************************************************ALL_PCA（tumor and  normal）*****************************************************************************
*********************************************************************************************************************************

############################
#move to bamfiles path  #####
#move to bamfiles path  #####
#############################
cd /mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/star_mm10_out


vi sample_sampletable_all_normal.csv

ids,sample,deal,source,order
1,ctrl_1.Aligned.sortedByCoord.out.bam,Pri_TMP_1,ctrl,primary,1
2,ctrl_2.Aligned.sortedByCoord.out.bam,Pri_TMP_2,ctrl,primary,2
3,ctrl_3.Aligned.sortedByCoord.out.bam,Pri_TMP_3,ctrl,primary,3
4,LMP1_2.Aligned.sortedByCoord.out.bam,Pri_TMPL_2,LMP1,primary,4
5,LMP1_3.Aligned.sortedByCoord.out.bam,Pri_TMPL_3,LMP1,primary,5
6,Lung_ctrl_1.Aligned.sortedByCoord.out.bam,Met_TMP_1,ctrl,metastasis,6
7,Lung_ctrl_2.Aligned.sortedByCoord.out.bam,Met_TMP_2,ctrl,metastasis,7
8,Lung_LMP1_1.Aligned.sortedByCoord.out.bam,Met_TMPL_1,LMP1,metastasis,8
9,Lung_LMP1_2.Aligned.sortedByCoord.out.bam,Met_TMPL_2,LMP1,metastasis,9
10,Lung_LMP1_3.Aligned.sortedByCoord.out.bam,Met_TMPL_3,LMP1,metastasis,10
11,normal_organoid_1.Aligned.sortedByCoord.out.bam,normal_1,normal,normal,11
12,normal_organoid_2.Aligned.sortedByCoord.out.bam,normal_2,normal,normal,12
13,normal_organoid_3.Aligned.sortedByCoord.out.bam,normal_3,normal,normal,13




############################
#work in R environment  #####
#work in R environment  #####
#############################


project <- c("all_PCA")
deal_design <- c("LMP1","ctrl")
significant_cutoff <- c(1)
organism <- "mouse"
##物种为人就是hsa，是老鼠就是mouse

file_path <- "/mnt/data/user_data/zlu/01_job/WXD_LMP1_rnaseq/all_PCA/"
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
csvfile <- file.path(indir, "sample_sampletable_all_normal.csv")
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
getwd()

colData(se) <- DataFrame(sampleTable)

#发现收得较早的样LMP1_1和ctrl聚到一起，故决定踢掉LMP1_1这个样。

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
count <- assay(se)
colnames(tpm) <- sampleTable$sample
colnames(count) <- sampleTable$sample
count_and_tpm <- cbind(count,tpm)

id <- c("Epha2","Tmprss2","Tmprss4","Slc3a2")
ZQ_gene <- count_and_tpm[rownames(count_and_tpm) %in% id,]
write.csv(ZQ_gene, "ZQ_gene_TPM.csv")
