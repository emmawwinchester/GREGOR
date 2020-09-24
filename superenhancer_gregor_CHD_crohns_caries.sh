#!/bin/bash
#SBATCH --job-name=consensusmouse
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=END
#SBATCH --mem=64G
#SBATCH --mail-user=wentworth@uchc.edu
#SBATCH -o consensus.out
#SBATCH -e consensus.err
export DIR=/home/CAM/ewentworth/cotney/analysis/gwas/superenhancers
ls /home/CAM/ewentworth/cotney/rawdata/chromatinsegmentations/superenhancers/human/*.bed > index_human.txt

for sample in *human/  
do
cd $sample
echo $sample
export NAME=`echo $sample | sed 's/\//-out.txt/g'`
awk 'NR==1{print;next}{if ($3=="NA") {print $0, 0} else {print $0, (log($2/$3)/log(2))}}' Statistic*.txt | tail -n +2 | sed 's/state_2/state2/g' | sed 's/_rconsensus//g' | sed 's/.bed//g' | sed 's/hg19-25pct-//g' | awk '{if ($1 ~ /CS|facial|cranio|reproducible|ooth|snout|palate/) {print $0, "craniofacial"} else if ($1 ~ /CD|Liver|Thymus/) {print $0, "immune"} else if ($1 ~ /heart|Ventricle|Atrium/) {print $0, "heart"} else {print $0, "noncraniofacial"}}' | sed 's/ /\t/g' | sed 's/_Superenhancer_sorted//g' | sed 's/-superenhancer//g' > $DIR/$NAME
echo $DIR/$NAME
cd $DIR
done

#col1: bed file col2:#observed Snps col3: expected Snps col4: pvalue (unadjusted) col5: log2(observed/expected) col6: tissue 
#change the NAs to 1s in R
#for consensus do p value correction in R! add a new column with corrected bonferroni p value
##need to replace all these spaces with tabs in each file
R
library(ggplot2)
library(ggrepel)
library(stringr)
library(writexl)
library(patchwork)
path = "~/cotney/analysis/gwas/superenhancers"
out.file<-""
file.names<-dir(path, pattern=(glob2rx("*human*out.txt")))
for (val in 1:length(file.names)){
x<-read.table(file.names[val], sep='\t', header=FALSE)
x[is.na(x)]<-1
w<-cbind(x, (p.adjust(x$V4, "BH")), (p.adjust(x$V4, "bonferroni")), -log10(p.adjust(x$V4, "BH")), -log10(p.adjust(x$V4, "bonferroni")))
x<-w
colnames(x)<-c("Sample", "ObservedSNPs", 'ExpectedSNPs', 'PValue', 'Log2SNPFoldEnrichment', 'Lineage', 'BHCorrectedPValue', 'BonferroniCorrectedPValue', 'BHCorrectedPValueLog', 'BonferroniCorrectedPValueLog')
write.table(x, str_remove(paste(file.names[val], '-table.txt', sep=''), "-out.txt"), sep='\t')
write_xlsx(x, str_remove(paste(file.names[val], '-table.xlsx', sep=''), "-out.txt"))
n<-head(x[order(-x$Log2SNPFoldEnrichment),], 15)
q<-head(x[order(-x$BHCorrectedPValueLog),], 15)
pdf(file=str_remove(paste(file.names[val], '.pdf', sep=''), "-out.txt"), height=8.5, width=11)
print(ggplot(x, aes(x=Log2SNPFoldEnrichment, y=BHCorrectedPValueLog, colour=Lineage, shape=Lineage, label=Sample)) + geom_point(aes(size=ObservedSNPs))+ geom_hline(yintercept=1.3) + ggtitle(str_replace_all(str_replace(str_remove(paste(file.names[val], ' SNPs in Human SuperEnhancers', sep=''), "-out.txt"), "-human", " "),"_", " ")) + geom_label_repel(aes(label=Sample), data=subset(x, Log2SNPFoldEnrichment >= n$Log2SNPFoldEnrichment[10] | BHCorrectedPValueLog >= q$BHCorrectedPValueLog[10])) + scale_color_manual(values=c("orange", "red", "grey", "grey")) + labs(x="Fold Enrichment (Log2)", y="BH Corrected P Value (-Log10)", colour="Lineage", size="Number of SNPs in SuperEnhancers", shape="Lineage"))
dev.off()
}
file.names<-dir(path, pattern=(glob2rx("*mouse*out.txt")))
for (val in 1:length(file.names)){
x<-read.table(file.names[val], sep='\t', header=FALSE)
x[is.na(x)]<-1
w<-cbind(x, (p.adjust(x$V4, "BH")), (p.adjust(x$V4, "bonferroni")), -log10(p.adjust(x$V4, "BH")), -log10(p.adjust(x$V4, "bonferroni")))
x<-w
colnames(x)<-c("Sample", "ObservedSNPs", 'ExpectedSNPs', 'PValue', 'Log2SNPFoldEnrichment', 'Lineage', 'BHCorrectedPValue', 'BonferroniCorrectedPValue', 'BHCorrectedPValueLog', 'BonferroniCorrectedPValueLog')
write.table(x, str_remove(paste(file.names[val], '-table.txt', sep=''), "-out.txt"), sep='\t')
write_xlsx(x, str_remove(paste(file.names[val], '-table.xlsx', sep=''), "-out.txt"))
n<-head(x[order(-x$Log2SNPFoldEnrichment),], 15)
q<-head(x[order(-x$BHCorrectedPValueLog),], 15)
pdf(file=str_remove(paste(file.names[val], '.pdf', sep=''), "-out.txt"), height=8.5, width=11)
print(ggplot(x, aes(x=Log2SNPFoldEnrichment, y=BHCorrectedPValueLog, colour=Lineage, shape=Lineage, label=Sample)) + geom_point(aes(size=ObservedSNPs))+ geom_hline(yintercept=1.3) + ggtitle(str_replace_all(str_replace(str_remove(paste(file.names[val], ' SNPs in Mouse SuperEnhancers', sep=''), "-out.txt"), "-mouse", " "),"_", " ")) + geom_label_repel(aes(label=Sample), data=subset(x, Log2SNPFoldEnrichment >= n$Log2SNPFoldEnrichment[10] | BHCorrectedPValueLog >= q$BHCorrectedPValueLog[10])) + scale_color_manual(values=c("orange", "red", "grey", "grey")) + labs(x="Fold Enrichment (Log2)", y="BH Corrected P Value (-Log10)", colour="Lineage", size="Number of SNPs in SuperEnhancers", shape="Lineage"))
dev.off()
}





#this last section is to export all the SNPs in each bed for each sample, which is for downstream specific investigations of specific SNPs
for sample in *human/ *mouse/
do
cd $sample
echo $sample
export NAME=`echo $sample | sed 's/-human\//-SNPsperbed-human.txt/g' | sed 's/-mouse\//-SNPsperbed-mouse.txt/g'`
for bed in *.bed
do
cd $bed
echo $bed | sed 's/hg19-25pct-//g' >> $DIR/$NAME
ls index*.txt | xargs cat | egrep -v "index" | sed 's/:/\t/g' | sed 's/^/chr/' | sed 's/ /\t/g' | awk '{print $0, FILENAME}' >> $DIR/$NAME
cd ../
done
cd $DIR
done

