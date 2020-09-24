#here we lumped together all active enhancer states into one file, instead of running them separately. 
#superenhancers were NOT included in this analysis. 
#samples were 25 state chromatin segmentation files from human embryonic craniofacial and all roadmap samples. 
#human embryonic heart samples were not included in this analysis. 
#included bed files are seen in the facial_variation_gwas_beds.txt file in this repository. 
cd /home/CAM/ewentworth/cotney/analysis/gwas/facialprofile2/
#this particular segment only took SNPs which passed the significance level of 5e-8 set by the GWAS field. 
for sample in agg*.txt
do
awk '{if ($3<=5E-8) {print $0}}' $sample | grep -v "E-07" | cut -f1,2 | tail -n +1 | sed 's|^|chr|' | sed 's/\t/:/g' > temp.txt
mv temp.txt $sample
done
for sample in hg19*-gregor.txt 
do
echo $sample
mkdir /home/CAM/ewentworth/cotney/analysis/gwas/facialprofile2/"${sample%%.*}"_out
echo -e "INDEX_SNP_FILE = $sample\nBED_FILE_INDEX = facial_variation_gwas_beds.txt\nREF_DIR = /home/FCAM/jcotney/TOOLS/GREGOR/ref\nR2THRESHOLD = 0.8\nLDWINDOWSIZE = 1000000\nOUT_DIR = /home/CAM/ewentworth/cotney/analysis/gwas/facialprofile2/"${sample%%.*}"_out\nMIN_NEIGHBOR_NUM = 500\nBEDFILE_IS_SORTED = False\nPOPULATION = EUR\nTOPNBEDFILES = 2\nJOBNUMBER = 10\n###############################################################################\nBATCHTYPE = slurm\nBATCHOPTS = --partition=general --qos=general --mem=16G " > ~/cotney/analysis/gwas/facialprofile2/"${sample%%.*}"-conf.file
echo -e "#/bin/bash\n#SBATCH --job-name=gregor\n#SBATCH -N 1\n#SBATCH -n 1\n#SBATCH -c 4\n#SBATCH --partition=general\n#SBATCH --qos=general\n#SBATCH --mail-type=END\n#SBATCH --mem=128G\n#SBATCH --mail-user=wentworth@uchc.edu\n#SBATCH -o %j.out\n#SBATCH -e %j.err\nperl -MCPAN -e install DBD::SQLite\nperl -MCPAN -e install Switch.pm\ncd /home/CAM/ewentworth/cotney/analysis/gwas/facialprofile2\nperl /home/CAM/ewentworth/cotney/tools/gregor/GREGOR/script/GREGOR.pl --conf "${sample%%.*}"-conf.file" | sed 's/\/bin/!\/bin/g' > gregor-"${sample%%.*}".sh
sbatch gregor-"${sample%%.*}".sh
done

#after these were done, ran this
for sample in *_out/
do
cd $sample
mv StatisticSummaryFile.txt  ~/cotney/analysis/gwas/facialprofile2/"${sample%%allenhancers_out*}"StatsisticSummaryFile.txt
cd ../
done
#yes I am aware I spelled statistic wrong
for sample in *Statsistic*
do
export NAME=`echo $sample | sed 's/-gregor--StatsisticSummaryFile.txt/-allenhancers/g' | sed 's/hg19-//g'`
export DIR=/home/CAM/ewentworth/cotney/analysis/gwas/facialprofile2
awk 'NR==1{print;next}{if ($3=="NA") {print $0, 0} else {print $0, (log($2/$3)/log(2))}}' $sample | tail -n +2 | sed 's/_25state_dense.enhancer_states.bed//g' | sed 's/_25_imputed12marks_dense.enhancer_states.bed//g' | sed 's/reproducible_enhancer_states.25.bed/All_CS_Enhancers/g' | sed 's/reproducible_craniofacial_specific_enhancer_states.bed/CS_Specific_Enhancers/g' | sed 's/hg19-//g' | awk '{if ($1 ~ /CS|F2/) {print $0, "craniofacial"} else {print $0, "noncraniofacial"}}' | awk '{if ($1 ~ /E080|E081|E082|E083|E084|E085|E086|E087|E088|E089/) {print $0, "fetal"} else if ($1 ~ /CS|F2|E001|E002|E003|E004|E005|E006|E007|E008|E009|E010|E011|E012|E013|E014|E015|E016|E024|E018|E019|E020|E021|E022/) {print $0, "embryonic"} else {print $0, "adult"}}' | sed 's/ /\t/g' > $DIR/scatterplot-craniovsnon-$NAME-out.txt
done
module load R/3.6.1
R
library(ggplot2)
library(ggrepel)
library(stringr)
library(writexl)
library(patchwork)
path = "~/cotney/analysis/gwas/facialprofile2/"
out.file<-""
file.names<-dir(path, pattern=(glob2rx("scatterplot*")))
for (val in 1:length(file.names)){
x<-read.table(file.names[val], sep='\t', header=FALSE)
x[is.na(x)]<-1
w<-cbind(x, (p.adjust(x$V4, "BH")), (p.adjust(x$V4, "bonferroni")), -log10(p.adjust(x$V4, "BH")), -log10(p.adjust(x$V4, "bonferroni")))
x<-w
colnames(x)<-c("Sample", "ObservedSNPs", 'ExpectedSNPs', 'PValue', 'Log2SNPFoldEnrichment', 'Lineage', 'Timepoint', 'BHCorrectedPValue', 'BonferroniCorrectedPValue', 'BHCorrectedPValueLog', 'BonferroniCorrectedPValueLog')
write.table(x, str_remove(str_remove(paste(file.names[val], '-table.txt', sep=''), "scatterplot-craniovsnon-"), "-allenhancers-out.txt"), sep='\t')
write_xlsx(x, str_remove(str_remove(paste(file.names[val], '-table.xlsx', sep=''), "scatterplot-craniovsnon-"), "-allenhancers-out.txt"))
n<-head(x[order(-x$Log2SNPFoldEnrichment),], 15)
q<-head(x[order(-x$BonferroniCorrectedPValueLog),], 15)
pdf(file=str_remove(str_remove(paste(file.names[val], 'allenh.pdf', sep=''), "scatterplot-craniovsnon-"), "-allenhancers-out.txt"), height=8.5, width=11)
print(ggplot(x, aes(x=Log2SNPFoldEnrichment, y=BonferroniCorrectedPValueLog, colour=Lineage, shape=Lineage, label=Sample)) + geom_point(aes(size=ObservedSNPs))+ geom_hline(yintercept=1.3) + ggtitle(str_replace(str_remove(str_remove(paste(file.names[val], ' SNPs in All Enhancers', sep=''), "scatterplot-craniovsnon-"), "-allenhancers-out.txt"),"_", " ")) + geom_label_repel(aes(Log2SNPFoldEnrichment, BonferroniCorrectedPValueLog, label=Sample), data=subset(x, Lineage == "craniofacial" ), nudge_x=-2, direction="y")+ scale_color_manual(values=c("orange", "grey")) + xlim(-1,3.5) + labs(x="Fold Enrichment (Log2)", y="Bonferroni Corrected P Value (-Log10)", colour="Lineage", size="Number of SNPs in Enhancers"))
print(ggplot(x, aes(x=Log2SNPFoldEnrichment, y=BonferroniCorrectedPValueLog, colour=Lineage, shape=Lineage, label=Sample)) + geom_point(aes(size=ObservedSNPs))+ geom_hline(yintercept=1.3) + ggtitle(str_replace(str_remove(str_remove(paste(file.names[val], ' SNPs in All Enhancers', sep=''), "scatterplot-craniovsnon-"), "-allenhancers-out.txt"),"_", " ")) + geom_label_repel(aes(Log2SNPFoldEnrichment, BonferroniCorrectedPValueLog, label=Sample), data=subset(x, BonferroniCorrectedPValueLog>=q$BonferroniCorrectedPValueLog[15] | Log2SNPFoldEnrichment>=n$Log2SNPFoldEnrichment[15]), nudge_x=-2, direction="y")+ scale_color_manual(values=c("orange", "grey")) + xlim(-1,3.5) + labs(x="Fold Enrichment (Log2)", y="Bonferroni Corrected P Value (-Log10)", colour="Lineage", size="Number of SNPs in Enhancers"))
dev.off()
}

for sample in *allenhancers_out
do
cd $sample/reproducible_craniofacial_specific_enhancer_states.bed
pwd
ls index*.txt | xargs cat | egrep -v "index" > ~/cotney/analysis/gwas/facialprofile2/"${sample%%allenhancers_out*}"SNPSincfspecificenh.txt
cd ~/cotney/analysis/gwas/facialprofile2/$sample/reproducible_enhancer_states.25.bed
ls index*.txt | xargs cat | egrep -v "index" > ~/cotney/analysis/gwas/facialprofile2/"${sample%%allenhancers_out*}"SNPSincfreproducibleenh.txt
cd ~/cotney/analysis/gwas/facialprofile2/
done
