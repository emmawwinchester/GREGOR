source ~/.bash_profile
perl -MCPAN -e install DBD::SQLite
perl -MCPAN -e install Switch.pm
export DIR=/home/CAM/ewentworth/cotney/analysis/gwas/superenhancers
export DIR1=/home/CAM/ewentworth/cotney/rawdata/chromatinsegmentations/superenhancers/human
export DIR2=/home/CAM/ewentworth/cotney/rawdata/chromatinsegmentations/superenhancers/mouse
export DIR3=/home/CAM/ewentworth/cotney/analysis/gwas/superenhancers/scripts
ls $DIR2/hg19*.bed > $DIR/index_mouse.txt
ls $DIR1/*.bed > $DIR/index_human.txt
cd $DIR
##these first steps are to generate the files required to run GREGOR. each *SNPs.txt file is the list of SNPs for each trait (one file per trait), in rsID format. See GREGOR webpage for details of setup. Note this requires a ton of memory and will take a long time to run, especially when SNP number > 300.
##the bed_file_index.txt contains the list of the path to all the bed files you want to analyze for enrichment. NOTE THAT YOU ARE NOT ONLY INCLUDING BEDS FOR THE EXPECTED TISSUE OF INTEREST. To get anything of statistical significance, include a lot of tissues. 
for sample in *SNPs.txt
do
export NAME=`echo $sample | sed 's/SNPs.txt//g'`
echo $sample
mkdir $DIR/$NAME-mouse
echo -e "INDEX_SNP_FILE = $sample\nBED_FILE_INDEX = index_mouse.txt\nREF_DIR = /home/FCAM/jcotney/TOOLS/GREGOR/ref\nR2THRESHOLD = 0.8\nLDWINDOWSIZE = 1000000\nOUT_DIR = $DIR/$NAME-mouse\nMIN_NEIGHBOR_NUM = 500\nBEDFILE_IS_SORTED = True\nPOPULATION = EUR\nTOPNBEDFILES = 2\nJOBNUMBER = 10\n###############################################################################\nBATCHTYPE = slurm\nBATCHOPTS = --partition=general --qos=general -n 2 -c 32 --mem=250G " > $DIR/$NAME-mouse-conf.file
echo -e "#/bin/bash\n#SBATCH --job-name=mouse\n#SBATCH -N 1\n#SBATCH -n 2\n#SBATCH -c 32\n#SBATCH --partition=general\n#SBATCH --qos=general\n#SBATCH --mail-type=END\n#SBATCH --mem=250G\n#SBATCH --mail-user=wentworth@uchc.edu\n#SBATCH -o %j.out\n#SBATCH -e %j.err\ncd $DIR/\nperl /home/CAM/ewentworth/cotney/tools/gregor/GREGOR/script/GREGOR.pl --conf $DIR/$NAME-mouse-conf.file" | sed 's/\/bin/!\/bin/g' > $DIR3/gregor-$NAME.sh
sbatch $DIR3/gregor-$NAME.sh
done
for sample in *SNPs.txt
do
export NAME=`echo $sample | sed 's/SNPs.txt//g'`
echo $sample
mkdir $DIR/$NAME-human
echo -e "INDEX_SNP_FILE = $sample\nBED_FILE_INDEX = index_human.txt\nREF_DIR = /home/FCAM/jcotney/TOOLS/GREGOR/ref\nR2THRESHOLD = 0.8\nLDWINDOWSIZE = 1000000\nOUT_DIR = $DIR/$NAME-human\nMIN_NEIGHBOR_NUM = 500\nBEDFILE_IS_SORTED = False\nPOPULATION = EUR\nTOPNBEDFILES = 2\nJOBNUMBER = 10\n###############################################################################\nBATCHTYPE = slurm\nBATCHOPTS = --partition=general --qos=general -c 8 --mem=200G " >$DIR/$NAME-human-conf.file
echo -e "#/bin/bash\n#SBATCH --job-name=human\n#SBATCH -N 1\n#SBATCH -n 1\n#SBATCH -c 8\n#SBATCH --partition=general\n#SBATCH --qos=general\n#SBATCH --mail-type=END\n#SBATCH --mem=200G\n#SBATCH --mail-user=wentworth@uchc.edu\n#SBATCH -o %j.out\n#SBATCH -e %j.err\nperl -MCPAN -e install DBD::SQLite\nperl -MCPAN -e install Switch.pm\ncd $DIR\nperl /home/CAM/ewentworth/cotney/tools/gregor/GREGOR/script/GREGOR.pl --conf $DIR/$NAME-human-conf.file" | sed 's/\/bin/!\/bin/g' > $DIR3/gregor-human-$NAME.sh
sbatch $DIR3/gregor-human-$NAME.sh
done
#this was run as a separate script after the others had finished. it generates output files containing all information about enrichment of SNPs within each input bed. 
for sample in *human/  
do
cd $sample
echo $sample
export NAME=`echo $sample | sed 's/\//-out.txt/g'`
awk 'NR==1{print;next}{if ($3=="NA") {print $0, 0} else {print $0, (log($2/$3)/log(2))}}' Statistic*.txt | tail -n +2 | sed 's/state_2/state2/g' | sed 's/_rconsensus//g' | sed 's/.bed//g' | sed 's/hg19-25pct-//g' | awk '{if ($1 ~ /CS|facial|cranio|reproducible|ooth|snout|palate/) {print $0, "craniofacial"} else if ($1 ~ /CD|Liver|Thymus/) {print $0, "immune"} else if ($1 ~ /heart|Ventricle|Atrium/) {print $0, "heart"} else {print $0, "noncraniofacial"}}' | sed 's/ /\t/g' | sed 's/_Superenhancer_sorted//g' | sed 's/-superenhancer//g' > $DIR/$NAME
echo $DIR/$NAME
cd $DIR
done
#output: col1: bed file col2:#observed Snps col3: expected Snps col4: pvalue (unadjusted) col5: log2(observed/expected) col6: tissue 
#change the NAs to 1s in R
#do p value correction in R! add a new column with corrected bonferroni p value and another with BH corrected p value. 
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
#another important thing to note here is that it contains all the relevant information for each analysis. AFTER THIS YOU CAN DELETE THE DIRECTORIES CONTAINING YOUR GREGOR ANALYSIS. These directories take up a lot of space, so if you're limited then make sure not to forget the rm step. 
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

