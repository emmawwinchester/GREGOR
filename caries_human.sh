source ~/.bash_profile
perl -MCPAN -e install DBD::SQLite
perl -MCPAN -e install Switch.pm
export DIR=/home/CAM/ewentworth/cotney/analysis/gwas/carieshuman
export DIR1=/home/CAM/ewentworth/cotney/rawdata/chromatinsegmentations/superenhancers/human
export DIR3=/home/CAM/ewentworth/cotney/analysis/gwas/carieshuman/scripts
export DIR2=/home/CAM/ewentworth/cotney/rawdata/chromatinsegmentations/human_18state
cd $DIR
ls $DIR2/*.bed $DIR1/*.bed> $DIR/index_human_caries.txt
egrep -hi "state9|state10|state8|uper" $DIR/index_human_caries.txt > $DIR/index_human_activeenhancers_caries.txt
egrep -hi  "state9|state10|state8|uper|state7|state15|state11" $DIR/index_human_caries.txt > $DIR/index_human_allenhancers_caries.txt


#Here is the script only to start running GREGOR on active enhancers and all states separated, as we previously noted that all enhancers was not an interesting find. we kept the index file in case we want to run it later. 


#these sections generate the gregor scripts and files for analyses, sbatch command starts running your scripts. 
for sample in *SNPs.txt
do
export NAME=`echo $sample | sed 's/SNPs.txt//g'`
echo $sample
mkdir $DIR/$NAME
echo -e "INDEX_SNP_FILE = $sample\nBED_FILE_INDEX = index_human_caries.txt\nREF_DIR = /home/FCAM/jcotney/TOOLS/GREGOR/ref\nR2THRESHOLD = 0.8\nLDWINDOWSIZE = 1000000\nOUT_DIR = $DIR/$NAME\nMIN_NEIGHBOR_NUM = 500\nBEDFILE_IS_SORTED = False\nPOPULATION = EUR\nTOPNBEDFILES = 2\nJOBNUMBER = 10\n###############################################################################\nBATCHTYPE = slurm\nBATCHOPTS = --partition=general --qos=general -c 8 --mem=200G " >$DIR/$NAME-conf.file
echo -e "#/bin/bash\n#SBATCH --job-name=$NAME-human\n#SBATCH -N 1\n#SBATCH -n 1\n#SBATCH -c 8\n#SBATCH --partition=general\n#SBATCH --qos=general\n#SBATCH --mail-type=END\n#SBATCH --mem=200G\n#SBATCH --mail-user=wentworth@uchc.edu\n#SBATCH -o %j.out\n#SBATCH -e %j.err\nperl -MCPAN -e install DBD::SQLite\nperl -MCPAN -e install Switch.pm\ncd $DIR\nperl /home/CAM/ewentworth/cotney/tools/gregor/GREGOR/script/GREGOR.pl --conf $DIR/$NAME-conf.file" | sed 's/\/bin/!\/bin/g' > $DIR3/gregor-$NAME.sh
sbatch $DIR3/gregor-$NAME.sh
done
for sample in *SNPs.txt
do
export NAME=`echo $sample | sed 's/SNPs.txt/-activeenhancers/g'`
echo $sample
mkdir $DIR/$NAME
echo -e "INDEX_SNP_FILE = $sample\nBED_FILE_INDEX = index_human_activeenhancers_caries.txt\nREF_DIR = /home/FCAM/jcotney/TOOLS/GREGOR/ref\nR2THRESHOLD = 0.8\nLDWINDOWSIZE = 1000000\nOUT_DIR = $DIR/$NAME\nMIN_NEIGHBOR_NUM = 500\nBEDFILE_IS_SORTED = False\nPOPULATION = EUR\nTOPNBEDFILES = 2\nJOBNUMBER = 10\n###############################################################################\nBATCHTYPE = slurm\nBATCHOPTS = --partition=general --qos=general -c 8 --mem=200G " >$DIR/$NAME-conf.file
echo -e "#/bin/bash\n#SBATCH --job-name=$NAME-human\n#SBATCH -N 1\n#SBATCH -n 1\n#SBATCH -c 8\n#SBATCH --partition=general\n#SBATCH --qos=general\n#SBATCH --mail-type=END\n#SBATCH --mem=200G\n#SBATCH --mail-user=wentworth@uchc.edu\n#SBATCH -o %j.out\n#SBATCH -e %j.err\nperl -MCPAN -e install DBD::SQLite\nperl -MCPAN -e install Switch.pm\ncd $DIR\nperl /home/CAM/ewentworth/cotney/tools/gregor/GREGOR/script/GREGOR.pl --conf $DIR/$NAME-conf.file" | sed 's/\/bin/!\/bin/g' > $DIR3/gregor-$NAME.sh
sbatch $DIR3/gregor-$NAME.sh
done

#here we begin to generate output files for downstream analysis. 
export DIR=/home/CAM/ewentworth/cotney/analysis/gwas/carieshuman
for sample in caries*/ crohns*
do
cd $sample
echo $sample
export NAME=`echo $sample | sed 's/\//-out.txt/g'`
awk 'NR==1{print;next}{if ($3=="NA") {print $0, 0} else {print $0, (log($2/$3)/log(2))}}' Statistic*.txt | tail -n +2 | sed 's/state_2/state2/g' | sed 's/_rconsensus//g' | sed 's/.bed//g' | sed 's/hg19-25pct-//g' | awk '{if ($1 ~ /CS|facial|cranio|reproducible|ooth|snout|palate/) {print $0, "craniofacial"} else if ($1 ~ /CD|Liver|Thymus/) {print $0, "immune"} else {print $0, "noncraniofacial"}}' | sed 's/ /\t/g' | sed 's/_Superenhancer_sorted//g' | sed 's/-superenhancer//g' > $DIR/$NAME
echo $DIR/$NAME
cd $DIR
done
#makefigures.Rscript is uploaded as a separate file; it generates .pdf files and .xlsx and .txt files of the tables generated in R. 
sbatch $DIR/makefigures.Rscript

#here we create files with all the information about which enhancers were in LD with a SNP or contained a SNP; this file makes it possible to go ahead and delete the directories that take up a ton of room. 
for sample in caries*/ crohns*/
do
cd $sample
echo $sample
export NAME=`echo $sample | sed 's/\//-SNPsperbed-human.txt/g'/g'`
for bed in *.bed
do
cd $bed
echo $bed | sed 's/hg19-25pct-//g' >> $DIR/$NAME
ls index*.txt | xargs cat | egrep -v "index" | sed 's/:/\t/g' | sed 's/^/chr/' | sed 's/ /\t/g' | awk '{print $0, FILENAME}' >> $DIR/$NAME
cd ../
done
cd $DIR
done
for sample in *.pdf
do
mv $sample figures/$sample
done
for sample in *table.xlsx *table.txt
do
mv $sample tables/$sample
done
for sample in *SNPsperbed*
do
mv $sample SNPs_per_bed/$sample
done
rm -r caries*/ crohns*/





