source ~/.bash_profile
perl -MCPAN -e install DBD::SQLite
perl -MCPAN -e install Switch.pm
module load bedtools
module load kent-tools
export DIR=/home/CAM/ewentworth/cotney/rawdata/chromatinsegmentations/mouse_18state
cd $DIR
#this is where the mouse encode samples are all lifted over to the hg19 genome. 
rm *unknown* *brain* 
for sample in *.bed
do
liftOver -minMatch=0.25 -bedPlus=9 $sample ~/cotney/genome/mm10/mm10ToHg19.over.chain.gz temp.bed temp.unmapped_0.10.bed
sort -k1,1 -k2,2n temp.bed > hg19-25pct-$sample
done
cd /home/CAM/ewentworth/cotney/rawdata/chromatinsegmentations/superenhancers/mouse
for sample in *.bed
do
liftOver -minMatch=0.25 -bedPlus=9 $sample ~/cotney/genome/mm10/mm10ToHg19.over.chain.gz temp.bed temp.unmapped_0.10.bed
sort -k1,1 -k2,2n temp.bed > hg19-25pct-$sample
done
export DIR=/home/CAM/ewentworth/cotney/rawdata/chromatinsegmentations/mouse_18state
export DIR2=/home/CAM/ewentworth/cotney/analysis/gwas/cariesmouse
export DIR3=/home/CAM/ewentworth/cotney/rawdata/chromatinsegmentations/superenhancers/mouse
#active enhancers contained states 8, 9, 10 and superenhancers from the 18 state model. 
ls $DIR/hg19*.bed $DIR3/hg19*| egrep -v "consensus" > $DIR2/index_mouse_caries.txt
ls $DIR/hg19*.bed $DIR3/hg19*| egrep -hi "state8|state9|state10|uper" > $DIR2/index_mouse_caries_activeenhancers.txt
ls $DIR/hg19*.bed $DIR3/hg19*| egrep -hi "state7|state8|state9|state10|state11|state15|uper" > $DIR2/index_mouse_caries_allenhancers.txt
#included in index snp files is only those rsIDs from GWAS catalog for "caries" which have a p value of less than 5e-8 (or crohns)
cd $DIR2
for sample in index*.txt
do
export NAME=`echo $sample | sed 's/index_mouse_//g' | sed 's/.txt//g'`
echo $sample
mkdir $DIR2/$NAME-out
echo -e "INDEX_SNP_FILE = cariesSNPs.txt\nBED_FILE_INDEX = $sample\nREF_DIR = /home/FCAM/jcotney/TOOLS/GREGOR/ref\nR2THRESHOLD = 0.8\nLDWINDOWSIZE = 1000000\nOUT_DIR = /home/CAM/ewentworth/cotney/analysis/gwas/cariesmouse/$NAME-out\nMIN_NEIGHBOR_NUM = 500\nBEDFILE_IS_SORTED = False\nPOPULATION = EUR\nTOPNBEDFILES = 2\nJOBNUMBER = 10\n###############################################################################\nBATCHTYPE = slurm\nBATCHOPTS = --partition=general --qos=general -c 8 --mem=200G " > ~/cotney/analysis/gwas/cariesmouse/caries-$NAME-conf.file
echo -e "#/bin/bash\n#SBATCH --job-name=gregormouse\n#SBATCH -N 1\n#SBATCH -n 1\n#SBATCH -c 8\n#SBATCH --partition=general\n#SBATCH --qos=general\n#SBATCH --mail-type=END\n#SBATCH --mem=200G\n#SBATCH --mail-user=wentworth@uchc.edu\n#SBATCH -o caries-$NAME.out\n#SBATCH -e caries-$NAME.err\nperl -MCPAN -e install DBD::SQLite\nperl -MCPAN -e install Switch.pm\ncd /home/CAM/ewentworth/cotney/analysis/gwas/cariesmouse\nperl /home/CAM/ewentworth/cotney/tools/gregor/GREGOR/script/GREGOR.pl --conf ~/cotney/analysis/gwas/cariesmouse/caries-$NAME-conf.file" | sed 's/\/bin/!\/bin/g' > gregor-caries-$NAME.sh
sbatch gregor-caries-$NAME.sh
done
for sample in index*.txt
do
export NAME=`echo $sample | sed 's/index_mouse_caries/crohns/g' | sed 's/.txt/-out/g'`
echo $sample
mkdir /home/CAM/ewentworth/cotney/analysis/gwas/cariesmouse/$NAME
echo -e "INDEX_SNP_FILE = crohnsSNPs.txt\nBED_FILE_INDEX = $sample\nREF_DIR = /home/FCAM/jcotney/TOOLS/GREGOR/ref\nR2THRESHOLD = 0.8\nLDWINDOWSIZE = 1000000\nOUT_DIR = /home/CAM/ewentworth/cotney/analysis/gwas/cariesmouse/$NAME\nMIN_NEIGHBOR_NUM = 500\nBEDFILE_IS_SORTED = False\nPOPULATION = EUR\nTOPNBEDFILES = 2\nJOBNUMBER = 10\n###############################################################################\nBATCHTYPE = slurm\nBATCHOPTS = --partition=general --qos=general -c 8 --mem=200G " > ~/cotney/analysis/gwas/cariesmouse/$NAME-conf.file
echo -e "#/bin/bash\n#SBATCH --job-name=gregormouse\n#SBATCH -N 1\n#SBATCH -n 1\n#SBATCH -c 8\n#SBATCH --partition=general\n#SBATCH --qos=general\n#SBATCH --mail-type=END\n#SBATCH --mem=200G\n#SBATCH --mail-user=wentworth@uchc.edu\n#SBATCH -o crohns-%j.out\n#SBATCH -e crohns-%j.err\nperl -MCPAN -e install DBD::SQLite\nperl -MCPAN -e install Switch.pm\ncd /home/CAM/ewentworth/cotney/analysis/gwas/cariesmouse\nperl /home/CAM/ewentworth/cotney/tools/gregor/GREGOR/script/GREGOR.pl --conf ~/cotney/analysis/gwas/cariesmouse/$NAME-conf.file" | sed 's/\/bin/!\/bin/g' > gregor-$NAME.sh
sbatch gregor-crohns-$NAME.sh
done
export DIR=/home/CAM/ewentworth/cotney/analysis/gwas/cariesmouse
cd $DIR
for sample in *-out/
do
cd $sample
echo $sample
export NAME=`echo $sample | sed 's/-out\//-out.txt/g'`
awk 'NR==1{print;next}{if ($3=="NA") {print $0, 0} else {print $0, (log($2/$3)/log(2))}}' Statistic*.txt | tail -n +2 | sed 's/state_2/state2/g' | sed 's/_rconsensus//g' | sed 's/.bed//g' | sed 's/hg19-25pct-//g' | awk '{if ($1 ~ /CS|facial|cranio|reproducible|ooth|snout|palate/) {print $0, "craniofacial"} else if ($1 ~ /CD|iver|hymus/) {print $0, "immune"} else {print $0, "noncraniofacial"}}' | sed 's/ /\t/g' | sed 's/_Superenhancer_sorted//g' | sed 's/-superenhancer//g' > $DIR/$NAME
echo $DIR/$NAME
cd $DIR
done
for x in *-out/
do
export NAME=`pwd | sed 's/-out//g'`
cd $x
for sample in *.bed
do
cd $sample
echo $sample | sed 's/hg19-25pct-//g' >> $NAME-SNPs.bed
ls index*.txt | xargs cat | egrep -v "index" | sed 's/:/\t/g' | sed 's/^/chr/' | sed 's/ /\t/g' | awk '{print $0, FILENAME}' >> $NAME-SNPs.bed
cd ../
done
cd ../
done
#from here follow the same makefigures.Rscript, swapping working directory and title to "mouse enhancers" rather than "human"
