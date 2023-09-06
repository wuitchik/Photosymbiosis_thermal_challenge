## Pipeline to call snips from A. poculata Tagseq data for the purpose of assessing clonality of samples

# Working directory

/projectnb/coral/MPCC_2018/oct_2020/SNP_analysis

# Modules

module load python3
module load bowtie2
module load htslib/1.9
module load samtools/1.9
module load angsd
module load picard

 
#converting .sam to .bam

>sortConvert
for file in *.sam
do echo "samtools sort -O bam -o ${file/.sam/}.bam $file" >> sortConvert
done 

/projectnb/coral/MPCC_2018/scc6_qsub_launcher.py -N bamming -P coral -M wuitchik@bu.edu -j y -h_rt 24:00:00 -jobsfile sortConvert

# picard index files

picard CreateSequenceDictionary R=../reference_files/coral_refgenome.fasta O=coral.dict

# samtool index files 

samtools faidx coral_refgenome.fasta 

# Picard dictionary 
picard CreateSequenceDictionary R=coral_refgenome.fasta O=coral.dict

# running angsd

FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -baq 1 -ref coral_refgenome.fasta -maxDepth 470" 
TODO="-doQsDist 1 -doDepth 1 -doCounts 1 -dumpCounts 2"

for file in *.bam
do echo "$file" >> bam_list
done 

samtools view AI2.bam|head

angsd -b bam_list -GL 1 $FILTERS $TODO -P 12 -out dd 

#Make new file called plotQC.R and copy script from Misha's website (https://github.com/z0on/2bRAD_denovo/blob/master/plotQC.R) into it

module load R
Rscript plotQC.R dd

AG1.bam     0.000000000
AK2.bam     0.000000000
AS1.bam     0.000000000
AL3.bam     0.009059787
AF2.bam     0.010305257
AD5.bam     0.053056261
AN3.bam     0.057123893
AI2.bam     0.057136653
AB3.bam     0.058161757
AH2.bam     0.059688424
AS3.bam     0.060481349
AL2.bam     0.067322085
AM1.bam     0.073883858
AE5.bam     0.074924838
AP2.bam     0.076458210
AJ3.bam     0.077264128
AF3.bam     0.078014900
AA4.bam     0.078148060
AB2.bam     0.084514864
AS5.bam     0.085747582
AH1.bam     0.086272251
AJ4.bam     0.087965904
AI5.bam     0.094025359
AL1.bam     0.095579943
AM3.bam     0.095892520
AH3.bam     0.096641378
AJ2.bam     0.096922016
AK3.bam     0.097237659
AE6.bam     0.097464228
AP4.bam     0.099112500
AI1.bam     0.102585199
AD4.bam     0.104234417
AE2.bam     0.105428060
AM2.bam     0.105811009
AA2.bam     0.107799754
AD6.bam     0.108655758
AK5.bam     0.109939100
AB1.bam     0.110433592
AG2.bam     0.112276291
AN1.bam     0.112664564
AA3.bam     0.112970147
AN2.bam     0.116806330
AP1.bam     0.117880535
AF5.bam     0.119398926
AC1.bam     0.130789303
AC3.bam     0.131211147
AG3.bam     0.132366454

# So we remove d the following Bams and renamed to bam_list2
AG1.bam     0.000000000
AK2.bam     0.000000000
AS1.bam     0.000000000

# Same filters as Oculina
FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 25 -dosnpstat 1 -doHWE 1 -sb_pval 1e-5 -hetbias_pval 1e-5 -skipTriallelic 1 -minInd 36 -snp_pval 1e-5 -minMaf 0.05"
TODO="-doMajorMinor 1 -ref coral_refgenome.fasta -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 8 -doBcf 1 -doPost 1 -doGlf 2"

angsd -b bam_list2 -GL 1 $FILTERS $TODO -P 12 -out coral_ang

# ADMIXTURE

module load admixture
NSITES=`zcat coral_ang.beagle.gz | wc -l`
echo $NSITES

## 17017 SNPs

for K in 2 3 4; 
do 
NGSadmix -likes coral_ang.beagle.gz -K $K -P 10 -o admix;
done
