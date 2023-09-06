## October 2020, mapping to A.poculata genome and B. psygmophilum transcriptome 
# working directory
/projectnb/coral/MPCC_2018/oct_2020

# local directory
/Users/danielwuitchik/Documents/Experiments/MPCC_2018/Astrangia/oct_2020

# raw backup
/projectnb/coral/raw_backup/Astrangia_GE_sequences

#First, I need to conver Hanny's transcriptome to be usable (ie, gff format). So, I turn the transcriptome into a dummy genome making fake reference genome (of 97 chromosomes) out ot cd-hit cluster representatives
# perl script can be found: https://github.com/z0on/2bRAD_denovo/blob/master/concatFasta.pl
../concatFasta.pl fasta=B_psygmophilum_transcriptome.fasta num=97 

# I used 97 chromosomes based on LaJeunesse, T.C., Lambert, G., Andersen, R.A., Coffroth, M.A. & Galbraith, D.W.
#Symbiodinium (Pyrrhophyta) genome sizes (DNA content) are smallest among
#dinoflagellates. J. Phycol. 41, 880-886 (2005).

# Rename files to make more readable
mv apoculata.assembly.scaffolds_chromosome_level.fasta coral_refgenome.fasta
mv B_psygmophilum_transcriptome_cc.fasta algae_reftranscriptome.fasta
mv B_psygmophilum_transcriptome_cc.tab algae.tab
mv apoculata_v2.0.gff3 coral.gff

# create 'gff' from scratch -  NB: 'BED' coords are actually not in proper coord format for a BED 
awk 'BEGIN {FS="\t"; OFS="\t"} {print $2}' algae.tab > algae.bed
sed -i -z 's/\n/\t'tech'\n/g' algae.bed 
sed -i -z 's/\n/\t'gene'\n/g' algae.bed 
awk 'BEGIN {FS=OFS="\t"} {print $3, $4}' algae.tab > algae.cols
paste algae.bed algae.cols > algae.v2.bed
sed -i -z 's/\n/\t\.\n/g' algae.v2.bed
sed -i -z 's/\n/\t\.\n/g' algae.v2.bed
sed -i -z 's/\n/\t\.\n/g' algae.v2.bed
awk 'BEGIN {FS=OFS="\t"} {print $1}' algae.tab > algae.ids
sed 's/^/'ID='/' algae.ids > algae.ids.v2
paste algae.v2.bed algae.ids.v2 > algae_ref.gff

# combine fastas and .gff3
cat algae_reftranscriptome.fasta coral_refgenome.fasta > holobiont.fasta
cat algae_ref.gff coral.gff > holobiont.gff3

# copy files into working directory
cp /projectnb/coral/raw_backup/Astrangia_GE_sequences/*.gz /projectnb/coral/MPCC_2018/oct_2020

# Load modules
module load python3
module load fastx-toolkit
module load bowtie2

# these files look like Astrangia_Pdam_AA2_S3_R1_001.fastq.gz, so we need to strip the characters into something meaningful, here I strip the first 15 characters
for i in *.gz
	do mv "$i" "`echo $i | sed 's/Astrangia_Pdam_//'`"
done

# Now I strip _R1_001

for i in *.gz
	do mv "$i" "`echo $i | sed 's/_R1_001//'`"
done

# Now I strip everything after the _

for i in *.gz
	do mv "$i" "`echo $i | sed 's/_.*//'`"
done

# that removed the extension, so we should fix that
for f in *;
	do mv "$f" "$f.fastq.gz";
done

# Let's unzip in parallel
for file in *.gz
	do echo "gunzip $file">> gunzip
done

# tidying
mv zippy* parallel_outputs

# concatenate the two genomes
cat apoculata.assembly.scaffolds_chromosome_level.fasta B_psygmophilum_transcriptome_cc.fasta  > holobiont.fasta
cat apoculata_v2.0.gff3 B_psg.gff3 > holobiont.gff3

# concatenate the two AC1 files
cat AC1.fastq AC1-2.fastq > AC1.fastq
rm -rf AC1-2.fastq


# trimming
for file in *.fastq
	do echo  "../tagseq_files/tagseq_clipper.pl $file | fastx_clipper -a AAAAAAAA -l 20 -Q33 | fastx_clipper -a AGATCGGAAG -l 20 -Q33 | fastq_quality_filter -Q33 -q 20 -p 90 > ${file/.fastq/}.trim" >> trimming
done 

../scc6_qsub_launcher.py -N trim -P coral -M wuitchik@bu.edu -j y -h_rt 24:00:00 -jobsfile trimming
qsub trim_array.qsub

mv trim* parallel_outputs


# creating bowtie2 index for your genome
bowtie2-build reference_files/holobiont.fasta reference_files/holobiont.fasta

# mapping 

for file in *.trim
	do echo "bowtie2 -x reference_files/holobiont.fasta -U $file --local -p 4 -S ${file/.trim/}.sam">> maps
done

../scc6_qsub_launcher.py -N mapping -P coral -M wuitchik@bu.edu -j y -h_rt 24:00:00 -jobsfile maps
qsub mapping_array.qsub

mv mapping* parallel_outputs

# sort the sam alignment files and convert to bams
module load samtools

for file in *.sam
	do echo "samtools sort -n -O bam -o ${file/.sam/}.bam $file" >> sortConvert
done 

../scc6_qsub_launcher.py -N sortconverting -P coral -M wuitchik@bu.edu -j y -h_rt 24:00:00 -jobsfile sortConvert
qsub sortconverting_array.qsub

mv sortconvert* parallel_outputs

#### Counting

module load subread

# Choose the GFF
MY_GFF="reference_files/holobiont.gff3"; GENE_ID="ID"
featureCounts -a $MY_GFF -p -t gene -g $GENE_ID -o feature_counts_out.txt -T 32 --primary *.bam


#get alignment counts
for file in *.bam
	do echo "samtools flagstat $file > ${file/.bam/}_flagStats.txt" >> getInitialAlignment
done

module load samtools
../scc6_qsub_launcher.py -N samtools_alignment_counts -P coral -M wuitchik@bu.edu -j y -h_rt 24:00:00 -jobsfile getInitialAlignment
qsub samtools_alignment_counts_array.qsub

mv samtools* parallel_outputs


# making raw read count table 
module purge
module load python2
module load htseq/0.11.0

MY_GFF="reference_files/holobiont.gff3"; GENE_ID="ID"

for file in *.sam
	do echo "htseq-count -t gene -i ID -m intersection-nonempty --stranded=no $file $MY_GFF > ${file/sam/counts.txt}" >> doCounts
done

../scc6_qsub_launcher.py -N counts -P coral -M wuitchik@bu.edu -j y -h_rt 24:00:00 -jobsfile doCounts
qsub counts_array.qsub 

# Compile to make big expression table

../tagseq_files/expression_compiler.pl *.counts.txt > holobiont_counts.txt

# Making a big o'l table compiling the number of counts in each step

>mapped_count.tsv
for file in *_flagStats.txt
do pp=$(grep "mapped" $file | head -n 1)
 echo -e "$file\t$pp" |\
 awk '{split($1, a, "_flagStats.txt")
 print a[1]"\t"$2"\tmapped"}' >> mapped_count.tsv
 done

# raw
wc -l *.fastq |\
	awk '{split($2, a, ".fastq")
	print a[1]"\t"$1/4"\rawCounts"}' |\
grep -v total > raw_read_counts.tsv 

# trimmed
wc -l *.trim |\
	awk '{split($2, a, ".trim")
	print a[1]"\t"$1/4"\ttrimmedCounts"}' |\
grep -v total > trimmed_read_counts.tsv 

# Compile pipeline counts for supps
cat *.tsv > final_pipeline_counts.txt