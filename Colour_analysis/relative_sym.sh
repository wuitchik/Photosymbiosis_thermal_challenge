# Counting relative symbiont counts from .sam files


module load samtools

MINQ=40
>countSymb
for file in *.sam
do echo "samtools view -F 256 -q $MINQ $file | sym_vs_host.awk > ${file/.sam/}_relative_Counts.tsv" >>countSymb
done

../scc6_qsub_launcher.py -N rel_counting -P coral -M wuitchik@bu.edu -j y -h_rt 24:00:00 -jobsfile countSymb
qsub rel_counting_array.qsub	

awk '{print$0, FILENAME}' *.tsv > rel_counts.txt