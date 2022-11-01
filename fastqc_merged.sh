#!/bin

# Run FastQC on all files, output into fastqc_merged folder

find ./data_merged -maxdepth 1 -name "*.fastq.gz"  -print|sort|uniq > filelist.txt

for F in $(cat filelist.txt) ; do
FULLSTRING=$F
fastqc ${FULLSTRING} -t 20 -o ./fastqc_merged
done

