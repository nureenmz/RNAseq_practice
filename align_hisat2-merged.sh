#!/bin
find ./data_merged/* -maxdepth 1 -name "GSM*" -print|cut -d "/" -f3| sort|uniq|cut -d "." -f1|uniq|sort >filelist.txt

for F in $(cat filelist.txt) ; do

        FULLSTRING=$F
        hisat2 -p 20 -x ./omes/chr/hisat_index -U ./data_merged/${F}.fastq.gz | samtools view -@ 20 -bS | samtools sort -@ 20 -O BAM -o ./aligned_merged/${F}.sorted.bam
done

