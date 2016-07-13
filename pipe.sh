#!/bin/bash

GRCh37="/home/idedios/reference/GRCh37.fa"
picard="/home/idedios/dev/picard-tools-2.4.1/picard.jar"
cores=32
picardMem="Xmx2g"
samMem="682M"
exitstatus=0
date="2016-07-11"
align_prep="$PWD/work/align_prep/tx/tmp"
align="$PWD/work/align"

for sample in $PWD/*_R1_001.fastq.gz; do
    sampleidpath="${sample%_R1_001.fastq.gz}"
    sampleid=${sampleidpath##*/}
    sampleidR1="${sampleid}_R1_001"
    samplename=${sampleid%%_*}
    echo $samplename
    fastqR1gzpath=$sample
    fastqR2gzpath="$PWD/${sampleid}_R2_001.fastq.gz"
    fastqR1gz="${sampleid}_R1_001.fastq.gz"
    fastqR2gz="${sampleid}_R2_001.fastq.gz"

    # Convert to bgzip format and rezip
    #mkdir -p $PWD/work/align_prep/tx/tmp/
    #unpigz -c $fastqR1gzpath | bgzip -@ $cores -c /dev/stdin > $align_prep/$fastqR1gz
    #unpigz -c $fastqR2gzpath | bgzip -@ $cores -c /dev/stdin > $align_prep/$fastqR2gz

    # Index the fastqs with grabix
    #grabix index $align_prep/$fastqR1gz
    #grabix index $align_prep/$fastqR2gz

    # Align with bwa mem, add read groups, mark duplicates with samblaster, sort with samtools
    bwa mem -c 250 -M -t $cores -R "@RG\tID:${sampleid}\tPL:illumina\tPU:1_${date}_${samplename}\tSM:${sampleidR1}" -v 1 $GRCh37 $align_prep/$fastqR1gz $align_prep/$fastqR2gz > $align_prep/${sampleid}.sam #| samblaster --addMateTags -M --splitterFile >(samtools sort -@ $cores -m $samMem -T $align/$sampleidR1/tx/tmp/${sampleidR1}-sort-sorttmp-sp1 -o $align/$sampleidR1/tx/tmp/${sampleidR1}-sort-sr.bam /dev/stdin) --discordantFile >(samtools sort -@ $cores -m $samMem -T $align/$sampleidR1/tx/tmp/${sampleidR1}-sort-sorttmp-disc -o $align/$sampleidR1/tx/tmp/${sampleidR1}-sort-disc.bam /dev/stdin) | samtools view -b -S -u - | sambamba sort -t $cores -m $samMem --tmpdir $align/$sampleidR1/tx/tmp/${sampleidR1}-sort-sorttmp-full -o $align/$sampleidR1/tx/tmp/${sampleidR1}-sort.bam /dev/stdin

done

exit $exitstatus
