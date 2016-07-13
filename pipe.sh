#!/bin/bash

GRCh37="/home/idedios/reference/GRCh37.fa"
picard="/home/idedios/dev/picard-tools-2.4.1/picard.jar"
cores=32
mem="Xmx2g"
exitstatus=0

for sample in $PWD/*_R1_001.fastq.gz; do
    sampleid="${sample%_R1_001.fastq.gz}"
    #sampleid=${sampleidpath##*/}
    echo $sampleid
    fastqR1gz=$sample
    fastqR2gz="${sampleid}_R2_001.fastq.gz"
    fastqR1=${fastqR1gz%.gz}
    fastqR2=${fastqR2gz%.gz}

    # Convert to pbgzip format
    exec unpigz -c $fastqR1gz | pbgzip -n $cores -c /dev/stdin > $PWD/work/align_prep/tx/tmp/$fastqR1gz
    exec unpigz -c $fastqR2gz | pbgzip -n $cores -c /dev/stdin > $PWD/work/align_prep/tx/tmp/$fastqR2gz

    # Grabix the fastqs
    exec grabix index $PWD/work/align_prep/tx/tmp/$fastqR1gz
    exec grabix index $PWD/work/align_prep/tx/tmp/$fastqR2gz

    # bwa mem
    #exec bwa mem -c 250 -M -t $cores -R '@RG\tID:'$sampleid'\tPL:illumina\tPU:'
done

exit $exitstatus
