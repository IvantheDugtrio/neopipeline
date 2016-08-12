#!/bin/bash

allelefreq=0.05
cosmic="/home/idedios/reference/cosmic-v68-GRCh37.vcf.gz"
dbsnp="/home/idedios/reference/dbsnp_138.vcf.gz"
refgenome="/home/idedios/reference/GRCh37.fa"
gatkexec="/home/idedios/dev/gatk-3.6/GenomeAnalysisTK.jar"
picardexec="/home/idedios/dev/picard-tools-2.4.1/picard.jar"
cores=32
picardMem="Xmx2g"
samMem="682M"
exitstatus=0
date="2016-07-11"
align_prep="$PWD/work/align_prep/tx/tmp"
align="$PWD/work/align"
bedfile="/home/idedios/bcbio-manifests/F4.bed"
bedprep="$PWD/work/bedprep/tx/tmp"
bed="$PWD/work/bed/tx/tmp"
freebayes="$PWD/work/freebayes"
gatk="$PWD/work/gatk"
mutect="$PWD/work/mutect"
mutect2="$PWD/work/mutect2"
platypus="$PWD/work/platypus"
vardict="$PWD/work/vardict"

# Create align_prep tmp directory
if [ ! -f $PWD/work/align_prep/tx/tmp ]; then
    mkdir -p $PWD/work/align_prep/tx/tmp/
fi

for sample in $PWD/*1_001.fastq.gz; do
    sampleidpath="${sample%.fastq.gz}" # Remove extension
    sampleid=${sampleidpath##*/}    # Strip path from sampleidpath
    bedfilepath="${bedfile%.bed}"
    bedfilename="${bedfilepath##*/}"
    echo "sampleid is $sampleid"

    # Convert fastqs to bgzip format
    if [ ! -f $align_prep/${sampleid%1_001}2_001.fastq.gz ]; then
        unpigz -c $sample | bgzip -@ $cores -c /dev/stdin > $align_prep/${sampleid}.fastq.gz
        unpigz -c $PWD/${sampleid%1_001}2_001.fastq.gz | bgzip -@ $cores -c /dev/stdin > $align_prep/${sampleid%1_001}2_001.fastq.gz
    fi

    # Index the fastqs with grabix
    if [ ! -f $align/${sampleid%1_001}2_001.fastq.gz.gbi ]; then
        grabix index $align_prep/${sampleid}.fastq.gz
        grabix index $align_prep/${sampleid%1_001}2_001.fastq.gz
    fi

    sampleAlignTMP="$align/$sampleid/tx/tmp"
    if [ ! -f $sampleAlignTMP ]; then
        mkdir -p $sampleAlignTMP
    fi

    # Align with bwa mem, add read groups, mark duplicates with samblaster, sort with samtools
    if [ ! -f $align/${sampleid}/tx/tmp/${sampleid}-sort.bam ]; then
        bwa mem -c 250 -M -t $cores -R "@RG\tID:${sampleid%%_*}\tPL:illumina\tPU:1_${date}_${sample##*/}\tSM:${sampleid%%_*}" -v 1 $refgenome $align_prep/${sampleid}.fastq.gz $align_prep/${sampleid%1_001}2_001.fastq.gz | samblaster --addMateTags -M --splitterFile >(samtools sort -@ $cores -m $samMem -T $sampleAlignTMP/${sampleid}-sort-sorttmp-sp1 -o $sampleAlignTMP/${sampleid}-sort-sr.bam /dev/stdin) --discordantFile >(samtools sort -@ $cores -m $samMem -T $sampleAlignTMP/${sampleid}-sort-sorttmp-disc -o $sampleAlignTMP/${sampleid}-sort-disc.bam /dev/stdin) | samtools sort -@ $cores -m $samMem -T $sampleAlignTMP/${sampleid}-sort-sorttmp -o $sampleAlignTMP/${sampleid}-sort.bam #| samtools view -b -S -u - | sambamba sort -t $cores -m $samMem --tmpdir $sampleAlignTMP/${sampleid}-sort-sorttmp-full -o $sampleAlignTMP/${sampleid}-sort.bam /dev/stdin
    fi

    # Index bams produced
    if [ ! -f $sampleAlignTMP/${sampleid}-sort.bam.bai ]; then
        samtools index $sampleAlignTMP/${sampleid}-sort-sr.bam
        samtools index $sampleAlignTMP/${sampleid}-sort-disc.bam
        samtools index $sampleAlignTMP/${sampleid}-sort.bam
        #sambamba index -t $cores $sampleAlignTMP/${sampleid}-sort-sr.bam
        #sambamba index -t $cores $sampleAlignTMP/${sampleid}-sort-disc.bam
        #sambamba index -t $cores $sampleAlignTMP/${sampleid}-sort.bam
    fi

    # Cleanup and sort input bed file by chromosome
    if [ ! -f $bedprep/${bedfilename}.bed.gz.tbi ]; then
        `cat $bedfile | grep -v ^track | grep -v ^browser | grep -v ^# | /home/idedios/dev/ngsutils/bin/bedutils clean | sort -V -k1,1 -k2,2n > $bedprep/${bedfile##*/}`
        #`cat $bedfile | grep -v ^track | grep -v ^browser | grep -v ^# | py -x 'bcbio.variation.bedutils.remove_bad(x)' | sort -V -k1,1 -k2,2n > $bedprep/${bedfile##*/}`
        cat $bedprep/${bedfile##*/} | bgzip -@ $cores -c > $bedprep/${bedfilename}.bed.gz
        tabix -f -p bed $bedprep/${bedfilename}.bed.gz
    fi

    # Combine overlapping intervals
    if [ ! -f $bedprep/${bedfilename}-merged.bed.gz.tbi ]; then
        bedtools merge -i $bedprep/${bedfile##*/} > $bedprep/${bedfilename}-merged.bed
        cat $bedprep/${bedfilename}-merged.bed | bgzip -@ $cores -c > $bedprep/${bedfilename}-merged.bed.gz
        tabix -f -p bed $bedprep/${bedfilename}-merged.bed.gz
    fi

    sampleCallableRegionsTMP="$align/$sampleid/${sampleid}-sort-callable-split/tx/tmp"
    if [ ! -f $sampleCallableRegionsTMP ]; then
        mkdir -p $sampleCallableRegionsTMP
    fi

    # Remove reads with mapping qualities =< 0 and split bams for improved parallelism
    if [ ! -f $sampleCallableRegionsTMP/${sampleid}-sort-${cores}-callable-regions.bed ]; then
        for chr in {1..24}; do
            if [ $chr == 23 ]; then
                chr="X"
            elif [ $chr == 24 ]; then
                chr="Y"
            fi
            samtools view -F 'mapping_quality > 0' -L $sampleCallableRegionsTMP/${sampleid}-sort-${chr}-callable-regions.bed -f bam -l 1 $sampleAlignTMP/${sampleid}-sort.bam | bedtools genomecov -ibam -stdin -bga -g ${refgenome}.fai > $sampleCallableRegionsTMP/${sampleid}-sort-${n}-callable-genomecov.bed
            #sambamba view -F 'mapping_quality > 0' -L $sampleCallableRegionsTMP/${sampleid}-sort-${chr}-callable-regions.bed -f bam -l 1 $sampleAlignTMP/${sampleid}-sort.bam | bedtools genomecov -ibam -stdin -bga -g ${refgenome}.fai > $sampleCallableRegionsTMP/${sampleid}-sort-${n}-callable-genomecov.bed
            cat $sampleCallableRegionsTMP/${sampleid}-sort-${chr}-callable-genomecov.bed | awk $({if ($4 == 0) {print $0`\tNO_COVERAGE`} else if ($4 < 4) {print $0`\tLOW_COVERAGE`} else {print $0`\tCALLABLE`}}) | bedtools groupby -prec 21 -g 1,5 -c 1,2,3,5 -o first,first,max,first | cut -f 3-6 | bedtools intersect -nonamecheck -a - -b $sampleCallableRegionsTMP/${sampleid}-sort-${chr}-callable-regions.bed | sort -V -k1,1 -k2,2n > $sampleCallableRegionsTMP/${sampleid}-sort-${chr}-callable.bed

            # Run Platypus and tabix output vcf
            if [ ! -f $platypus/$chr/tmp/${sampleid%%_*}-${chr}.vcf.gz ]; then
                platypus callVariants --regions=$platypus/$chr/${sampleid%%_*}-${chr}-nolcr.bed --bamFiles=$sampleAlignTMP/${sampleid}-sort.bam --refFile=$refgenome --output=- --logFileName /dev/null --verbosity=1 | awk -F$'\t' -v OFS='\t' '{if ($0 !~ /^#/) gsub(/[KMRUSWBVHDX]/, "N", $4) } {print}' | vcfallelicprimitives -t DECOMPOSED --keep-geno | vcffixup - | vcfstreamsort | bgzip -c > $platypus/$chr/tmp/${sampleid%%_*}-${chr}.vcf.gz
                tabix -f -p vcf $platypus/$chr/tmp/${sampleid%%_*}-${chr}.vcf.gz
            fi

            # Run mutect2 and tabix output vcf
            if [ ! -f $mutect2/$chr/tmp/${sampleid%%_*}-${chr}.vcf.gz ]; then
                java -Xms681m -Xmx1818m -XX:+UseSerialGC -Djava.io.tmpdir=$mutect2/$chr/tmp -jar $gatkexec -T MuTect2 -R $refgenome --annotation ClippingRankSumTest --annotation DepthPerSampleHC --annotation BaseQualityRankSumTest --annotation FisherStrand --annotation GCContent --annotation HaplotypeScore --annotation HomopolymerRun --annotation MappingQualityRankSumTest --annotation MappingQualityZero --annotation ReadPosRankSumTest --annotation RMSMappingQuality --annotation DepthPerAlleleBySample --annotation Coverage -I:tumor $sampleAlignTMP/${sampleid}-sort.bam -L $mutect2/$chr/${sampleid%%_*}-${chr}-regions-nolcr.bed --interval_set_rule INTERSECTION --dbsnp $dbsnp --cosmic $cosmic -ploidy 2 -U LENIENT_VCF_PROCESSING --read_filter BadCigar --read_filter NotPrimaryAlignment | py 'bcbio.variation.mutect2.fix_mutect2_output(x, "$sampleid", "")' | bgzip -c > $mutect2/$chr/tmp/${sampleid%%_*}-${chr}.vcf.gz
                tabix -f -p vcf $mutect2/$chr/tmp/${sampleid%%_*}-${chr}.vcf.gz
            fi

            # Run vardict-java, tabix output vcf, and annotate with bcftools
            if [ ! -f $vardict/$chr/tmp/${sampleid%%_*}-${chr}.vcf.gz ]; then
                unset R_HOME && export VAR_DICT_OPTS='-Xms750m -Xmx3000m -XX:+UseSerialGC -Djava.io.tmpdir=$vardict/$chr/tmp' && VarDictJava -G $refgenome -f $allelefreq -N $sampleid $sampleAlignTMP/${sampleid}-sort.bam -c 1 -S 2 -E 3 -g 4 -Q 10 $vardict/$chr/${sampleid%%_*}-${chr}-unmerged-regions-regionlimit.bed | teststrandbias.R | var2vcf_valid.pl -N $sampleid -E -f $allelefreq | awk -F$'\t' -v OFS='\t' '{if ($0 !~ /^#/) gsub(/ [KMRYSWBVFDX]/, "N", $4) } {print}' | awk -F$'\t' -v OFS='\t' '$1!~/^#/ && $4 == $5 {next} {print}' | vcfstreamsort | bgzip -c $vardict/$chr/tmp/${sampleid%%_*}.vcf.gz
                tabix -f -p vcf $vardict/$chr/tmp/${sampleid%%_*}-${chr}.vcf.gz
                bcftools annotate -c ID -a $dbsnp -o $vardict/$chr/tmp/${sampleid%%_*}-${chr}-wdbsnp.vcf.gz -O z $vardict/$chr/tmp/${sampleid%%_*}-${chr}.vcf.gz
            fi
        done
    fi

    # Concatenate vcfs from all variant callers used
    for variantcaller in $PWD/work/*; do
        gatk-framework org.broadinstitute.gatk.tools.CatVariants -R $refgenome -V $variantcaller/${sampleid%%_*}-files.list -out $variantcaller/tmp/${sampleid%%_*}.vcf.gz -assumeSorted -Xms750m -Xmx2000m -XX:+UseSerialGC

        # Annotate concatenated vcfs with snpEff
        snpEff -Xms750m -Xmx6g -Djava.io.tmpdir=$variantcaller/tmp eff -dataDir $refgenome/snpeff -cancer -noLog -i vcf -o vcf -s $variantcaller/${sampleid%%_*}-effects-stats.html ${refgenome##*/} $variantcaller/${sampleid%%_*}.vcf.gz | bgzip -@ 5 -c > $variantcaller/tmp/${sampleid%%_*}-effects.vcf.gz
        tabix -f -p vcf $variantcaller/tmp/${sampleid%%_*}-effects.vcf.gz

        # Filter annotated vcfs
        bcftools filter -O v -T $bedprep/$bedfile --soft-filter 'PlatQualDepth' -e '(FR[0] <= 0.5 && TC < 4 && %QUAL < 20) || (TC < 13 && %QUAL < 10) || (FR[0] > 0.5 && TC < 4 && %QUAL < 50)' -m '+' $variantcaller/${sampleid%%_*}-effects.vcf.gz | sed 's/\tQ20\t/\tPASS\t/' | bgzip -c > $variantcaller/tmp/${sampleid%%_*}-effects-filter.vcf.gz
        tabix -f -p vcf $variantcaller/tmp/${sampleid%%_*}-effects-filter.vcf.gz

        # Get germline vcfs
        cat $variantcaller/${sampleid%%_*}-effects-filter-germline.vcf | bgzip -@ 5 -c > $variantcaller/tmp/${sampleid%%_*}-effects-filter-germline.vcf.gz
        tabix -f -p vcf $variantcaller/tmp/${sampleid%%_*}-effects-filter-germline.vcf.gz

    done


done

exit $exitstatus
