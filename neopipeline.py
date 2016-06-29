#!/bin/python

import call
import ipyparallel as ipp
import os

class Sample:
    # Global parameters
    str sampleid
    str runid

    # Standard options for Picard tools. More info at http://broadinstitute.github.io/picard/command-line-overview.html
    def __picard__(sampleid):
        str TMP_DIR
        bool QUIET
        str VALIDATION_STRINGENCY
        int COMPRESSION_LEVEL
        int MAX_RECORDS_IN_RAM
        bool CREATE_INDEX
        bool CREATE_MD5_FILE
        str REFERENCE_SEQUENCE
        str GA4GH_CLIENT_SECRETS

    # Parameters for Picard IlluminaBasecallsToFastq tool
    def __picardIlluminaBasecallsToFastq__(sampleid):
        str BASECALLS_DIR
        str BARCODES_DIR
        int LANE
        str OUTPUT_PREFIX
        str RUN_BARCODE
        str MACHINE_NAME
        str FLOWCELL_BARCODE
        str READ_STRUCTURE
        str MULTIPLEX_PARAMS
        str ADAPTERS_TO_CHECK
        int NUM_PROCESSORS
        int FIRST_TILE
        int TILE_LIMIT
        bool APPLY_EAMSS_FILTER
        bool FORCE_GC
        int MAX_READS_IN_RAM_PER_TILE
        int MINIMUM_QUALITY
        bool INCLUDE_NON_PF_READS
        str READ_NAME_FORMAT
        bool COMPRESS_OUTPUTS

    # Parameters for Picard MarkDuplicates tool
    def __picardMarkDuplicates__(sampleid):

    # bwa global parameters. See http://bio-bwa.sourceforge.net/bwa.shtml
    def __bwa__(sampleid):
        str reference
        str sampleid.fastq
        str sampleid.sam

    # Parameters for bwa mem
    def __bwamem__(sampleid):
        int nThreads
        int minSeedLen
        int bandWidth
        int zDropoff
        float seedSplitRatio
        int maxOcc
        bool pairedEndMode
        int matchScore
        int mmPenalty
        int gapOpenPen
        int gapExtPen
        int clipPen
        int unpairPen
        bool pairedEndInterleaved
        str RGline
        str dbPrefix
        str readsFQ
        str matesFQ
        bool allAlign
        bool appendComment
        bool hardClipping
        bool markShortSplit
        int verboseLevel

    # Samtools global options. See http://www.htslib.org/doc/samtools.html
    def __samtools__(sampleid)
        int nthreads
        str reference
        str sampleid.sam
        str sampleid.bam


def main:
    # todo
    # Pipeline is as follows:
    # Alignment and filtering
    # * picard IlluminaBasecallsToFastq
    # * bwa mem to align samples
    # * samtools view to compress samples to .bam
    # * samtools sort
    # * samtools index
    # * picard MarkDuplicates
    # Variant calling
    # * Vardict
    # * GATK - Mutect2
    # Variant annotation
    # * Annovar
    call(["java", picardMem, "-jar", "picard", "IlluminaBasecallsToFastq", picardOptions, "|", "bwa", "mem", bwaOptions, "|", ])
