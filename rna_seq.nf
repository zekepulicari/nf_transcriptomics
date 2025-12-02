#!/usr/bin/env nextflow

// Output params
params.outdir = "${projectDir}/results"

// Input params
params.reads = "${projectDir}/data/reads/*fastq.gz"
params.hisat2_index = "${projectDir}/data/aligned/genome_index.tar.gz"

// Import modules
include { FASTQC } from "${projectDir}/modules/fastqc.nf"
include { TRIM_GALORE } from "${projectDir}/modules/trimming.nf"
include { HISAT2_ALIGN } from "${projectDir}/modules/hisat2_align.nf"

// Main workflow
workflow {

    // Create a channel from input reads
    read_ch = channel.fromPath(params.reads)

    // Call processes
    FASTQC(read_ch)

    // Adapter trimming and post-trimming QC
    TRIM_GALORE(read_ch)

    // Alignment to a reference genome
    HISAT2_ALIGN(TRIM_GALORE.out.trimmed_reads, file(params.hisat2_index))

}
