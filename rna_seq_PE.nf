#!/usr/bin/env nextflow

// Output params
params.outdir = "${projectDir}/results"

// Input params
params.reads = "${projectDir}/data/reads/*fastq.gz"
params.hisat2_index = "${projectDir}/data/aligned/genome_index.tar.gz"
params.paired_end = "${projectDir}/data/paired-end.csv"
params.report_id = "paired-end_reports"

// Import modules
include { FASTQC } from "${projectDir}/modules/fastqc_pe.nf"
include { TRIM_GALORE } from "${projectDir}/modules/trimming_pe.nf"
include { HISAT2_ALIGN } from "${projectDir}/modules/hisat2_align_pe.nf"
include { MULTI_QC } from "${projectDir}/modules/multiqc.nf"


// Main workflow
workflow {
 
    // Create input channel from the contents of a CSV file
    read_ch = channel.fromPath(params.paired_end)
        .splitCsv(header:true)
        .map { row -> [file(row.fastq_1), file(row.fastq_2)] }
    
    // Call processes
    FASTQC(read_ch)

    // Adapter trimming and post-trimming QC
    TRIM_GALORE(read_ch)

    // Alignment to a reference genome
    HISAT2_ALIGN(TRIM_GALORE.out.trimmed_reads, file(params.hisat2_index))

    ch_all_outputs = FASTQC.out.zip.mix(
        FASTQC.out.html,
        TRIM_GALORE.out.trimmed_reads,
        TRIM_GALORE.out.trimming_reports,
        TRIM_GALORE.out.fastqc_reports_1,
        TRIM_GALORE.out.fastqc_reports_2,
        HISAT2_ALIGN.out.log
    )

    MULTI_QC(ch_all_outputs.collect(), params.report_id)
}
