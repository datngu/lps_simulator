#!/usr/bin/env nextflow
/*
========================================================================================
                                LPS-simulator
========================================================================================
 Author: Dat T Nguyen
 Contact: ndat<at>utexas.edu
----------------------------------------------------------------------------------------
*/



/*
 Define the default parameters
*/ 
params.metadata        = "$baseDir/data/10_samples.csv"
params.genome          = "$baseDir/data/hg38.fa"

params.trace_dir       = "trace_dir"
params.outdir          = "results"

nextflow.enable.dsl=2

workflow {

    // Create a channel of samples from the metadata
    samplesChannel = Channel.fromPath(params.metadata) \
        | splitCsv(header:true) \
        | map { row-> tuple(row.SAMPLE_NAME, row.ENA_FILE_PATH, row.READ_COUNT, row.POPULATION, row.MD5SUM)}

    Sample2BAM(samplesChannel, params.genome)
}

process Sample2BAM {
    input:
    tuple val(sample_id), val(url_path), val(read_count), val(population), val(md5)
    path genome

    cpus 1
    memory '8GB'

    output:
    path "${sample_id}.bam"

    script:
    """
    wget ${url_path} -O ${sample_id}.cram

    md5sum ${sample_id}.cram
    
    samtools view -b -T ${genome} -o ${sample_id}.bam ${sample_id}.cram

    rm ${sample_id}.cram
    """
}