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
params.metadata        = "$baseDir/data/60_samples.csv"
params.genome          = "$baseDir/data/hg38.fa"

params.trace_dir       = "trace_dir"
params.outdir          = "results"

nextflow.enable.dsl=2

workflow {

    // Create a channel of samples from the metadata
    samplesChannel = file(params.metadata).readCSV().map { row ->
        tuple(row.sample_id, row.url_path, row.read_count, row.population, row.md5)
    }

    Sample2BAM(samplesChannel, params.genome)
}

process Sample2BAM {
    input:
    tuple sample_id, url_path, read_count, population, md5
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

    """
}