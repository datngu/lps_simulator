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

    DownSampling(Sample2BAM.out)
}

process Sample2BAM {

    publishDir "${params.trace_dir}/hc_bam", mode: 'symlink', overwrite: true

    input:
    tuple val(sample_id), val(url_path), val(read_count), val(population), val(md5)
    path genome

    cpus 1
    memory '8GB'

    output:

    tuple val("${sample_id}"), val("${read_count}"), path ("${sample_id}.bam")

    script:
    """

    ## try 5 times

    for ((i=1; i<=5; i++)); do
        wget ${url_path} -O ${sample_id}.cram

        actual_md5="\$(md5sum ${sample_id}.cram | awk '{print \$1}')"

        if [ "\$actual_md5" == "$md5" ]; then
            echo "MD5 sum matches! File downloaded successfully."
            break
        else
            echo "MD5 sum does not match. Retrying..."
            rm -f ${sample_id}.cram
            wget ${url_path} -O ${sample_id}.cram
        fi
    done

    samtools view -b -T ${genome} -o ${sample_id}.bam ${sample_id}.cram
    
    ## release space
    rm -f ${sample_id}.cram

    """
}



process DownSampling {

    publishDir "${params.trace_dir}/lowpass_bam_files", mode: 'symlink', overwrite: true

    input:
    tuple val(sample_id), val(read_count), path(bam_file)

    cpus 1
    memory '16GB'

    output:
    path "*_lps.bam"

    script:
    """

    bam_sampling.py --bam ${sample_id}.bam \
        --depth 0.3 0.5 0.67 0.8 1.0 1.25 1.5 2.0 \
        --out ${sample_id}_0.3_lps.bam ${sample_id}_0.5_lps.bam ${sample_id}_0.67_lps.bam ${sample_id}_0.8_lps.bam ${sample_id}_1.0_lps.bam ${sample_id}_1.25_lps.bam ${sample_id}_1.5_lps.bam ${sample_id}_2.0_lps.bam \
        --bam_size ${read_count}


    """
}

