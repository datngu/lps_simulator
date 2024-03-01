#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --nodes=1                
#SBATCH --job-name=lps   
#SBATCH --mem=4G                
#SBATCH --partition=gpu
#SBATCH --mail-user=nguyen.thanh.dat@nmbu.no
#SBATCH --mail-type=ALL


module load BCFtools/1.10.2-GCC-8.3.0
module load git/2.23.0-GCCcore-9.3.0-nodocs
module load Nextflow/21.03
module load singularity/rpm


export NXF_SINGULARITY_CACHEDIR=/mnt/users/ngda/sofware/singularity


nextflow run main.nf \
    -resume \
    -profile cluster \