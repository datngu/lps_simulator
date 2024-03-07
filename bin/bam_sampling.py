#!/usr/bin/env python

import argparse
import pysam
import subprocess
import random


## set similar to https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8015847/
GENOME_SIZE = 3.3e9
## Standard deviation when doing lowpass sequencing
## Emperically computed from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8015847/
SD_COV_0_5X = 0.2314565
SD_COV_1X = 0.2138574
SD_COV = (SD_COV_0_5X + SD_COV_1X) /2
## read length, 1000 NYGC dataset is 150bp
READ_LENGTH = 150

def depth2count(depth):
    # avoid zero sampling
    ## 1.28 is z-score to exclude values < quantile(0.1) and > quantile(0.9) of the distribution
    min_cut = max(0.1, depth - 1.28 * SD_COV)
    max_cut = depth + 1.28 * SD_COV
    ran_depth = random.gauss(depth, SD_COV)
    
    ## repeat sampling until have the expected sampling
    while not (min_cut <= ran_depth <= max_cut):
        ran_depth = random.gauss(depth, SD_COV)

    n_read = (ran_depth * GENOME_SIZE) // READ_LENGTH
    return int(n_read)




def count_total_reads(bam_path):
    command = ["samtools", "index", bam_path]
    subprocess.run(command, check=True)
    bam = pysam.AlignmentFile(bam_path, "rb")
    index_stats = bam.get_index_statistics()
    return sum(i[3] for i in index_stats)


def write_read(read, idx, bam_obj, set_idxs):
    if idx in set_idxs:
        bam_obj.write(read) 


def parse_args():
    parser = argparse.ArgumentParser(description="Process BAM files with depth values.")
    # Required arguments
    parser.add_argument("--bam", help="Path to the BAM file.")
    parser.add_argument("--depth", nargs="+", type=float, help="List of float values of target depths.")
    parser.add_argument("--out", nargs="+", help="List of output file paths.")
    # Optional argument
    parser.add_argument("--bam_size", type=int, help="Optional parameter for BAM size.")
    return parser.parse_args()


if __name__ == "__main__":
    
    args = parse_args()
    # Accessing the arguments
    print("BAM file:", args.bam)
    print("Depth values:", args.depth)
    print("Output files:", args.out)
    assert len(args.depth) == len(args.out), 'number of given --depth must equal number of outputs'
    

    if args.bam_size:
        BAM_SIZE = args.bam_size
        print('uses prodived bam_size!')
    else:
        print('BAM size will be computed!')
        BAM_SIZE = count_total_reads(args.bam)

    print("BAM size:", BAM_SIZE)


    ## sampling reads to pick for each depth
    pick_reads = [set(random.sample(range(0, BAM_SIZE), depth2count(d))) for d in args.depth]

    ## open input bam file
    main_bam = pysam.AlignmentFile(args.bam, "rb")

    ## prepare output bam files
    out_bams = []
    for fo in args.out:
        tem_bam = pysam.AlignmentFile(fo, "wb", header=main_bam.header)
        out_bams.append(tem_bam)

    #main_bam.close()
    #tem_bam.close()
    idx = 0
    for read in main_bam:
        
        for bam_obj, set_idxs in zip(out_bams, pick_reads):
            write_read(read, idx, bam_obj, set_idxs)

        idx += 1

    ## close all_files
    main_bam.close()
    for fo in out_bams:
        fo.close()




# python bin/bam_sampling.py --bam data/test_NA12878.bam --depth 0.01 0.002 --out output1.bam output2.bam

# python bam_sampling.py --bam NA12842.bam --depth 0.67 1.25 --out '0.67.bam' '1.25.bam' --bam_size 729859386