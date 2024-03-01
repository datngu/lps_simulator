#!/bin/bash

echo "Downloading"

wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz -O data/hg38.fa.gz

echo "Unzip...!"

gunzip data/hg38.fa.gz

echo "DONE!"