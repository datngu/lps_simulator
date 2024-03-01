require(data.table)
require(tidyverse)
setwd('/Users/datn/DATA_ANALYSES/KTEST_codes/nf_lps_simulator/data')
HUMAN_HG_38_SIZE = 3209286105
cov = fread('real_data_coverage.csv')
src = fread('1000G_2504_high_coverage.sequence.index')
samples = fread('sellect_ids-0-Illumina_HiSeq_X.txt', header = FALSE)

cov$n_read = cov$bases_all/150

d1 = cov[cov$assay_type == 'ilmn_0.5x',]
d2 = cov[cov$assay_type == 'ilmn_1x',]

# exclude extreme values
d3 = d1[d1$coverage_all < quantile(d1$coverage_all, 0.99),]
d4 = d2[d2$coverage_all < quantile(d2$coverage_all, 0.99),]

print(median(d1$coverage_all))
print(sd(d3$coverage_all))

print(median(d2$coverage_all))
print(sd(d4$coverage_all))


## process src
pick = c("SAMPLE_NAME", "#ENA_FILE_PATH", "READ_COUNT", "POPULATION", "MD5SUM")
src = src[,..pick]
colnames(src) = pick = c("SAMPLE_NAME", "ENA_FILE_PATH", "READ_COUNT", "POPULATION", "MD5SUM")
src = src[src$SAMPLE_NAME %in% samples$V1]
src$ENA_FILE_PATH = sub('ftp://', 'https://', src$ENA_FILE_PATH)
fwrite(src, file = '60_samples.csv', sep = ',', row.names = F)

