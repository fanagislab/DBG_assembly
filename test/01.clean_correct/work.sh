
##remove the low quality sequences in reads
../../clean_illumina/clean_lowqual -e 0.01 -r 75 ../00.raw_reads/Ecoli_readlen250_insert400_20X_250_400_1.fq.gz ./Ecoli_readlen250_insert400_20X_250_400_1.fq.gz.nonLowQual.gz ./Ecoli_readlen250_insert400_20X_250_400_1.fq.gz.nonLowQual.stat &
../../clean_illumina/clean_lowqual -e 0.01 -r 75 ../00.raw_reads/Ecoli_readlen250_insert400_20X_250_400_2.fq.gz ./Ecoli_readlen250_insert400_20X_250_400_2.fq.gz.nonLowQual.gz ./Ecoli_readlen250_insert400_20X_250_400_2.fq.gz.nonLowQual.stat &
../../clean_illumina/clean_lowqual -e 0.01 -r 75 ../00.raw_reads/Ecoli_readlen250_insert800_20X_250_800_1.fq.gz ./Ecoli_readlen250_insert800_20X_250_800_1.fq.gz.nonLowQual.gz ./Ecoli_readlen250_insert800_20X_250_800_1.fq.gz.nonLowQual.stat &
../../clean_illumina/clean_lowqual -e 0.01 -r 75 ../00.raw_reads/Ecoli_readlen250_insert800_20X_250_800_2.fq.gz ./Ecoli_readlen250_insert800_20X_250_800_2.fq.gz.nonLowQual.gz ./Ecoli_readlen250_insert800_20X_250_800_2.fq.gz.nonLowQual.stat &

## remove the adapter contaminations in reads
../../clean_illumina/clean_adapter -a ../../clean_illumina/illumina_NEB_adapter.fa -r 75 -s 12 Ecoli_readlen250_insert400_20X_250_400_1.fq.gz.nonLowQual.gz Ecoli_readlen250_insert400_20X_250_400_1.fq.gz.nonLowQual.gz.nonAdapter.gz Ecoli_readlen250_insert400_20X_250_400_1.fq.gz.nonLowQual.gz.nonAdapter.stat &
../../clean_illumina/clean_adapter -a ../../clean_illumina/illumina_NEB_adapter.fa -r 75 -s 12 Ecoli_readlen250_insert400_20X_250_400_2.fq.gz.nonLowQual.gz Ecoli_readlen250_insert400_20X_250_400_2.fq.gz.nonLowQual.gz.nonAdapter.gz Ecoli_readlen250_insert400_20X_250_400_2.fq.gz.nonLowQual.gz.nonAdapter.stat &
../../clean_illumina/clean_adapter -a ../../clean_illumina/illumina_NEB_adapter.fa -r 75 -s 12 Ecoli_readlen250_insert800_20X_250_800_1.fq.gz.nonLowQual.gz Ecoli_readlen250_insert800_20X_250_800_1.fq.gz.nonLowQual.gz.nonAdapter.gz Ecoli_readlen250_insert800_20X_250_800_1.fq.gz.nonLowQual.gz.nonAdapter.stat &
../../clean_illumina/clean_adapter -a ../../clean_illumina/illumina_NEB_adapter.fa -r 75 -s 12 Ecoli_readlen250_insert800_20X_250_800_2.fq.gz.nonLowQual.gz Ecoli_readlen250_insert800_20X_250_800_2.fq.gz.nonLowQual.gz.nonAdapter.gz Ecoli_readlen250_insert800_20X_250_800_2.fq.gz.nonLowQual.gz.nonAdapter.stat &



##error correction by k-mer frequency
ls ./*.nonAdapter.gz > clean_reads.lib
../../kmerfreq/kmerfreq -k 17 -m 1 -q 10 ./clean_reads.lib
../../correct_error/correct_error_reads -k 17 -c 2 ./clean_reads.lib.kmer.freq.cz  ./clean_reads.lib &

ls ../00.raw_reads/*.fq.gz > raw_reads.lib
../../kmerfreq/kmerfreq -k 17 raw_reads.lib &

ls ./*.correct.fa.gz > corrected_reads.lib
../../kmerfreq/kmerfreq -k 17 -f 2 corrected_reads.lib &
