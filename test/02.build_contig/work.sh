
ls ../01.clean_correct/*.correct.fa.gz
ls ../01.clean_correct/*.correct.fa.gz > corrected_reads.lib


../../DBG_contig/debruijn_contig -k 31 -r 250  -t 10  -i 0.1 -M 125  -o Ecoli_raw_reads ./raw_reads.lib  2> Ecoli_raw_reads.contig.log &

../../DBG_contig/debruijn_contig -f 2 -k 31 -r 250 -t 10 -i 0.1 -M 125 -o Ecoli_corrected_reads ./corrected_reads.lib 2> Ecoli_corrected_reads.contig.log &


fastaDeal.pl -attr id:len Ecoli_corrected_reads.contig.seq.fa > Ecoli_corrected_reads.contig.seq.fa.len
fastaDeal.pl -attr id:len Ecoli_raw_reads.contig.seq.fa > Ecoli_raw_reads.contig.seq.fa.len


seqlen_stat.pl -col 2 Ecoli_corrected_reads.contig.seq.fa.len > Ecoli_corrected_reads.contig.seq.fa.len.stat
seqlen_stat.pl -col 2 Ecoli_raw_reads.contig.seq.fa.len > Ecoli_raw_reads.contig.seq.fa.len.stat
