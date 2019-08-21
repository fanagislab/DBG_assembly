../../../link_scaffold/map_pair -l 125 -r 250 -o ./maping_results/ ./Ecoli_corrected_reads.contig.seq.fa  ./raw_reads.lib
../../../link_scaffold/link_scaffold -i 400 -o Ecoli_corrected_reads.contig  Ecoli_corrected_reads.contig.seq.fa ./raw_reads.lib.map_pair.2ctg.lib
fastaDeal.pl -attr id:len Ecoli_corrected_reads.contig.insert400.scaffold.seq.fa > Ecoli_corrected_reads.contig.insert400.scaffold.seq.fa.len
seqlen_stat.pl -col 2 Ecoli_corrected_reads.contig.insert400.scaffold.seq.fa.len > seqlen_stat.pl -col 2 Ecoli_corrected_reads.contig.insert400.scaffold.seq.fa.len.stat
