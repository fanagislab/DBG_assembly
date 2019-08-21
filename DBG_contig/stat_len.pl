use strict;

my $contig_file = shift;

 `perl /ifs1/GAG/assemble/fanw/bin/fastaDeal.pl  -attr id:len $contig_file > $contig_file.len`;
 `perl /ifs1/GAG/assemble/fanw/Assembly/bin/seqlen_stat.pl -column 2 -cutoff 100 $contig_file.len > $contig_file.len.cut100.stat`;
`perl /ifs1/GAG/assemble/fanw/Assembly/bin/seqlen_stat.pl -column 2 -cutoff 0 $contig_file.len > $contig_file.len.all.stat`;
