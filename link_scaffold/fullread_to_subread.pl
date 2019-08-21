#!/usr/bin/perl

=head1 Name



=head1 Description



=head1 Version

  Author: Fan Wei, fanw@genomics.org.cn
  Version: 1.0,  Date: 2006-12-6
  Note:

=head1 Usage

  --verbose   output running progress information to screen  
  --help      output help information to screen  

=head1 Exmple

perl fullread_to_subread.pl  Ecoli_pacbio_30X_0001.fastq  m151107_052705_42221_c100913462550000001823202104291630_s1_p0

=cut

use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname); 
use Data::Dumper;


##get options from command line into variables and set default values
my ($Verbose,$Help);
GetOptions(
	"verbose"=>\$Verbose,
	"help"=>\$Help
);
die `pod2text $0` if (@ARGV == 0 || $Help);


my $pbsim_reads_file = shift;
my $smart_cell_id = shift;

my $zmw_hole_id = 1;

open IN, $pbsim_reads_file || die "fail $pbsim_reads_file";
while (<IN>) {
		if (/^@/) {
			my $read_head = $_;
			chomp $read_head;
			my $read_seq = <IN>;
			chomp $read_seq;
			my $read_len = length($read_seq);
			$read_head = "\@$smart_cell_id/$zmw_hole_id/0_$read_len RQ=0.84";
			my $read_head2 = <IN>;
			$read_head2 = "+\n";
			my $read_qual = <IN>;
			print "$read_head\n$read_seq\n$read_head2$read_qual";

			$zmw_hole_id ++;
		}
}
close IN;



