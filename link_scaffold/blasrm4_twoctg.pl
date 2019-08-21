#!/usr/bin/perl

=head1 Name

blasrm4_twoctg.pl  convert the blasr map file into 2ctg format which can be used by link_scafftig

=head1 Description

The first input file is blasr_besthit_map file
The second input file is optional, in table list format, column 1 must be repeat scafftig IDs, the alignments of which should be filtered

=head1 Version

  Author: Fan Wei, fanw@genomics.org.cn
  Version: 1.0,  Date: 2006-12-6
  Note:

=head1 Usage

  blasrm4_twoctg.pl   <*.blasrm4.best.map>  [*.repeatctg.list]
  --verbose   output running progress information to screen  
  --help      output help information to screen  

=head1 Exmple


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

my $input_file = shift;
my $repeatctg_file = shift;



my %RepeatCtg;

if (-f $repeatctg_file){
	open IN, $repeatctg_file || die "fail";
	while(<IN>)
	{	if(/^(\S+)/){
			$RepeatCtg{$1} = 1;
		}
	}
	close IN;
}

print STDERR "Repeat Nodes:\n";
print STDERR Dumper \%RepeatCtg;


##S1_1    1       11      7431    7440    0.997311827956989       sct_27  68844   75937   86966   0.0815606098935216      1       87.5964
open IN, $input_file || die "fail $input_file ";
while (<IN>)
{	chomp;
	my @t = split /\s+/;
	my $reads_id = shift @t;
	my $hit_num = shift @t;
	next if($hit_num < 2);
	
	my @FilteredRecords;
	for (my $i = 0; $i < @t; $i += 11)
	{	
		my $curr_i = $i ;
		my $read_start = $t[$curr_i];
		my $read_end = $t[$curr_i+1];
		my $read_len  = $t[$curr_i+2];

		my $ctg_id = $t[$curr_i+4];
		my $strand = $t[$curr_i+9];
		my $identity = $t[$curr_i+10];
		my $ctg_len = $t[$curr_i+7];
		my $ctg_start = ($strand == 0) ? $t[$curr_i+5] : ($ctg_len - $t[$curr_i+6]);
		my $ctg_end = ($strand == 0) ? $t[$curr_i+6] : ($ctg_len - $t[$curr_i+5]);
		my $strand_char = ($strand == 0) ? "F" : "R";
		
		##ignore the aligments whose scafftig ids existed in the second input file [repeat nodes]
		if (! exists $RepeatCtg{$ctg_id}){
			push @FilteredRecords, "$reads_id\t$read_len\t$read_start\t$read_end\t$ctg_id\t$ctg_len\t$ctg_start\t$ctg_end\t$strand_char\t$identity%";
		}else{
			print STDERR "$ctg_id records deleted\n";
		}
	}

	next if(@FilteredRecords < 2);
	
	for (my $i = 1; $i < @FilteredRecords; $i ++)
	{	print $FilteredRecords[$i-1]."\t".$FilteredRecords[$i]."\n";
	
	}

}
close IN;


