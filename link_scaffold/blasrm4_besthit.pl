#!/usr/bin/perl

=head1 Name

blasrm4_besthit.pl  --  choose the best aligning hit for each pacbio-reads and scafftig pairs

=head1 Description

The best aligning hit will have the highest value of overlap_length * identity

=head1 Version

  Author: Fan Wei, fanw@genomics.org.cn
  Version: 1.0,  Date: 2007-2-7

=head1 Usage

  --fileformat   set input file format, paf, blasrm4 are supported
  --cutoff       set a cut off to filter low quality alignments
  --verbose      output verbose information to screen  
  --help         output help information to screen  

=head1 Exmple

blasrm4_besthit.pl  pacbio30X_vs_scafftig.blasrm4 > pacbio30X_vs_scafftig.blasrm4.best

=cut

use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname); 
use Data::Dumper;

my ($Fileformat,$Cutoff,$Verbose,$Help);
GetOptions(
	"fileformat:s"=>\$Fileformat,
	"cutoff:s"=>\$Cutoff,
	"verbose"=>\$Verbose,
	"help"=>\$Help
);
die `pod2text $0` if (@ARGV == 0 || $Help);

my %Data;

read_paf() if($Fileformat eq 'paf' || $ARGV[0] =~ /\.paf$/);
read_blasrm4() if($Fileformat eq 'blasrm4' || $ARGV[0] =~ /\.blasrm4$/);


foreach my $qname (sort keys %Data) {
	print $Data{$qname}{line},"\n";
}

####################################################
################### Sub Routines ###################
####################################################





##read the paf format by minimap
sub read_paf{
	while (<>) {
		chomp;
		my @temp = split /\t/;
		my $qname = $temp[0]."-".$temp[5];
		my $value = ($temp[3] - $temp[2])/$temp[1]; ##overlap length / read length, i.e.  overlap ratio
		next if(defined $Cutoff && $value > $Cutoff);
		if (! exists $Data{$qname} || $Data{$qname}{value} < $value) {
			$Data{$qname}{value} = $value;
			$Data{$qname}{line} = $_;
		}
	}
}




##read the blasrm4 format by blasr with -m 4 option
#qName tName score percentSimilarity qStrand qStart qEnd qLength tStrand tStart tEnd tLength mapQV
sub read_blasrm4{
	while (<>) {
		next if(/^qName/);
		chomp;
		my @temp = split /\s+/;
		my $qname = $temp[0]."-".$temp[1];
		my $value = ($temp[6] - $temp[5]) * $temp[3];  ##overlap length * identity
		my $len = ($temp[6] - $temp[5]);
		##print STDERR $qname."\t".$value."\t".$temp[3]."\t".$len."\n";
		next if(defined $Cutoff && $value > $Cutoff);
		if (! exists $Data{$qname} || $Data{$qname}{value} < $value) {
			$Data{$qname}{value} = $value;
			$Data{$qname}{line} = $_;
		}
	}
}


