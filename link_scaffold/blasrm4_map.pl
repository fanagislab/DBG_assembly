#!/usr/bin/perl

=head1 Name

blasrm4_map.pl  --  convert the blasrm4 best hits to the mapping scafftigs for each pacbio_reads

=head1 Description

Input is the best hit file generated by blasrm4_besthit.pl
Output is to STDOUT, and the statistics to STEDERR

=head1 Version

  Author: Fan Wei, fanw@genomics.org.cn
  Version: 1.0,  Date: 2006-12-6
  Note:

=head1 Usage
  
  blasrm4_map.pl [options] <blasrm4_besthit_file>
  --endlencut <int>      max unaligned end length allowed for semiglobal alignment, default=100
  --alignlencut <int>    min overlap length required for confident alignment, default=1000
  --identitycut <float>  min sequence identity required for confident alignment, default=0.7
  --verbose   output running progress information to screen  
  --help      output help information to screen  

=head1 Exmple

    blasrm4_map.pl pacbio30X_vs_scafftig.blasrm4.best > pacbio30X_vs_scafftig.blasrm4.best.map 2> pacbio30X_vs_scafftig.blasrm4.best.map.stat

=cut

use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname); 
use Data::Dumper;


##get options from command line into variables and set default values
my ($AlignEndLenCut, $AlignLenCut, $AlignIdentityCut);
my ($Verbose,$Help);
GetOptions(
	"endlencut:i"=>\$AlignEndLenCut, 
	"alignlencut:i"=>\$AlignLenCut,
	"identitycut:f"=>\$AlignIdentityCut,
	"verbose"=>\$Verbose,
	"help"=>\$Help
);
$AlignEndLenCut ||= 100;
$AlignLenCut ||= 1000;
$AlignIdentityCut ||= 0.7;
die `pod2text $0` if (@ARGV == 0 || $Help);


my $input_blasrm4_file = shift;


my $total_alignment_num = 0;
my $unconfident_alignment_num = 0;
my $nonsemiglobal_alignment_num = 0;
my $multiple_alignment_num = 0;
my $fine_alignment_num = 0;


my %AlignData;

##qName tName score percentSimilarity qStrand qStart qEnd qLength tStrand tStart tEnd tLength mapQV
open IN, $input_blasrm4_file || die "fail $input_blasrm4_file";
while (<IN>)
{	my ($qid,$tid,$score, $identity, $qstrand, $qstart,$qend,$qlen, $tstrand,$tstart,$tend, $tlen) = split /\s+/;
	#$qid = $1 if($qid =~ /^(\w+)/);
	my $qrate = ($qend - $qstart) / $qlen;
	my $trate = ($tend - $tstart) / $tlen;
	

	if ($qend - $qstart >= $AlignLenCut && $tend - $tstart >= $AlignLenCut && $identity > $AlignIdentityCut ){
		push @{$AlignData{$qid}},  [$qstart,$qend, $qlen, $qrate, $tid, $tstart,$tend, $tlen, $trate, $tstrand, $identity];
	}else{  ##filter the unconfident alignment firstly
		$unconfident_alignment_num ++;
	}
	$total_alignment_num ++;
}
close IN;


my $output_line_num = 0;
my $output_align_num = 0;

print "#pacbio_Id\tscafftig_num\tpacbio_start\tpacbio_end\tpacbio_length\tpacbio_coverage\tscafftig_id\tscafftig_start\tscafftig_end\tscafftig_length\tscafftig_coverage\talign_strand\talign_identity\n";
foreach my $qid (sort keys %AlignData)
{	my $qid_p = $AlignData{$qid};
	my @ary = sort {$a->[0] <=> $b->[0]} @$qid_p;
	my $num = @ary;
	my @ary2;
	my @ary3;
	
	##get alignment that meets the requirement of semi-global alignment, i.e. query and target are fragments from the same genome
	foreach my $p (@ary){
		my ($qstart,$qend, $qlen, $qrate, $tid, $tstart,$tend, $tlen, $trate, $tstrand, $identity) = @$p;
		my $qleft = ($qstart - 0) ;
		my $qright = ($qlen - $qend) ;
		my $tleft = ($tstart - 0) ;
		my $tright = ($tlen - $tend) ;
		
		if ( ($qleft > $AlignEndLenCut && $tleft > $AlignEndLenCut) || ($qright > $AlignEndLenCut && $tright > $AlignEndLenCut)  ){
			$nonsemiglobal_alignment_num ++;
			next;
		}
		push @ary2, $p;
	}

	##get rid of the overlapped fragments, the one with lower mapping identity is thought to be wrong mapping and removed
	next if(@ary2 < 1);

	##put overlapped alignments into each groups
	my @group;
	push @group, [ $ary2[0] ];
	for (my $i = 1; $i < @ary2; $i ++)
	{	if ( $ary2[$i][0] < $ary2[$i-1][1] )   ##this alignment has overlap to the last alignment
		{	push @{$group[-1]}, $ary2[$i];
		}else
		{	push @group, [ $ary2[$i] ];
		}
	}
	
	##find and keep the alignment with max identity in each group
	foreach my $group_p (@group)
	{	
		my $max_p;

		if (@$group_p > 1)  ##has overlapped fragments
		{	my $max_identity = 0.0;
			foreach my $p (@$group_p)
			{	if ($p->[-1] > $max_identity )
				{	$max_identity  = $p->[-1];
					$max_p = $p;
				}
			}
			$multiple_alignment_num += @$group_p - 1;
		}else
		{	$max_p = $group_p->[0];
		}

		push @ary3, $max_p;
	}
	
	##output the clean scafftig alignment results for each pacbio-reads
	my $num = @ary3;
	my $out_head = "$qid\t$num";
	my $out_content;
	$output_line_num ++;
	foreach my $p ( @ary3)
	{	$output_align_num ++;
		foreach (@$p){
			$out_content .= "\t".$_;
		}
		#$out_content .= "\n";
	}

	print $out_head.$out_content."\n";

}

$fine_alignment_num = $total_alignment_num - $unconfident_alignment_num - $nonsemiglobal_alignment_num - $multiple_alignment_num;

print STDERR "Total alignment number: ".$total_alignment_num."\n";
print STDERR "Deleted Unconfident alignment number: ".$unconfident_alignment_num."\n";
print STDERR "Deleted Nonsemiglobal alignment number: ".$nonsemiglobal_alignment_num."\n";
print STDERR "Deleted multiple alignment number: ".$multiple_alignment_num."\n";
print STDERR "Usable fine alignment number: ".$fine_alignment_num."\n";
print STDERR "Output lines (reads) number: ".$output_line_num."\n";

####################################################
################### Sub Routines ###################
####################################################

