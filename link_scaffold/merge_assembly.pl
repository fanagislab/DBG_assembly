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

=head1 Example

perl merge_assembly.p --output_prefix Tmulticeps  scafftig_vs_smartdenovo.psl.best   Tmulticeps.scafftig.seq.all.fa    wt.zmo.cns.fa

=cut

use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname); 
use Data::Dumper;


##get options from command line into variables and set default values
my ($Output_prefix, $SeqIdPrefix, $Verbose,$Help);
GetOptions(
	"output_prefix:s"=>\$Output_prefix,
	"seqidprefix:s"=>\$SeqIdPrefix,
	"verbose"=>\$Verbose,
	"help"=>\$Help
);
$Output_prefix ||= "Output";
$SeqIdPrefix ||= "TMC_";
die `pod2text $0` if (@ARGV == 0 || $Help);

my $input_psl = shift;
my $scafftig_fa = shift;
my $utg_fa = shift;

##require these 3 conditions to define a confident alignment
my $AlignLenCut = 125;
my $IdentityCut = 0.9;
my $QueryRateCut = 0.5;

##store the scafftig and utg sequences
my %SCAFFTIG;
my %UTG;

my %AlignedScafftig;  ##store scafftigs that aligned to utg
my %AlignedUTG;  ##store utgs that align to scafftig
my %AlignData;  ##store all the alignment information

Read_fasta($scafftig_fa, \%SCAFFTIG);

Read_fasta($utg_fa, \%UTG);


my $lines_num = 0;

##Load the blat psl best result
open IN, $input_psl  || die "fail $input_psl";
while (<IN>){
	chomp;
	my ($match_bases, $mismatch_bases, $repeat_match, $Ns, $QgapCount, $QgapBases, $TgapCount, $TgapBases, $strand, $Qname, $Qsize, $Qstart, $Qend, $Tname, $Tsize, $Tstart, $Tend, $block_count, $blockSizes, $qStarts, $tStarts) = split /\s+/;
	
	
	my $Qrate = ($Qend - $Qstart) / $Qsize;
	my $Trate = ($Tend - $Tstart) / $Tsize;
	my $identity = $match_bases / ($match_bases + $mismatch_bases + $QgapBases + $TgapBases);
	
	if ( ($identity > $IdentityCut) && ($Qend - $Qstart > $AlignLenCut) && ($Qrate > $QueryRateCut)  ){
		push @{$AlignData{$Tname}},  [$Tstart,$Tend, $Tsize, $Trate, $Qname, $Qstart,$Qend, $Qsize, $Qrate, $strand, $identity];
		$AlignedScafftig{$Qname} = 1;
		$AlignedUTG{$Tname} = 1;
		
	}

	$lines_num ++;
	
}
close IN;



##print Dumper \%AlignData;

open SEQFILE, ">$Output_prefix.merged_assembly.seq.fa";
open POSFILE, ">$Output_prefix.merged_assembly.pos.tab";

print POSFILE "#Seq_Id\tblock_start\tblock_end\tblock_size\ttarget_block_start\ttarget_block_end\ttarget_block_size\ttarget_name\t+\toutput_block_length\toutput_block_sequence\n";


my $Aligned_utg_num = 0;
my $Unaligned_utg_num = 0;
my $Unaligned_scafftig_num = 0;

my $SeqIdNum = 0;

##replace the pacbio-cns with illumina-contig for each utg that has alignment to scafftigs
foreach my $Tname (sort keys %AlignData)
{	
	$Aligned_utg_num ++;
	$SeqIdNum ++;
	
	my $Tname_p = $AlignData{$Tname};
	my @ary2 = sort {$a->[0] <=> $b->[0]} @$Tname_p;   ##sorted two-dimension array
	my @ary3;  ##store the non-redundant two-dimension array
	
	##layout out to get a non-redundant contig set [remove those contigs which are fully included in others]
	push @ary3, $ary2[0];
	for (my $i = 1; $i < @ary2; $i ++){
		if ( $ary2[$i][1] > $ary3[-1][1] + 100 ) 
		{	push @ary3, $ary2[$i];
		}
	}
	
	my @Out; ##store the results for preparing output
	for (my $i = 0; $i < @ary3; $i ++){
				
		my $sct_p = $ary3[$i];
		my ($Tstart,$Tend, $Tsize, $Trate, $Qname, $Qstart,$Qend, $Qsize, $Qrate, $strand, $identity) = @$sct_p;
		
		my $gap_size;
		if ($i == 0){
			$gap_size = $Tstart - 0;
		}else{
			my $last_p = $ary3[$i-1];
			my ($last_Tstart, $last_Tend) = @$last_p;
			$gap_size = $Tstart - $last_Tend;
		}

		my $Qseq = $SCAFFTIG{$Qname};
		
		if ($strand eq "-"){
			Complement_Reverse(\$Qseq);
			my $qstart = $Qsize - $Qend;
			my $qend = $Qsize - $Qstart;
			$Qstart = $qstart;
			$Qend = $qend;
		}##if minus alignment, the sequence should be reverse-complement, the coordinates should be reversed
		
		if ($gap_size < 0){
			$Qseq = substr($Qseq, $Qstart+abs($gap_size), $Qend-$Qstart-abs($gap_size));
		}else{
			$Qseq = substr($Qseq, $Qstart, $Qend-$Qstart);
		}

		push @Out, [$Tstart,$Tend, $Qname, $Qsize, $Qstart, $Qend, $strand, $gap_size, $Qseq];

	}
	
	my $Tseq = $UTG{$Tname};
	my $final_seq;
	my $final_seq_len = 0;

	my $pos_output;
	my $ctg_all_names;
	for (my $i=0; $i<@Out; $i++){
		
		my ($Tstart,$Tend, $Qname, $Qsize, $Qstart, $Qend, $strand, $gap_size, $ctg_seq) = @{$Out[$i]};
		
		my $gap_seq = "";   ##pacbio-cns
		my $gap_end = $Tstart;
		my $gap_start = $Tstart;
		if($gap_size > 0){
			$gap_end = $Tstart;
			$gap_start = $gap_end - $gap_size;
			$gap_seq = substr($Tseq, $gap_start, $gap_size);
		}
		
		$final_seq .= $gap_seq.$ctg_seq;
		
		
		my $block_start = $final_seq_len + 1;
		my $block_size = ($gap_size > 0) ? $gap_size : 0;
		$final_seq_len += $block_size;
		$pos_output .= "$SeqIdPrefix$SeqIdNum\t$block_start\t$final_seq_len\t$block_size\t$gap_start\t$gap_end\t$gap_size\t$Tname\t+\t".length($gap_seq)."\t".$gap_seq."\n";
		
		
		my $ctg_start = $Tstart;  ##illumina-scafftig
		my $ctg_end = $Tend;
		my $ctg_size = $ctg_end  - $ctg_start;
		my $ctg_name = $Qname;
		
		$ctg_all_names .= ",".$Qname;
		$block_start = $final_seq_len + 1;
		$block_size = length($ctg_seq);
		$final_seq_len += $block_size;
		
		my $QalignLen = $Qend - $Qstart;
		my $output_ctg_len = length($ctg_seq);
		$pos_output .= "$SeqIdPrefix$SeqIdNum\t$block_start\t$final_seq_len\t$block_size\t$ctg_start\t$ctg_end\t$ctg_size\t$ctg_name\t$strand\t$output_ctg_len\t$ctg_seq\n";

	}
	
	##output the last gap sequence at the utg end
	if ($Out[-1][1] < length($Tseq))
	{	
		my $gap_start = $Out[-1][1];
		my $gap_end = length($Tseq);
		my $gap_size = $gap_end - $gap_start;
		my $gap_seq = substr($Tseq, $gap_start, $gap_size);

		my $block_start = $final_seq_len + 1;
		my $block_size = ($gap_size > 0) ? $gap_size : 0;
		$final_seq_len += $block_size;

		$pos_output .= "$SeqIdPrefix$SeqIdNum\t$block_start\t$final_seq_len\t$block_size\t$gap_start\t$gap_end\t$gap_size\t$Tname\t+\t".length($gap_seq)."\t".$gap_seq."\n";
		$final_seq .= $gap_seq;
		
	}
	
	my $final_seq_len = length($final_seq);
	Display_seq(\$final_seq, 100);
	print SEQFILE ">$SeqIdPrefix$SeqIdNum     Length: $final_seq_len     Category: Merged_illumina_pacbio      Source: $Tname$ctg_all_names\n$final_seq\n";
	print POSFILE $pos_output;
	
}


##output the unaligned utgs
foreach my $utg_id (keys %UTG){
		if (! exists $AlignedUTG{$utg_id}){
		$Unaligned_utg_num ++;
		$SeqIdNum ++;
		my $utg_seq = $UTG{$utg_id};
		my $utg_len = length($utg_seq);
		
		my $pos_output = "$SeqIdPrefix$SeqIdNum\t1\t$utg_len\t$utg_len\t0\t$utg_len\t$utg_len\t$utg_id\t+\t$utg_len\t$utg_seq\n";
		Display_seq(\$utg_seq, 100);
		print SEQFILE ">$SeqIdPrefix$SeqIdNum    Length: $utg_len    Category: Derived_from_pacbio     Source: $utg_id\n$utg_seq\n";
		print POSFILE $pos_output;

	}
}

##output the unaligned scafftigs
foreach my $scafftig_id (keys %SCAFFTIG){
	
	my $scafftig_seq = $SCAFFTIG{$scafftig_id};
	my $scafftig_len = length($scafftig_seq);

	if (! exists $AlignedScafftig{$scafftig_id} && $scafftig_len >= 250)
	{	$Unaligned_scafftig_num ++;
		$SeqIdNum ++;
		my $pos_output = "$SeqIdPrefix$SeqIdNum\t1\t$scafftig_len\t$scafftig_len\t0\t$scafftig_len\t$scafftig_len\t$scafftig_id\t+\t$scafftig_len\t$scafftig_seq\n";
		Display_seq(\$scafftig_seq, 100);
		print SEQFILE ">$SeqIdPrefix$SeqIdNum     Length: $scafftig_len    Category: Derived_from_illumina     Source: $scafftig_id\n$scafftig_seq\n";
		print POSFILE $pos_output;

	}
}

close SEQFILE;
close POSFILE;




print STDERR "Aligned_utg_num:  $Aligned_utg_num\n";
print STDERR "Unaligned_utg_num: $Unaligned_utg_num\n";
print STDERR "Unaligned_scafftig_num: $Unaligned_scafftig_num\n";


#read fasta file
#usage: Read_fasta($file,\%hash);
#############################################
sub Read_fasta{
	my $file=shift;
	my $hash_p=shift;
	
	my $total_num;
	open(IN, $file) || die ("can not open $file\n");
	$/=">"; <IN>; $/="\n";
	while (<IN>) {
		chomp;
		my $head = $_;
		my $name = $1 if($head =~ /^(\S+)/);
		
		$/=">";
		my $seq = <IN>;
		chomp $seq;
		$seq=~s/\s//g;
		$/="\n";
		
		if (exists $hash_p->{$name}) {
			warn "name $name is not uniq";
		}

		$hash_p->{$name} = $seq;

		$total_num++;
	}
	close(IN);
	
	return $total_num;
}



##complement and reverse the given sequence
#usage: Complement_Reverse(\$seq);
#############################################
sub Complement_Reverse{
	my $seq_p=shift;
	if (ref($seq_p) eq 'SCALAR') { ##deal with sequence
		$$seq_p=~tr/AGCTagct/TCGAtcga/;
		$$seq_p=reverse($$seq_p);  
	}
}
#############################################



#display a sequence in specified number on each line
#usage: disp_seq(\$string,$num_line);
#		disp_seq(\$string);
#############################################
sub Display_seq{
	my $seq_p=shift;
	my $num_line=(@_) ? shift : 50; ##set the number of charcters in each line
	my $disp;

	$$seq_p =~ s/\s//g;
	for (my $i=0; $i<length($$seq_p); $i+=$num_line) {
		$disp .= substr($$seq_p,$i,$num_line)."\n";
	}
	$$seq_p = ($disp) ?  $disp : "\n";
}
#############################################
