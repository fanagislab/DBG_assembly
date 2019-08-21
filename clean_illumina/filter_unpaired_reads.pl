#!/usr/bin/perl	-w
#
#Author: wanghengchao
#Date: 2017-08-05
#the program can deal with sam.noHost.fq.gz form.Input should be reads1 and reads2 files.Output will \
##be .pe.gz form,which means two outfile have reads pired.
use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Data::Dumper;
use File::Path;

my $reads1 = shift;
my $reads2 = shift;
open FILE1,"gzip -dc $reads1 | ";#input file be opened.
open FILE2,"gzip -dc $reads2 | ";

$reads1 = basename($reads1);
$reads2 = basename($reads2);

my $out_reads1 = "$reads1.pe1.gz";
my $out_reads2 = "$reads2.pe2.gz";

open OUT1,"| gzip - > $out_reads1";
open OUT2,"| gzip - > $out_reads2";
#@ST-E00144:238:HNHY3CCXX:3:1101:11383:1344  Filtered_reads
##@ST-E00144:238:HNHY3CCXX:3:1101:11383:1344  Filtered_reads

my $both_end_number = 0;
my $single_end_number = 0;
my $both_empty_number = 0;
my $read1_base_number = 0;
my $read2_base_number = 0;

my ($line1,$line2);
while(defined($line1=<FILE1>) && defined($line2=<FILE2>)) {

	if(($line1 =~ /^@/) && ($line2 =~ /^@/)) {
		my $read1_name = $1 if($line1 =~ /^(\S+)/);
		my $base1 = <FILE1>;#MAY NEED CHOMP
		my $add_char = <FILE1>;
		my $base_q1 = <FILE1>;

		my $read2_name = $1 if($line2 =~ /^(\S+)/);
		my $base2 = <FILE2>;
		$add_char = <FILE2>;
		my $base_q2 = <FILE2>;
			
		chomp($base1);
		chomp($base2);
		if ($read1_name ne $read2_name){
			print STDERR "Error: Read1_name does not equal Read2_name\n";
			print STDERR "$read1_name\n$base1\n$add_char$base_q1";
			print STDERR "$read2_name\n$base2\n$add_char$base_q2";
			exit();
		}

		if($base1 && $base2) {
			print OUT1 "$read1_name\n$base1\n$add_char$base_q1";
			print OUT2 "$read2_name\n$base2\n$add_char$base_q2";
			$both_end_number ++;
			$read1_base_number += length($base1);
			$read2_base_number += length($base2);
		}else
		{	if($base1 || $base2){
				$single_end_number ++;
			}else{
				$both_empty_number ++;
			}		
		}		
	}
}
close FILE1;
close FILE2;
close OUT1;
close OUT2;

my $total_end_number = $both_end_number + $single_end_number + $both_empty_number;

print STDERR "Both end with reads: ".$both_end_number."\t".($both_end_number/$total_end_number)."\n";
print STDERR "Single end with reads: ".$single_end_number."\t".($single_end_number/$total_end_number)."\n";
print STDERR "Both end without reads: ".$both_empty_number."\t".($both_empty_number/$total_end_number)."\n";

print STDERR "\n\nboth_end_number / (both_end_number + single_end_number) : " . ($both_end_number / ($both_end_number + $single_end_number)) . "\n\n";

print STDERR "Read1 output: $both_end_number  $read1_base_number\n";
print STDERR "Read2 output: $both_end_number  $read2_base_number\n";


