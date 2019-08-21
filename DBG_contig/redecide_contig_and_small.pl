use strict;

my $contig_file = shift;
my $small_file = shift;
my $len_cutoff = shift || 100;


my %hash;


Read_fasta($contig_file, \%hash);

Read_fasta($small_file, \%hash);


my $ctg_id = -1;
my $sma_id = 0;


open CTG, ">$contig_file.len$len_cutoff.fa" || die "fail";
open SMA, ">$small_file.len$len_cutoff.fa" || die "fail";

foreach  my $id (sort keys %hash){
	my $p = $hash{$id};
	if ($p->{"len"} > $len_cutoff){
		$ctg_id += 2;
		my $head = $p->{"head"};
		$head =~ s/^\S+/$ctg_id/;
		my $seq = $p->{"seq"};
		print CTG ">ctg_".$head."\n".$seq."\n";
		
	}else{
		$sma_id ++;
		my $head = $p->{"head"};
		$head =~ s/^\S+/$sma_id/;
		my $seq = $p->{"seq"};
		print SMA ">small_".$head."\n".$seq."\n";
	}

}


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

		$hash_p->{$name}{"head"} =  $head;
		$hash_p->{$name}{"len"} = length($seq);
		$hash_p->{$name}{"seq"} = $seq;

		$total_num++;
	}
	close(IN);
	
	return $total_num;
}


