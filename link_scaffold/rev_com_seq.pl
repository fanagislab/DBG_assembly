use strict;

my $seq_file = shift;

my %Data;

Read_fasta_simple($seq_file, \%Data);

foreach my $id (sort keys %Data) {
	my $seq = $Data{$id};
	Complement_Reverse(\$seq);
	Display_seq(\$seq);
	print ">".$id."_rc\n".$seq;
}


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


#read fasta file
#usage: Read_fasta($file,\%hash);
#############################################
sub Read_fasta_simple{
	my $file=shift;
	my $hash_p=shift;
	
	my $total_num;
	open(IN, $file) || die ("can not open $file\n");
	$/=">"; <IN>; $/="\n";
	while (<IN>) {
		chomp;
		my $head = $_;
		my $name = $1 if($head =~ /^(\S+)/);
		
		$/="\n>";
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



