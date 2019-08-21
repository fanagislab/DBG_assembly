my $lib_file = shift;

my $number = 1;
open IN, $lib_file || die "fail";
while(<IN>){
	open OUT, ">$lib_file.$number";
	$number ++;
	print OUT $_;
	close OUT;
}
close IN;

