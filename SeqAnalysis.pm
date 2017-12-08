#!_usr_bin/perl
# use strict;
package SeqAnalysis;

my $myseq;
my $name;
my $desc;

my %dinucleotide_hash = (
	"AA"=>0,	"AT"=>0,	"AG"=>0,	"AC"=>0,
	"TA"=>0,	"TT"=>0,	"TG"=>0,	"TC"=>0,
	"GA"=>0,	"GT"=>0,	"GG"=>0,	"GC"=>0,
	"CA"=>0,	"CT"=>0,	"CG"=>0,	"CC"=>0
);

my %codon_distribution = ();
my @dinucleotide = ('A','T','G','C');
foreach my $first (@dinucleotide){
	foreach my $second (@dinucleotide){
		foreach my $third (@dinucleotide){
			$codon_distribution{$first.$second.$third} = 0;
		}
	}
}
sub new
{
    my $class = shift;
	my $omit;
	($omit,$name,$omit,$myseq,$omit,$desc) = @_;
	my $self = {};
    bless $self, $class;
	# print join(', ',$myseq,$name,$desc);
    return $self;
}

sub printWithSpacer {
	# Put your core code for Assignment No.1 here, which will print a scale, label and spacer.
	# You can use print statements here to print nucleotides with spacers.
	$myseq =~ s/\s//g;
	print $desc;
	print '              1          2          3          4          5          6          7          8          9         10'."\n";
	print 'line 1234567890 1234567890 1234567890 1234567890 1234567890 1234567890 1234567890 1234567890 1234567890 1234567890'."\n";
	my $line_number=0;
	for (my $i=0; $i<length($myseq); $i=$i+100){
		my $seq100=substr($myseq,$i,100);
		my $space=' ';
		if ($line_number<9){                             #for line1-line9
			$line_number=$line_number+1;
			print $space.$space.$space.$line_number.$space;
		}
		else{                                            #for line10-line99
			$line_number=$line_number+1;
			print $space.$space.$line_number.$space;
		}
		   for (my $j=0; $j<length($seq100); $j=$j+10) {
				my $seq10=substr($seq100,$j,10);
				print $seq10.$space;
				}
		print "\n";
	}
}


sub codonUsage { # For CS students only
	# my %codon_distribution_ref = \%codon_distribution;
	# my ($myseq,$codon_distribution_ref)=shift;
	# You need to examine all codons existing in the given sequence
	# You also need to consider the reverse complementary sequence
	# The ratio for a particular codon will be its number divided by the total codon
	# In main program, you need to define %codon_distribution, which can be updated here.
	# This function is for computer science students only. The relevant print statements in main
	# program will need to generate the output similar to [6] shown in the previous page
	print "[7] Codon Usage:\n";
	$myseq =~ s/\s//g;
	@reverse_array = $myseq =~ /\w{1}/g;
	for(my $index =0 ;$index < length($myseq) - 2; $index ++ ){
		$codon_distribution{substr($myseq,$index,3)} ++;
	}
	# reverse the gene code
	$myseq = join "" ,reverse(@reverse_array);

	for(my $index =0 ;$index < length($myseq) - 2; $index ++ ){
		$codon_distribution{substr($myseq,$index,3)} ++;
	}
	foreach my $var(keys %codon_distribution){
		$codon_distribution{$var}  = $codon_distribution{$var} / length($myseq) / 2;
	}
	my $counts = 0;
	foreach (sort keys %codon_distribution){
		printf "[%s] = %.3f\t",$_,$codon_distribution{$_};
		$counts ++;
		if($counts % 4 == 0){
			print "\n";
		}
	}
	return \%codon_distribution;

}

sub detectPolyaSignal { # CS and BIO students have different tasks
	# my $myseq=shift;
	# The print statements for the detected poly(A) signals can be here in this subroutine.
	# Therefore, you should not have a return statement here.
	# For biology students, you can use exact match.
	# For computer science students, both exact and fuzzy matches (i.e., one base difference, but not
	# in first, third and last positions of [A]A[T]AA[A]) should be allowed.
	print "[6] Detection of poly(A) signal (AATAAA):\n";
	$myseq =~ s/\s//g;
	my $no = 0;
	print "No.\tStart\tEnd\tSignal\n";
	while ($myseq =~ m/A[ACGT]T[ACGT]{2,2}A/g){
		my $len = length($&);
		my $end = pos($myseq);
		my $start = $end - $len + 1;
		$no ++;
		print "$no\t$start\t\t$end\t$&\n";
	}
}

sub dinucleotideFrequency {
	# my ($myseq,%dinucleotide_hash_ref)=shift;
	# You need to examine all dinucleotides existing in the given sequence
	# You do not need to consider the reverse complementary sequence
	# In your main program, you need to set up a hash for 16 dinucleotide with 0 as initial value.
	# Then, you need to have print statement in main program to print this hash to get output similar
	# to what is shown as [5] in the previous page.
	print "[5] Dinucleotide Frequency (%):\n";
	$myseq =~ s/\s//g;
	for(my $index =0 ;$index < length($myseq) - 1; $index ++ ){
		$dinucleotide_hash{substr($myseq,$index,2)}++;
	}
	foreach my $var(keys %dinucleotide_hash){
		$dinucleotide_hash{$var} = $dinucleotide_hash{$var} / length($myseq);
	}

	my $counts = 0;
	foreach (sort keys %dinucleotide_hash){
		printf "[%s] = %.2f\t",$_,$dinucleotide_hash{$_};
		$counts ++;
		if($counts % 4 == 0){
			print "\n";
		}
	}

	return \%dinucleotide_hash;
}

sub detectEnzyme {

	$myseq =~ s/\s//g;
	print "[4] Restriction Sites:\n";
	print "Name\t Pos\t Seq\t\t IUPAC\t ALT\n";
	$myseq =~ m/ACAAGGG/g;
	my $pos = pos($myseq);
	if(defined($pos)){
		print "BclI\t $pos\t $&\t ACAAGGG\t ACAAGGG","\n";
	}else{
		print "BclI\t -1\t  \t ACAAGGG\t ACAAGGG","\n";
	}
	$myseq =~ m/GG[AG]GCA[CT]T/g;
	$pos = pos($myseq);
	if(defined($pos)){
		print "BfmI\t $pos\t $&\t GGRGCAYT\t GG[AG]GCA[CT]T\n";
	}else{
		print "BfmI\t -1\t  \t GGRGCAYT\t GG[AG]GCA[CT]T\n";
	}
	$myseq =~ m/AGG[ACGT]TTTA/g;
	$pos = pos($myseq);
	if(defined($pos)){
		print "EcoRI\t $pos\t $&\t AGGNTTTA\t AGG[ACGT]TTTA\n";
	}else{
		print "EcoRI\t -1\t  \t AGGNTTTA\t AGG[ACGT]TTTA\n";
	}
	$myseq =~ m/TGGCC[CT]/g;
	$pos = pos($myseq);
	if(defined($pos)){
		print "Cac8I\t $pos\t $&\t TGGCCY\t\t TGGCC[CT]\n";
	}else{
		print "Cac8I\t -1\t     \t\t TGGCCY\t\t TGGCC[CT]\n";
	}

}

sub gcContentSeqLength {
	# my $myseq=shift;
	my ($GCcontent,$SeqLength)=(0,0);
	# You need to invoke nucleotideCounter() get relevant nucleotide account
	# Then, you will calculaate total sequence length and GC content value
	$myseq =~ s/\s//g;
	$SeqLength = length($myseq);
	while ($myseq =~ m/GC/g){
		$GCcontent ++;
	}
	printf "[2] GC Content: %.6f\n",$GCcontent / $SeqLength;
	print "[3] Sequence length: $SeqLength\n";
	return ($GCcontent,$SeqLength);
}

sub nucleotideCounter {
	# my $myseq=shift;
	my ($a,$t,$g,$c,$n)=(0,0,0,0,0);
    for(my $i=0;$i<length($myseq);$i++){
    	my $nucleotide = substr($myseq,$i,1);
    	if($nucleotide eq 'A'){
    		$a=$a+1;
    	}
    	elsif($nucleotide eq 'T'){
    		$t=$t+1;
    	}
    	elsif($nucleotide eq 'G'){
    		$g=$g+1;
    	}
    	elsif($nucleotide eq 'C'){
    		$c=$c+1;
    	}
    	else{
    		$n=$n+1;
    	}
    }
	print "[1] Nucleotide Counts: A=$a T=$a G=$g C=$c Other=$n\n";
}

sub detectMotif{
	shift;
	my (@seq) = @_;
	print "[8] Detection of motifs without labels:\n";
	print "(",join(", ",@seq),")\n";

	print  "Motif\tStart\tEnd\n";
	foreach(@seq){
		$myseq =~ m/$_/g;
		my $pos = pos($myseq);
		my $len = length($&);
		my $start = $pos - $len;
		$pos = defined($pos)? $pos:0;
		$start = defined($start)? $start:0;
		print "$_\t$start\t\t$pos\n";
	}

}
sub detectMotifsWithLables{
	shift;
	my (%seq) = @_;

	print "[9] Detection of motifs with labels: ";
	print "(",join(", ",@_),")\n\n";

	print  "Label\tMotif\tStart\tEnd\n";
	foreach(sort keys %seq){
		my $value = $seq{$_};
		# print $value;
		$myseq =~ m/$value/g;
		my $pos = pos($myseq);
		my $len = length($&);
		my $start = $pos - $len;
		$pos = defined($pos)? $pos:0;
		$start = defined($start)? $start:0;
		print "$_\t$value\t$start\t\t$pos\n";

	}
	# for(my $var = 0;$var * 2 < scalar(@seq); $var ++){
		# $myseq =~ m/$seq[$var+1]/g;
		# my $pos = pos($myseq);
		# my $len = length($&);
		# my $start = $pos - $len;
		# $pos = defined($pos)? $pos:0;
		# $start = defined($start)? $start:0;
		# print "$seq[$var]\t$seq[$var+1]\t$start\t\t$pos\n";
	# }

}


1;
