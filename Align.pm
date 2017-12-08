#!_usr_bin/perl
# use strict;
package Align;
my $myseq;
my $seq2;
my $match;
my $mismatch;
my $gap;
my $min_map_len;
my $max_error;
my $mode;

sub new
{
  my $class = shift;
	my $omit;
	($omit,$myseq,$omit,$seq2,$omit,$match,$omit,$mismatch,$omit,$gap,$omit,
      $min_map_len,$omit,$max_error,$omit,$mode) = @_;
	my $self = {};
    bless $self, $class;
    return $self;
}

sub printSeqWithSpacer {
	# Put your core code for Assignment No.1 here, which will print a scale, label and spacer.
	# You can use print statements here to print nucleotides with spacers.
	$myseq =~ s/\s//g;
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

sub getAlignment(){
  # print "\n----------------------------------------------------------\n";
  # print $myseq;
  # print $seq2;
  # print $match;
  # print $mismatch;
  # print $gap;
  # print $min_map_len;
  # print $max_error;
  # print $mode;
  # print "\n-----------------------------------------------------------\n";
  $maxlen = length($seq2)+$max_error;
  $minlen = length($seq2)-$max_error;

  print "[Scoring schema]: match=$match, mismatch=$mismatch, gap=$gap\n";
  print "[Search Target]: $seq2\n";
  print "[Maximum target length]: $maxlen\n";
  print "[Minimum target length]: $minlen\n";
  print "[Maximum allowed error bases]: $max_error\n";


  my $e = $max_error;
  my @target_arr=();
  $target_arr[0] = $seq2;
  $target_arr[1] = substr($seq2,0,9);
  $target_arr[2] = substr($seq2,0,8);
  $target_arr[3] = substr($seq2,1,8);
  $target_arr[4] = substr($seq2,1,9);
  $target_arr[5] = substr($seq2,2,8);

  # print $max_error."\n";
  $myseq =~ s/\s//g;
  my $seq = $myseq;
  # print $seq;
  my @dic = ();
  my $seq_len = length($seq);
  my $max_score = -100000;
  my $pre_end = -1;
  for($mis1 = -$e;$mis1<=$e; $mis1++){
  # $mis1 = 0;
    for(my $i =0; $i< $seq_len-9+$mis1;$i++){
        $s1 = substr($seq,$i,10-$mis1);
        my $start_index = $i+1;
        my $end_index = $start_index + (10 -$mis1);
        for(my $j=0; $j<@target_arr; $j++){
        $s2 = $target_arr[$j];
        @value = needleman($s1,$s2,$gap,$match,$mismatch,$start_index,$end_index);
        if(@value){
          if($value[0]>$max_score && $pre_end != $end_index){
            $max_score = $value[0];
            # shift @value;
            unshift @dic,@value;
            $pre_end = $end_index;
          }
        }else{
          next;
        }
        }#for
      }#for
  }#for
  my $len1 = @dic/5;

  my $cnt1 = 0;
  my $num_highest_score = 0;
  foreach(@dic){
    $cnt1++;
    if($cnt1 % 5 == 1){
      if($_!=$max_score){
        last;
      }
      else{
        $num_highest_score++;
      }
    }
   }

  print "[The highest alignment score]: $max_score\n";
  print "[The alignments with the highest score]: $num_highest_score\n";
  # print $mode."\n";
  # for($k=0;$k<@dic;$k++){
  #   @val = $dic[$k];
  #   foreach(@val){
  #     print "$_\n";
  #   }
  # }
  my $cnt = 0;
  foreach(@dic){
    $cnt++;
    if($cnt % 5 == 1){
      if($_!=$max_score){
        last;
      }
      else{
        next;
      }
    }
    print "$_\n";
    if($mode=="one" && $cnt==5){
      last;
    }
   }


  # print "the best score is $max_score\t$len1\n";




}
sub needleman($$$$$$$) {
  # print "ok";
  my $seq1 = $_[0];
  my $seq2 = $_[1];
  my $gap_penalty = $_[2];
  my $ma = $_[3];
  my $mism  = $_[4];
  my $ss = $_[5];
  my $ee = $_[6];

  if ( not $seq1 or not $seq2 ){
     die "Must present two strings on command line";
  }
  my $scoreMatrix;

  my @lettersA = split //, $seq1;
  my @lettersB = split //, $seq2;

  for ( my $i = 0 ; $i < length($seq1) + 1 ; $i++ ) {
     $scoreMatrix->[$i][$j] = 0;
  }

  for ( my $i = 0 ; $i < length($seq1) ; $i++ ) {
     $scoreMatrix->[$i+1][0] = ($i+1) * $gap_penalty;
  }

  for (my $j=0; $j < length($seq2); $j++) {
  	$scoreMatrix->[0][$j+1] = ($j+1) * $gap_penalty;
  }

  my $tracebackMatrix;
  $tracebackMatrix->[0][0] = 'Diagonal';

  for ( my $i = 1 ; $i < length($seq1) + 1 ; $i++ ) {
     for ( my $j = 1 ; $j < length($seq2) + 1 ; $j++ ) {
        no warnings;

        # Diagonal
        # my $temp   = score( $lettersA[ $i - 1 ], $lettersB[ $j - 1 ],$ma, $mism);
        my $match_mismatch = $lettersA[$i-1] eq $lettersB[$j-1]?$ma:$mism;
        my $match  =  $scoreMatrix->[ $i - 1 ][ $j - 1 ] + $match_mismatch;
        # up
        my $delete = $scoreMatrix->[ $i - 1 ][$j] + $gap_penalty;

        # left
        my $insert = $scoreMatrix->[$i][ $j - 1 ] + $gap_penalty;

        my $max = max( $match, $delete, $insert );

        my $whatMatched = 'Diagonal';
        for ($max) {
           $whatMatched = 'Left'     if ( $_ eq $insert );
           $whatMatched = 'Up'       if ( $_ eq $delete );
          #  $whatMatched = 'Diagonal' if ( $_ eq $match );
        }
        $scoreMatrix->[$i][$j]     = $max;
        $tracebackMatrix->[$i][$j] = $whatMatched;

     }
  }


  my $alignA = '';
  my $alignB = '';
  my $line = "";
  my $operator = "";
  my $i = length($seq1);
  my $j = length($seq2);
  my $num_I = 0;
  my $num_D = 0;
  my $num_M = 0;
  my $num_MM = 0;
  while ( $i != 0 && $j != 0 ) {
     no warnings;
     if ( $i > 0 && $j > 0 && $tracebackMatrix->[$i][$j] eq 'Diagonal' ) {
        $alignA = $lettersA[ $i - 1 ] . $alignA;
        $alignB = $lettersB[ $j - 1 ] . $alignB;
        if($lettersA[ $i - 1 ] eq $lettersB[ $j - 1 ]){
          $line = "|".$line;
          $operator ="m".$operator;
          $num_M++;
        }else{
          $line = " ".$line;
          $operator = "M".$operator;
          $num_MM++;
        }
        $i--;
        $j--;
     }
     elsif ( $i > 0 && $tracebackMatrix->[$i][$j] eq 'Up' ) {
        $alignA = $lettersA[ $i - 1 ] . $alignA;
        $alignB = '-' . $alignB;
        $i--;
        $line = " ".$line;
        $operator ="I".$operator;
        $num_I++;
     }
     elsif ($j > 0 && $tracebackMatrix->[$i][$j] eq 'Left') {
        $alignA = '-' . $alignA;
        $alignB = $lettersB[ $j - 1 ] . $alignB;
        $j--;
        $line = " ".$line;
        $operator ="D".$operator;
        $num_D++;
     }
  }

  if( $num_I + $num_D+$num_MM >2){
    @t=();
    return @t;
  }
  my $n = length($seq1);
  my $m = length($seq2);

  my $mapped_length = length($alignA);

  my $tmpstr = "[Mapped-length=$mapped_length m=$num_M I=$num_I D=$num_D M=$num_MM]";

  my $line1 = $ss." ".$alignA." ".$ee." ".$tmpstr;
  my $num_of_space = length($ss." ");
  for(my $k=0; $k<$num_of_space; $k++){
    $line = " ".$line;
    $alignB = " ".$alignB;
    $operator = " ".$operator;
  }
  my $myspace = $num_of_space * " ";

  @arr = ($scoreMatrix->[$n][$m],$line1,$line,$alignB,$operator);

  return @arr;


}
sub score {
   my ( $letterA, $letterB,$ma,$mism ) = @_;

   return $mism if ( scalar @_ < 2 );           #Only received one letter
   return $ma  if ( $letterA eq $letterB );    #Match
   return $mism if ( $letterA ne $letterB );
     # Mismatch
}

sub max {
   my ( $max, @vars ) = @_;
   for (@vars) {
      $max = $_ if ( $_ > $max );
   }
   return $max;
}

1;
