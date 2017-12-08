sub needleman($$$$$) {
  # print "ok";
  my $seq1 = $_[0];
  my $seq2 = $_[1];
  my $gap_penalty = $_[2];
  my $ma = $_[3];
  my $mism  = $_[4];
  # foreach(@_){
  #   print $_."\n";
  # }
  if ( not $seq1 or not $seq2 ){
     die "Must present two strings on command line";
  }
  my $scoreMatrix;

  my @lettersA = split //, $seq1;
  my @lettersB = split //, $seq2;
  for ( my $i = 0 ; $i < length($seq1) + 1 ; $i++ ) {
     $scoreMatrix->[$i][0] = $i * $gap_penalty;
  }

  for (my $j=0; $j < length($seq2) + 1; $j++) {
  	$scoreMatrix->[0][$j] = $j * $gap_penalty;
  }

  my $tracebackMatrix;
  $tracebackMatrix->[0][0] = 'Diagonal';

  for ( my $i = 1 ; $i < length($seq1) + 1 ; $i++ ) {
     for ( my $j = 1 ; $j < length($seq2) + 1 ; $j++ ) {
        no warnings;

        # Diagonal
        my $temp   = score( $lettersA[ $i - 1 ], $lettersB[ $j - 1 ],$ma, $mism);

        my $match  =  $scoreMatrix->[ $i - 1 ][ $j - 1 ] + $temp;
        # up
        my $delete = $scoreMatrix->[ $i - 1 ][$j] + $gap_penalty;

        # left
        my $insert = $scoreMatrix->[$i][ $j - 1 ] + $gap_penalty;

        my $max = max( $match, $delete, $insert );

        my $whatMatched = '';
        for ($max) {
           $whatMatched = 'Left'     if ( $_ eq $insert );
           $whatMatched = 'Up'       if ( $_ eq $delete );
           $whatMatched = 'Diagonal' if ( $_ eq $match );
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
  while ( $i > 0 or $j > 0 ) {
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
     else {
        $alignA = '-' . $alignA;
        $alignB = $lettersB[ $j - 1 ] . $alignB;
        $j--;
        $line = " ".$line;
        $operator ="D".$operator;
        $num_D++;
     }
  }

  # if( $num_I + $num_D+$num_MM >2){
  #   @t=();
  #   return @t;
  # }
  my $n = length($seq1);
  my $m = length($seq2);

  $mapped_length = length($alignA);

  $tmpstr = "[Mapped-length=$mapped_length m=$num_M I=$num_I D=$num_D M=$num_MM]";

  @arr = ($scoreMatrix->[$n][$m],$alignA." ".$tmpstr,$line,$alignB,$operator);

  # foreach(@arr){
  #   print $_."\n";
  # }
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
