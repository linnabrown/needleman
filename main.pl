#!/usr/bin/env perl
# use strict;
# use warnings;
use Bio::SeqIO;
use Bio::Seq;
use SeqAnalysis;
use Align;
# $gp = -1;
# $ma = 1;
# $mm = -1;
# $s1 = qw/ATGCATGCT/;
# $s2 = qw/ATGACAGCT/;
# $filename = GeneSeq.fa


# $file = "GeneSeq.fa";
$file = $ARGV[0];
$target = $ARGV[1];
$ma = $ARGV[2];
$mm = $ARGV[3];
$gp = $ARGV[4];
$e = $ARGV[5];
$mode = $ARGV[6];


$filein ="GeneSeq.fa"; #$ARGV[1];
open f,"<$filein";
$tt = 0;
while($line=<f>){
  if(substr($line,0,1) eq '>'){

    chomp $line;
    $name[$tt++] = $line;
  }
}



# $target = "ATGCATGCTT";

# $ma = 1;
# $mm = -1;
# $gp = -1;
# $e = 2;
# $mode = "all"; #all or one

$maxlen = length($target)+$e;
$minlen = length($target)-$e;

$tt=0;
my $file = "GeneSeq.fa";
my $seqio_object = Bio::SeqIO->new(-file => $file);


while(my $seq_object = $seqio_object->next_seq){
    $proname = $name[$tt++];
    $seq = $seq_object->seq;
    $seq_len = length($seq);
    print $proname."\n";
    # print $mm."\n";
    my $ao=Align->new(-seq1=>$seq,-seq2=>$target,
                                      -match=>$ma,
                                      -mismatch=>$mm,
                                      -gap=>$gp,
                                      -min_map_len=>$minlen,
                                      -max_error=>$e,
                                      -mode=>$mode # mode=all or one
                                      );
    $ao->printSeqWithSpacer();
    $ao->getAlignment();


    # last;

}#while
$s1 = qw/ATGCATGCT/;
$s2 = qw/ATGACAGCT/;
