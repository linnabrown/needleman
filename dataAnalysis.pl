use Bio::SeqIO;
use Bio::Seq;
use SeqAnalysis;


$filein ="GeneSeq.fa"; #$ARGV[1];
open f,"<$filein";
$tt = 0;
while($line=<f>){
  if(substr($line,0,1) eq '>'){

    chomp $line;
    $name[$tt++] = $line;
  }
}

$tt= 0;
my $file = "GeneSeq.fa";

my $seqio_object = Bio::SeqIO->new(-file => $file);
while(my $seq_object = $seqio_object->next_seq){
  $proname = $name[$tt++];

  $seq = $seq_object->seq;

  my $object2=SeqAnalysis->new(-seq_name=>$proname,
                 -sequence=>$seq,
                   -desc=>"");
 print $proname."\n";
 $object2->nucleotideCounter();
 $object2->detectEnzyme();
 $object2->gcContentSeqLength();
 $object2->detectPolyaSignal();
 print "\n";

}
