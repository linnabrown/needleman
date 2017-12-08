#!/usr/bin/env perl
#pls edit parameter in this file
system("perl main.pl GeneSeq.fa ATGCATGCTT 1 -1 -1 2 all > alignment.txt");
system("perl dataAnalysis.pl > dataAnalysis_geneseq.txt");
