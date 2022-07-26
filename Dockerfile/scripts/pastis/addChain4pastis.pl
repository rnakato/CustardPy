#!/usr/bin/env perl
use strict;
use warnings;

my $filename=$ARGV[0];
my $output=$filename . ".addCONECT.pdb";

open(PDB, $filename) || die "error: can't open $filename\n";
my @array=();
my @linearray=();
while(<PDB>) {
    next if($_ eq "\n" || $_ =~ /nan/);
    push (@linearray, $_);
    
    chomp;
    my @clm = split(/\s+/, $_);
    push (@array, $clm[1]);
#    print "$clm[1]\n";
#    print "$clm[5]\t$clm[6]\t$clm[7]\n";
}
close(PDB);

#my $num = `wc -l $filename | cut -f1 -d' '`;
#system("cp $filename $output");

open(FILE, ">", $output) or die;
foreach my $line (@linearray){
  print FILE $line;
}

my $num=-100;
foreach my $id (@array){
    if($id-$num==1) {
	printf FILE "CONECT%5d%5d\n", $num, $id;
    }
    $num = $id;
}
close(FILE)
