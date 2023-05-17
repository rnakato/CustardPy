#!/usr/bin/env perl

=head1 DESCRIPTION

    Parse stats outputted by Juicer.

=head1 SYNOPSIS

    parseJuicerstats.pl <file>

=cut

use strict;
use warnings;
use autodie;
use Path::Class;
use Getopt::Long;
use Pod::Usage qw/pod2usage/;
my $header=0;
GetOptions('header' => \$header);

my $filename=shift;
pod2usage unless $filename;
my $file = file($filename);

my $num_total="";
my $num_paired = "";
my $num_chimpaired = "";
my $num_chimambi = "";
my $num_unaligned = "";
my $num_ligationmotif = "";

my $fh = $file->open('r') or die $!;
while(<$fh>){
    next if($_ eq "\n");
    chomp;

    if($_ =~ /Experiment description:/) {
        if ($header) {
            print "Experiment description";
        }else{
            my @clm = split(/;/, $_);
            print "$clm[1],$clm[2]";
        }
    }elsif($_ =~ /(.+): (.+)/) {
        if ($header) {
            print "\t$1";
        }else{
            print "\t$2";
        }
    }

}
$fh->close;

print "\n";
