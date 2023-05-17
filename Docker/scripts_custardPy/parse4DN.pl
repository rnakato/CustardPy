#!/usr/bin/env perl

=head1 DESCRIPTION

    Parse stats outputted by 4DN pipeline

=head1 SYNOPSIS

    parse4DN.pl <file>

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

    if($_ =~ /(.+)\t(.+)/) {
        if ($header) {
            print "$1\t";
        }else{
            print "$2\t";
        }
    }

}
$fh->close;

print "\n";
