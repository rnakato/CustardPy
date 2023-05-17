#!/usr/bin/env perl

=head1 DESCRIPTION

    split multifasta into single fastas in <dir>

=head1 SYNOPSIS

    % splitmulitfasta.pl <input.fa> [--dir <dir>] [--underbar|-u]

=cut

use strict;
use warnings;
use autodie;
use Path::Class;
use Getopt::Long qw/:config no_ignore_case bundling auto_help/;
use Pod::Usage qw/pod2usage/;

my $underbar=0;
my $dir="./";
GetOptions(
    "dir|d=s" => \$dir,
    "underbar|u" => \$underbar
    ) or pod2usage(1);

my $filename = $ARGV[0];
pod2usage(2) unless $filename;
my $head ="";
my $outfile = "";
my $seq = "";
my $file = file($filename);
my $fh = $file->open('r') or die $!;
while(<$fh>){
    next if($_ eq "\n");
    chomp;
    if($_ =~ />(.+)/){
	if($seq ne ""){
	    my $out = file($dir . "/" . $outfile . ".fa");
	    my $writer = $out->open('w') or die "Can't read $out: $!";
	    $writer->print(">$head\n");
	    $writer->print("$seq\n");
	    $writer->close;
	    $seq = "";
	}
	$head =$1;
	if($underbar) {
	    $outfile = $head;
	    $outfile =~ s/\s+/_/g;
	}
	else {
	    my @clm = split(/ /, $head);    
	    $outfile = $clm[0];
	}
#	print "$head\n"
    }else{
	$seq .= $_;
    }
} 
$fh->close;

if($seq ne ""){
    my $out = file($dir . "/" . $outfile . ".fa");
    my $writer = $out->open('w') or die "Can't read $out: $!";
    $writer->print(">$head\n");
    $writer->print("$seq\n");
    $writer->close;
}
