#!/usr/bin/perl -w

$filename = $ARGV[0];
$len=0;

open(InputFile, $filename) ||die "error: can't open $filename.\n";
while(<InputFile>) {
    next if ($_ eq "\n");
    chomp;
    if ($_ =~ ">"){
	print "$len\n" if ($len);
	$len=0;
	if ($' =~ /([\S]+)\s(.+)/){ $name = $1; }
	else { $name = $'; }
	print "$name\t";
	next;
    } else {
	chomp;
	$len += length($_);
    }
}

close (InputFile);

print "$len\n";
