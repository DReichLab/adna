#!/usr/bin/env perl

use strict;
use warnings;

die "Usage: gen-bc.pl <barcode.txt> <sam.hdr>\n" if @ARGV < 2;

my %h;

open(FH, $ARGV[0]) || die;
while (<FH>) {
	chomp;
	my @t = split("\t");
	$h{$t[0]} = [$t[1], $t[2]];
}
close(FH);

my $sample;
open(FH, $ARGV[1]) || die;
while (<FH>) {
	$sample = $1 if (/^\@RG.*\tSM:([^\t]+)/);
}
close(FH);

if (defined $h{$sample}) {
	print join("\n", @{$h{$sample}}), "\n";
} else {
	warn("ERROR: Failed to find sample $sample in file $ARGV[1]\n");
}
