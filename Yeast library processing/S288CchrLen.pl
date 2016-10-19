#Combine location of each read in genome. The location of a read aligned to Crick strand is represented by -(distance from right end of chr - 1)
#OUTPUT: Column 1: read name
#        Column 2: alignment report: field6_field19 of SAM file
#        Column 3: genome location: chr_location
#! usr/bin/perl
use strict;
use warnings;
require hiseq;

my ($chr, $len, %chrLen) = ();
my $config = $ARGV[0];
my %cf;
open(CF, $config) || die "$config\n";
while(<CF>){
	if($_ =~ /(\S+)\t(\S+)/){
		$cf{$1} = $2;
	}
}
my $fa = $cf{"genome"};
open(I, "$fa.fa") || die "$fa.fa";
open(O, ">$fa.Len.txt") || die;
while(<I>) {
    if(/^\>([\w\d\.]+)/) {
	if($chr) {
	    $len || die;
	    exists($chrLen{$chr}) ? die : ($chrLen{$chr} = $len);
	    print O "$chr\t$chrLen{$chr}\n";
	}

	$chr = $1;
	$len = 0;
    } elsif(/^([ACTGNRYMKSWHDBVactgnrymkswhdbv]+)\s*$/) {
	$len += length($1);
    } else {
	die "$_\n";
    }
}
$chr || die;
$len || die;
exists($chrLen{$chr}) ? die : ($chrLen{$chr} = $len);
print O "$chr\t$chrLen{$chr}\n";
