#OUTPUT: original genome and reverse complemented genome where $chr is labelled as $chr_rc
#        bowtie index
#! usr/bin/perl
use strict;
use warnings;
require hiseq;

my $config = $ARGV[0];
my %cf;
open(CF, $config) || die "$config\n";
while(<CF>){
	if($_ =~ /(\S+)\t(\S+)/){
		$cf{$1} = $2;
	}
}
my $genome = $cf{'genome'};
my ($chr, $t) = ();

system "cp $genome.fa $genome"."WC.fa";
open(G, "$genome.fa") || die "$genome.fa";
$genome .= 'WC';
open(RC, ">>$genome.fa") || die;
$chr = '';
while(<G>) {
    if(/^(\>[^\s]+)/) {
	$t = $1;
	$chr && print RC &hiseq::revcomDNA($chr)."\n";
	$chr = '';
	print RC $t."_rc\n";
    } elsif(/^([ACGTNRYMKSWHBVD]+)\s?$/i) {
	$chr .= $1;
    } else {
	/^\s*$/ || die "$_*";
    }
}
$chr ? print RC &hiseq::revcomDNA($chr)."\n" : die;
system "bowtie2-build ./$genome".".fa ./$genome";
