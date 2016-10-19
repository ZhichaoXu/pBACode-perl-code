#remove vector sequences at both ends of merged read pairs
#Assuming HindIII site is (1) perfect or (2) HindIII has 1bp mismatch but barcode is exactly 20 bp
#1. 10 bp flanking and between two barcodes have no more than 1bp mismatch to vector sequence
#2. If 10 bp flanking two barcodes have more than 1bp mismatch to vector sequence, check 10bp covering BamHI or HindIII site. If these 10 bp regions perfectly match vector sequence, take 20 bp flanking these regionsand bridges as barcodes.
#2. If 10 bp vector sequence flanking two barcodes have no more than 1bp mismatch, but sequence between two barcodes have more than 1bp mismatch, as long as there is perfect HindIII site, take 20 bp juxtaposed to flanking vector sequence as barcodes
#Priority: at first, consider mismatch type, 1 substitution > 1 deletion > 1 insertion. Then consider location, pick the sequence with 1 mismatch closest to barcode
#! usr/bin/perl
use strict;
use warnings;
use IO::Handle;
require hiseq;

my $config = $ARGV[0];
my $lr = $ARGV[1];
my %cf;
open(CF, $config) || die "$config\n";
while(<CF>){
	if($_ =~ /(\S+)\t(\S+)/){
		$cf{$1} = $2;
	}
}
#my $fq = $cf{$lr."postpool_reads"};
my $fq = $ARGV[2]; 
my (@barcode, @barcode1bpMid, @convert, @lines, @readNumShort, @readNumMid, @readNumLong, @t, %barcode, %count, %doubt, %mismatchL, %mismatchR, %t, %tl, %tr, %ttl, %ttr, $barcode, $bh, $bl, $br, $i, $j, $k, $l1, $l2, $num, $readNum, $rmL, $rmR, $t) = ();
my ($BARCODELEN, $BARCODELENMIN, $BARCODELENMAX, $DS, $DSinternal) = ($cf{"barcode_length_standard"}, $cf{"barcode_length_min"}, $cf{"barcode_length_max"}, $cf{"pool_distance"}, $cf{"pool_distance_internal_postpool"});
my $CONVERTLEN = $BARCODELEN/2 - 2;
my $INPUT = 'barcodeRMvectorSite'.$lr.$fq.'all';
my $OUTPUT = 'barcode'.$fq;
my $tmpFile = $OUTPUT.'N';

open(TMP, "$INPUT.txt") || die;
while(<TMP>) {
     /^([ACGTN]+)\t(\d+)/ || die "$_*$INPUT.txt";
     exists($tl{$1}) ? die : ($tl{$1}=$2);
}
#system "rm $tmpFile.txt";
$num = 0;
@barcode = sort {$tl{$b} <=> $tl{$a}} keys %tl;
for($i = 0; $i < @barcode; $i++) {
	$tl{$barcode[$i]} > 1 || last;
	for($j = $#barcode; $j > $i; $j--) {
		if(&hiseq::nomorethan1indelxS($barcode[$i], $barcode[$j], $DSinternal)) {
			$tl{$barcode[$i]} += $tl{$barcode[$j]};
			$num += $tl{$barcode[$j]};
			delete($tl{$barcode[$j]});
			@barcode = grep {$_ ne $barcode[$j]} @barcode;
		}
	}
}
open(O, ">$OUTPUT.txt") || die;
@barcode = sort {$tl{$b} <=> $tl{$a}} keys %tl;
for($i = 0; $i < @barcode; $i++) {
	print O	"$barcode[$i]	$tl{$barcode[$i]}\n";
}
close O;
