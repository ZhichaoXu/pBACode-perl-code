#remove vector sequences at both ends of merged read pairs
#Assuming HindIII site is (1) perfect or (2) HindIII has 1bp mismatch but barcode is exactly 20 bp
#1. 10 bp flanking and between two barcodes have no more than 1bp mismatch to vector sequence
#2. If 10 bp flanking two barcodes have more than 1bp mismatch to vector sequence, check 10bp covering BamHI or HindIII site. If these 10 bp regions perfectly match vector sequence, take 20 bp flanking these regionsand bridges as barcodes.
#3. If 10 bp vector sequence flanking two barcodes have no more than 1bp mismatch, but sequence between two barcodes have more than 1bp mismatch, as long as there is perfect HindIII site, take 20 bp juxtaposed to flanking vector sequence as barcodes
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
my (@barcode, @convert, @lines, @t, %barcode, %count, %SQ, %t, %tl, %tr, $barcode, $bh, $br, $i, $j, $k, $l1, $l2, $lowQ, $num, $Q, $readNum, $rmL, $rmR, $t) = ();
my ($BARCODELEN, $BARCODELENMIN, $BARCODELENMAX, $MERGE, $QUALITY, $VAR) = ($cf{"barcode_length_standard"}, $cf{"barcode_length_min"}, $cf{"barcode_length_max"}, '', $cf{"SQ_threshold"}, $cf{"pool_variance"});
my $CONVERTLEN = $BARCODELEN/2 - 2;
my $OUTPUT = 'barcodeRMvectorSite'.$lr.$fq;
my $up = $cf{$lr."_upstream"};
my $down = $cf{$lr."_downstream"};
my $BRIDGELEN = 1; #nucleotide juxtaposed to barcodes has high mutation rate
my $BARCODERMSITE = $BARCODELENMIN + $BRIDGELEN;
$QUALITY !~ /^\d+$/ && die "SQ threshold as a number!";
#my ($LEFT, $SITELEFT, $SITERIGHT, $RIGHT) = ($cf{"left_upstream"}, $cf{"left_downstream"}, $cf{"right_upstream"}, $cf{"right_downstream"});
$SQ{$fq}->[0] = $cf{$lr."_anchor"};
$SQ{$fq}->[0] > length($up) + $VAR || die;
my ($sl, $slLen) = ($SQ{$fq}->[0]-length($up)-$VAR-$BRIDGELEN, length($up) + 2 * $VAR);
my($LEFTMAXLEN) = ($SQ{$fq}->[0] + $VAR);
my $miniLen = $SQ{$fq}->[0] + $BARCODELEN + length($down) + $VAR;
my $tmpFile = $OUTPUT.'all';
exists($SQ{$fq}) || die;

$t = 0;
$lowQ = 0;
$readNum = 0;
open(FQ, "$fq$MERGE") || die "$fq$MERGE";
while(<FQ>){
    push @lines,$_;

    if(@lines == 4) {
	$lines[2] =~ /^\+\s*$/ || die $lines[2];

	if(length($lines[1]) < $miniLen) {
	    $t++;
	} elsif(&hiseq::SQmean(substr($lines[3],$SQ{$fq}->[0]-$VAR, $BARCODELEN+2*$VAR)) < $QUALITY) {
	    $lowQ++;
	} else {
	    chomp($lines[1]);
	    $t{$lines[1]}++;
	    $readNum++;
	}

	@lines = ();
    } elsif(@lines > 4) {
	die;
    }
}
print "$readNum reads.\n";
print "low quality reads $lowQ\tshort reads $t\n";

$t = 0;
$readNum = 0;
%tr = ();
while(($barcode, $num) = each (%t)) {
    $rmL = &hiseq::rmLendBridge($barcode, $up, $sl, $slLen, 1, $BRIDGELEN);
    if($rmL) {
	($rmL < $sl + $BRIDGELEN || $rmL > ($sl + $slLen + $BRIDGELEN)) && die;
	$i = substr($barcode, $rmL);
    } elsif($barcode =~ s/^[ACGTN]{0,$LEFTMAXLEN}([ACGTN]{$BARCODELEN}[ACGTN]{$BRIDGELEN}$down)/$1/) {
            #If vector sequence is wrong but BamHI region is perfect, still take it assuming barcode is 20 bp long
	$i = $barcode;
    } else {
	next;
    }
	$tr{$i} += $num;
	$readNum += $num;
#    }
 
}
print scalar(keys %tr), " raw barcode pairs by $readNum reads have expected vector sequence at left.\n";
print "short reads $t\n";
%t = %tr;
%tr = ();
$readNum = 0;
while(($barcode, $num) = each %t) {
	$rmR = &hiseq::rmRend($barcode, $down, length($barcode) - $BARCODELENMIN, length($barcode) - $BARCODELENMIN, 1);
	if($rmR) {
	    length($barcode) > $rmR + $BRIDGELEN || die;
		if($barcode =~ /^([ACGTN]{1,$BARCODELENMAX})[ATCGN][ACGTN]{$rmR}/){
			$i = $1;
		}else{
			die length($barcode),"	$rmR\n";
		}
	}else {
	    next;
	}
	$tr{$i} += $num;
	$readNum += $num;
}
print scalar(keys %tr), " barcodes by $readNum reads.\n";

open(TMP, ">$tmpFile.txt") || die;
while(($barcode, $num) = each (%tr)) {
	print TMP "$barcode\t$num\n";
}
close(TMP);
