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
my %cf;
open(CF, $config) || die "$config\n";
while(<CF>){
	if($_ =~ /(\S+)\t(\S+)/){
		$cf{$1} = $2;
	}
}
my $fq = $cf{"prepool_reads"};
my (@barcode, @convert, @lines, @readNumShort, @readNumMid, @readNumLong, @t, %barcode, %count, %SQ, %t, %tl, %tr, $barcode, $bh, $br, $i, $j, $k, $l1, $l2, $lowQ, $num, $Q, $readNum, $rmL, $rmR, $t) = ();


my ($BARCODELEN, $BARCODELENMIN, $BARCODELENMAX, $DS, $MERGE, $QUALITY, $VAR,$middle) = ($cf{"barcode_length_standard"}, $cf{"barcode_length_min"}, $cf{"barcode_length_max"}, $cf{"pool_distance"}, '.extendedFrags', $cf{"SQ_threshold"}, $cf{"pool_variance"},$cf{"middle"});
my $CONVERTLEN = $BARCODELEN/2 - 2;
my $OUTPUT = 'barcodeRMvectorSite'.$fq;
my ($LEFT, $SITELEFT, $SITERIGHT, $RIGHT) = ($cf{"prepool_left_upstream"}, $cf{"prepool_left_downstream"}, $cf{"prepool_right_upstream"}, $cf{"prepool_right_downstream"});
my $BRIDGELEN = 1; #nucleotide juxtaposed to barcodes has high mutation rate
$QUALITY !~ /^\d+$/ && die "SQ threshold as a number!";

$SQ{$fq}->[0] = $cf{"prepool_left_anchor"};
$SQ{$fq}->[1] = $cf{"prepool_right_anchor"};
$SQ{$fq}->[0] > length($LEFT) + $VAR || die;
$SQ{$fq}->[1] > length($RIGHT) + $VAR || die;
my ($sl, $slLen, $sr, $srLen) = ($SQ{$fq}->[0]-length($LEFT)-$VAR-$BRIDGELEN, length($LEFT) + 2 * $VAR, $SQ{$fq}->[1] + $VAR + $BRIDGELEN, length($RIGHT) + 2 * $VAR);
my($LEFTMAXLEN, $RIGHTMAXLEN) = ($SQ{$fq}->[0] + $VAR, $SQ{$fq}->[1] + $VAR);
my $miniLen = $SQ{$fq}->[0] + $SQ{$fq}->[1] + $VAR;
my $tmpFile = $OUTPUT.'all';
exists($SQ{$fq}) || die;
open LOG, ">t.txt";
LOG->autoflush(1); 

$t = 0;
$lowQ = 0;
$readNum = 0;
open(FQ, "$fq$MERGE.fastq") || die "$fq$MERGE.fastq";
while(<FQ>){
    push @lines,$_;

    if(@lines == 4) {
	$lines[2] =~ /^\+\s*$/ || die $lines[2];

	if(length($lines[1]) < $miniLen) {
	    $t++;
	} elsif(&hiseq::SQmean(substr($lines[3],$SQ{$fq}->[0], length($lines[3])-$SQ{$fq}->[0]-$SQ{$fq}->[1])) < $QUALITY) {
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
print LOG scalar(keys %t), " raw barcode pairs by $readNum reads.\n";
print LOG "low quality reads $lowQ\tshort reads $t\n";

#$t = 0;
#while(($barcode, $num) = each (%t)) {
#    if($barcode =~ /TTTGAAACGTAATGTTG/) {
#	print "$barcode\t$num\n";
#	$t += $num;
#    }
#}
#print "$t\n";

$t = 0;
$readNum = 0;
%tr = ();
while(($barcode, $num) = each (%t)) {
    $rmL = &hiseq::rmLendBridge($barcode, $LEFT, $sl, $slLen, 1, $BRIDGELEN);

    if($rmL) {
	($rmL < $sl + $BRIDGELEN || $rmL > ($sl + $slLen + $BRIDGELEN)) && die;
	$i = substr($barcode, $rmL);
    } elsif($barcode =~ s/^[ACGTN]{0,$LEFTMAXLEN}([ACGTN]{$BARCODELEN}[ACGTN]{$BRIDGELEN}$SITELEFT)/$1/) {
            #If vector sequence is wrong but BamHI region is perfect, still take it assuming barcode is 20 bp long
	$i = $barcode;
    } else {
	next;
    }

    if(length($i) > $sr) {
	$tr{$i} += $num;
	$readNum += $num;
    } else {
	$t += $num;
    }
 
}
print LOG scalar(keys %tr), " raw barcode pairs by $readNum reads have expected vector sequence at left.\n";
print LOG "short reads $t\n";
%t = ();

#$t = 0;
#while(($barcode, $num) = each (%tr)) {
#    if($barcode =~ /TTTGAAACGTAATGTTG/) {
#	print "$barcode\n";
#	$t += $num;
#    }
#}
#print "$t\n";

$readNum = 0;
%tl = ();
while(($barcode, $num) = each (%tr)) {
    $rmR = &hiseq::rmRendBridge($barcode, $RIGHT, $sr, $srLen, 1, $BRIDGELEN);

    if($rmR) {
	($rmR > $sr + $BRIDGELEN || $rmR < ( $sr - $srLen)) && die "$rmR, $sr, $srLen\t$barcode";
	$i = substr($barcode, 0, length($barcode) - $rmR);
    } elsif($barcode =~ s/($SITERIGHT[ACGTN]{$BRIDGELEN}[ACGTN]{$BARCODELEN})[ACGTN]{0,$RIGHTMAXLEN}$/$1/) {
            #If vector sequence is wrong but HindIII region is perfect, still take it assuming barcode is 20 bp long
	$i = $barcode;
    } else {
	next;
    }

    if(length($i) > 2 * $BARCODELENMIN + length($middle)) {
	    $tl{$i} += $num;
	    $readNum += $num;
    } else {
	    $t += $num;
    }
}
print LOG scalar(keys %tl), " raw barcode pairs by $readNum reads have expected vector sequences at both sides.\n";
print LOG "short reads $t\n";
%tr = ();

#$t = 0;
#while(($barcode, $num) = each (%tl)) {
#    if($barcode =~ /CCGAGCAAGACGTACTATCA/) {
#	print "$barcode\t$num\n";
#	$t += $num;
#    }
#}
#print "$t\n";

$t = 0;
$readNum = 0;
while(($barcode, $num) = each (%tl)) {
    #print "$barcode, $SITERIGHT, ",length($barcode) - $BARCODELENMAX - length($SITERIGHT),",", $BARCODELENMAX - $BARCODELENMIN + length($SITERIGHT),", 1\n";
    $rmL = &hiseq::rmLend($barcode, $SITERIGHT, length($barcode) - $BARCODELENMAX - length($SITERIGHT), $BARCODELENMAX - $BARCODELENMIN + length($SITERIGHT), 1);
    #print "$rmL\n";
    if($rmL) {
	$barcode =~ s/^([ACGTN]{$rmL})[ACGTN]{$BRIDGELEN}// || next; #There is 1bp between tested vector sequence and BARCODE
	$i = $1;
#    } elsif($barcode =~ s/([ACGTN]+$site[ACGTN])[ACGTN]{$BRIDGELEN}//) { #take it as long as there is HindIII site
#	$i = $1;
    } else {
	next;
    }

	if(length($i) > $BARCODELENMIN + length($middle) -1) {
	    #print "$barcode	$i\n";
	    $tr{$barcode}->{$i} += $num;
	    $readNum += $num;
	} else {
	    $t += $num;
	}
}
print LOG scalar(keys %tr), " right barcodes by $readNum reads have expected HindIII-like site.\n";
print LOG "short reads $t\n";

#$t = 0;
#exists($tr{'CCGAGCAAGACGTACTATCA'}) || die;
#while(($barcode, $num) = each (%{$tr{'CCGAGCAAGACGTACTATCA'}})) {
#	print "$barcode\t$num\n";
#	$t += $num;
#}
#print "*$t*\n";

$readNum = 0;
%tl = ();
while(($br, $bh) = each (%tr)) {
    while(($barcode, $num) = each (%{$bh})) {
	$rmR = &hiseq::rmRend($barcode, $SITELEFT, length($barcode) - $BARCODELENMIN, length($barcode) - $BARCODELENMIN, 1);

	if($rmR) {
	    length($barcode) > $rmR + $BRIDGELEN || die;
	    $i = substr($barcode, 0, length($barcode) - $rmR - $BRIDGELEN); #There is 1bp between tested vector sequence and BARCODE
	} elsif(length($barcode) > $BARCODELEN) {
	    $i = substr($barcode, 0, $BARCODELEN);
	} else {
	    next;
	}

	unless(length($i) < $BARCODELENMIN || length($br) < $BARCODELENMIN) {
	    $tl{$i}->{$br} += $num;
	    $readNum += $num;
	}
    }
}
print LOG "$readNum reads have expected HindIII-like site and BamHI site.\n";
%tr = ();

#$t = 0;
#while(($br, $bh) = each (%tl)) {
#    while(($barcode, $num) = each (%{$bh})) {
#	if($barcode =~ /TTTGAAACGTAATGTTG/) {
#	    print "$br\t$barcode\t$num\n";
#	    $t += $num;
#	}
#    }
#}
#print "*$t*\n";

my %bl = ();
my %br = ();
open(TMP, ">$tmpFile.txt") || die;
while(($br, $bh) = each (%tl)) {
	while(($barcode, $num) = each (%{$bh})) {
		print TMP "$br:$barcode\t$num\n";
		$bl{$br} += $num;
		$br{$barcode} += $num;
	}
}
close(TMP);
#system "rm $tmpFile.txt";
print LOG scalar(keys %bl), " left raw barcodes\t", scalar(keys %br), " right raw barcodes\n";
