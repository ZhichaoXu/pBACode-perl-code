#remove vector sequences at both ends of merged read pairs
#Assuming HindIII site is (1) perfect or (2) HindIII has 1bp mismatch but barcode is exactly 20 bp
#1. 10 bp flanking and between two barcodes have no more than 1bp mismatch to vector sequence, separately
#2. If 10 bp flanking two barcodes have more than 1bp mismatch to vector sequence, check 10bp covering BamHI or HindIII site. If these 10 bp regions perfectly match vector sequence, take 20 bp flanking these regionsand bridges as barcodes.
#3. If 10 bp vector sequence flanking two barcodes have no more than 1bp mismatch, but sequence between two barcodes have more than 1bp mismatch, as long as there is perfect HindIII site, take 20 bp juxtaposed to flanking vector sequence as barcodes
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
my (@barcode, @barcode1bpMid, @convert, @lines, @readNumShort, @readNumMid, @readNumLong, @t, %barcode, %bl, %br, %count, %mismatchL, %mismatchR, %SQ, %t, %tl, %tr, %ttl, %ttr, $barcode, $bh, $bl, $br, $i, $j, $k, $l1, $l2, $lowQ, $num, $Q, $readNum, $rmL, $rmR, $t) = ();
my ($BARCODELEN, $BARCODELENMIN, $BARCODELENMAX, $DS, $DSinternal, $MERGE, $QUALITY, $VAR) = ($cf{"barcode_length_standard"}, $cf{"barcode_length_min"}, $cf{"barcode_length_max"}, $cf{"pool_distance"}, $cf{"pool_distance_internal_prepool"}, '.extendedFrags', $cf{"SQ_threshold"}, $cf{"pool_variance"});
my $CONVERTLEN = $BARCODELEN/2 - 2;
my $INPUT = 'barcodeRMvectorSite'.$fq.'all';
my $OUTPUT = 'barcodePrepool'.$fq;
my ($LEFT, $SITELEFT, $SITERIGHT, $RIGHT) = ($cf{"left_upstream"}, $cf{"left_downstream"}, $cf{"right_upstream"}, $cf{"right_downstream"});
my $BRIDGELEN = 1; #nucleotide juxtaposed to barcodes has high mutation rate
$QUALITY !~ /^\d+$/ && die "SQ threshold as a number!";

$SQ{$fq}->[0] = $cf{"prepool_left_anchor"};
$SQ{$fq}->[1] = $cf{"prepool_right_anchor"};
$SQ{$fq}->[0] > length($LEFT) + $VAR || die;
$SQ{$fq}->[1] > length($RIGHT) + $VAR || die;
my ($sl, $slLen, $sr, $srLen) = ($SQ{$fq}->[0]-length($LEFT)-$VAR-$BRIDGELEN, length($LEFT) + 2 * $VAR, $SQ{$fq}->[1] + $VAR + $BRIDGELEN, length($RIGHT) + 2 * $VAR);
my($LEFTMAXLEN, $RIGHTMAXLEN) = ($SQ{$fq}->[0] + $VAR, $SQ{$fq}->[1] + $VAR);
my $miniLen = $SQ{$fq}->[0] + $SQ{$fq}->[1] + $VAR;
my $tmpFile = $OUTPUT.'N';
exists($SQ{$fq}) || die;
open LOG, ">tb.txt";
LOG->autoflush(1); 

open(TMP, "$INPUT.txt") || die;
while(<TMP>) {
     /^([ACGTN]+)\:([ACGTN]+)\t(\d+)/ || die "$_*$INPUT.txt";
     exists($tl{$1}->{$2}) ? die : ($tl{$1}->{$2} = $3);
}
#system "rm $tmpFile.txt";

$num = 0;
%tl = &hiseq::rmHybrid1(%tl);

while(($br, $bh) = each (%tl)) {
    %t = %{$tl{$br}};
    @barcode = sort {$t{$b} <=> $t{$a}} keys %t;
    for($i = 0; $i < @barcode; $i++) {
	$tl{$br}->{$barcode[$i]} > 1 || last;

	for($j = $#barcode; $j > $i; $j--) {
	    if(&hiseq::nomorethan1indelxS($barcode[$i], $barcode[$j], $DSinternal)) {
		$tl{$br}->{$barcode[$i]} += $tl{$br}->{$barcode[$j]};
		$num += $tl{$br}->{$barcode[$j]};
		delete($tl{$br}->{$barcode[$j]});
		@barcode = grep {$_ ne $barcode[$j]} @barcode;
	    }
	}
    }
}
print LOG "processed 1l:Nr barcode-pairs by correctifying $num reads whose right barcodes differ by no more than $DSinternal bp from other multi-read right barcodes in the 1l:Nr group\n";

%tr = ();
while(($br, $bh) = each (%tl)) {
    while(($barcode, $num) = each (%{$bh})) {
	$tr{$barcode}->{$br} = $num;
    }
}
%tl = ();
print scalar(keys %tr), " right barcodes\n";

$num = 0;
while(($br, $bh) = each (%tr)) {
    %t = %{$tr{$br}};
    @barcode = sort {$t{$b} <=> $t{$a}} keys %t;
    for($i = 0; $i < @barcode; $i++) {
	$tr{$br}->{$barcode[$i]} > 1 || last;

	for($j = $#barcode; $j > $i; $j--) {
	    if(&hiseq::nomorethan1indelxS($barcode[$i], $barcode[$j], $DSinternal)) {
		$tr{$br}->{$barcode[$i]} += $tr{$br}->{$barcode[$j]};
		$num += $tr{$br}->{$barcode[$j]};
		delete($tr{$br}->{$barcode[$j]});
		@barcode = grep {$_ ne $barcode[$j]} @barcode;
	    }
	}
    }
}
print LOG "processed Nl:1r barcode-pairs by correctifying $num reads whose left barcodes differ by no more than $DSinternal bp from other multi-read left barcodes in the Nl:1r group\n";

open(TMP, ">$tmpFile.txt") || die;
while(($br, $bh) = each (%tr)) {
        while(($barcode, $num) = each (%{$bh})) {
                print TMP "$barcode:$br\t$num\n";
                $bl{$barcode} += $num;
                $br{$br} += $num;
        }
}
close(TMP);
print LOG scalar(keys %bl), " left raw barcodes\t", scalar(keys %br), " right raw barcodes\n";

$readNum = 0;
@readNumLong = ({}, {});
@readNumMid = ({}, {});
@readNumShort = ({}, {});
while(($br, $bh) = each (%tr)) {
    while(($barcode, $num) = each (%{$bh})) {
	$readNum += $num;
	$l1 = length($br);
	$l2 = length($barcode);
	$t = $CONVERTLEN*2;

	if($l1 < $t - 1) {
	    $readNumShort[1]->{$br} += $num;
	} elsif($l1 < $t) {
	    $readNumShort[1]->{$br} += $num;
	    $readNumMid[1]->{$br} += $num;
	} elsif($l1 < $t + 1) {
	    $readNumMid[1]->{$br} += $num;
	    $readNumLong[1]->{$br} += $num;
	} else {
	    $readNumLong[1]->{$br} += $num;
	}

	if($l2 < $t - 1) {
	    $readNumShort[0]->{$barcode} += $num;
	} elsif($l2 < $t) {
	    $readNumShort[0]->{$barcode} += $num;
	    $readNumMid[0]->{$barcode} += $num;
	} elsif($l2 < $t + 1) {
	    $readNumMid[0]->{$barcode} += $num;
	    $readNumLong[0]->{$barcode} += $num;
	} else {
	    $readNumLong[0]->{$barcode} += $num;
	}
    }
}

@convert = ();
@readNumLong != 2 && die scalar @readNumLong;
for($k = 0; $k < @readNumLong; $k++) {
    $convert[$k] = {&hiseq::mismatch1($CONVERTLEN, %{$readNumLong[$k]})};
    $k ? ($t = 'right') : ($t = 'left');
    print LOG scalar(keys %{$convert[$k]}), " $t long barcodes have 1 mutation from more enriched ones\n";

    %t = %{$readNumShort[$k]};
    @barcode = sort {$t{$b} <=> $t{$a}} keys %t;

    for($i = 0; $i < @barcode; $i++) {
	$t{$barcode[$i]} > 1 || last;
#	$i =~ /0000$/ && print "$i right barcodes have been paired processed\n";

	for($j = $#barcode; $j > $i; $j--) {
	    if(&hiseq::nomorethan1indelxS($barcode[$i], $barcode[$j], $DS)) {
		if(exists($convert[$k]->{$barcode[$j]})) {
		    die "$barcode[$j]:", $convert[$k]->{$barcode[$j]}, "\t$barcode[$i]";
		} else {
		    $convert[$k]->{$barcode[$j]} = $barcode[$i];
		    @barcode = grep {$_ ne $barcode[$j]} @barcode;

		}
	    }
	}
    }

    $k ? ($t = 'right') : ($t = 'left');
    print LOG scalar(keys %{$convert[$k]}), " $t long-mid barcodes have 1 mutation from more enriched ones\n";

    %t = %{$readNumMid[$k]};
    @barcode = sort {$t{$b} <=> $t{$a}} keys %t;

    for($i = 0; $i < @barcode; $i++) {
	$t{$barcode[$i]} > 1 || last;
#	$i =~ /0000$/ && print "$i right barcodes have been paired processed\n";

	for($j = $#barcode; $j > $i; $j--) {
	    abs(length($barcode[$i]) - length($barcode[$j])) > 1 && die;
	    if(&hiseq::nomorethan1indelxS($barcode[$i], $barcode[$j], 0)) {#only check indel because substitution has been removed in last 2 steps
		if(exists($convert[$k]->{$barcode[$j]})) {
		    print "warning: indel $barcode[$j]:", $convert[$k]->{$barcode[$j]}, "\t$barcode[$i]\n";
		} else {
		    $convert[$k]->{$barcode[$j]} = $barcode[$i];
		    @barcode = grep {$_ ne $barcode[$j]} @barcode;
		}
	    }
	}
    }

    $k ? ($t = 'right') : ($t = 'left');
    print LOG scalar(keys %{$convert[$k]}), " $t barcodes have 1 mutation from more enriched ones\n";

	$t = 0;
	while(($barcode, $bh) = each(%{$convert[$k]})) {
	    if(exists($convert[$k]->{$bh})) {
#		&hiseq::nomorethan1indelxS($barcode, $convert[$k]->{$bh}, 2*$DS) || die "$barcode\t$bh\t", $convert[$k]->{$bh};
		$barcode1bpMid[$k]->{$bh} = 1;
#		$convert[$k]->{$barcode} = $convert[$k]->{$bh};
#		print "$barcode is $DS bp from $bh, which is $DS bp from ", $convert[$k]->{$bh}, ". Now consider the 1st barcode $DS bp from the 3rd one\n";
		$t++;
	    }
	}
}
%tl = ();

while(($br, $bh) = each (%tr)) {
    while(($bl, $num) = each (%{$bh})) {
	if(exists($convert[0]->{$bl})) {
	    exists($mismatchL{$bl}->{$br}) ? die : ($mismatchL{$bl}->{$br} = $num);
	}
	if(exists($convert[1]->{$br})) {
	    exists($mismatchR{$br}->{$bl}) ? die : ($mismatchR{$br}->{$bl} = $num);
	}
    }
}
scalar(keys %mismatchL) != scalar(keys %{$convert[0]}) && die;
scalar(keys %mismatchR) != scalar(keys %{$convert[1]}) && die;

$j = 0;
BR: while(($br, $bh) = each (%mismatchR)) {
    exists($convert[1]->{$br}) || die;
    exists($tr{$br}) || die;

    if(exists($tr{$convert[1]->{$br}})){
        %t = %{$tr{$convert[1]->{$br}}};
        @barcode = sort {$t{$b} <=> $t{$a}} keys %t;

        while(($barcode, $num) = each (%{$bh})) {
            $tr{$br}->{$barcode} < $num && die "$br:$barcode should have $num reads";

            for($i = 0; $i < @barcode; $i++) {
                $t{$barcode[$i]} > 1 || last;

                if(&hiseq::nomorethan1indelxS($barcode[$i], $barcode, $DSinternal)){
                    $t{$barcode[$i]} > $num || print "warning: ", $convert[1]->{$br}, ":$barcode[$i] should have more reads than $br:$barcode, $t{$barcode[$i]} $num\n";
                    $tr{$convert[1]->{$br}}->{$barcode[$i]} += $tr{$br}->{$barcode};
		        $j++;

                    if(scalar(keys %{$tr{$br}}) > 1) {
                        delete($tr{$br}->{$barcode});
                        last;
                    } else {
                        delete($tr{$br});
                        next BR;
                    }
                }
            }
        }
    } else {
        exists($barcode1bpMid[1]->{$convert[1]->{$br}}) || die ("$br\t", $convert[1]->{$br});
    }
}
print LOG scalar(keys %tr), " right barcodes\n";
print LOG "$j left barcodes differ from more eriched left barcodes by <=$DSinternal bp in case of only 1 mismatch between their right barcodes. ";
print LOG scalar(keys %tr), " right barcodes left after removing them\n";

$t = 0;
while(($br, $bh) = each (%tr)) {
    while(($bl, $num) = each (%{$bh})) {
        $tl{$bl}->{$br} += $num;
        $t += $num;
    }
}
$t != $readNum && die;
print LOG scalar(keys %tl), " left barcodes\n";
%tr = ();

$j = 0;
BL: while(($bl, $bh) = each(%mismatchL)) {
    exists($convert[0]->{$bl}) || die;
    exists($tl{$bl}) || next;
    exists($tl{$convert[0]->{$bl}}) || next;
    %t = %{$tl{$convert[0]->{$bl}}};
    @barcode = sort {$t{$b} <=> $t{$a}} keys %t;

    while(($barcode, $num) = each (%{$bh})) {
	exists($tl{$bl}->{$barcode}) || next;
	$tl{$bl}->{$barcode} < $num && die "$bl:$barcode should have $num reads";

	for($i = 0; $i < @barcode; $i++) {
	    $t{$barcode[$i]} > 1 || last;

	    if(&hiseq::nomorethan1indelxS($barcode[$i], $barcode, $DSinternal)) {
		$t{$barcode[$i]} > $num || print "warning: ", $convert[0]->{$bl}, ":$barcode[$i] should have more reads than $bl:$barcode, $t{$barcode[$i]} $num\n";
		$tl{$convert[0]->{$bl}}->{$barcode[$i]} += $tl{$bl}->{$barcode};
		$j++;

                if(scalar(keys %{$tl{$bl}}) > 1) {
                    delete($tl{$bl}->{$barcode});
                    last;
                } else {
                    delete($tl{$bl});
                    next BL;
                }
	    }
	}
    }

    scalar(keys %{$tl{$bl}}) || delete($tl{$bl});
}
print LOG scalar(keys %tl), " left barcodes\n";
print LOG "$j right barcodes differ from more eriched right barcodes by <=$DSinternal bp in case of only 1 mismatch between their left barcodes. ";
print LOG scalar(keys %tl), " left barcodes left after removing them\n";
%tl = &hiseq::rmHybrid1(%tl);
print LOG scalar(keys %tl), " left barcodes after 2nd rmoval of hybrid\n";

while(($bl, $bh) = each (%tl)) {
    while(($br, $num) = each (%{$bh})) {
        $ttl{$bl} += $num;
        $ttr{$br} += $num;
    }
}

while(($bl, $bh) = each (%tl)) {
    while(($br, $num) = each (%{$bh})) {
        if(exists($convert[0]->{$bl}) && exists($convert[1]->{$br})) {
	    #print "$bl	$br	$num\n";
	    #print "$convert[0]->{$bl}	$convert[1]->{$br}\n";
            exists($ttl{$convert[0]->{$bl}}) || die;
            exists($ttr{$convert[1]->{$br}}) || die;
	    #print "$ttl{$convert[0]->{$bl}}	$ttr{$convert[1]->{$br}}\n";
            if($num > 1) {
                $ttl{$convert[0]->{$bl}}/$num < 10 && next;
                $ttr{$convert[1]->{$br}}/$num < 10 && next;
            }

            print LOG "$bl\t", $convert[0]->{$bl}, "\t", scalar(keys %{$tl{$bl}}), "\n";

            if(scalar(keys %{$tl{$bl}}) > 1) {
                delete($tl{$bl}->{$br});
            } elsif(scalar(keys %{$tl{$bl}}) < 1) {
                die;
            } else {
                delete($tl{$bl});
                last;
            }
        }
    }
}
print LOG scalar(keys %tl), " left barcodes after removal of barcode pairs, each of whose barcode differs by 1bp from more enriched ones\n";

%count = ();
open(O, ">$OUTPUT.txt") || die;
open(L, ">$OUTPUT"."similarbarcodeL.txt") || die;
open(R, ">$OUTPUT"."similarbarcodeR.txt") || die;
while(($br, $bh) = each (%tl)) {
    if(exists($convert[0]->{$br})) {
	exists($tl{$convert[0]->{$br}}) || print LOG "warning: $br should be 1bp from a real left barcode\n";
	print L "$br\t", $convert[0]->{$br}, "\n";
	delete($convert[0]->{$br});
    }

    %t = %{$tl{$br}};
    @barcode = sort {$t{$b} <=> $t{$a}} keys %t;

    for($i = 0; $i < @barcode; $i++) {
	print O "$br:$barcode[$i]\t", $tl{$br}->{$barcode[$i]}, "\n";

	if(exists($convert[1]->{$barcode[$i]})) {
	    print R "$barcode[$i]\t", $convert[1]->{$barcode[$i]}, "\n";
	    delete($convert[1]->{$barcode[$i]});
	}
    }

    if(exists($tl{$br}->{$barcode[0]})) {
	$tl{$br}->{$barcode[0]} =~ /^\d+$/ || print LOG "warning: ", $tl{$br}->{$barcode[0]}, " should be number $br:$barcode[0]\n";
	$count{$tl{$br}->{$barcode[0]}}++;
    } else {
	die "should exist $br:$barcode[0]";
    }
}

@t = sort {$b <=> $a} keys %count;
for($i = 0; $i < @t; $i++) {
    $t[$i] < 10 && print LOG "$t[$i]\t$count{$t[$i]}\n";
}
