#remove vector and barcode sequences of right barcode-linked genome inserts
#Assuming HindIII site is (1) perfect or (2) HindIII has 1bp mismatch but barcode is exactly 20 bp
#INPUT: barcodes are in file1 of mated reads
#! usr/bin/perl
use strict;
use warnings;
use List::Util qw(max);
require hiseq;
sub SQend($$$$);

my $config = $ARGV[0];
my $jump = $ARGV[1];
my $lr;
if(defined($jump)){
	$lr = $jump;
	$lr eq 'jumpleft' || die;
}else{
	$lr = 'left';
}
my %cf;
open(CF, $config) || die "$config\n";
while(<CF>){
	if($_ =~ /(\S+)\t(\S+)/){
		$cf{$1} = $2;
	}
}
my $fq = $cf{$lr."_reads"};
my (@lines, @rm, @t, %barcode, %barcodeBrokenBAMHI, %barcodeBrokenHINDIII, %BRIDGELEN, %noHindIII, %RESTRICTIONSITE, %SITE, %SQ, %t, %UPSTREAM);
my ($bridgeLen, $count, $i, $lowQ, $MINILEN, $Q, $readNum, $readNum1, $readNum2, $readNum3, $readNum4, $rmL, $rmR, $site, $sl, $slLen, $ss, $t, $upstream) = ();
my ($BARCODELEN, $BARCODELENMIN, $BARCODELENMAX) = ($cf{"barcode_length_standard"}, $cf{"barcode_length_min"}, $cf{"barcode_length_max"});
my ($MERGE, $NOHINDIII, $QUALITY,  $QUALITYTRIM, $TRIMLEN, $VAR) = ('32.notCombined_', 'noHindIII', $cf{"SQ_threshold"}, $cf{"SQG_threshold"}, $cf{"flag_length"}, $cf{"BAC_variance"});
my @MINILEN = (10, 20);
my $OUTPUT = 'rmVectorBarcode1mate'.$fq;
my $BARCODEBRIDGELEN = 1; #$UPSTREAMBARCODE is 1bp from barcode becasue that bp has high mutation rate
                          #PvuII half-site tend to have 1bp mutation during ligation
$BRIDGELEN{$fq} = [1, 1]; #1bp juxtaposed to barcodes has high mutation rate
$UPSTREAM{$fq} = [$cf{$lr."_upstream1"}, $cf{$lr."_upstream2_rc"}];
$QUALITY !~ /^\d+$/ && die "SQ threshold as a number!";
$SITE{$fq} = [$cf{$lr."_cutsite3"}, '']; #HINDIII and BamHI with 1bp from barcode
my %RESTRICTIONSITELEFT = ($fq, $cf{$lr."_cutsite2"}); #BamHI with 1bp from barcode left end
my %SITEBRIDGELEN = ($fq, 12); #barcdoe is 12bp from HindIII site
my $CHOPLEN = length($SITE{$fq}->[0]) + $SITEBRIDGELEN{$fq};
$SQ{$fq} = [$cf{$lr."_anchor1"}, $cf{$lr."_anchor2"}];
$SQ{$fq}->[0] > length($UPSTREAM{$fq}->[0]) + $VAR || die;
$SQ{$fq}->[1] > length($UPSTREAM{$fq}->[1]) + $VAR || die;
my $miniLen = max($SQ{$fq}->[0], $SQ{$fq}->[1]) + $VAR;
my $tmpFile = $OUTPUT.'tmp';
exists($SQ{$fq}) || die;

for($i = 0; $i < @{$SQ{$fq}}; $i++) {
    $lowQ = 0;
    $readNum = 0;
    $readNum1 = 0;
    $readNum2 = 0;
    $readNum3 = 0;
    $readNum4 = 0;
    $MINILEN = $MINILEN[$i];
    exists($UPSTREAM{$fq}->[$i]) ? ($upstream = $UPSTREAM{$fq}->[$i]) : die;
    $site = $SITE{$fq}->[$i];
    $ss = $SQ{$fq}->[$i];
    $sl = $ss - length($upstream) - $VAR;
    $slLen = length($upstream) + 2 * $VAR;
    $bridgeLen = $BRIDGELEN{$fq}->[$i];
    open(FQ, "$fq$MERGE".($i+1).".fastq") || die "$fq$MERGE".($i+1).".fastq";
    open(TMP, ">$tmpFile".($i+1).".fastq") || die;
    $i || (open(NO, ">$tmpFile".($i+1)."$NOHINDIII.fastq") || die);
    while(<FQ>){
	push @lines,$_;

	if(@lines == 4){
	    $lines[2] =~ /^\+\s*$/ || die $lines[2];

	    if(length($lines[1]) < $miniLen) {
		$readNum4++;
		@lines = ();
		next;
	    }

	    $Q = &hiseq::SQmean(substr($lines[3], $ss, 20));
	    ++$readNum =~ /000000$/ && print "$readNum\t$readNum3\n";
	    @rm = ();

	    if($Q < $QUALITY) {
		$lowQ++;
		@lines = ();
		next;
	    }

	    chomp($lines[0]);
	    chomp($lines[1]);
	    chomp($lines[2]);
	    chomp($lines[3]);
	    $rmL = &hiseq::rmLendBridge($lines[1], $upstream, $sl, $slLen, 1, $bridgeLen);

	    if($rmL) {
		$rmL >= $sl && $rmL <= ( $sl + $slLen + $bridgeLen) || die "$lines[1]\t$rmL\t$bridgeLen";
		    $lines[1] =~ s/^[ACGTN]{$rmL}([ACGTN]+)$/$1/ || die;
		    $lines[3] = substr($lines[3], 0 - length($lines[1]));

		if(length($lines[1]) > length($site) + $SITEBRIDGELEN{$fq} + $BARCODELEN) {
		    @rm = &hiseq::rmSite($lines[1], $site, 0, $BARCODELEN+$SITEBRIDGELEN{$fq}, $BARCODELENMAX+$SITEBRIDGELEN{$fq}, 'l');
		} else {
		    $rmL = 0;
		    $readNum4++;
		}
	    } elsif($lines[1] =~ s/^([ACGTN]+$upstream)//) {
		    if(length($lines[1]) > length($site) + $SITEBRIDGELEN{$fq} + $BARCODELEN) {
			$lines[3] = substr($lines[3], 0 - length($lines[1]));
			@rm = &hiseq::rmSite($lines[1], $site, 0, $BARCODELEN+$SITEBRIDGELEN{$fq}, $BARCODELENMAX+$SITEBRIDGELEN{$fq}, 'l');
			$rmL = 1;
		    } else {
			$readNum4++;
		    }
	    }

	    if(@rm == 2) {
		$lines[1] = $site.$rm[0];
		$lines[3] = substr($lines[3], 0 - length($lines[1]));
		$lines[3] = &hiseq::SQtrimEnd3($lines[3], $TRIMLEN, $QUALITYTRIM);

		if(length($lines[3]) > $MINILEN) {
		    $lines[1] = substr($lines[1], 0, length($lines[3]));

                    @t = &hiseq::rmSite($rm[1], $RESTRICTIONSITELEFT{$fq}, $BARCODEBRIDGELEN, $BARCODELEN, $BARCODELENMAX, 'l');

                    if(@t == 2) {
                        $t = $t[1];
                    } elsif($rm[1] =~ /^([ACGTN]{$BARCODELEN})/) {
                        $t = $1;
                    } else {
                        $t = $rm[1];
			$barcodeBrokenBAMHI{$t}++;
                    }


		    $lines[0] =~ s/^([^\s]+)[\/\s].*$/$1BARCODE$t/ || die;
		    $barcode{$1}->[$i] = $t;
		    $lines[1] =~ /\s/ && die;
		    $lines[2] =~ /\s/ && die;
		    $lines[3] =~ /\s/ && die;
		    print TMP join("\n", @lines), "\n";
		    $readNum3++;
		} else {
		    $readNum1++;
		}
	    } elsif($rmL) {
		if($site) {
		    $i && die;
		    $lines[1] =~ s/^([ACGTN]{$BARCODELEN})[ACGTN]{$CHOPLEN}// || die $lines[1];
		    $t = $1;
		    $lines[3] = substr($lines[3], 0 - length($lines[1]));
		    $lines[3] = &hiseq::SQtrimEnd3($lines[3], $TRIMLEN, $QUALITYTRIM);

		    if(length($lines[3]) > $MINILEN) {
			$lines[1] = substr($lines[1], 0, length($lines[3]));
			$barcodeBrokenHINDIII{$t}++;
			$lines[0] =~ s/^([^\s]+)[\/\s].*$/$1BARCODE$t/ || die;
			$noHindIII{$1} = $t;
			$lines[1] =~ /\s/ && die;
			$lines[2] =~ /\s/ && die;
			$lines[3] =~ /\s/ && die;
			$readNum2++;
			print NO join("\n", @lines), "\n";
		    } else {
			$readNum1++;
		    }
		} else {
		    $i || die;
		    $lines[3] = &hiseq::SQtrimEnd3($lines[3], $TRIMLEN, $QUALITYTRIM);

		    if(length($lines[3]) > $MINILEN && &SQend($lines[3], $TRIMLEN, $QUALITYTRIM, 5)) {
			$lines[1] = substr($lines[1], 0, length($lines[3]));
			$lines[0] =~ s/^([^\s]+)[\/\s].*$/$1/ || die;
			$barcode{$1}->[$i] = 1;
			$lines[1] =~ /\s/ && die;
			$lines[2] =~ /\s/ && die;
			$lines[3] =~ /\s/ && die;
			$readNum3++;
			print TMP join("\n", @lines), "\n";
		    } else {
			$readNum1++;
		    }
		}
	    }

	    @lines = ();
	} elsif(@lines > 4) {
	    die @lines;
	}
    }

    print "$lowQ low quality reads, $readNum4 too short and $readNum1 ones with low quality genome among $readNum total reads\n";
    $i ? print "reads with expected non-barcode-side vector sequence\t$readNum3\n" :  print "reads with expected barcode-side vector sequence and restriction site\t$readNum3\nreads with expected barcode-side vector sequence but no restriction site\t$readNum2\n";
}

$t = 0;
foreach $count (values %barcodeBrokenBAMHI) {
    $count > 2 && $t++;
}
print "$t barcodes have broken BamHI and NheI which are likely picked up in post-pool.\n";

open(ONLY, ">$OUTPUT"."barcodesideOnly.fastq") || die;
open(NO, ">$OUTPUT$NOHINDIII"."2.fastq") || die;
for($i = 0; $i < @{$SQ{$fq}}; $i++) {
    open(FQ, ">$OUTPUT".($i+1).".fastq") || die;
    open(TMP, "$tmpFile".($i+1).".fastq") || die;
 
    while(<TMP>) {
	push @lines,$_;

	if(@lines == 4){
	    if($lines[0] =~ /^([^\s]+)BARCODE([ACGTN]+)/) {
		$i && die;
		exists($barcode{$1}->[0]) || die;
		exists($barcode{$1}->[1]) ? print FQ (@lines) : print ONLY (@lines);
	    } elsif($lines[0] =~ /^([^\s]+)$/) {
		$i || die;
		exists($barcode{$1}->[1]) || die;

		if(exists($barcode{$1}->[0])) {
		    exists($noHindIII{$1}) && die;
		    $t = $barcode{$1}->[0];
		    $lines[0] =~ s/^([^\s]+)/$1BARCODE$t/ || die;
		    print FQ @lines;
		} elsif(exists($noHindIII{$1})) {
		    $t = $noHindIII{$1};
		    $lines[0] =~ s/^([^\s]+)/$1BARCODE$t/ || die;
		    print NO @lines;
		}
	    } else {
		die;
	    }

	    @lines = ();
	} elsif(@lines > 4) {
	    die scalar(@lines);
	}
    }
    close(FQ);
    close(TMP);
    unlink "$tmpFile".($i+1).".fastq";
}
close(NO);

open(ONLY, ">$OUTPUT"."barcodesideOnly$NOHINDIII.fastq") || die;
    open(FQ, ">$OUTPUT$NOHINDIII"."1.fastq") || die;
    open(TMP, "$tmpFile"."1$NOHINDIII.fastq") || die;
 
    while(<TMP>) {
	push @lines,$_;

	if(@lines == 4){
	    if($lines[0] =~ /^([^\s]+)BARCODE([ACGTN]+)/) {
		exists($noHindIII{$1}) || die;
		exists($barcode{$1}->[1]) ? print FQ (@lines) : print ONLY (@lines);
	    } else {
		die;
	    }

	    @lines = ();
	} elsif(@lines > 4) {
	    die scalar(@lines);
	}
    }
    close(FQ);
    close(TMP);
    unlink "$tmpFile"."1$NOHINDIII.fastq";
close(ONLY);
#open(O, ">$OUTPUT"."BarcodeWithoutHindIII.txt") || die;
#@t = sort {$barcodeBrokenHINDIII{$b} <=> $barcodeBrokenHINDIII{$a}} keys %barcodeBrokenHINDIII;
#for($i = 0; $i < @t; $i++) {
#    print O "$t[$i]\t", $barcodeBrokenHINDIII{$t[$i]}, "\n";
#}
#close(O);

sub SQend($$$$) {
    my ($seq, $len, $threshold, $end) = @_;
    $len > 0 || die;
    $threshold < 20 && die; #must be wrong parameter
    my ($seqLen, $i) = (length($seq), );
    $seqLen < $len && return 0;

    if($end eq '3') {
        for($i = $seqLen - 1; $i >= $seqLen - $len; $i--) {
            ord(substr($seq, $i, 1)) < $threshold+33 && return 0;
        }
    } elsif($end eq '5') {
        for($i = 0; $i < $len; $i++) {
            ord(substr($seq, $i, 1)) < $threshold+33 && return 0;
        }
    }

    return 1;
}
