#remove vector sequences at both ends of merged read pairs of Left Barcode-mates which are separated from genome inserts by HindIII site and other two restriction sites
#barcode is at 5' side of merged reads, while the non-barcode end $END3-bp of genome insert needs to reach Quality Score $QUALITYG
#Assuming HindIII site is (1) perfect or (2) HindIII has 1bp mismatch but barcode is exactly 20 bp
#! usr/bin/perl
use strict;
use warnings;
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

my (@lines, @rm, @t, %barcode, %barcodeBrokenHINDIII, %t, $count, $i, $file, $lowQ, $Q, $readNum, $readNum2, $readNum3, $readNum4, $readNum5, $rmL, $rmR, $site, $t) = ();
my ($BARCODELEN, $BARCODELENMIN, $BARCODELENMAX, $END3, $MERGE, $QUALITY, $QUALITYG, $VAR) = ($cf{"barcode_length_standard"}, $cf{"barcode_length_min"}, $cf{"barcode_length_max"}, $cf{"flag_length"}, '.extendedFrags', $cf{"SQ_threshold"}, $cf{"SQG_threshold"}, $cf{"BAC_variance"});
my $readlength = $cf{$lr."_readlength"};
my $OUTPUT = 'rmVectorBarcode'.$fq;
my %VECTOR;
$VECTOR{$fq}->[0] = $cf{$lr."_upstream1"};
$VECTOR{$fq}->[1] = $cf{$lr."_upstream2"};
#my %VECTOR = ('Yeast_BAC_pvuii_L' => ['GGGGCTCCCAC', 'CTGGCGAAAGG']);
exists($VECTOR{$fq}) || die;
my ($UPSTREAMBARCODE, $UPSTREAMSITE) = @{$VECTOR{$fq}};
my $BARCODEBRIDGELEN = 1; #$UPSTREAMBARCODE is 1bp from barcode becasue that bp has high mutation rate
                          #PvuII half-site tend to have 1bp mutation during ligation
$QUALITY !~ /^\d+$/ && die "SQ threshold as a number!";

my %RESTRICTIONSITE;
$RESTRICTIONSITE{$fq}->[0] = $cf{$lr."_cutsite1"};
$RESTRICTIONSITE{$fq}->[1] = $cf{$lr."_cutsite2"};
#my %RESTRICTIONSITE = ('Yeast_BAC_pvuii_L' => ['AAGCTT', 'GATCCG']); #HINDIII and BamHI with 1bp from barcode
my %BRIDGELEN = ($fq => [12, 1]);
exists($RESTRICTIONSITE{$fq}) || die;
exists($BRIDGELEN{$fq}) || die;
my %SQ;
$SQ{$fq}->[0] = $cf{$lr."_anchor1"};
$SQ{$fq}->[1] = $cf{$lr."_anchor2"};
#my %SQ = ('Yeast_BAC_pvuii_L' => [47, 32]);
$SQ{$fq}->[1] > length($VECTOR{$fq}->[1]) + $VAR || die;
$SQ{$fq}->[0] > length($VECTOR{$fq}->[0]) + $VAR || die;
my %FILE;
if($readlength <150){
	%FILE = ($fq => ["$fq$MERGE.fastq", $fq."2$MERGE.fastq"]);
}else{
	%FILE = ($fq => ["$fq$MERGE.fastq", $fq."2$MERGE.fastq", $fq."3$MERGE.fastq"]);
}
exists($FILE{$fq}) || die;
my ($sl, $slLen, $sr, $srLen) = ($SQ{$fq}->[0]-length($VECTOR{$fq}->[0])-$VAR, length($VECTOR{$fq}->[0]) + 2 * $VAR, $SQ{$fq}->[1] + $VAR, length($VECTOR{$fq}->[1]) + 2 * $VAR);
my $miniLen = $SQ{$fq}->[0] + $SQ{$fq}->[1] + $VAR;
my $CHOPLEN = length($RESTRICTIONSITE{$fq}->[0]) + $BRIDGELEN{$fq}->[0];
my $tmpFile = $OUTPUT.'noHindIII';
#my ($regx1) = (qr/^\s*$/);
exists($SQ{$fq}) || die;

open(NO, ">$tmpFile$MERGE.fastq") || die;
open(TMP, ">$OUTPUT$MERGE.fastq") || die;
foreach $file (@{$FILE{$fq}}) {
    $lowQ = 0;
$readNum = 0;
$readNum2 = 0;
$readNum3 = 0;
$readNum4 = 0;
$readNum5 = 0;
    open(FQ, $file) || die $file;
    while(<FQ>){
	push @lines,$_;
	#$_ =~ $regx1 && die;

	if(@lines == 4){
	    $lines[1] =~ /^[ACGTN]+\s*$/ || die $lines[1];
	    $lines[2] =~ /^\+\s*$/ || die $lines[2];

	    if(length($lines[1]) < $miniLen) {
		@lines = ();
		$readNum5++;
		next;
	    }

	    $Q = &hiseq::SQmean(substr($lines[3],$SQ{$fq}->[0], length($lines[3])-$SQ{$fq}->[0]-$SQ{$fq}->[1]));
	    ++$readNum =~ /00000$/ && print "$readNum\t$readNum3\n";
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
	    $file eq $fq."2$MERGE.fastq" ? ($rmL = &hiseq::rmLendBridge($lines[1], $VECTOR{$fq}->[0], $sl+length($lines[1])-$readlength, $slLen, 1, $BARCODEBRIDGELEN)) : ($rmL = &hiseq::rmLendBridge($lines[1], $VECTOR{$fq}->[0], $sl, $slLen, 1, $BARCODEBRIDGELEN));
	    if($rmL) {
		$file eq $fq."2$MERGE.fastq" ?($rmR = &hiseq::rmRendBridge($lines[1], $VECTOR{$fq}->[1], $sr+length($lines[1])-$readlength, $srLen, 1, $BARCODEBRIDGELEN)):($rmR = &hiseq::rmRendBridge($lines[1], $VECTOR{$fq}->[1], $sr, $srLen, 1, $BARCODEBRIDGELEN));

		if($rmR) {
		    if($lines[1] =~ s/^[ACGTN]{$rmL}([ACGTN]+)[ACGTN]{$rmR}$/$1/){
			    $lines[3] =~ s/^.{$rmL}(.+).{$rmR}$/$1/ || die "$lines[1]\n$rmL\n$rmR";
		    }else{
			    @lines=();
			    next;
		    }
		    #print "1:\n$lines[1], $RESTRICTIONSITE{$fq}->[0], 0, $BARCODELEN + $BRIDGELEN{$fq}->[0], $BARCODELENMAX + $BRIDGELEN{$fq}->[0], l\n";
		    length($lines[3]) > length($RESTRICTIONSITE{$fq}->[0])+$BRIDGELEN{$fq}->[0]+$BARCODELEN ? (@rm = &hiseq::rmSite($lines[1], $RESTRICTIONSITE{$fq}->[0], 0, $BARCODELEN + $BRIDGELEN{$fq}->[0], $BARCODELENMAX + $BRIDGELEN{$fq}->[0], 'l')) : ($rmR = 0);
		    #defined($rm[0])? print "$rm[0]\n$rm[1]\n\n":print"undefined\n\n";
		    $lines[1] =~ /^[ACGTN]+$/ || die;
		}
	    }else {
		$rmR = 0;
	    }

	    if(@rm == 2) {
		$lines[1] = $RESTRICTIONSITE{$fq}->[0].$rm[0];
		$lines[3] = substr($lines[3], 0 - length($lines[1]));

		if(length($lines[1]) > 20 && &SQend($lines[3], $END3,$QUALITYG, 3)) {
		    #print "3:\n$rm[1], $RESTRICTIONSITE{$fq}->[1], $BRIDGELEN{$fq}->[1], $BARCODELEN, $BARCODELENMAX, l\n";
		    @t = &hiseq::rmSite($rm[1], $RESTRICTIONSITE{$fq}->[1], $BRIDGELEN{$fq}->[1], $BARCODELEN, $BARCODELENMAX, 'l');
		    #defined($t[0])? print "$t[0]\n$t[1]\n\n":print"undefined\n\n";
		    if(@t == 2) {
			$t = $t[1];
		    } elsif($rm[1] =~ /^([ACGTN]{1,$BARCODELENMAX})/) {
			$t = $1;
		    } else {
			die "$lines[1]\n$rm[1]\n$file";
		    }

		    $barcode{$t}++;
		    $lines[0] =~ s/^([^\s]+)/$1BARCODE$t/ || die;
		    $lines[1] =~ /\s/ && die;
		    $lines[2] =~ /\s/ && die;
		    $lines[3] =~ /\s/ && die;

		    $readNum3++;
		    print TMP join("\n", @lines), "\n";
		} else {
#		    print "$lines[3]\n";
		    $readNum4++;
		}
	    } elsif($rmL && $rmR) {
		$lines[1] =~ s/^([ACGTN]{$BARCODELEN})[ACGTN]{$CHOPLEN}// || die;
		$t = $1;
		$lines[3] = substr($lines[3], 0 - length($lines[1]));

		if(length($lines[1]) > 20 && &SQend($lines[3], $END3,$QUALITYG, 3)) {
		    $barcodeBrokenHINDIII{$t}++;
		    $lines[0] =~ s/^([^\s]+)/$1BARCODE$t/ || die;
		    $lines[1] =~ /\s/ && die $lines[0], "*", $lines[1], "*$rmL*$rmR";
		    $lines[2] =~ /\s/ && die;
		    $lines[3] =~ /\s/ && die;
		    $readNum2++;
		    print NO join("\n", @lines), "\n";
		} else {
#		    print "$lines[3]\n";
		    $readNum4++;
		}
	    }

	    @lines = ();
	} elsif(@lines > 4) {
	    die @lines;
	}
    }

    print "$file\n\t$lowQ low quality barcodes, $readNum5 reads too short and $readNum4 ones with low quality non-barcode end among $readNum total reads\n\t$readNum3 reads with expected two ends and HindIII site\n\t$readNum2 reads with expected two ends but no HindIII site\n";
}

open(O, ">$OUTPUT$MERGE"."BarcodeWithoutHindIII.txt") || die;
@t = sort {$barcodeBrokenHINDIII{$b} <=> $barcodeBrokenHINDIII{$a}} keys %barcodeBrokenHINDIII;
for($i = 0; $i < @t; $i++) {
    print O "$t[$i]\t", $barcodeBrokenHINDIII{$t[$i]}, "\n";
}
close(O);

sub SQend($$$$) {
    my ($seq, $len, $threshold, $end) = @_;
    $len > 0 || die;
    $threshold || return 0;
    $threshold < 0 && die; #must be wrong parameter
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
