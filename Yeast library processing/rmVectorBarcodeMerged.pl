#remove vector sequences at both ends of merged read pairs of Right Barcode-mate which are separated from genome inserts only by HindIII site. and its bridging sequnce
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
	$lr eq 'jumpright' || die;
}else{
	$lr = 'right';
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
#print $VECTOR{$fq}->[0],"	",$VECTOR{$fq}->[1],"\n";
exists($VECTOR{$fq}) || die;
my ($UPSTREAMBARCODE, $UPSTREAMPVUII) = @{$VECTOR{$fq}};
my $BRIDGELEN = 2;
my $BARCODEBRIDGELEN = 1; #$UPSTREAMBARCODE is 1bp from barcode becasue that bp has high mutation rate
                          #PvuII half-site tend to have 1bp mutation during ligation
$QUALITY !~ /^\d+$/ && die "SQ threshold as a number!";
my %RESTRICTIONSITE = ($fq => $cf{$lr."_cutsite"}); #HINDIII
exists($RESTRICTIONSITE{$fq}) || die;
my %SQ;
$SQ{$fq}->[0] = $cf{$lr."_anchor1"};
$SQ{$fq}->[1] = $cf{$lr."_anchor2"};
#print $SQ{$fq}->[0],"	",$SQ{$fq}->[1],"\n";
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
my $SITE = $RESTRICTIONSITE{$fq};
my $SITELEN = length($SITE);
my $CHOPLEN = length($SITE) + $BRIDGELEN;
my $tmpFile = $OUTPUT.'noHindIII';
#my ($regx1) = (qr/^\s*$/);
exists($SQ{$fq}) || die;

open(NO, ">$tmpFile$MERGE.fastq") || die;
open(TMP, ">$OUTPUT$MERGE.fastq") || die;
foreach $file (@{$FILE{$fq}}) {
    print "$file\n";
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
	    $lines[2] =~ /^\+\s*$/ || die $lines[2];

	    if(length($lines[1]) < $miniLen) {
#		$file eq $fq."3$MERGE.fastq" || die "$miniLen\n$file\n", @lines;
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
	    #print "$lines[1], $VECTOR{$fq}->[0], ",length($lines[1])+$sl-$readlength,", $slLen, 1, $BARCODEBRIDGELEN\n" if($file eq $fq."2$MERGE.fastq");
	    $file eq $fq."2$MERGE.fastq" ? ($rmL = &hiseq::rmLendBridge($lines[1], $VECTOR{$fq}->[0], length($lines[1])+$sl-$readlength, $slLen, 1, $BARCODEBRIDGELEN)) : ($rmL = &hiseq::rmLendBridge($lines[1], $VECTOR{$fq}->[0], $sl, $slLen, 1, $BARCODEBRIDGELEN));
	    #print "$rmL\n" if($file eq $fq."2$MERGE.fastq");
	    if($rmL) {
		#$rmL >= $sl && $rmL <= ( $sl + $slLen + $BARCODEBRIDGELEN) || die "$lines[1]\t$rmL";
		#print "$lines[1], $VECTOR{$fq}->[1], ",$sr+length($lines[1])-$readlength,", $srLen, 1, $BARCODEBRIDGELEN\n"if($file eq $fq."2$MERGE.fastq");
		$file eq $fq."2$MERGE.fastq" ? ($rmR = &hiseq::rmRendBridge($lines[1], $VECTOR{$fq}->[1], $sr+length($lines[1])-$readlength, $srLen, 1, $BARCODEBRIDGELEN)):($rmR = &hiseq::rmRendBridge($lines[1], $VECTOR{$fq}->[1], $sr, $srLen, 1, $BARCODEBRIDGELEN));
		#print length($lines[1]),"	$rmL	$rmR\n"if($file eq $fq."2$MERGE.fastq");
		#print length($lines[1])-$rmL-$rmR,"\n"if($file eq $fq."2$MERGE.fastq");
		if($rmR) {
		    #$rmR <= $sr + $BARCODEBRIDGELEN && $rmR >= ( $sr - $srLen) || die "$lines[1]\t$rmR";
		    if($lines[1] =~ s/^[ACGTN]{$rmL}([ACGTN]+)[ACGTN]{$rmR}$/$1/){
			    $lines[3] =~ s/^.{$rmL}(.+).{$rmR}$/$1/ || die "$lines[1]\n$rmL\n$rmR";
		    }else{
			    @lines=();
			    next;
		    }
		    #print "$lines[1]\n"if($file eq $fq."2$MERGE.fastq");
		    length($lines[3]) > length($SITE)+$BRIDGELEN+$BARCODELEN ? (@rm = &hiseq::rmSite($lines[1], $SITE, $BRIDGELEN, $BARCODELEN, $BARCODELENMAX, 'l')) : ($rmR = 0);
		    #print "$lines[1]\n$rm[0]	$rm[1]\n\n"if($file eq $fq."2$MERGE.fastq" && defined($rm[0]));
		    $lines[1] =~ /^[ACGTN]+$/ || die;
		}
	    }else{
		    $rmR = 0;
	    }
# elsif($lines[1] =~ s/^([ACGTN]+$UPSTREAMBARCODE)([ACGTN]+)($UPSTREAMPVUII[ACGTN]+)$/$2/) {
#		$rmL = length($1);
#		$rmR = length($3);
#		$lines[3] =~ s/^.{$rmL}(.+).{$rmR}$/$1/ || die $lines[1];
#		length($lines[3]) > length($SITE)+$BRIDGELEN+$BARCODELEN ? (@rm = &hiseq::rmSite($lines[1], $SITE, $BRIDGELEN, $BARCODELEN, $BARCODELENMAX, 'l')) : ($rmR = 0);
#		$lines[1] =~ /^[ACGTN]+$/ || die;
#	    } else {
#		$rmR = 0;
#	    }

	    if(@rm == 2) {
		$lines[1] = $SITE.$rm[0];
		$lines[3] = substr($lines[3], 0 - length($lines[1]));

		if(length($lines[1]) > 20 && &SQend($lines[3], $END3,$QUALITYG, 3)) {
		    $barcode{$rm[1]}++;
		    $lines[0] =~ s/^([^\s]+)/$1BARCODE$rm[1]/ || die;
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

    print "$lowQ low quality barcodes, $readNum5 reads too short and $readNum4 ones with low quality non-barcode end among $readNum total reads\n$readNum3 reads with expected two ends and HindIII site\n$readNum2 reads with expected two ends but no HindIII site\n";
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
