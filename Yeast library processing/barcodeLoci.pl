#! usr/bin/perl
use strict;
use warnings;
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
my $FQ = $cf{$lr."_reads"};
my $POOL = $cf{'pool_file'};
my ($INPUT, $OUTPUT, $RANGE, $RANGELOCAL, $RANGESITE, $count, $barcode, $barcodeNum, $bt, $chr, $cluster, $i, $l, $loc, $location, $mate2, $num, $read, $t, %barcode, %barcodeHindIII, %barcodeLocation, %barcodeNum, %bl, %l, %pool, %t,%tt,  @cluster, @t) = ('groupBarcode', 'barcodeLoci', 10000, 1000, 20, 0);
my %INPUT = ($FQ => ['']);
exists($INPUT{$FQ}) || die;
my %BRIDGE = ($FQ => [$cf{$lr."_bridge"}, $cf{$lr."_bridge"}]);
exists($BRIDGE{$FQ}) || die;
my $SITE = $cf{'inverse_PCR_cutsite'};
$RANGE > $RANGELOCAL || die;

open(TAG, "$POOL.txt") || die;
while(<TAG>) {
    /^([ACGTN]+)\:([ACGTN]+)\s+(\d+)\s*$/ || die "$_*";
    $num = $3;

    if($lr eq 'left') {
	$bt = $1;
	$t = $BRIDGE{$FQ}->[0].$bt.$BRIDGE{$FQ}->[1];
	$t =~ /$SITE/ && next;
	$pool{$bt} = $num;
    } elsif($lr eq 'right') {
	$bt = $2;
	&hiseq::revcomDNA($BRIDGE{$FQ}->[0].$bt.$BRIDGE{$FQ}->[1]) =~ /$SITE/ && next;
	$pool{&hiseq::revcomDNA($bt)} = $num;
    }else{die;}
}
#print scalar(keys %pool),"\n";
open(I, $INPUT.$INPUT{$FQ}->[0].$FQ.'.txt') || die;
while(<I>) {
    chomp;
    s/^([ACGTN]+)// || die;
    $barcode = $1;
    %barcode = ();
    $count = 0;
    exists($barcodeNum{$barcode}) && die;
    exists($barcodeHindIII{$barcode}) && die;

    while(s/\t([a-z\d\.]+)(_rc)?_(\d+)LOCNUM(\d+)READNUM(\d+)//i) {
	$chr = $1;
	$2 && ($chr .= $2);
	$location = $3;
	$num = $5;
	exists($barcodeHindIII{$barcode}->{$chr}->{$location}) ? die : ($barcodeHindIII{$barcode}->{$chr}->{$location} = $num);
	$barcodeNum{$barcode} += $num;
	$count < $num && ($count = $num);
    }
    $_ && die "$_*";

    while(($chr, $location) = each(%{$barcodeHindIII{$barcode}})) {
	@t = sort {$a <=> $b} keys %{$location};
	$t = 0;

	    for($i = 1; $i < @t; $i++) {
		if($t[$i] - $t[$i-1] < $RANGESITE) {
		    $t++;
		    $location->{$t[$i]} < $location->{$t[$i-1]} && (($t[$i-1], $t[$i]) = ($t[$i], $t[$i-1]));
		    $barcodeHindIII{$barcode}->{$chr}->{$t[$i]} += $barcodeHindIII{$barcode}->{$chr}->{$t[$i-1]};
		    delete($barcodeHindIII{$barcode}->{$chr}->{$t[$i-1]});
		}
	    }

	$t > 5 && warn "HindIII-linked barcode $barcode should not spread $t loci along chromosome $chr local region\n"; #mostly 6bp separated 2 hits, but errors in repeats can result in many hits 
	@t = sort {$a <=> $b} keys %{$barcodeHindIII{$barcode}->{$chr}};
	$t = 0;

	for($i = 1; $i < @t; $i++) {
	    if($t[$i] - $t[$i-1] < $RANGE) {
		$t++;
		$barcodeHindIII{$barcode}->{$chr}->{$t[$i]} < $barcodeHindIII{$barcode}->{$chr}->{$t[$i-1]} && (($t[$i-1], $t[$i]) = ($t[$i], $t[$i-1]));
		$chr =~ /pcc/i || (exists($pool{$barcode}) && $barcodeHindIII{$barcode}->{$chr}->{$t[$i]} > $count/2 && $barcodeHindIII{$barcode}->{$chr}->{$t[$i-1]} > $count/2 && warn "should be rare that HindIII-linked inserts are less than $RANGE bp apart $barcode $chr $t[$i] ", $t[$i-1], "\t$count*", $barcodeHindIII{$barcode}->{$chr}->{$t[$i]}, "*\n");
		$barcodeHindIII{$barcode}->{$chr}->{$t[$i]} += $barcodeHindIII{$barcode}->{$chr}->{$t[$i-1]};
		delete($barcodeHindIII{$barcode}->{$chr}->{$t[$i-1]});
	    }
	}
	$chr =~ /^pcc/i || ($t > 2 && warn "HindIII-linked barcode $barcode should not spread $t loci along chromosome $chr\n");
    }
}

$count = 0;
open(NP, ">$OUTPUT$FQ".'NotInPool.txt') || die;
open(P, ">$OUTPUT$FQ".'InPoolunique.txt') || die;
open(N, ">$OUTPUT$FQ".'InPoolNonunique.txt') || die;

while(($barcode, $l) = each(%barcodeHindIII)) {
    exists($barcodeNum{$barcode}) ? ($barcodeNum = $barcodeNum{$barcode}) : die;
    @t = ();

    BL: while(($chr, $loc) = each(%{$l})) {
	    while(($location, $num) = each(%{$loc})) {
		if($num < $barcodeNum) {
		    unless($num < 3*($barcodeNum-$num)) {
			@t = ($chr.'_'.$location, $num);
			last BL;
		    }
		} elsif($num > $barcodeNum) {
		    die "$barcode\t$chr\t$location\t$num\t$barcodeNum";
		} elsif($num > 2) {
		    @t = ($chr.'_'.$location, $num);
		}
	    }
    }

    if(scalar(@t)) {
	if(exists($pool{$barcode})) {
	    print P "$barcode\t$t[0]NUM$t[1]\n";
	    $count++;
	} else {
	    print NP "$barcode\t$t[0]NUM$t[1]\n";
	}
    } else {
	%l = ();
	while(($chr, $location) = each(%{$l})) {
	    @t = sort {$a <=> $b} keys %{$location};
	    for($i = 0; $i < @t; $i++) {
		exists($l{$chr.'_'.$t[$i]}) ? die : ($l{$chr.'_'.$t[$i]} = $location->{$t[$i]});
	    }
	}

	@t = sort {$l{$b} <=> $l{$a}} keys %l;
	@t || die;

	if(exists($pool{$barcode})) {
	    if(@t > 1) {
		if($l{$t[1]} > 10 + $barcodeNum/100) {
		    print N "$barcode";

		    for($i = 0; $i < @t; $i++) {
			print N "\t$t[$i]NUM$l{$t[$i]}";
		    }
		    print N "\n";
		} elsif($l{$t[0]} > $barcodeNum-$l{$t[0]}) {
		    if($l{$t[0]}/$l{$t[1]} > $l{$t[1]}+1) {
			print P "$barcode\t$t[0]NUM$l{$t[0]}\n";
			$count++;
		    } elsif($l{$t[1]} > 2) {
			if(@t > 2 && $l{$t[2]} > 2) {
			    print N "$barcode\t$t[0]NUM$l{$t[0]}\t$t[1]NUM$l{$t[1]}\t$t[2]NUM$l{$t[2]}\n";
			} else {
			    print N "$barcode\t$t[0]NUM$l{$t[0]}\t$t[1]NUM$l{$t[1]}\n";
			}
		    } elsif($l{$t[0]} > 1) {
			print  "$barcode\t$t[0]NUM$l{$t[0]}\t$t[1]NUM$l{$t[1]} too few to judge among total $barcodeNum\n";
		    }
		} elsif($l{$t[0]} > 2) {
		    print N "$barcode";

		    for($i = 0; $i < @t; $i++) {
			print N "\t$t[$i]NUM$l{$t[$i]}";
		    }
		    print N "\n";
		}
	    } elsif($l{$t[0]} > 2) {
		print P "$barcode\t$t[0]NUM$l{$t[0]}\n";
		$count++;
	    }
	} else {
	    if(@t > 1) {
		if(2*$l{$t[0]} > $barcodeNum) {
		    if($l{$t[0]}/$l{$t[1]} > $l{$t[1]}+1) {
			print NP "$barcode\t$t[0]NUM$l{$t[0]}\n";
		    } elsif($l{$t[1]} > 2) {
			if(@t > 2 && $l{$t[2]} > 2) {
			    print NP "$barcode\t$t[0]NUM$l{$t[0]}\t$t[1]NUM$l{$t[1]}\t$t[2]NUM$l{$t[2]}\n";
			} else {
			    print NP "$barcode\t$t[0]NUM$l{$t[0]}\t$t[1]NUM$l{$t[1]}\n";
			}
		    }
		}
	    } elsif($l{$t[0]} > 2) {
		print NP "$barcode\t$l{$t[0]}\n";
	    }
	}
    }
}
print "$count pooled barcodes have unique genome inserts\n";
close(P);
close(NP);

if($FQ eq 'Library129_yeast_pvuii_digestion_first') {
%t = ();
open(P, "$OUTPUT$FQ".'InPoolunique.txt') || die;
while(<P>) {
    /^([ACGTN]+)\s+([a-z\d]+)(_rc)?_(\d+)NUM(\d+)\s*$/i || die "$_*";
    $barcode = $1;
    $chr = uc $2;
    $3 ? ($location = 0 - $4) : ($location = $4);
    $num = $5;
    exists($t{$barcode}) ? die : ($t{$barcode} = $chr);
}
print scalar(keys %t), " barcodes in $OUTPUT$FQ", "InPoolunique.txt\n";

open(P, "barcodeLocation".'.txt') || die;
while(<P>) {
    /^([ACGTN]+)\s+([a-z\d]+)/i || die;
    $barcode = $1;
    $chr = uc $2;

    if(exists($t{$barcode})) {
	if($t{$barcode} eq $chr) {
	    delete $t{$barcode};
	} else {
	    die "old $_", $t{$barcode};
	}
    } else {
	print "old $_";
    }
}

while(($barcode, $chr) = each(%t)) {
    print "$barcode\t$chr\n";
}
}
my (%detected,%m,$key,$value)=();
open(NP, "<$OUTPUT$FQ".'NotInPool.txt') || die;
open(P, "<$OUTPUT$FQ".'InPoolunique.txt') || die;
open(N, "<$OUTPUT$FQ".'InPoolNonunique.txt') || die;
open(M, $INPUT.$INPUT{$FQ}->[0].$FQ.'multiple.txt') || die;
open(MO, ">$OUTPUT$FQ".'InPoolMultiple.txt') || die;
while(<P>){
	$_ =~ /^([ACGT]+)\t\w+/ || die;
	defined($detected{$1}) ? die:($detected{$1} = 1);
}
while(<N>){
	$_ =~ /^([ACGT]+)\t\w+/ || die;
	defined($detected{$1}) ? die:($detected{$1} = 1);
}
while(<M>){
	$_ =~ /BARCODE([ACGT]+)/ || die $_;
	if(defined($pool{$1})){
		if(defined($detected{$1})){
		}else{
			$m{$1}++;
		}
	}
}
close NP;
close P;
close N;
close I;
close M;
while(($key,$value)=each %m){
	print MO "$key	$value\n";
}
close MO;
