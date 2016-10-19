#For each set of barcodes whose insertion loci are close, identify barcodes differ from pre- or post-pool barcodes by no more than 3bp and convert them to the pooled ones
#1. Slide a 10kb window along chromosome 5kb/step, consider reads from _0.sam and _67.sam aligned to a window as hitting to same locus.
#2. Identify barcodes in each read set of same locus. If there are multiple barcodes occurred in pre- or post-pools, but differ by only 1 bp (i.e., in "barcodePrepoolLibrary-129similarbarcode[LR].txt"), the less enriched barcode is likely due to sequencing error so that convert it to the more enriched one.
#3. For barcodes in  each read set of same locus, barcodes differ by no more than 3bp from those in pools are likely due to sequencing error, correctify them in each chromosome locus
#predefine %BRIDGE, %INPUT, %SITE, %LEFTBARCODE
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
my ($DS2, $INPUT, $OUTPUT, $RANGE, $barcode, $barcodeSet, $bt, $chr, $count, $fastq, $i, $l, $location, $n, $num, $read, $r, $t, @barcode, @lines, %barcode, %barcodeInLibrary, %barcodeInLibraryLocation, %barcodeNotInLibraryLocation, %convert, %empty, %location, %locationRead, %repeat, %similar, %t, %tt) = (3, 'groupLocation', 'poolizeBarcode', 10000);
my %BRIDGE = ($FQ => [$cf{$lr."_bridge"}, $cf{$lr."_bridge"}]);
exists($BRIDGE{$FQ}) || die;
my %INPUT = ($FQ => ['']);
exists($INPUT{$FQ}) || die;
my $SITE = $cf{'inverse_PCR_cutsite'};
my %LEFTBARCODE = ();

foreach $t (@{$INPUT{$FQ}}) {
	print "$INPUT$t$FQ"."unique.txt\n";
    open(UNIQUE, "$INPUT$t$FQ".'unique.txt') || die "$INPUT$t$FQ".'unique.txt';
    while(<UNIQUE>) {
	/^([^\s]+)BARCODE([ACGTN]+).+\s([_\w\.]+)_(\d+)_\d+\s*$/ || die "$_\n";
	$read = $1;
	$barcode = $2;
	$chr = $3;
	$l = $4;
	%tt = ();

	if($l < 0) {
	    foreach $location ((int($l/$RANGE), int($l/$RANGE-0.5))) {
		$tt{$location} = 1;
	    }
	} else {
	    foreach $location ((int($l/$RANGE), int($l/$RANGE+0.5))) {
		$tt{$location} = 1;
	    }
	}
#	$barcode eq 'TTGAATCGGTTGTGTTGTGG' && print; #This barcode differs from much more enriched TTGTATCGGTTGTGTTGTGT by 2bp, but its left barcode mate is 2bp from standard one too
	foreach $location (keys %tt) {
	    $location{$chr.'_'.$location}->{$barcode}++;
	    $barcode{$barcode}->{$chr.'_'.$location}++;
	    $locationRead{$chr.'_'.$location}->{$barcode}->[@{$locationRead{$chr.'_'.$location}->{$barcode}}] = $read;
	}
    }
}
print scalar(keys %location), " genome loci\n";

open(TAG, "$POOL.txt") || die;
while(<TAG>) {
    /^([ACGTN]+)\:([ACGTN]+)\s+(\d+)\s*$/ || die "$_*";

    if($lr eq 'left') {
        $bt = $1;
        $t = $BRIDGE{$FQ}->[0].$bt.$BRIDGE{$FQ}->[1];
        $t =~ /$SITE/ && next;
    } else {
        $bt = $2;
        &hiseq::revcomDNA($BRIDGE{$FQ}->[0].$bt.$BRIDGE{$FQ}->[1]) =~ /$SITE/ && next;
        $bt = &hiseq::revcomDNA($bt);
    }

    exists($barcode{$bt}) && ($barcodeInLibrary{$bt} = $barcode{$bt});
}


%convert = ();
while(($location, $barcodeSet) = each(%location)) {
	#print "$location, $barcodeSet\n";
    %barcodeInLibraryLocation = ();
    %barcodeNotInLibraryLocation = ();

    while(($barcode, $num)  = each(%{$barcodeSet})) {
	exists($barcodeInLibrary{$barcode}) ? ($barcodeInLibraryLocation{$barcode} = $num) : ($barcodeNotInLibraryLocation{$barcode} = $num);
    }

    @barcode = keys %barcodeInLibraryLocation;
    while(($barcode, $num)  = each(%barcodeInLibraryLocation)) {
	while(($bt, $n)  = each(%barcodeNotInLibraryLocation)) {
            if(&hiseq::nomorethan1indelxS($barcode, $bt, $DS2)) {
#		if($num > $n || &hiseq::nomorethan1indelxS($barcode, $bt, 1)) {
		if($num > $n) {
		    exists($locationRead{$location}->{$bt}) || die;
		    foreach $read (@{$locationRead{$location}->{$bt}}) {
			if(exists($convert{$bt}->{$read})) {
			    if($convert{$bt}->{$read} ne $barcode) {
				if($num > 10 * $barcodeInLibraryLocation{$convert{$bt}->{$read}}) {
				    $convert{$bt}->{$read} = $barcode;
				} else {
				   # $barcodeInLibraryLocation{$convert{$bt}->{$read}} > 10 * $num || warn "$bt\t$barcode\t$read\t", $convert{$bt}->{$read},"\n";
				}
			    }
			} else {
			    $convert{$bt}->{$read} = $barcode;
			}
		    }
		}
	    }
	}
    }
}

%t = ();
while(($bt, $barcodeSet) = each(%convert)) {
    while(($read, $barcode) = each(%{$barcodeSet})) {
	if(exists($t{$bt})) {
	   # $t{$bt} eq $barcode || warn "$bt\n$t{$bt}\n$barcode\t$read";
	} else {
	    $t{$bt} = $barcode;
	}

	exists($barcodeInLibrary{$barcode}) || die;
    }
}
print scalar(keys %convert), " bracodes differ from pooled ones by up to $DS2 bp\n";

foreach $t (@{$INPUT{$FQ}}) {
    %repeat = ();
    open(UNIQUE, "$INPUT$t$FQ".'multiple.txt') || die;
    while(<UNIQUE>) {
	chomp;
	/^([^\s]+)(BARCODE)([ACGTN]+)(\s.+)$/ || die;
	$read = $1;
	$repeat{$read} = 1;
    }

    open(UNIQUE, "$INPUT$t$FQ".'unique.txt') || die;
    while(<UNIQUE>) {
	chomp;
	/^([^\s]+)(BARCODE)([ACGTN]+)(\s.+)$/ || die;
	$read = $1;
	exists($repeat{$read}) && delete($repeat{$read});
    }
    scalar(keys %repeat) && die scalar(keys %repeat);

    foreach $l (('unique', 'multiple')) {
	open(O, ">$OUTPUT$t$FQ$l".'.txt') || die;
	open(UNIQUE, "$INPUT$t$FQ$l".'.txt') || die "$OUTPUT$FQ$l".'.txt';
	    while(<UNIQUE>) {
		chomp;
		/^([^\s]+)(BARCODE)([ACGTN]+)(\s.+)$/ || die;
		$read = $1;
		$i = $2;
		$barcode = $3;
		$chr = $4;

		if(exists($convert{$barcode}->{$read})) {
		    exists($barcodeInLibrary{$convert{$barcode}->{$read}}) || die;
		    print O $read, $i, $convert{$barcode}->{$read}, &hiseq::strDistance($barcode, $convert{$barcode}->{$read}), $barcode, $chr, "\n";
		} else {
		    print O "$_\n";
		}
	    }
	close(O);
	close(UNIQUE);
    }
}
