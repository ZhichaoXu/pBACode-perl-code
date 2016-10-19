#Calculate the distance between unique inserts of each left barcode and its corresponding right barcode
#OUTPUT: barcodepairDistance*txt: distance between concordant left and right barcodes along genome
#        barcodepairWrong: location of discordant barcodepairs, i.g. either in different chromsomes or in wrong direction
#! usr/bin/perl
use strict;
use warnings;
use List::Util qw(min);
use List::Util qw(max);
require hiseq;

my $config = $ARGV[0];
my %cf;
open(CF, $config) || die "$config\n";
while(<CF>){
	if($_ =~ /(\S+)\t(\S+)/){
		$cf{$1} = $2;
	}
}
my $FILELEFT = $cf{"left_reads"};
my $FILERIGHT = $cf{"right_reads"};
my $POOL = $cf{'pool_file'};
my $fa = $cf{"genome"};
my $CHRLEN = $fa.'.Len';
my($INPUTPREFIX, $INPUTSURFIX, $OUTPUT, $OUTPUTWRONG, $barcode, $chr, $file, $i, $j, $leftBarcode, $location, $num, $output, $rightBarcode, $t, %jump, %pair, %wrongPair, %CHRLEN, @FILE, @barcode, @exist, @pooledbarcode, @t,$left,$right,$tmp) = ('barcodeLoci', 'InPoolunique', 'barcodepairDistance', 'barcodepairWrong');

$FILE[0] = [split(/\+/, $FILELEFT)];
$FILE[1] = [split(/\+/, $FILERIGHT)];
$output = '';
for($i = 0; $i < @FILE; $i++) {
    for($j = 0; $j < @{$FILE[$i]}; $j++) {
	$output .= $FILE[$i]->[$j];
    }
}
#print "$output\n";
open(I, "$CHRLEN.txt") || die;
while(<I>) {
    /^([\w_\d\.]+)\s(\d+)\s*$/ || die "$_\n";
    exists($CHRLEN{uc $1}) ? die : ($CHRLEN{uc $1} = $2);
}
#print %CHRLEN;
for($i = 0; $i < @FILE; $i++) {
    foreach $file (@{$FILE[$i]}) {
	open(I, "$INPUTPREFIX$file$INPUTSURFIX.txt") || die;
	print "$INPUTPREFIX$file$INPUTSURFIX.txt\n";
	while(<I>) {
	    /^([ACGT]+)\s+(\w+_?\d+\.?\d*)(_rc)?_(\d+)NUM(\d+)$/i || die "$_\n";
	    $barcode = $1;
#	    exists($pooledbarcode[1]->{$barcode}) || die;
	    $chr = uc $2;
	    $location = $4;
	    $num = $5;
	    #print "\n$chr\n";
	    if($3) {
		if($chr =~ /PCC/) {
		    $location = $location - 9000;
		} elsif(exists($CHRLEN{$chr})) {
		    $location = $location - $CHRLEN{$chr} - 1;
		} else {
		    die $chr;
		}
	    }

	    if(exists($barcode[$i]->{$barcode})) {
		$barcode[$i]->{$barcode} =~ /^(\w+\.)_(\-?\d+)NUM(\d+)$/ || die;

		if($chr ne $1) {
		    die "$barcode\t$chr\t$1";
		} elsif(abs($location - $2) > 1000) {#some barcodes in 'Yeast_BAC_pvuii_R' have no HinIII-side read-mate
		    die "$_$location\t", $barcode[$i]->{$barcode}, "\n$i";
		} elsif($2 < 0) {
		    $t = max($location, $2);
		} elsif($2 > 0) {
		    $t = min($location, $2);
		} else {
		    die;
		}

		$barcode[$i]->{$barcode} = $chr.'_'.$t.'NUM'.($num+$3);
	    } else {
		$barcode[$i]->{$barcode} = $chr.'_'.$location.'NUM'.$num;
	    }
	}
    }
}

print scalar(keys %{$barcode[0]}), " left barcodes\n";
print scalar(keys %{$barcode[1]}), " right barcodes\n";
#my ($key,$value)=();
#while(($key,$value)= each %{$barcode[1]}){
#	print "$key	$value\n";
#}
%wrongPair = ();
%pair = ();
open(TAG, "$POOL.txt") || die;
open SINGLE,">barcodepairSingle$output$POOL.txt"||die;
while(<TAG>) {
    /^([ACGTN]+)\:([ACGTN]+)\s+(\d+)\s*$/ || die "$_*";
    $leftBarcode = $1;
    $rightBarcode = &hiseq::revcomDNA($2);
    $pooledbarcode[1]->{$leftBarcode}->[@{$pooledbarcode[1]->{$leftBarcode}}] = $rightBarcode;
    $pooledbarcode[0]->{$rightBarcode}->[@{$pooledbarcode[0]->{$rightBarcode}}] = $leftBarcode;

    if(exists($barcode[0]->{$leftBarcode}) && exists($barcode[1]->{$rightBarcode})) {
#	print "$barcode[0]->{$leftBarcode}	$barcode[1]->{$rightBarcode}\n";
#	exists($exist[0]->{$leftBarcode}) ? die "$_*" : ($exist[0]->{$leftBarcode} = 1);
#	exists($exist[1]->{$rightBarcode}) ? die "$_*" : ($exist[1]->{$rightBarcode} = 1);
	   # print "$barcode[0]->{$leftBarcode}\n";
	    $barcode[0]->{$leftBarcode} =~ /^([\w_\d\.]+)_(\-?\d+)NUM(\d+)$/ || die;
	    $chr = uc $1;
	    $location = $2;
	    $num = $3;
	    $barcode[1]->{$rightBarcode} =~ /^([\w_\d\.]+)_(\-?\d+)NUM(\d+)$/ || die;
	    if($chr ne uc $1 || ($location * $2 > 0)) {
#		print "wrong\n";
		$wrongPair{"$leftBarcode $rightBarcode"} = $barcode[0]->{$leftBarcode}.' '.$barcode[1]->{$rightBarcode};
	    } else {
		$pair{"$leftBarcode $rightBarcode"} = $location + $2;
#		print $pair{"$leftBarcode $rightBarcode"},"\n";
		if($pair{"$leftBarcode $rightBarcode"} > 0){
			delete $pair{"$leftBarcode $rightBarcode"};
			$wrongPair{"$leftBarcode $rightBarcode"} = $barcode[0]->{$leftBarcode}.' '.$barcode[1]->{$rightBarcode};
			next;
		}
		($left,$right) = sort {$a<=>$b} (abs($location),abs($2));
		#($left,$right) = ($location,$2);
		$jump{"$leftBarcode $rightBarcode"} = $chr."\t".$left."\t".$right;
	    }
    }elsif(exists($barcode[0]->{$leftBarcode}) || exists($barcode[1]->{$rightBarcode})){
	if(exists($barcode[0]->{$leftBarcode})){
		$tmp = $barcode[0]->{$leftBarcode};
		$tmp =~ s/NUM\d+//;
		$tmp =~ s/_/\t/;
		print SINGLE "$leftBarcode	$tmp\n";
	}elsif(exists($barcode[1]->{$rightBarcode})){
		$tmp = $barcode[1]->{$rightBarcode};
		$tmp =~ s/NUM\d+//;
		$tmp =~ s/_/\t/;
		print SINGLE "$rightBarcode	$tmp\n";
	}else{die;}
    }
}
close SINGLE;
print scalar keys %wrongPair, " barcode-pairs are in different chromosomes or discordant directions\n";
print scalar keys %pair, " barcode-pairs are in same chromsomes and concordant\n";


open(O, ">$OUTPUT$output$POOL.txt") || die;
@t = sort {$pair{$a} <=> $pair{$b}} keys %pair;
for($i = 0; $i < @t; $i++) {
    print O "$t[$i]\t",$jump{$t[$i]},"\t",0-$pair{$t[$i]}, "\n";
}
close(O);

open(O, ">$OUTPUTWRONG$output$POOL.txt") || die;
while(($barcode, $location) = each(%wrongPair)) {
    print O "$barcode\t$location\n";
}
close(O);

for($i = 0; $i < @pooledbarcode; $i++) {
    while(($rightBarcode, $t) = each(%{$pooledbarcode[$i]})) {
	if(@{$t} > 1) {
	    $num = 0;

	    foreach $leftBarcode (@{$t}) {
		exists($barcode[$i]->{$leftBarcode}) && $num++;
	    }

            $num > 1 && warn "misleading $i barcode mate $rightBarcode because it is 1:N and multiple mate partnes exist in BAC\n"; 
	}
    }
}
