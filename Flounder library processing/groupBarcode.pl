#Try to identify reads aligned to a repeat unit with unique polymorphism so that it has a unique genome locus
#If a read hits to multiple loci (_256.sam, _323.sam), find their Mapping Quality scores. If the MQ score in _0.sam/_67.sam is 7 more than all of those in _256.sam/_323.sam, suggesting the hit in _0.sam/_67.sam has at least one high QC nucleotide matching to genome so that consider it is significantly better than other hits so as to discard those in  _256.sam/_323.sam
#INPUT: poolizeBarcode*.txt
#OUTPUT: groupbarcode*uniqueORmultiple.txt: read\tlocations with indistinguishable mapping quality
#        groupbarcode*.txt: barcode\tlocations of independent reads, 'LOCNUM', number of independent reads, i.e., different startpoint or length (if two reads start at same location and their lengths differ only by 1, considering them as non-independent reads), 'READNUM', numbers of reads mapped to this location
#! usr/bin/perl
use strict;
use warnings;
require hiseq;
sub bestMappingQuality($$$); #INPUT: $1: mapping score of the best hit(_0.sam or _67.sam)
                             #       $2 mapping score of addtional hit(_256.sam or _323.sam)
                             #       $3 threshold of Mapping Score differece that read1 has significant mapping quality
                             #die if $1 < $2
                             #OUTPUT: 1 if $1 > $2 + 7 and neither mate of read1 is significantly worse than that of read2
                             #        0 if read1 is not signifcantly better than read2
sub betterMappingQuality($$$); #INPUT: $1 and $2 mapping scores of read1 and read2, $3 threshold of Mapping Score differece that read1 has significantly better mapping quality
                               #OUTPUT: 1 if $1 > $2 + 7
                               #        0 if read1 is not signifcantly better than read2
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
my ($DS, $INPUT, $OUTPUT, $barcode, $count, $len, $loc, $location, $mq, $num, $read, $t, $unique, @t, %barcode, %mq,%multiple,%mqbarcode,%t,%repeat) = (7, 'poolizeBarcode', 'groupBarcode');

open(MULTIPLE, "$INPUT$FQ".'multiple.txt') || die "$OUTPUT$FQ".'multiple.txt';
while(<MULTIPLE>) {
    chomp;
    s/^([^\s]+BARCODE[ACGTN]+\d?[ACGTN]*)// || die;
    $read = $1;
    $mq = '';
    exists($repeat{$read}) && die;

    while(s/\s+(\d+_?\d*)\s+([_\w\.]+)//) {
	$t = $1;
	$location = $2;

	if($mq) {
	    &betterMappingQuality($t, $mq, $DS) && die;
	    &betterMappingQuality($mq, $t, $DS) && last;
	} else {
	    $mq = $t;
	}

	$repeat{$read}->[@{$repeat{$read}}] = $location;
    }

    exists($mq{$read}) ? die : ($mq{$read} = $mq);
}
print scalar keys %repeat, " reads have multiple hits\n";

%multiple = ();
%mqbarcode = ();
$count = 0;
open(UNIQUE, "$INPUT$FQ".'unique.txt') || die "$OUTPUT$FQ".'unique.txt';
open(O, ">$OUTPUT$FQ".'multiple.txt') || die;
while(<UNIQUE>) {
    chomp;
    /^([^\s]+BARCODE[ACGTN]+\d?[ACGTN]*)\s(\d+_?\d*)\s([_\w\.]+)*$/ || die "$_";
    $read = $1;
    exists($t{$read}) && die;
    $mq = $2;
    $location = $3;
    $unique = 1;

    if(exists($repeat{$read})) {
	exists($mq{$read}) || die;
	if($mq =~ /^\d+$/) {
	    $mq{$read} =~ /^\d+$/ || die;
	} elsif($mq =~ /^\d+_\d+$/) {
	    $mq{$read} =~ /^\d+_\d+$/ || die;
	} else {
	    die;
	}

	unless(&bestMappingQuality($mq, $mq{$read}, $DS)) {
	    $unique = 0;
	    $count++;
	    $read =~ /BARCODE([ACGTN]+)/ || die;
	    exists($multiple{$1}) ? ($multiple{$1} .= "\n$_") : ($multiple{$1} = $_);

	    foreach $location (@{$repeat{$read}}) {
		$multiple{$1} .= "\t$location";
	    }
	}

	delete($repeat{$read});
	delete($mq{$read});
    }

    if($unique) {
	$read =~ /BARCODE([ACGTN]+)/ || die;
	$barcode{$1}->{$location}++;
    }
}
scalar keys %repeat && die scalar keys %repeat;
scalar keys %mq && die;
print "$count reads have significantly better MQ in hit than in other hits, suggesting it hits to a copy of reapts with unique polymorphism\n";
$t  =0;
foreach $barcode (keys %barcode) {
    if(exists($multiple{$barcode})) {
	delete($multiple{$barcode});
	$t++;
    }
}
print "$t barcodes have both unique-genome-locus and multiple-loci reads\n";
print scalar(keys %barcode), " barcodes have unique hits while ", scalar(keys %multiple), " do multiple\n";
foreach (values %multiple) {
    print O "$_\n";
}
close(O);

open(OO, ">$OUTPUT$FQ".'.txt') || die;
while(($barcode, $location) = each(%barcode)) {
    print OO $barcode;
    %t = ();
    while(($t, $num) = each(%{$location})) {
	$t =~ s/(\d+)_(\d+)$/$1/ || die;
	exists($t{$t}->{$2}) ? die : ($t{$t}->{$2} = $num);
    }

    while(($loc, $len) = each(%t)) {
	@t = keys %{$len};

	if(@t == 2) { #if two reads start at same location and their lengths differ only by 1, considering them as non-independent reads
	    abs($len->{$t[0]} - $len->{$t[1]}) > 1 ? ($t = scalar(@t)) : ($t = 1);
	} else {
	    $t = scalar(@t);
	}

	print OO "\t$loc", 'LOCNUM', $t;

	$t = 0;
	foreach $num (values %{$len}) {
	    $t += $num;
	}

	print OO "READNUM$t";
    }
    print OO "\n";
}

sub bestMappingQuality($$$) {
    my ($mq1, $mq2, $threshold) = @_;
    my ($m11, $m12, $m21, $m22);

    if($mq1 =~ /^\d+$/ && $mq2 =~ /^\d+$/) {
	$mq1 < $mq2 && die "$mq1,$mq2";
	$mq1 - $mq2 < $threshold ? return 0 : return 1;
    } elsif($mq1 =~ /^(\d+)_(\d+)$/) {
	$m11 = $1;
	$m12 = $2;
	$mq1 eq $mq2 && return 0;
	$mq2 =~ /^(\d+)_(\d+)$/ || die;
	$m21 = $1;
	$m22 = $2;
	$m11+$m12 < $m21+$m22 && die "$mq1,$mq2";
        ($m11+$m12 - $m21-$m22 < $threshold) ? return 0 : return 1;
    }
}

sub betterMappingQuality($$$) {
    my ($mq1, $mq2, $threshold) = @_;
    my ($m11, $m12, $m21, $m22);

    if($mq1 =~ /^\d+$/ && $mq2 =~ /^\d+$/) {
	$mq1 - $mq2 < $threshold ? return 0 : return 1;
    } elsif($mq1 =~ /^(\d+)_(\d+)$/) {
	$m11 = $1;
	$m12 = $2;
	$mq1 eq $mq2 && return 0;
	$mq2 =~ /^(\d+)_(\d+)$/ || die "$mq1,$mq2";
	$m21 = $1;
	$m22 = $2;
        ($m11+$m12 - $m21-$m22 < $threshold) ? return 0 : return 1;
    }
}
