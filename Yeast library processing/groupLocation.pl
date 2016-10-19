#Location inserts linked to barcodes in genome. discard low mapping quality reads, especially those with multiple hits
#'--norc' and 'genome in both W and C directions' are used in bowtie because, if a query is aligned to genome sequence with a few mismatches at end, W and C directions can give different aligment patterns and scores, i.e., \dS in one direction while discarding end in the other direction.
#mate2 is temperally  reverse complemented before bowtie because '--ff' has to be used due to '--norc'
#INPUT: $1: fastqfile, predefined %SINGLE: merged or non-barcode-side read mate
#                      predefined %PAIRED: uncombined read pair
#                      @NOHINDIII: 1bp mismatch at 5' end of non-unique-match reads with HindIII is tolerable
#                                  7bp mismatch at 5' end of non-unique-match reads without HindIII is tolerable
#                      %NOHINDIII: 2bp mismatch at 5' end of non-unique-match non-HindIII-site read mates is tolerable
#       $2: genome combining the orginal genome sequence and the reverse complemented sequence
#       $3: vector name
#OUTPUT: "*unique.txt": hits from _0.sam or _67.sam
#        "*multiple.txt": hits from _256.sam or _323.sam
#! usr/bin/perl
use strict;
use warnings;
require hiseq;
sub cutawk($);
sub chrLocationCW($%);
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
my $GENOME = $cf{'genome'}.'WC'; #S288CpccFOS2WC
my $VECTOR = $cf{'vector'}; #pccFOS2
my ($OUTPUT, $barcode, $count, $fastq, $i, $j, $location, $t, @lines, @t, %genome, %location, %locationRead, %vector,$sam) = ('groupLocation');
my %SINGLE;
if($lr =~ /jump/){
	%SINGLE = ($FQ => ['.extendedFrags']);
}else{
	%SINGLE = ($FQ => ['.extendedFrags', 'noHindIII.extendedFrags']);
}
exists($SINGLE{$FQ}) || die;
my %PAIRED;
if($lr =~ /jump/){
	%PAIRED = ($FQ => ['']);
}else{
	%PAIRED = ($FQ => ['', 'noHindIII']);
}
my %SINGLEINPUT = ($FQ => 'rmVectorBarcode');
my %PAIREDINPUT = ($FQ => 'rmVectorBarcode1mate');
my %NOHINDIII;
if($lr =~ /jump/){
	%NOHINDIII = ($FQ => [7]);
}else{
	%NOHINDIII = ($FQ => [1, 7]);
}
@{$SINGLE{$FQ}} != @{$NOHINDIII{$FQ}} && die;

#my %CHRLEN = ();
#open(I, $CHRLEN.'.txt') || die;
#while(<I>) {
#    /^(chr\w+)\s(\d+)\s*$/ || die;
#    exists($CHRLEN{uc $1}) ? die : ($CHRLEN{uc $1} = $2);
#}

open(UNIQUE, ">$OUTPUT$FQ".'unique.txt') || die;
open(MULTI, ">$OUTPUT$FQ".'multiple.txt') || die;
for($j = 0; $j < @{$NOHINDIII{$FQ}}; $j++) {
    $fastq = $SINGLE{$FQ}->[$j];

	print "bowtie2-align -x ./$GENOME $SINGLEINPUT{$FQ}$FQ$fastq.fastq --local --norc -k 10 --quiet --no-hd -S $OUTPUT$FQ$fastq$GENOME.sam\n";
	system "bowtie2-align -x ./$GENOME $SINGLEINPUT{$FQ}$FQ$fastq.fastq --local --norc -k 10 --quiet --no-hd -S $OUTPUT$FQ$fastq$GENOME.sam";
	&cutawk("$OUTPUT$FQ$fastq$GENOME");
	#unlink "$OUTPUT$FQ$fastq$GENOME.sam";

#   open(I, "$OUTPUT$FQ$fastq$GENOME".'_0.sam') || die;
#	while(<I>) {
#	    @lines = split;
#	    $lines[0] =~ /BARCODE([ACGTN]+)$/ || die;
#	    $barcode = $1;
#
#	    if($lines[2] =~ /^$VECTOR/) {
#		$vector{$barcode}++;
#	    } else {
#		$genome{$barcode}++;
#	    }
#	}
#    close(I);
    %locationRead = &hiseq::locationReadMulti("$OUTPUT$FQ$fastq$GENOME", 0, $NOHINDIII{$FQ}->[$j], 1);
    print "$OUTPUT$FQ$fastq$GENOME\t", scalar keys %locationRead, " merged or single reads\n";

    while(($barcode, $location) = each(%locationRead)) {
#	print $barcode,"\t";
#	$i = 0;
#	while(defined($location->[$i])){
#		print $location->[$i],"\n";
#		$i++;
#	}
	#print scalar(@{$location}),"\n";
	if(scalar(@{$location})> 1) {
	    print MULTI $barcode;

	    for($i = 0; $i < @{$location} - 1; $i++) {
#		if($location->[$i] =~ s/^([\d_]+\t)(c[a-z\d]+_rc)/$2/i) {
#		    print MULTI "\t$1", &chrLocationCW($location->[$i], %CHRLEN);
#		} else {
#		    $location->[$i] =~ /^([\d_]+\t)(p|c[a-z\d]+_\d+)/i || die $location->[$i];
		    print MULTI "\t", $location->[$i];
#		}
	    }

	    print MULTI "\n";
	}

	$i = @{$location}-1;
	#print "$i\n";
#	if($location->[$i] =~ s/^([\d_]+\t)(c[a-z\d]+_rc)/$2/i) {
#	    print UNIQUE "$barcode\t$1",  &chrLocationCW($location->[$i], %CHRLEN), "\n";
#	} else {
#	    $location->[$i] =~ /^([\d_]+\t)(p|c[a-z\d]+_\d+)/i || die $location->[$i];
	    print UNIQUE "$barcode\t", $location->[$i], "\n";
#	}
    }

exists($PAIRED{$FQ}) || next;
    exists($PAIREDINPUT{$FQ}) || die;
@{$PAIRED{$FQ}} != @{$NOHINDIII{$FQ}} && die;
    $fastq = $PAIRED{$FQ}->[$j];
    open(FQ, "$PAIREDINPUT{$FQ}$FQ$fastq"."2.fastq") || die;
    open(TMP, ">$OUTPUT$PAIREDINPUT{$FQ}$FQ$fastq.fastq") || die;

    @lines = ();
    while(<FQ>) {
	push @lines,$_;

	if(@lines == 4){
		chomp($lines[1]);
		chomp($lines[3]);
		$lines[1] = &hiseq::revcomDNA($lines[1])."\n";
		$lines[3] = reverse($lines[3])."\n";
		print TMP @lines;
	    $lines[0] =~ /BARCODE([ACGTN]+)/ || die $lines[0];
	    @lines = ();
	} elsif(@lines > 4) {
	    die scalar @lines;
	}
    }
    close(TMP);
    print "bowtie2-align -x ./$GENOME -1 $PAIREDINPUT{$FQ}$FQ$fastq"."1.fastq -2 $OUTPUT$PAIREDINPUT{$FQ}$FQ$fastq.fastq --local --ff --norc -k 10 --quiet --no-hd -S $OUTPUT$FQ$fastq$GENOME.sam\n";
    system "bowtie2-align -x ./$GENOME -1 $PAIREDINPUT{$FQ}$FQ$fastq"."1.fastq -2 $OUTPUT$PAIREDINPUT{$FQ}$FQ$fastq.fastq --local --ff --norc -k 10 --quiet --no-hd -S $OUTPUT$FQ$fastq$GENOME.sam";
    #unlink "$OUTPUT$PAIREDINPUT{$FQ}$FQ$fastq.fastq";
    &cutawk("$OUTPUT$FQ$fastq$GENOME");
    #unlink "$OUTPUT$FQ$fastq$GENOME.sam";
#    print "../Yeast_BAC_library/bowtie2-2.1.0/bowtie2-align -x ./$GENOME -1 $PAIREDINPUT{$FQ}$FQ$fastq"."1.fastq -2 $OUTPUT$PAIREDINPUT{$FQ}$FQ$fastq.fastq --local --ff --norc -k 10 --no-hd -S $OUTPUT$FQ$fastq$GENOME.sam\n";
#    open(I, "$OUTPUT$FQ$fastq$GENOME".'_67.sam') || die "$OUTPUT$FQ$fastq$GENOME".'_67.sam';
#	while(<I>) {
#	    @lines = split;
#	    $lines[0] =~ /BARCODE([ACGTN]+)$/ || die;
#	    $barcode = $1;
#
#	    if($lines[2] =~ /^$VECTOR/) {
#		$vector{$barcode}++;
#	    } else {
#		$genome{$barcode}++;
#	    }
#	}

    foreach $sam ((65, 67)) {
	%locationRead = &hiseq::locationReadMulti("$OUTPUT$FQ$fastq$GENOME", $sam, $NOHINDIII{$FQ}->[$j], 1000);
	print "$OUTPUT$FQ$fastq$GENOME\t", scalar keys %locationRead, " mated reads in $sam.sam\n";

	while(($barcode, $location) = each(%locationRead)) {
	    if(@{$location} > 1) {
		print MULTI $barcode;

		for($i = 0; $i < @{$location} - 1; $i++) {
		    print MULTI "\t", $location->[$i];
		}

		print MULTI "\n";
	    }

	    $i = @{$location}-1;
	    print UNIQUE "$barcode\t", $location->[$i], "\n";
	}
     }
}
close(UNIQUE);
close(MULTI);


#open(O, ">barcodeVectorGenomeHitNum.txt") || die;
#@t = sort {$vector{$b} <=> $vector{$a}} keys %vector;
#for($i = 0; $i < @t; $i++) {
#    print O "$t[$i]\t", $vector{$t[$i]};
#    exists($genome{$t[$i]}) && print O "\t$genome{$t[$i]}";
#    print O "\n";
#}
#close(O);

#sub chrLocationCW($%) {
#    my ($location, %chrLen) = @_;
#    my ($chr, $chrt, $loc);

#    $location =~ s/^([a-z\d]+)_rc//i || die $location;
#    $chrt = $1;
#    $chr = uc $1;
#    exists($chrLen{$chr}) || die;
#    $location =~ s/^_(\d+)// || die;
#    $chrLen{$chr} < $1 ? die : ($loc = $1 - $chrLen{$chr} - 1);
#    return $chrt.'_'."$loc$location";
#}

sub cutawk($) {
    my ($output) = @_;
    my $tmpFile = $output.'tmp';

    system "cut -f 2 $output".'.sam | sort | uniq -c > '.$tmpFile.'.txt';
    open(TMP, $tmpFile.'.txt') || die;
    while(<TMP>) {
	/^\s+(\d+)\s+(\d+)\s*$/ || die "$_$output";
	#print;
	system 'awk \'$2=='.$2.'\' '."$output.sam > $output"."_$2.sam";
#	    ($2 eq '0' || $2 eq '4' || $2 eq '16') && ($count += $1);
    }
    unlink "$tmpFile".'.txt';
}
