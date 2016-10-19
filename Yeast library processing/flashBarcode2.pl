#Trim low quality 3' end and then merge mated reads
#INPUT:($fileNameKey): 
#predefined %fi{$fileNameKey} = [$file1, $file2]: $file1 has mates with barcodes
#predefined %OVERLAP{$fileNameKey} = [250, 25]: maximum and minimum overlapping region for FLASH
#trim 3' end of each read mate so that 3' 5bp all have quality score > 20.
#FLASH 1st: merge paired reads whose 5' ends are not complement with the other read of the pair
#           OUTPUT: INPUT.extendedFrags.fastq
#FLASH 2st: To merge paires whose 5' ends are complemented by their mate reads, uncombined reads from FLASH 1 will be reverse complemented and merged. OUTPUT reads will be reversed back
#           OUTPUT: INPUT2.extendedFrags.fastq INPUT2.notCombined_1.fastq INPUT2.notCombined_2.fastq
#OUTPUT: Barcodes are at 5' side of merged mates and in "_1.fastq" of uncombined read files.
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
my $fq = $cf{$lr."_reads"}."2.notCombined";
#my $fq = $ARGV[0]; #Library129_yeast_pvuii_digestion_first2.notCombined
my $TRIM = $cf{$lr."_reads"}."3";
#my $TRIM = $ARGV[1]; #Library129_yeast_pvuii_digestion_first3
my (@lines, %f1, %f2, %lr, $file, $i, $j, ,$k, $t,%fi) = ();
$fi{'right'} = ["_1.fastq", "_2.fastq"];
$fi{'left'} = ["_1.fastq", "_2.fastq"];
$fi{'jumpright'} = ["_1.fastq", "_2.fastq"];
$fi{'jumpleft'} = ["_1.fastq", "_2.fastq"];
#my @fi = ("_1.fastq", "_2.fastq");
my ($fq1, $tmpFile, $rc) = ("$TRIM.notCombined", "flashBarcode$fq"."tmp", 'rc');
my $readlength = $cf{$lr."_readlength"};
my $overlapmin = int($readlength / 10);
#my(%OVERLAP) = ('Library129_yeast_pvuii_digestion_first2.notCombined' => [300, 30]);
#exists($OVERLAP{$fq});

open(FQ, $fq.$fi{$lr}->[0]) || die $fq.$fi{$lr}->[0];
open(TMP, ">$tmpFile") || die;

%f1 = ();
$i = 0;
$j = 0;
$k = 0;
while(<FQ>){
    push @lines,$_;

    if(@lines == 4) {
	chomp($lines[0]);
	chomp($lines[1]);
	chomp($lines[2]);
	chomp($lines[3]);
	$lines[0] =~ /^([^\s]+)/ || die;
	$lines[1] =~ /^([ACGTN]+)$/ || die;
	$t = length($1);
	$lines[3] =~ /^[^\s]{$t}$/ || die "@lines$fq.$fi{$lr}->[0]";
	$lines[2] =~ /^\+$/ || die;
	$lines[3] = &hiseq::SQtrimEnd3($lines[3], 5, 20);
	$k++ if length($lines[3])<$t;
	$lines[1] = substr($lines[1], 0, length($lines[3]));
	$lines[1] =~ /\s/ && die;
	$lines[2] =~ /\s/ && die;
	$lines[3] =~ /\s/ && die;

	if(length($lines[1]) > 50) {
	    $lines[1] =~ /^([ACGTN]+)$/ || die "$lines[1]*";
	    $t = length($1);
	    $lines[3] =~ /^[^\s]{$t}$/ || die "@lines$fq.$fi{$lr}->[0]";
	    $lines[0] =~ /^([^\s]+)\s*[^\s]*$/ || die $lines[0];
	    length($lines[1]) != length($lines[3]) && die;
	    exists($f1{$1}) ? die : ($f1{$1} = 1);
	    print TMP join("\n", @lines), "\n";
	    $j++;
	} else {
	    $i++;
	}

	@lines = ();
    } elsif(@lines > 4) {
	die scalar(@lines);
    }
}
print "$k mate1 trimmed shorter\n";
print "$i mate1 trimmed off while ", scalar(keys %f1), " left\t$j\n";
close(TMP);
close(FQ);

$j = 0;
open(FQ, $tmpFile) || die;
while(<FQ>){
    push @lines,$_;

    if(@lines == 4) {
	chomp($lines[0]);
	chomp($lines[1]);
	chomp($lines[2]);
	chomp($lines[3]);
	$lines[1] =~ /\s/ && die;
	$lines[2] =~ /\s/ && die;
	$lines[3] =~ /\s/ && die;
	$lines[0] =~ /^([^\s]+)/ || die;
	$lines[1] =~ /^([ACGTN]+)$/ || die;
	$t = length($1);
	$lines[3] =~ /^[^\s]{$t}$/ || die;
	$lines[2] =~ /^\+$/ || die;

	$j++;
	@lines = ();
    } elsif(@lines > 4) {
	die scalar(@lines);
    }
}
print "$j\n";
close(FQ);

open(FQ, $fq.$fi{$lr}->[1]) || die $fq.$fi{$lr}->[1];
open(F, ">$TRIM".$fi{$lr}->[1]) || die;

%f2 = ();
$i = 0;
$j = 0;
$k = 0;
while(<FQ>){
    push @lines,$_;

    if(@lines == 4) {
	$lines[0] =~ /^([^\s]+)/ || die;

	if(exists($f1{$1})) {
	    chomp($lines[1]);
	    chomp($lines[3]);
	    $t = length($lines[3]);
	    $lines[3] = &hiseq::SQtrimEnd3($lines[3], 5, 20);
	    $k++ if length($lines[3])<$t;
	    $lines[1] = substr($lines[1], 0, length($lines[3]));
	    $lines[1] .= "\n";
	    $lines[3] .= "\n";

	    if(length($lines[1]) > 50) {
		print F @lines;
		$lines[0] =~ /^([^\s]+)/ || die;
		exists($f2{$1}) ? die : ($f2{$1} = 1);
		$j++;
		exists($f1{$1}) || die;
	    } else {
		$i++;
	    }
	}

	@lines = ();
    } elsif(@lines > 4) {
	die scalar(@lines);
    }
}
print "$k mate2 trimmed shorter\n";
print "$i mate2 trimmed off while ", scalar(keys %f2), " pairs left\t$j\n";
close(FQ);

$j = 0;
open(FQ, $tmpFile) || die $fq.$fi{$lr}->[1];
open(F, ">$TRIM".$fi{$lr}->[0]) || die;
while(<FQ>){
    push @lines,$_;

    if(@lines == 4) {
	$lines[0] =~ /^([^\s]+)/ || die;

	if(exists($f2{$1})) {
	    print F @lines;
	    delete($f2{$1});
	}

	$j++;
	@lines = ();
    } elsif(@lines > 4) {
	die scalar(@lines);
    }
}
print "$j\n";
scalar(keys %f2) > 0 && die scalar(keys %f2);
#unlink $tmpFile;
system "flash -M ".$readlength." -m ".$overlapmin." -x 0.10 $TRIM"."$fi{$lr}->[0] $TRIM"."$fi{$lr}->[1] -o $TRIM";

#system "perl SEflash.pl ".$readlength.' '.$overlapmin." $TRIM".$fi{$lr}->[0]." $TRIM".$fi{$lr}->[1]." $TRIM";
print "1\n";
#unlink glob "$TRIM*";


for($i = 1; $i <= 2; $i++) {
    open(FQ, $fq1."_$i.fastq") || die $fq1."_$i.fastq";
    open(TMP, ">$TRIM$rc"."_$i.fastq") || die;

    while(<FQ>){
	push @lines,$_;

	if(@lines == 4) {
	    chomp($lines[1]);
	    chomp($lines[3]);
	    $lines[1] = &hiseq::revcomDNA($lines[1])."\n";
	    $lines[3] = reverse($lines[3])."\n";
	    print TMP @lines;
	    @lines = ();
	} elsif(@lines > 4) {
	    die scalar(@lines);
	}
    }
}
unlink glob "$fq1*";
system "flash -M ".$readlength." -m ".$overlapmin." -x 0.10 $TRIM$rc"."$fi{$lr}->[0] $TRIM$rc"."$fi{$lr}->[1] -o $TRIM$rc";
#system "perl SEflash.pl ".$readlength.' '.$overlapmin." $TRIM$rc"."_1.fastq $TRIM$rc"."_2.fastq $TRIM$rc";
unlink glob "$TRIM$rc"."_*";

while(defined($file = <$TRIM$rc*fastq>)) {
    open(FQ, $file) || die;
    open(TMP, ">$tmpFile.fastq") || die;
    while(<FQ>){
	push @lines,$_;

	if(@lines == 4) {
	    chomp($lines[1]);
	    chomp($lines[3]);
	    $lines[1] = &hiseq::revcomDNA($lines[1])."\n";
	    $lines[3] = reverse($lines[3])."\n";
	    print TMP @lines;
	    @lines = ();
	} elsif(@lines > 4) {
	    die scalar(@lines);
	}
    }

    $t = $file;
    $file =~ s/^$TRIM$rc/2/ || die;
    $file = $TRIM.$file;
    unlink $t;
    system "mv $tmpFile.fastq $file";
}
