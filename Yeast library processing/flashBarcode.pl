#Use FLASH to merge paired reads.
#INPUT:($fileNameKey): 
#predefined %fi{$fileNameKey} = [$file1, $file2]: $file1 has mates with barcodes
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
my $fq = $cf{$lr."_reads"};
my $readlength = $cf{$lr."_readlength"};
my $overlapmin = int($readlength / 10);

my (@lines, %fi, $file, $t) = ();
my ($fq1, $tmpFile, $rc) = ("$fq.notCombined", "flashBarcode$fq"."tmp", 'rc');

$fi{'right'} = ["_2.fastq", "_1.fastq"];
$fi{'left'} = ["_1.fastq", "_2.fastq"];
$fi{'jumpright'} = ["_2.fastq", "_1.fastq"];
$fi{'jumpleft'} = ["_1.fastq", "_2.fastq"];
system "flash -M ".$readlength." -m ".$overlapmin." -x 0.10 $fq".$fi{$lr}->[0]." $fq".$fi{$lr}->[1]." -o $fq";


&hiseq::revcomFASTQ($fq1."_1.fastq", "$fq$rc"."_1.fastq");
&hiseq::revcomFASTQ($fq1."_2.fastq", "$fq$rc"."_2.fastq");
unlink glob "$fq1*";
system "flash -M ".$readlength." -m ".$overlapmin." -x 0.10 $fq$rc"."_1.fastq $fq$rc"."_2.fastq -o $fq$rc";
unlink glob "$fq$rc"."_*";

while(defined($file = <$fq$rc*fastq>)) {
    $file =~ /^$fq$rc(.+)$/ || die;
    $t = $1;
    &hiseq::revcomFASTQ($file, $fq."2$t");
    unlink $file;
}
