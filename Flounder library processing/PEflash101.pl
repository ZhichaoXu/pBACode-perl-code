#Use FLASH to merge paired reads. Merged reads use file1's 5' end as 5' and file2's 5' end as 3'. If direction of merged reads is reversed, print the read
#! usr/bin/perl
use strict;
use warnings;
require hiseq;
my $config = $ARGV[0];
my $prepost = $ARGV[1];
my %cf;
open(CF, $config) || die "$config\n";
while(<CF>){
	if($_ =~ /(\S+)\t(\S+)/){
		$cf{$1} = $2;
	}
}
my $fq = $cf{$prepost."_reads"};
my $readlength = $cf{$prepost."_readlength"};
my $fragmentlength = $cf{$prepost."_fragmentlength"};
my $overlapmin = int((2*$readlength - $fragmentlength)/2-1);
my ($fq1, $tmpFile, $LEN) = ("$fq.notCombined", "SEflash$fq"."tmp", 20);
my @fi = ("_1.fastq", "_2.fastq");
my (@lines, %lr, $i, $num, $t) = ();

system "flash -M $readlength -m $overlapmin -x 0.10 $fq"."_1.fastq $fq"."_2.fastq -o $fq";
