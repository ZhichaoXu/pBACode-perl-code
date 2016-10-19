#! usr/bin/perl
use strict;
use warnings;

my $config = $ARGV[0];
my %cf;
open(CF, $config) || die "$config\n";
while(<CF>){
	if($_ =~ /(\S+)\t(\S+)/){
		$cf{$1} = $2;
	}
}
my $leftfq = $cf{"left_reads"};
my $rightfq = $cf{"right_reads"};
my $nohind3 = 'nosite';

system "cat rmVectorBarcode1mate$leftfq"."1.fastq > tmpfile";
system "mv rmVectorBarcode1mate$leftfq"."2.fastq rmVectorBarcode1mate$leftfq"."1.fastq";
system "mv tmpfile rmVectorBarcode1mate$leftfq"."2.fastq";

system "cat rmVectorBarcode1mate$leftfq$nohind3"."1.fastq > tmpfile";
system "mv rmVectorBarcode1mate$leftfq$nohind3"."2.fastq rmVectorBarcode1mate$leftfq$nohind3"."1.fastq";
system "mv tmpfile rmVectorBarcode1mate$leftfq$nohind3"."2.fastq";

system "cat rmVectorBarcode1mate$rightfq"."1.fastq > tmpfile";
system "mv rmVectorBarcode1mate$rightfq"."2.fastq rmVectorBarcode1mate$rightfq"."1.fastq";
system "mv tmpfile rmVectorBarcode1mate$rightfq"."2.fastq";

system "cat rmVectorBarcode1mate$rightfq$nohind3"."1.fastq > tmpfile";
system "mv rmVectorBarcode1mate$rightfq$nohind3"."2.fastq rmVectorBarcode1mate$rightfq$nohind3"."1.fastq";
system "mv tmpfile rmVectorBarcode1mate$rightfq$nohind3"."2.fastq";


