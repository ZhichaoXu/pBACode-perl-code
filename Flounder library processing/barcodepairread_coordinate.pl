#! usr/bin/perl
use strict;
use warnings;
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
my $genome = $cf{'genome'};
my (@lines,%readl,%readr,$i,$barcode,$seq,%chrlen,$chr,$start,$end)=();
open CHRLEN,"<$genome.Len.txt"||die;
while(<CHRLEN>){
	chomp;
	@lines = split;
	defined($chrlen{$lines[0]}) && die;
	$chrlen{$lines[0]} = $lines[1];
}
close CHRLEN;
print "barcoderead$FILELEFT.txt","\n";
open LF,"barcoderead$FILELEFT.txt"||die;
while(<LF>){
	chomp;
	@lines = split;
	$barcode = $lines[0];
	$chr = $lines[1];
	$start = $lines[2];
	$end = $lines[3];
	if($chr =~ s/_rc$//){
		defined($chrlen{$chr}) || die;
		$start = $chrlen{$chr} + 1 - $start;
		$end = $chrlen{$chr} + 1 - $end;
	}
	#print "$barcode	$chr	$start	$end\n";
	$seq = $chr."	".$start."	".$end;
	$readl{$barcode} = $seq;
}
close LF;
print "barcoderead$FILERIGHT.txt","\n";
open RF,"barcoderead$FILERIGHT.txt"||die;
while(<RF>){
	chomp;
	@lines = split;
	$barcode = $lines[0];
	$barcode = reverse $barcode;
	$barcode =~ tr/ATCGNatcgn/TAGCNtagcn/;
	$chr = $lines[1];
	$start = $lines[2];
	$end = $lines[3];
	if($chr =~ s/_rc$//){
		defined($chrlen{$chr}) || die;
		$start = $chrlen{$chr} + 1 - $start;
		$end = $chrlen{$chr} + 1 - $end;
	}
	#print "$barcode	$chr	$start	$end\n";
	$seq = $chr."	".$start."	".$end;
	$readr{$barcode} = $seq;
}
close LF;
print scalar(keys %readl),"	",scalar(keys %readr),"\n";
open TAG, "$POOL.txt" || die;
open OUTPUT, ">barcodepairread$FILELEFT$FILERIGHT.tab.txt" || die;
while(<TAG>) {
	/^([ACGTN]+)\:([ACGTN]+)\s+(\d+)\s*$/ || die "$_*";
	if(defined($readl{$1})&&defined($readr{$2})){
		print OUTPUT "$readl{$1}	$readr{$2}\n";
		#print READ1 ">LeftBarcode$1RightBarcode$2\n$readl{$1}\n";
		#for($i=0;$i<length($readl{$1});$i++){
		#	print READ1 "I";
		#}
		#print READ1 "\n";
		#print READ2 ">LeftBarcode$1RightBarcode$2\n$readr{$2}\n";
		#for($i=0;$i<length($readr{$2});$i++){
		#	print READ2 "I";
		#}
		#print READ2 "\n";
	}
}
close TAG;
close OUTPUT;
