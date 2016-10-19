#! usr/bin/perl
use strict;
use warnings;

my $config = $ARGV[0];
my %cf;
my (%nl,%nr,$seq,@lines) = ();
open(CF, $config) || die "$config\n";
while(<CF>){
	if($_ =~ /(\S+)\t(\S+)/){
		$cf{$1} = $2;
	}
}
my $lr = 'left';
my $FQ = $cf{$lr."_reads"};
my $POOL = $cf{"pool_file"};
open LN,"<barcodeNonunique$FQ.txt"||die;
while(<LN>){
	chomp;
	defined($nl{$_}) && die;
	$nl{$_} = 1;
}
close LN;
open LociL,"<barcodeLocus$FQ"."InPoolunique.txt"||die;
open LL,">barcodeLoci$FQ"."InPoolunique.txt" || die;
while(<LociL>){
	@lines = split;
	defined($nl{$lines[0]}) || print LL $_;
}
close LociL;
close LL;

$lr = 'right';
$FQ = $cf{$lr."_reads"};
open RN,"<barcodeNonunique$FQ.txt"||die;
while(<RN>){
	chomp;
	$seq = reverse($_);
	$seq =~ tr/ATCGNatcgn/TAGCNtagcn/;
	defined($nr{$seq}) && die;
	$nr{$seq} = 1;
}
close RN;
open LociR,"<barcodeLocus$FQ"."InPoolunique.txt"||die;
open LR,">barcodeLoci$FQ"."InPoolunique.txt" || die;
while(<LociR>){
	@lines = split;
	$lines[0] = reverse($lines[0]);
	$lines[0] =~ tr/ATCGNatcgn/TAGCNtagcn/;
	defined($nr{$lines[0]}) || print LR $_;
}
close LociR;
close LR;
#open POOL,"<$POOL.txt" || die;
#while(<POOL>){
#	/^([ACGTN]+)\:([ACGTN]+)\s+(\d+)\s*$/ || die "$_*";
#	print if(defined($nl{$1})&&defined($nr{$2}));
#}
#close POOL;
