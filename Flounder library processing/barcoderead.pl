#! usr/bin/perl
use strict;
use warnings;

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
my $GENOME = $cf{'genome'}.'WC';
my ($output,$seq,%chromosome,$chr,@lines,%read,%count,$key,$i,@array,$start,$end,$temp)=('barcoderead');
print "groupLocation$FQ$GENOME.sam\n";
open SAM, "groupLocation$FQ$GENOME.sam"||die;
while(<SAM>){
	$_ =~ /BARCODE([ATCGN]+)/ || die;
	if(defined($count{$1})){
		$count{$1}++;
	}else{
		$count{$1}=0;
	}
	$read{$1}->[$count{$1}] = $_;
}
close SAM;
print "groupLocation$FQ.extendedFrags$GENOME.sam\n";
open SAM, "groupLocation$FQ.extendedFrags$GENOME.sam"||die;
while(<SAM>){
	$_ =~ /BARCODE([ATCGN]+)/ || die;
	if(defined($count{$1})){
		$count{$1}++;
	}else{
		$count{$1}=0;
	}
	$read{$1}->[$count{$1}] = $_;
}
close SAM;
print "groupLocation$FQ"."nosite$GENOME.sam\n";
open SAM, "groupLocation$FQ"."nosite$GENOME.sam"||die;
while(<SAM>){
	$_ =~ /BARCODE([ATCGN]+)/ || die;
	if(defined($count{$1})){
		$count{$1}++;
	}else{
		$count{$1}=0;
	}
	$read{$1}->[$count{$1}] = $_;
}
close SAM;
print "groupLocation$FQ"."nosite.extendedFrags$GENOME.sam\n";
open SAM, "groupLocation$FQ"."nosite.extendedFrags$GENOME.sam"||die;
while(<SAM>){
	$_ =~ /BARCODE([ATCGN]+)/ || die;
	if(defined($count{$1})){
		$count{$1}++;
	}else{
		$count{$1}=0;
	}
	$read{$1}->[$count{$1}] = $_;
}
close SAM;
print "$GENOME.fa\n";
open FA,"$GENOME.fa"||die;
while(<FA>){
	chomp;
	if($_ =~ /^>(\S+)/){
		if(defined($seq)&&defined($chr)){
			$chromosome{$chr}=$seq;
		}
		$chr = $1;
		$seq = undef;
	}elsif($_ =~ /^[ACTGNRYMKSWHDBVactgnrymkswhdbv]+$/){
		$seq .= $_;
	}else{die;}
}
if(defined($seq)&&defined($chr)){
	$chromosome{$chr}=$seq;
}
#foreach $key(keys%chromosome){
#	print $key,"	",length($chromosome{$key}),"\n";
#}
close FA;
open LOCI,"barcodeLoci$FQ"."InPoolunique.txt"||die;
open READ,">$output$FQ.txt";
while(<LOCI>){
	@lines = split /\t/,$_;
	$lines[1] =~ /^(\S+)_(\d+)NUM\d+$/ ||die;
	$chr = $1;
	$start = $2;
	$end = $2;
	#$chr =~ /scaffold/ || next;
	print READ "$lines[0]	$chr	$start	";
	defined($read{$lines[0]}) || die;
	defined($chromosome{$chr}) || die;
	$i = 0;
	while(defined($read{$lines[0]}->[$i])){
		@array = split /\t/,$read{$lines[0]}->[$i];
		if($array[2] eq $chr && $array[3]>=$start && $array[3]-$start < 1000){
			if($array[5] =~ /(\d+)M/){
				$temp = $array[3] + $1 -1;
				$end = $temp if($end<$temp);
			}
		}
		$i++;
	}
	#$end = $end - $start +1;
	#print length($chromosome{$chr}),"	$chr	$start	$end\n";
	print READ  "$end","\n";
}
close LOCI;
close READ;
