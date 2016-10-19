#! usr/bin/perl
use strict;
use warnings;

my %tagpair_count;
my $key;
my $value;
my @array;
my $i;
my $HL;
my $HR;
my $count;
my $j;
my @sample;
my $Hkey;
my $Hvalue;
my $sum;
my $sample;
my $k;
my %cf;
my $input = $ARGV[0];
open CF,"<$input"||die;
while(<CF>){
        if($_ =~ /(\S+)\t(\S+)/){
                $cf{$1} = $2;
        }
}
close CF;
my $BC='barcodePrepool'.$cf{"prepool_reads"}.'.txt';
open BARCODE,"<$BC"||die;
while(<BARCODE>){
	/(\w+)(\:+)(\w+)\s+(\d+)/ || die "ERROR: wrong input format!\n";
	$key = $1.$2.$3;
	$value = $4;
	$tagpair_count{$key} = $value;
	#print "$key	$value\n$1	$3\n\n";
}
while(($key,$value)=each %tagpair_count){
	for($i=1;$i<=$value;$i++){
	 push @array,$key;
	}
}
#print $#array;
for($sample=2;$sample<=6;$sample+=0.1){
	$k = int(10**$sample);
	for($j=0;$j<100;$j++){
		$count = 0;
		for($i=0;$i<$k;$i++){
			if($array[int(rand($#array+1))] =~ /(\w+)\:+(\w+)/){
				if(defined($HL->{$1}->{$2})){$HL->{$1}->{$2}++;}else{$HL->{$1}->{$2} = 1;}
				if(defined($HR->{$1}->{$2})){$HR->{$1}->{$2}++;}else{$HR->{$1}->{$2} = 1;}
			}
		}
		while(($Hkey,$Hvalue)=each %$HL){
			if(scalar(keys %$Hvalue) > 1){
				#while(($key,$value)=each %$Hvalue){
					#delete $HL->{$Hkey}->{$key};
				#}
				next;
			}else{
				while(($key,$value)=each %$Hvalue){
					if(scalar (keys %{$HR->{$key}}) > 1){
						#delete $HL->{$Hkey}->{$key};
					}else{
						if($value > 1){
							#delete $HL->{$Hkey}->{$key};
						}else{
							$count++;
						}
					}
				}
			}
		}
		$sum += $count/($k);
		$HL = ();
		$HR = ();
		#printf "%6.4f\n",$sum;
		#$sum = 0;
	}
	printf "%d	%6.4f\n",$k,$sum/100;
	$sum = 0;
}
close BARCODE;
