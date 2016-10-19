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
my $SINGLEINPUT = 'rmVectorBarcode';
my $PAIREDINPUT = 'rmVectorBarcode1mate';
my (@lines,%reads,$barcode,%count,$seq,$key,$value,%hash,$count,@seq,$sum,$key_hash,$value_hash) = ();
#print "$SINGLEINPUT$FQ.extendedFrags.fastq\n";
open FQ,"<$SINGLEINPUT$FQ.extendedFrags.fastq" || die;
while(<FQ>){
	push @lines,$_;
	if(@lines == 4){
		if($lines[0] =~ /^(\S+)BARCODE([ACGTN]+)/){
			$barcode = $2;
			if(length($lines[1])<28){
				@lines = ();
				next;
			}
			if(defined($count{$barcode})){
				$count{$barcode}++;
			}else{
				$count{$barcode}=0;
			}
			$lines[0] = $1."\n";
			$lines[1] = substr($lines[1],6,21)."\n";
			$lines[3] = substr($lines[3],6,21)."\n";
			$seq = $lines[0].$lines[1].$lines[2].$lines[3];
			$reads{$barcode}->[$count{$barcode}] = $seq;
			#print "$seq";
		}else{die;}
		@lines = ();
	} elsif(@lines > 4) {
	    die;
	}
}
close FQ;
#print "$PAIREDINPUT$FQ"."1.fastq\n";
open FQ,"<$PAIREDINPUT$FQ"."1.fastq"||die;
while(<FQ>){
	push @lines,$_;
	if(@lines == 4){
		if($lines[0] =~ /^(\S+)BARCODE([ACGTN]+)/){
			$barcode = $2;
			if(length($lines[1])<28){
				@lines = ();
				next;
			}
			if(defined($count{$barcode})){
				$count{$barcode}++;
			}else{
				$count{$barcode}=0;
			}
			$lines[0] = $1."\n";
			$lines[1] = substr($lines[1],6,21)."\n";
			$lines[3] = substr($lines[3],6,21)."\n";
			$seq = $lines[0].$lines[1].$lines[2].$lines[3];
			$reads{$barcode}->[$count{$barcode}] = $seq;
			#print "$seq";
		}else{die;}
		@lines = ();
	} elsif(@lines > 4) {
	    die;
	}
}
close FQ;
#print scalar(keys%reads);
open TEMP, ">barcodeNonunique$FQ.txt"||die;
while(($key,$value) = each %reads){
	#print "$key	";
	scalar(@{$value}) <= 5 && next;
	open TEMPINPUT, ">./SEEDtemp/temp$lr.fastq";
	print TEMPINPUT @{$value};
	close TEMPINPUT;
	system "SEED --input ./SEEDtemp/temp$lr.fastq --output ./SEEDtemp/tempout$lr --short";
	open TEMPOUTput, "./SEEDtemp/tempout$lr";
	%hash = ();
	$barcode = undef;
	$sum = 0;
	while(<TEMPOUTput>){
		#print;
		if(/^([ATCGN]+)\n$/){
			if(defined($barcode)){
				$hash{$barcode} = $count;
				$barcode = $1;
				$count = 0;
			}else{
				$barcode = $1;
				$count = 0;
			}
		}elsif(/\d+\s+\d+/){
			$count++;
			$sum++;
		}
	}
	if(defined($barcode)){
		#$count == 0 && die;
		$hash{$barcode} = $count;
		$barcode = $1;
		$count = 0;
	}
	#print $sum,"\n";
	close TEMPOUTput;
	@seq = sort{$hash{$b} <=> $hash{$a}} keys %hash;
	if($hash{$seq[0]}/$sum >0.9){
	}elsif($hash{$seq[1]}>1){
		print TEMP "$key\n";
	}
	#last;
}
close TEMP;















