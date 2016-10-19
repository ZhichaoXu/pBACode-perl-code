#! usr/bin/perl
use strict;
use warnings;
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
my $POOL = $cf{"pool_file"};
my ($OUTPUT,@lines,%reads,%qual,%count,%pair,%cov,$key,$value,$i,$contig,$len,$seq,%Contig,$key_contig,$value_contig,@array,%raw_assembly,%corrected_assembly,$last) = ('groupAssemble');
my %SINGLE = ($FQ => ['.extendedFrags']);
exists($SINGLE{$FQ}) || die;
my $PAIRED = 'noHindIII';
my $SINGLEINPUT = 'rmVectorBarcode';
my $PAIREDINPUT = 'rmVectorBarcode1mate';
sub revdnacmp{
	my $dna = shift;
	my $revcomp = reverse($dna);
	$revcomp =~ tr/ACTGNactgn/TGACNtgacn/;
	return $revcomp;
}
print "$SINGLEINPUT$FQ$SINGLE{$FQ}->[0].fastq\n";
open FQ,"<$SINGLEINPUT$FQ$SINGLE{$FQ}->[0].fastq"||die;
while(<FQ>){
	push @lines,$_;
	if(@lines == 4){
		if($lines[0] =~ /^(\S+)BARCODE([ACGTN]+)/){
			if(defined($count{$2})){
				$count{$2}++;
			}else{
				$count{$2}=0;
			}
			$seq = $lines[1];
			#$seq = $lines[0].$lines[1].$lines[2].$lines[3];
			$reads{$2}->[$count{$2}] = $seq;
			$qual{$2}->[$count{$2}] = $lines[3];
			#print "$seq\n";
		}else{die;}
		@lines = ();
	} elsif(@lines > 4) {
	    die scalar @lines;
	}
}
close FQ;
print "$PAIREDINPUT$FQ"."2.fastq\n";
open FQ,"<$PAIREDINPUT$FQ"."2.fastq"||die;
while(<FQ>){
	chomp;
	push @lines,$_;
	if(@lines == 4){
		if($lines[0] =~ /^(\S+)BARCODE([ACGTN]+)/){
			if(defined($count{$2})){
				$count{$2}++;
			}else{
				$count{$2}=0;
			}
			$lines[1] = &revdnacmp($lines[1]);
			$lines[3] = reverse $lines[3];
			$seq = $lines[1]."\n";
			#$seq = $lines[0]."\n".$lines[1]."\n".$lines[2]."\n".$lines[3]."\n";
			$reads{$2}->[$count{$2}] = $seq;
			$qual{$2}->[$count{$2}] = $lines[3];
			#print "$seq\n";
		}else{die;}
		@lines = ();
	} elsif(@lines > 4) {
	    die scalar @lines;
	}
}
close FQ;
my %leftbarcode;
my %rightbarcode;
open(TAG, "$POOL.txt") || die;
while(<TAG>) {
	/^([ACGTN]+)\:([ACGTN]+)\s+(\d+)\s*$/ || die "$_*";
	$leftbarcode{$1}+= $3;
	$rightbarcode{&revdnacmp($2)} += $3;
}
my %hash;
if($lr eq 'left'){
	%hash = %leftbarcode;
}elsif($lr eq 'right'){
	%hash = %rightbarcode;
}else{die;}
my $count=0;
my $found=0;
my $notfound=0;
print scalar(keys %reads),"\n";
#$last = 0;
open LOG,">groupAssembletemp_phrap$FQ.log"||die;
while(($key,$value)=each %reads){
	if(defined($hash{$key})){
		$count++;
		$i=0;
		$found++;
		open TMP,">groupAssembletemp_phrap$FQ.fa"||die;
		open QUAL,">groupAssembletemp_phrap$FQ.fa.qual"||die;
		while(defined($reads{$key}->[$i])){
			print TMP ">$key.$i\n$reads{$key}->[$i]";
			print QUAL ">$key.$i\n";
			while($qual{$key}->[$i] =~ s/^(\S)//){
				print QUAL ord($1)-33," ";
			}
			print QUAL "\n";
			$i++;
		}
		close TMP;
		close QUAL;
		system "phrap groupAssembletemp_phrap$FQ.fa -new_ace 1>temp.log 2>temp.err";
		%cov = ();
		%Contig = ();
		open(ACE, "<groupAssembletemp_phrap$FQ.fa.ace")||die;
		while(<ACE>){
			if(/^CO (Contig\d+) \d+ (\d+) /){
				$cov{$1} = $2;				
			}
		}
		close ACE;
		open(CONTIG, "<groupAssembletemp_phrap$FQ.fa.contigs")||die;
		$contig = undef;
		$seq = undef;
		while(<CONTIG>){
			#print;
			chomp;
			if(/^\>(\S+)/){
				if(defined($contig) && defined($seq)){
					$Contig{$contig} = $seq;
				}
				$contig = $1;
				$seq = undef;
			}else{
				$seq .= $_;
			}
		}
		close CONTIG;
		$Contig{$contig} = $seq if(defined($contig) && defined($seq));
		#while(($key_contig,$value_contig) = each %Contig){print ">$key_contig\n$value_contig\n";}
		@array = sort{$b<=>$a} values %cov;
		print LOG "@array\n";
		if(scalar(@array)>1 && $array[1]*5>$array[0]){
			
		}else{
			open TMP,">groupAssembletemp_contig_phrap$FQ.fa"||die;
			foreach $key_contig (sort{$cov{$b}<=>$cov{$a}} keys %cov){
				$contig = "groupAssembletemp_phrap$FQ.fa.".$key_contig;
				print TMP ">$contig\n$Contig{$contig}\n";
				$raw_assembly{$key} = $Contig{$contig};
				last;
			}
			close TMP;
			system "bowtie2-build -q groupAssembletemp_contig_phrap$FQ.fa groupAssembletemp_contig_phrap$FQ";
			system "bowtie2-align -x groupAssembletemp_contig_phrap$FQ -f groupAssembletemp_phrap$FQ.fa -k 10 --local --quiet -S groupAssembletemp_phrap$FQ.sam";
			system "perl ./contig_correct.pl groupAssembletemp_contig_phrap$FQ.fa groupAssembletemp_phrap$FQ.sam groupAssembletemp_phrap$FQ.fa.qual > temp";
			open TMP,"<temp"||die;
			while(<TMP>){
				chomp;
				$contig = $_ if(/^[ATCGN]+$/);
				$corrected_assembly{$key} = $contig;
			}
			close TMP;
			print LOG "raw assembly: ",scalar(keys %raw_assembly),"\ncorrected assembly: ",scalar(keys %corrected_assembly),"\n";
		}
		#$last <0? $last++ : last;
	}else{$notfound++;}
}
close LOG;
open RAWUNIQUE,">groupAssembly$FQ"."RawUniqueInpool.fa"||die;
while(($key,$value) = each %raw_assembly){
	print RAWUNIQUE ">$key\n$value\n";
}
close RAWUNIQUE;
open CORRECTUNIQUE,">groupAssembly$FQ"."CorrectedUniqueInpool.fa"||die;
while(($key,$value) = each %corrected_assembly){
	print CORRECTUNIQUE ">$key\n$value\n";
}
close CORRECTUNIQUE;
#print "Found: $found\nNot found: $notfound\n";
