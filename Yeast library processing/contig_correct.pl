#! usr/bin/perl
use strict;
use warnings;

my $tempfile = $ARGV[0];
my $samfile = $ARGV[1];
my $qualfile = $ARGV[2];
my ($template,$seq,@temparray,@samline,%cov,%qualtemp,$i,$point,$pointread,@pos,%insertion,%mismatch,%deletion,$match,@snp,$deletion_count,$count,%qualread,$read,$key,$value,%qualsnp,$s,$ii)=();

open TP,"<$tempfile"||die;
while(<TP>){
	chomp;
	if($_ =~ /^[ATCGN]+$/){
		$seq .= $_;
	}
}
$template = $seq;
@temparray = split //,$seq;
#print scalar(@temparray),"\n";
close TP;
for($i=0;$i<=2000;$i++){
	$cov{$i} = 0;
	$qualtemp{$i} = 0;
}
open QUAL,"<$qualfile"||die;
while(<QUAL>){
	#print;
	if(/^\>(\S+)$/){
		$read = $1;
		$i = 0;
	}else{
		defined($read)||die;
		while($_ =~ s/(\d+)\s//){
			#print "$read	$i	$1\n";
			$qualread{$read}->[$i] = $1;
			$i++;
		}
	}
}
#while(($key,$value)=each%qualread){
#	print "$key	";
#	$i=0;
#	while(defined($qualread{$key}->[$i])){
#		print $qualread{$key}->[$i]," ";
#		$i++;
#	}
#	print "\n";
#}
open SAM,"<$samfile"||die;
while(<SAM>){
	#print;
	$_ =~ /^@/ && next;
	chomp;
	@samline = split /	/,$_;
	if($samline[1] == 0){
		$point = $samline[3] - 1;
		$pointread = 0;
		$deletion_count = 0;
		#print "$point\n";
		$samline[-2] =~ s/^MD\:Z\://;
		#print "$samline[5]	$samline[9]	$samline[-2]\n";# if($samline[5]=~/[SID]/);
		while($samline[5] =~ s/^(\d+)([SMID])//){
			#print "$1	$2\n";
			if($2 eq 'M'){
				for($i=0;$i<$1;$i++){
					$cov{$point}++;
					#print "$qualtemp{$point}	$qualread{$samline[0]}->[$pointread]\n";
					$qualtemp{$point} = $qualread{$samline[0]}->[$pointread] if($qualtemp{$point}<$qualread{$samline[0]}->[$pointread]);
					$point++;
					$pointread++;
				}
			}elsif($2 eq 'S'){
				#print "aaaaaaaaaaa$1\n";
				$s = $1;
				#print "$samline[9]	",$point-$samline[3]+1-$deletion_count,"	$s\n";
				substr($samline[9],$point-$samline[3]+1-$deletion_count,$s) =~ s/[ATCG]+//;
				for($i=0;$i<$s;$i++){
					$pointread++;
				}
			}elsif($2 eq 'I'){
				$ii = $1;
				#print "$point	$samline[9]	";
				$seq = substr($samline[9],$point-$samline[3]+1,$1);
				substr($samline[9],$point-$samline[3]+1,$1) =~ s/[ATCG]+//;
				$insertion{$point}->{$seq}++;
				#print "$seq\n";
				for($i=0;$i<$ii;$i++){
					$pointread++;
				}
			}elsif($2 eq 'D'){
				#print "$point	$1$2\n";
				for($i=0;$i<$1;$i++){
					$deletion_count += $1;
					$deletion{$point}++;
					$point++;
				}
			}else{die;}
		}
		$point = $samline[3] - 1;
		$pointread = 0;
		$deletion_count = 0;
		while($samline[-2] =~ s/([0-9|\^]+)([ATCG]*)//){
			$match = $1;
			$seq = $2;
			@snp = split //,$seq;
#			print "$match	@snp\n";
			if($match =~ s/\^$//){
				#print scalar(@snp),"\n";
				$point += $match+scalar(@snp);
				$pointread += $match;
				$deletion_count += scalar(@snp);
#				print "$match	@snp\n";
			}else{
				$point += $match;
				$pointread += $match;
				$i = 0;
				while(defined($snp[$i])){
#					print "$point:",$point-$samline[3]+1-$deletion_count,"	$snp[$i]->",substr($samline[9],$point-$samline[3]+1-$deletion_count,1),"\n";
					#print substr($samline[9],$point-$samline[3]+1-$deletion_count,1),"	$qualread{$samline[0]}->[$pointread]\n";
					$qualsnp{$point}->{substr($samline[9],$point-$samline[3]+1-$deletion_count,1)} = $qualread{$samline[0]}->[$pointread];
					$mismatch{$point}->{substr($samline[9],$point-$samline[3]+1-$deletion_count,1)}++;
					substr($template,$point,1) eq $snp[$i] || die;
					substr($samline[9],$point-$samline[3]+1-$deletion_count,1) eq $snp[$i] && die;
					$point++;
					$pointread++;
					$i++;
				}
			}
		}
	}else{next;}
}
#@pos = sort {$a<=>$b} keys %cov;
#for($i=0;$i<=$#pos;$i++){
#	print "$i	$cov{$i}\n" if($cov{$i}>0);
#}
my ($sup,$max,$con)=();
@pos = sort {$a<=>$b} keys %mismatch;
for($i=0;$i<=$#pos;$i++){
	#print "$pos[$i]	";
	$sup = $cov{$pos[$i]};
	$max = 0;
	while(($seq,$count) = each %{$mismatch{$pos[$i]}}){
		#print "$count$seq	";
		$sup -= $count;
		if($count>$max){
			$max = $count;
			$con = $seq;
		}
	}
	#print "$sup\n";
	if($max>$sup){
		#print "$pos[$i]:	$temparray[$pos[$i]]->$con\n";
		$temparray[$pos[$i]] = $con;
	}elsif($max == $sup && $max==1){
		#print "11 $pos[$i]:$temparray[$pos[$i]]->$con\n";
		#print $qualtemp{$pos[$i]},"	",$qualsnp{$pos[$i]}->{$con},"\n";
		if($qualsnp{$pos[$i]}->{$con} >= $qualtemp{$pos[$i]}){
			$temparray[$pos[$i]] = $con;
		}
	}
}
@pos = sort {$a<=>$b} keys %deletion;
for($i=0;$i<=$#pos;$i++){
	#print "$pos[$i]	$deletion{$pos[$i]}	$cov{$pos[$i]}\n";
	if($deletion{$pos[$i]} > $cov{$pos[$i]}){
		$temparray[$pos[$i]] = undef;
	}
}
@pos = sort {$a<=>$b} keys %insertion;
for($i=0;$i<=$#pos;$i++){
	#print "$pos[$i]	$cov{$pos[$i]}	";
	$sup = $cov{$pos[$i]};
	$max = 0;
	while(($seq,$count) = each %{$insertion{$pos[$i]}}){
		#print "$count$seq	";
		$sup -= $count;
		if($count>$max){
			$max = $count;
			$con = $seq;
		}
	}
	#print "$sup\n";
	if($max>$sup){
		$temparray[$pos[$i]] = $con.$temparray[$pos[$i]];
	}
}
#print scalar(@temparray),"\n";
#for($i=0;$i<=$#temparray;$i++){
	#print "$cov{$i}\n"; 
#}
for($i=$#temparray;$i>=0;$i--){
	if($cov{$i}<=1){
		delete $temparray[$i];
	}else{last;}
}
for($i=0;$i<=$#temparray;$i++){
	if($cov{$i}<=1){
		delete $temparray[$i];
	}else{last;}
}
#print scalar(@temparray),"\n";
for($i=0;$i<=$#temparray;$i++){
	print $temparray[$i] if(defined($temparray[$i]));
}
print "\n";
