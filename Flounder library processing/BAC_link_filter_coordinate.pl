#! usr/bin/perl
use strict;
use warnings;

my $input1 = $ARGV[0];
my $input2 = $ARGV[1];
my $input3 = $ARGV[2];
#my $input4 = $ARGV[3];
my $sh = 30000;
my (%cr,$chromosome,@lines,$count,$posleft,$posright,$key,%chr,$sum,$halfsum,$distance,%gaphash,$seq,%left,%right,%hash,$value,@arrayL,@arrayR,@key,$i,@nodes,$left,$right,%graph,$numofdisedges,$numofedges,$numofnodes,$seqL,$seqR,$leftend,$rightend,%lg,%lgcoor) = ();
#open ANCHOR,"<$input4"||die;
#while(<ANCHOR>){
#	#print;
#	@lines = split;
#	$lines[0] =~ tr/a-z/A-Z/;
#	#print "$lines[0]	$lines[1]\n";
#	$lg{$lines[0]} = $lines[1];
#	$lgcoor{$lines[0]} = $lines[2];
#}
#close ANCHOR;
open FA,"<$input3"||die;
while(<FA>){
	chomp;
	if(/^>(\S+)/){
		#print $1,"\n";
		if(defined($chromosome)){
			$cr{$chromosome} = $seq if(defined($seq));
			$chromosome = $1;
			$chromosome =~ tr/a-z/A-Z/;
			$seq = undef;
		}else{
			$chromosome = $1;
			$chromosome =~ tr/a-z/A-Z/;
			$seq = undef;
		}
	}else{
		$seq .= $_;
	}
}
$cr{$chromosome} = $seq if(defined($seq));
$seq = undef;
close FA;

open INPUT2,"<$input2"||die;
$sum = 0;
while(<INPUT2>){
	#print;
	@lines = split /\s/,$_;
	$lines[0] =~ tr/[a-z]/[A-Z]/;
	#print "$lines[0]	$lines[1]\n";
	$chr{$lines[0]} = $lines[1];
	$sum += $lines[1];
}
$halfsum = $sum / 2;
close INPUT2;
open INPUT1,"<$input1"||die;
while(<INPUT1>){
	#print;
	@lines = split /\s/,$_;
	$count = 0;
	#print "$lines[0]	$lines[1]	";
	$lines[2] =~ s/_([-]?\d+)(NUM\d+)$// || die;	
	$lines[2] =~ tr/[a-z]/[A-Z]/;
	$posleft = $1;
	if($posleft =~ /\-/){
		$count++;
		#print 0-$posleft+1,"	";
		next if(0-$posleft+1>300000);
	}else{
		#print $chr{$lines[2]} - $posleft+1,"	";
		next if($chr{$lines[2]} - $posleft+1 >300000);
	}
	$lines[3] =~ s/_([-]?\d+)(NUM\d+)$// || die;
	$lines[3] =~ tr/[a-z]/[A-Z]/;
	$posright = $1;
	if($posright =~ /\-/){
		$count++;
		#print 0-$posright,"\n";
		next if(0-$posright>300000);
	}else{
		#print $chr{$lines[3]} - $posright+1,"\n";
		next if($chr{$lines[3]} - $posright+1 >300000);
	}
	next if ($lines[2] eq $lines[3]);
	next if ($chr{$lines[2]} <= 10000);
	next if ($chr{$lines[3]} <= 10000);
	if($lines[2] gt $lines[3]){
		($lines[2],$lines[3]) = ($lines[3],$lines[2]);
		($posleft,$posright) = ($posright,$posleft);
		#print "$lines[2]	len$chr{$lines[2]}	$posleft	$lines[3]	len$chr{$lines[3]}	$posright	",$count,"\n\n";
	}elsif($lines[2] lt $lines[3]){
	}else{die;}
	$key = $lines[2].'-'.$lines[3];
	if($count%2==1){
		if($posleft<0){
			die if($posright<0);
			$distance = 0 - $posleft + $chr{$lines[3]} - $posright +2;
			if($distance>175000){
				#print"warn!\n";
				next;
			}
			$seq = 'R'.$distance.'F';
			#print "$key	$seq\n";
			push @{$gaphash{$key}},$seq;
		}else{
			die if($posleft<0);
			$distance = $chr{$lines[2]} - $posleft -$posright + 2;
			if($distance>175000){
				#print"warn!\n";
				next;
			}
			$seq = 'F'.$distance.'F';
			#print "$key	$seq\n";
			push @{$gaphash{$key}},$seq;
		}
	}elsif($count%2==0){
		if($count==2){
			$distance = 0-$posright-$posleft+2;
			if($distance>175000){
				#print"warn!\n";
				next;
			}
			$seq = 'F'.$distance.'R';
			#print "$key	$seq\n";
			push @{$gaphash{$key}},$seq;
		}elsif($count==0){
			$distance = $chr{$lines[3]}+$chr{$lines[2]}-$posleft-$posright+2;
			if($distance>175000){
				#print"warn!\n";
				next;
			}
			$seq = 'R'.$distance.'R';
			#print "$key	$seq\n";
			push @{$gaphash{$key}},$seq;
		}else{die;}
	}
	#print "$key	$posleft	$posright\n";
	$hash{$key}++;
	push @{$left{$key}},$posleft;
	push @{$right{$key}},$posright;
}
close INPUT1;
#print scalar(keys %hash),"\n";
while(($key,$value)=each %hash){
	if($value<=1){
		delete $hash{$key};
		next;
	}
	#print "$key	$value\n";
	($left,$right) = split '-',$key;
	@arrayL = @{$left{$key}};
	@arrayR = @{$right{$key}};
	@arrayL = sort{$a <=> $b} @arrayL;
	@arrayR = sort{$a <=> $b} @arrayR;
	#print $arrayL[$#arrayL]-$arrayL[0],"	",$arrayR[$#arrayR]-$arrayR[0],"\n";
	if($arrayL[$#arrayL]-$arrayL[0]>$sh && $arrayR[$#arrayR]-$arrayR[0]>$sh && $arrayL[$#arrayL]-$arrayL[0]<150000 && $arrayR[$#arrayR]-$arrayR[0]<150000){
#		if(defined($lg{$left}) && defined($lg{$right})){
#			print "$lg{$left}	$lg{$right}\n" if($lg{$left} != $lg{$right});
#		}
	}else{
		delete $hash{$key};
		next;
	}
}
#print scalar(keys %hash),"\n";
#$sum = 0;
#while(($key,$value) = each %hash){
#	$sum+=$value;
#}
#print $sum,"\n";
open INPUT1,"<$input1"||die;
while(<INPUT1>){
	#print;
	@lines = split /\s/,$_;
	$count = 0;
	#print "$lines[0]	$lines[1]	";
	$lines[2] =~ s/_([-]?\d+)(NUM\d+)$// || die;	
	$lines[2] =~ tr/[a-z]/[A-Z]/;
	$posleft = $1;
	if($posleft =~ /\-/){
		$count++;
		#print 0-$posleft+1,"	";
		next if(0-$posleft+1>300000);
	}else{
		#print $chr{$lines[2]} - $posleft+1,"	";
		next if($chr{$lines[2]} - $posleft+1 >300000);
	}
	$lines[3] =~ s/_([-]?\d+)(NUM\d+)$// || die;
	$lines[3] =~ tr/[a-z]/[A-Z]/;
	$posright = $1;
	if($posright =~ /\-/){
		$count++;
		#print 0-$posright,"\n";
		next if(0-$posright>300000);
	}else{
		#print $chr{$lines[3]} - $posright+1,"\n";
		next if($chr{$lines[3]} - $posright+1 >300000);
	}
	next if ($lines[2] eq $lines[3]);
	next if ($chr{$lines[2]} <= 10000);
	next if ($chr{$lines[3]} <= 10000);
	if($lines[2] gt $lines[3]){
		($lines[2],$lines[3]) = ($lines[3],$lines[2]);
		($posleft,$posright) = ($posright,$posleft);
		#print "$lines[2]	len$chr{$lines[2]}	$posleft	$lines[3]	len$chr{$lines[3]}	$posright	",$count,"\n\n";
	}elsif($lines[2] lt $lines[3]){
	}else{die;}
	$key = $lines[2].'-'.$lines[3];
	if($count%2==1){
		if($posleft<0){
			die if($posright<0);
			$distance = 0 - $posleft + $chr{$lines[3]} - $posright +2;
			if($distance>175000){
				#print"warn!\n";
				next;
			}
			$seq = 'R'.$distance.'F';
			#print "$key	$seq\n";
			push @{$gaphash{$key}},$seq;
		}else{
			die if($posleft<0);
			$distance = $chr{$lines[2]} - $posleft -$posright + 2;
			if($distance>175000){
				#print"warn!\n";
				next;
			}
			$seq = 'F'.$distance.'F';
			#print "$key	$seq\n";
			push @{$gaphash{$key}},$seq;
		}
	}elsif($count%2==0){
		if($count==2){
			$distance = 0-$posright-$posleft+2;
			if($distance>175000){
				#print"warn!\n";
				next;
			}
			$seq = 'F'.$distance.'R';
			#print "$key	$seq\n";
			push @{$gaphash{$key}},$seq;
		}elsif($count==0){
			$distance = $chr{$lines[3]}+$chr{$lines[2]}-$posleft-$posright+2;
			if($distance>175000){
				#print"warn!\n";
				next;
			}
			$seq = 'R'.$distance.'R';
			#print "$key	$seq\n";
			push @{$gaphash{$key}},$seq;
		}else{die;}
	}
	if(defined($hash{$key})){
		#print $_;
		$left = $lines[2];
		$right = $lines[3];
		$left =~ tr/A-Z/a-z/;
		$right =~ tr/A-Z/a-z/;
		if($posleft < 0){
			$leftend = abs($posleft) - 100;
			$leftend = 0 if($leftend<0);
		}else{
			#print "$posleft	",$chr{$lines[2]},"\n";
			if($posleft+100>$chr{$lines[2]}){
				$leftend = $chr{$lines[2]};
			}else{
				$leftend = $posleft+100;
			}
		}
		if($posright < 0){
			$rightend = abs($posright) - 100;
			$rightend = 0 if($rightend<0);
		}else{
			if($posright+100>$chr{$lines[3]}){
				$rightend = $chr{$lines[3]};
			}else{
				$rightend = $posright+100;
			}
		}
		#print "$_\n";
		print "$left	",abs($posleft),"	$leftend	$right	",abs($posright),"	$rightend\n";
	}
}
close INPUT1;
