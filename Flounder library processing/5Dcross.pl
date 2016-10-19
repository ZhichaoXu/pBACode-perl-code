#! usr/bin/perl
use strict;
use warnings;
use List::Util qw(min);
my $cutoff = $ARGV[0];
my (%countones,%countcol,%countrow,%counttens,%counthundred,$p,$i,$j,$k,$t,$key,$value,$barcode,$num,$temp,%delete,$key_delete,$value_delete,$count,@array,$temp2,%unique,$seq,%solved,$unsolved) = ();
my (%ones,%tens,%col,%row,%hundred)=();
for($i=0;$i<10;$i++){
	open INPUT,"<barcodeflounder_all_ones_".$i.".fastq.txt"||die;
	$countones{$i} = 0;
	while(<INPUT>){
		/^([ACGTN]+)\t(\d+)/ || die;
		defined($ones{$i}->{$1})? die :($ones{$i}->{$1}=$2);
		$countones{$i} += $2;
	}
	close INPUT;
}
while(($key,$value) = each %ones){
	while(($barcode,$num) = each %{$value}){
		$ones{$key}->{$barcode} = 10*96*$num/$countones{$key};
		#print "$ones{$key}->{$barcode}\n";
		delete $ones{$key}->{$barcode} if($ones{$key}->{$barcode}<=$cutoff);
	}
	#print $key,"	$countones{$key}	",scalar(keys %{$value}),"\n";
}
for($i=0;$i<10;$i++){
	open INPUT,"<barcodeflounder_all_tens_".$i.".fastq.txt"||die;
	$counttens{$i} = 0;
	while(<INPUT>){
		/^([ACGTN]+)\t(\d+)/ || die;
		defined($tens{$i}->{$1})? die :($tens{$i}->{$1}=$2);
		$counttens{$i} += $2;
	}
	close INPUT;
}
while(($key,$value) = each %tens){
	while(($barcode,$num) = each %{$value}){
		$tens{$key}->{$barcode} = 10*96*$num/$counttens{$key};		
		#print "$tens{$key}->{$barcode}\n";
		delete $tens{$key}->{$barcode} if($tens{$key}->{$barcode}<=$cutoff);
	}
	#print $key,"	$counttens{$key}	",scalar(keys %{$value}),"\n";
}
for($i=0;$i<10;$i++){
	open INPUT,"<barcodeflounder_all_hundred_".$i.".fastq.txt"||die;
	$counthundred{$i} = 0;
	while(<INPUT>){
		/^([ACGTN]+)\t(\d+)/ || die;
		defined($hundred{$i}->{$1})? die :($hundred{$i}->{$1}=$2);
		$counthundred{$i} += $2;
	}
	close INPUT;
}
while(($key,$value) = each %hundred){
	while(($barcode,$num) = each %{$value}){
		$hundred{$key}->{$barcode} = 10*96*$num/$counthundred{$key};		
		#print "$hundred{$key}->{$barcode}\n";
		delete $hundred{$key}->{$barcode} if($hundred{$key}->{$barcode}<=$cutoff);
	}
	#print $key,"	$counthundred{$key}	",scalar(keys %{$value}),"\n";
}

for($i=0;$i<12;$i++){
	$temp = $i+1;
	open INPUT,"<barcodeflounder_all_col_".$temp.".fastq.txt"||die;
	$countcol{$i} = 0;
	while(<INPUT>){
		/^([ACGTN]+)\t(\d+)/ || die;
		defined($col{$i}->{$1})? die :($col{$i}->{$1}=$2);
		$countcol{$i} += $2;
	}
	close INPUT;
}
while(($key,$value) = each %col){
	while(($barcode,$num) = each %{$value}){
		$col{$key}->{$barcode} = 100*8*$num/$countcol{$key};
		delete $col{$key}->{$barcode} if($col{$key}->{$barcode}<=$cutoff);
	}
	#print $key,"	$countcol{$key}	",scalar(keys %{$value}),"\n";
}
for($i=0;$i<8;$i++){
	#print chr($i+97),"\n";
	open INPUT,"<barcodeflounder_all_row_".chr($i+97).".fastq.txt"||die;
	$countrow{$i} = 0;
	while(<INPUT>){
		/^([ACGTN]+)\t(\d+)/ || die;
		defined($row{$i}->{$1})? die :($row{$i}->{$1}=$2);
		$countrow{$i} += $2;
	}
	close INPUT;
}
while(($key,$value) = each %row){
	while(($barcode,$num) = each %{$value}){
		$row{$key}->{$barcode} = 100*12*$num/$countrow{$key};
		delete $row{$key}->{$barcode} if($row{$key}->{$barcode}<=$cutoff);
	}
	#print $key,"	$countrow{$key}	",scalar(keys %{$value}),"\n";
}
for($p=0;$p<=9;$p++){#hundred
	for($j=0;$j<=9;$j++){#ones
		for($i=0;$i<=9;$i++){#tens
			for($k=0;$k<=7;$k++){#row
				for($t=0;$t<=11;$t++){#col
					$solved{$p}->{$i}->{$j}->{$k}->{$t} = 0;
				}
			}
		}
	}
}
for($i=0;$i<10;$i++){
	if(defined($ones{$i})){
		$count = 0;
		while(($key,$value) = each %{$ones{$i}}){
			%delete = ();
			$delete{$i} = $value;
			for($j=$i+1;$j<10;$j++){
				if(defined($ones{$j}->{$key})){
					defined($delete{$j})&&die;
					$delete{$j} = $ones{$j}->{$key};
				}
			}
			scalar(keys %delete) == 1 && next;
			scalar(keys %delete) < 1 &&die;
			@array = values %delete;
			#print "@array	";
			@array = sort {$b <=> $a} @array;
			if($array[0]>5*$array[1]){
				while(($key_delete,$value_delete) = each %delete){
					$value_delete == $array[0] && next;
					$key_delete == $i && $count++;
					delete $ones{$key_delete}->{$key};
				}
			}else{
				while(($key_delete,$value_delete) = each %delete){
					$key_delete == $i && $count++;
					delete $ones{$key_delete}->{$key};
				}
			}
		}
		#print $count,"	$i	",scalar keys($ones{$i}),"\n";
	}else{die;}
}
for($i=0;$i<10;$i++){
	if(defined($tens{$i})){
		$count = 0;
		while(($key,$value) = each %{$tens{$i}}){
			%delete = ();
			$delete{$i} = $value;
			for($j=$i+1;$j<10;$j++){
				if(defined($tens{$j}->{$key})){
					defined($delete{$j})&&die;
					$delete{$j} = $tens{$j}->{$key};
				}
			}
			scalar(keys %delete) == 1 && next;
			scalar(keys %delete) < 1 &&die;
			@array = values %delete;
			#print "@array	";
			@array = sort {$b <=> $a} @array;
			if($array[0]>5*$array[1]){
				while(($key_delete,$value_delete) = each %delete){
					$value_delete == $array[0] && next;
					$key_delete == $i && $count++;
					delete $tens{$key_delete}->{$key};
				}
			}else{
				while(($key_delete,$value_delete) = each %delete){
					$key_delete == $i && $count++;
					delete $tens{$key_delete}->{$key};
				}
			}
		}
		#print $count,"	$i	",scalar keys($tens{$i}),"\n";
	}else{die;}
}

for($i=1;$i<9;$i++){
	if(defined($hundred{$i})){
		$count = 0;
		while(($key,$value) = each %{$hundred{$i}}){
			%delete = ();
			$delete{$i} = $value;
			for($j=$i+1;$j<10;$j++){
				if(defined($hundred{$j}->{$key})){
					defined($delete{$j})&&die;
					$delete{$j} = $hundred{$j}->{$key};
				}
			}
			scalar(keys %delete) == 1 && next;
			scalar(keys %delete) < 1 &&die;
			@array = values %delete;
			#print "@array	";
			@array = sort {$b <=> $a} @array;
			if($array[0]>5*$array[1]){
				while(($key_delete,$value_delete) = each %delete){
					$value_delete == $array[0] && next;
					$key_delete == $i && $count++;
					delete $hundred{$key_delete}->{$key};
				}
			}else{
				while(($key_delete,$value_delete) = each %delete){
					$key_delete == $i && $count++;
					delete $hundred{$key_delete}->{$key};
				}
			}
		}
		#print $count,"	$i	",scalar keys($hundred{$i}),"\n";
	}else{die;}
}
for($i=0;$i<12;$i++){
	if(defined($col{$i})){
		$count = 0;
		while(($key,$value) = each %{$col{$i}}){
			%delete = ();
			$delete{$i} = $value;
			#print "$key	$value\n";
			for($j=$i+1;$j<12;$j++){
				if(defined($col{$j}->{$key})){
					defined($delete{$j})&&die;
					$delete{$j} = $col{$j}->{$key};
				}
			}
			scalar(keys %delete) == 1 && next;
			scalar(keys %delete) < 1 &&die;
			@array =values %delete;
			@array =  sort {$b <=> $a} @array;
			if($array[0]>5*$array[1]){
				while(($key_delete,$value_delete) = each %delete){
					$value_delete == $array[0] && next;
					$key_delete == $i && $count++;
					delete $col{$key_delete}->{$key};
				}
			}else{
				while(($key_delete,$value_delete) = each %delete){
					$key_delete == $i && $count++;
					delete $col{$key_delete}->{$key};
				}
			}
			
		}
		#print $count,"	$i	",scalar keys($col{$i}),"\n";
	}else{die;}
}
for($i=0;$i<8;$i++){
	if(defined($row{$i})){
		$count = 0;
		while(($key,$value) = each %{$row{$i}}){
			%delete = ();
			$delete{$i} = $value;
			for($j=$i+1;$j<8;$j++){
				if(defined($row{$j}->{$key})){
					defined($delete{$j})&&die;
					$delete{$j} = $row{$j}->{$key};
				}
			}
			scalar(keys %delete) == 1 && next;
			scalar(keys %delete) < 1 &&die;
			@array =values %delete;
			@array =  sort {$b <=> $a} @array;
			if($array[0]>5*$array[1]){
				while(($key_delete,$value_delete) = each %delete){
					$value_delete == $array[0] && next;
					$key_delete == $i && $count++;
					delete $row{$key_delete}->{$key};
				}
			}else{
				while(($key_delete,$value_delete) = each %delete){
					$key_delete == $i && $count++;
					delete $row{$key_delete}->{$key};
				}
			}
			
		}
		#print $count,"	$i	",scalar keys($row{$i}),"\n";
	}
}

my ($key_hundred,$key_ones,$key_tens,$key_row,$key_col,$value_hundred,$value_ones,$value_tens,$value_row,$value_col,%Flatfish,%Flatfish2)=();
open O,">Flatfish_BAC_library_3Dcross_cutoff_$cutoff"||die;
open O2,">Flatfish_BAC_library_3Dcross_cutoff_$cutoff".".log"||die;
while(($key_hundred,$value_hundred) = each %hundred){
	while(($key_ones,$value_ones) = each %ones){
		while(($key_tens,$value_tens) = each %tens){
			while(($key_row,$value_row) = each %row){
				while(($key_col,$value_col) = each %col){
					$count = 0;
					$temp2 = $key_col+1;
					$seq = 'Hundred'.$key_hundred.'Tens'.$key_tens.'Ones'.$key_ones.'Row'.chr($key_row+97).'Column'.$temp2;
					print O2 "Hundred",$key_hundred,"	Tens",$key_tens,"	Ones",$key_ones,"	Row",chr($key_row+97),"	Column",$key_col+1,"\n";
					print O "Hundred",$key_hundred,"	Tens",$key_tens,"	Ones",$key_ones,"	Row",chr($key_row+97),"	Column",$key_col+1,"\n";
					while(($barcode,$num) = each %{$value_col}){
						if(defined($hundred{$key_hundred}->{$barcode})&&defined($row{$key_row}->{$barcode})&&defined($ones{$key_ones}->{$barcode})&&defined($tens{$key_tens}->{$barcode})){
							if($num>$cutoff&&$hundred{$key_hundred}->{$barcode}>$cutoff&&$row{$key_row}->{$barcode}>$cutoff&&$ones{$key_ones}->{$barcode}>$cutoff&&$tens{$key_tens}->{$barcode}>$cutoff){
								$count++;
								print O $barcode,"	",min($hundred{$key_hundred}->{$barcode},$tens{$key_tens}->{$barcode},$ones{$key_ones}->{$barcode},$row{$key_row}->{$barcode},$num),"\n";
								#print O "$barcode	$plate{$key_plate}->{$barcode}\n";
								#print O "$barcode	$row{$key_row}->{$barcode}\n";
								#print O "$barcode	$num\n\n";
								delete $hundred{$key_hundred}->{$barcode};
								delete $ones{$key_ones}->{$barcode};
								delete $tens{$key_tens}->{$barcode};
								delete $row{$key_row}->{$barcode};
								delete $col{$key_col}->{$barcode};
								$temp = $barcode;
							}
						}
					}
					print O2 "	$count\n";
					if($count==1){
						$Flatfish{$seq} = $temp;
						$unique{$temp} = 1;
						$solved{$key_hundred}->{$key_ones}->{$key_tens}->{$key_row}->{$key_col} == 0 || die;
						$solved{$key_hundred}->{$key_ones}->{$key_tens}->{$key_row}->{$key_col} = 1;
					}else{
						$unsolved += $count;
					}
				}
			}
		}
	}
}
close O;
close O2;
open O3,">Flatfish_BAC_library_3Dcross_cutoff_$cutoff.out";
#print scalar(keys %Flatfish),"	",scalar(keys %unique),"\n",$unsolved,"\n";
while(($key,$value)= each %Flatfish){
	print O3 "$key	$value\n";
}
close O3;
