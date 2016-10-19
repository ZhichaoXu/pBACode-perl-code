#! usr/bin/perl
use List::Util qw(min);
my $cutoff = $ARGV[0];
my ($num_plate,$i,$j,$k,$plate384,$ab,$l,$temp,%countplate,%countcol,%countrow,$key,$value,$barcode,$num,$x,$y,$z,%count,$count,$temp2,$seq,%unique,%unique2,%delete,@array,$key_delete,$value_delete,$unsolved,%solved,%pool,$sum)=(16);
my (%plate,%col,%row)=();
for($i=0;$i<$num_plate;$i++){
	$plate384 = int($i/4)+1;
	$temp = $i%4;
	$ab = int($temp/2)+1;
	$l = $temp%2+1;
	#print $plate384,"	",chr($ab+96),"	",$l,"\n";
	open INPUT,"<barcodeYeast_BAC_library_".$plate384."-".chr($ab+96).$l.".extendedFrags.fastq.txt"||next;
	$countplate{$i} = 0;
	while(<INPUT>){
		/^([ACGTN]+)\t(\d+)/ || die;
		defined($plate{$i}->{$1})? die :($plate{$i}->{$1}=$2);
		$countplate{$i} += $2;
	}
	close INPUT;
}
#print scalar(keys %plate),"\n";
while(($key,$value) = each %plate){
	while(($barcode,$num) = each %{$value}){
		$plate{$key}->{$barcode} = 96*$num/$countplate{$key};
		delete $plate{$key}->{$barcode} if($plate{$key}->{$barcode}<=$cutoff);
	}
}
for($i=0;$i<12;$i++){
	$temp = $i+1;
	open INPUT,"<barcodeYeast_BAC_library_col-".$temp.".extendedFrags.fastq.txt"||next;
	$countcol{$i} = 0;
	while(<INPUT>){
		/^([ACGTN]+)\t(\d+)/ || die;
		defined($col{$i}->{$1})? die :($col{$i}->{$1}=$2);
		$countcol{$i} += $2;
	}
	close INPUT;
}
while(($key,$value) = each %col){
	#print "$key	$countcol{$key}\n";
	while(($barcode,$num) = each %{$value}){
		$col{$key}->{$barcode} = $num_plate*8*$num/$countcol{$key};
		#print "$barcode	$col{$key}->{$barcode}	$num\n";
		delete $col{$key}->{$barcode} if($col{$key}->{$barcode}<=$cutoff);
	}
}
for($i=0;$i<8;$i++){
	#print chr($i+97),"\n";
	open INPUT,"<barcodeYeast_BAC_library_row-".chr($i+97).".extendedFrags.fastq.txt"||next;
	$countrow{$i} = 0;
	while(<INPUT>){
		/^([ACGTN]+)\t(\d+)/ || die;
		defined($row{$i}->{$1})? die :($row{$i}->{$1}=$2);
		$countrow{$i} += $2;
	}
	close INPUT;
}
while(($key,$value) = each %row){
	#print "$key	$countrow{$key}\n";
	while(($barcode,$num) = each %{$value}){
		$row{$key}->{$barcode} = $num_plate*12*$num/$countrow{$key};
		#print "$barcode	$row{$key}->{$barcode}	$num\n";
		delete $row{$key}->{$barcode} if($row{$key}->{$barcode}<=$cutoff);
	}
}
for($i=0;$i<=15;$i++){
	for($j=0;$j<=7;$j++){
		for($k=0;$k<=11;$k++){
			#$k==1 && next;
			#$k==3 && next;
			$solved{$i}->{$j}->{$k} = 0;
		}
	}
}
%count = ();
#print "plate	delete	left\n";
for($i=0;$i<$num_plate;$i++){
	if(defined($plate{$i})){
		$plate384 = int($i/4)+1;
		$temp = $i%4;
		$ab = int($temp/2)+1;
		$l = $temp%2+1;
		#print $plate384."-".chr($ab+96).$l,"	";
		#print scalar keys($plate{$i}),"\n";
		while(($key,$value) = each %{$plate{$i}}){
			%delete = ();
			$delete{$i} = $value;
			#print "$key	$value\n";
			for($j=$i+1;$j<$num_plate;$j++){
				if(defined($plate{$j}->{$key})){
					#print "$key	$plate{$j}->{$key}\n";
					defined($delete{$j})&&die;
					$delete{$j} = $plate{$j}->{$key};
				}
			}
			if(scalar(keys %delete) == 1){
				#print "100\n";
				next;
			}
			scalar(keys %delete) < 1 &&die;
			#print scalar(keys %delete),"\n";
			@array = values %delete;
			#print "@array	";
			@array = sort {$b <=> $a} @array;
			$sum = 0;
			for($k=0;$k<=$#array;$k++){$sum += $array[$k];}
			#print 100*$array[0]/$sum,"\n";
			if(100*$array[0]/$sum>=60){
				#print "yes\n";
				while(($key_delete,$value_delete) = each %delete){
					$value_delete == $array[0] && next;
					$count{$key_delete}++;
					delete $plate{$key_delete}->{$key};
				}
			}else{
				#print "no\n";
				while(($key_delete,$value_delete) = each %delete){
					$count{$key_delete}++;
					delete $plate{$key_delete}->{$key};
				}
			}
		}
		#print $count,"	",scalar keys($plate{$i}),"\n";
	}
}
while(($key,$value) = each %count){
	#print $key,"	",$value,"	",scalar keys($plate{$key}),"\n";
}
#print scalar keys($plate{2}),"\n";
%count = ();
#print "column	delete	left\n";
for($i=0;$i<12;$i++){
	#print "col",$i+1,"	"if(defined($col{$i}));
	#print scalar(keys($col{$i})),"\n" if(defined($col{$i}));
	if(defined($col{$i})){
		while(($key,$value) = each %{$col{$i}}){
			%delete = ();
			$delete{$i} = $value;
			#print "$key	$value\n";
			for($j=$i+1;$j<12;$j++){
				if(defined($col{$j}->{$key})){
					#print "$key	$col{$j}->{$key}\n";
					defined($delete{$j})&&die;
					$delete{$j} = $col{$j}->{$key};
				}
			}
			if(scalar(keys %delete) == 1){
				#print "100\n";
				next;
			}
			scalar(keys %delete) < 1 &&die;
			@array =values %delete;
			#print "@array	";
			@array =  sort {$b <=> $a} @array;
			$sum = 0;
			for($k=0;$k<=$#array;$k++){$sum += $array[$k];}
			#print 100*$array[0]/$sum,"\n";
			if(100*$array[0]/$sum>=60){
				#print "yes\n";
				while(($key_delete,$value_delete) = each %delete){
					$value_delete == $array[0] && next;
					$count{$key_delete}++;
					delete $col{$key_delete}->{$key};
				}
			}else{
				#print "no\n";
				while(($key_delete,$value_delete) = each %delete){
					$count{$key_delete}++;
					delete $col{$key_delete}->{$key};
				}
			}
			
		}
		#print $count,"	",scalar keys($col{$i}),"\n";
	}
}
while(($key,$value) = each %count){
	#print $key,"	",$value,"	",scalar keys($col{$key}),"\n";
}
%count = ();
#print "row	delete	left\n";
for($i=0;$i<8;$i++){
	#print "row",$i+1,"	"if(defined($row{$i}));
	#print scalar(keys($row{$i})),"\n" if(defined($row{$i}));
	if(defined($row{$i})){
		while(($key,$value) = each %{$row{$i}}){
			%delete = ();
			$delete{$i} = $value;
			#print "$key	$value\n";
			for($j=$i+1;$j<8;$j++){
				if(defined($row{$j}->{$key})){
					#print "$key	$row{$j}->{$key}\n";
					defined($delete{$j})&&die;
					$delete{$j} = $row{$j}->{$key};
				}
			}
			if(scalar(keys %delete) == 1){
				#print "100\n";
				next;
			}
			scalar(keys %delete) < 1 &&die;
			@array =values %delete;
			#print "@array	";
			@array =  sort {$b <=> $a} @array;
			$sum = 0;
			for($k=0;$k<=$#array;$k++){$sum += $array[$k];}
			#print 100*$array[0]/$sum,"\n";
			if(100*$array[0]/$sum>=60){
				#print "yes\n";
				while(($key_delete,$value_delete) = each %delete){
					$value_delete == $array[0] && next;
					$count{$key_delete}++;
					delete $row{$key_delete}->{$key};
				}
			}else{
				#print "no\n";
				while(($key_delete,$value_delete) = each %delete){
					$count{$key_delete}++;
					delete $row{$key_delete}->{$key};
				}
			}
			
		}
		#print $count,"	",scalar keys($row{$i}),"\n";
	}
}
while(($key,$value) = each %count){
	#print $key,"	",$value,"	",scalar keys($row{$key}),"\n";
}
my ($key_plate,$key_row,$key_col,$value_plate,$value_row,$value_col,%yeast,%yeast2)=();
open O,">Yeast_BAC_library_3Dcross_cutoff_$cutoff"||die;
open O2,">Yeast_BAC_library_3Dcross_cutoff_$cutoff".".log"||die;
while(($key_plate,$value_plate) = each %plate){
	while(($key_row,$value_row) = each %row){
		while(($key_col,$value_col) = each %col){
			$count = 0;
			$plate384 = int($key_plate/4)+1;
			$temp = $key_plate%4;
			$ab = int($temp/2)+1;
			$l = $temp%2+1;
			$temp2 = $key_col+1;
			$seq = 'Plate'.$plate384.chr($ab+96).$l.'Row'.chr($key_row+97).'Column'.$temp2;
			print O2 "Plate",$plate384."-".chr($ab+96).$l,"	Row",chr($key_row+97),"	Column",$key_col+1;
			print O "Plate",$plate384."-".chr($ab+96).$l,"	Row",chr($key_row+97),"	Column",$key_col+1,"\n";
			while(($barcode,$num) = each %{$value_col}){
				if(defined($row{$key_row}->{$barcode})&&defined($plate{$key_plate}->{$barcode})){
					if($num>$cutoff&&$row{$key_row}->{$barcode}>$cutoff&&$plate{$key_plate}->{$barcode}>$cutoff){
					$count++;
					print O $barcode,"	",min($plate{$key_plate}->{$barcode},$row{$key_row}->{$barcode},$num),"\n";
					delete $plate{$key_plate}->{$barcode};
					delete $row{$key_row}->{$barcode};
					delete $col{$key_col}->{$barcode};
#					print O "$barcode	$plate{$key_plate}->{$barcode}\n";
#					print O "$barcode	$row{$key_row}->{$barcode}\n";
#					print O "$barcode	$num\n\n";
					$temp = $barcode;
					}
				}
			}
			print O2 "	$count\n";
			if($count==1){
				$yeast{$seq} = $temp;
				$unique{$temp} = 1;
				$solved{$key_plate}->{$key_row}->{$key_col} == 0 || die;
				$solved{$key_plate}->{$key_row}->{$key_col} = 1;
			}else{
				$unsolved += $count;
			}
		}
	}
}
close O;
close O2;
open O3,">Yeast_BAC_library_3Dcross_cutoff_$cutoff.out";
while(($key,$value)= each %yeast){
	print O3 "$key	$value\n";
}
close O3;

