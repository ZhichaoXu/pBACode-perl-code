#!/usr/bin/perl -w
#A package for constants and subroutines

package hiseq;
BEGIN{};
use List::Util qw(min);
use List::Util qw(max);

sub betterMatch($$$);
sub bowtieMatchLeft($$); #INPUT: sam output item 6
                         #       maximal mismatch at left end telerate, if more than $2, return 0
                         #OUTPUT: the number of bp left end matching reference. If there is deletion and insertion between two M and/or a few bp mistach at the left end, ignore these unmatch and report the length of two Ms
sub bowtieMatchRight($$);
sub checkSAM($$$);
sub checkSite($$$$$$); #Check whether a restriction site exists in expected location. The site can differ from expected one by 1bp mismatch if the length of barcode sequence is fixed $4
                       #INPUT: $1: sequence
                       #       $2: restrictione site
                       #       $3: length of bridge sequence between barcode and restriction site
                       #       $4: barcode length I will assume in case there is no perfect restriction site
                       #       $5: maximaum length of barcode in case there is a perfect restriction site $2
                       #       $6: barcode is at left or right
                       #OUTPUT: $_[0]: INPUT $1 after chopping barcode and bridge
                       #        $_[1]: barcode
sub convertTOfa($); #INPUT: name of a .txt file
                    #OUTPUT: a .fa file
sub convertfqTOfa($$); #INPUT: $1 name of a .fq file
                       #       $2 name of output fa file, $1 by default
                       #OUTPUT: print a fa file
sub convertHashTOHashOfArray(%);
sub differby1D($$); #INPUT: $1 is 1 bp shorter than $2
sub diffnomorethan($$$);
sub localDistance($$$); #INPUT: string1
                        #       string2
                        #       size-of-maximal-indel-considered, 2 by default
                        #OUTPUT: 0-based left location of the rightest best match
                        #        distance score
sub localDistanceR($$$);
#sub locationReadHindIII($$$);
sub locationRead($$$$);
sub locationReadMultiple($$$$); #OUTPUT: hash: key: $read; value: "mappingScore\tlocation_mappingLength"
                                #        in case a read aligned to multiple loci(i.e., in *_256.sam or *_323.sam), value: an array of "mappingScore\tlocation_mappingLength", whose last element is the best alignment(i.e., in *_0.sam or *_67.sam)
sub locnomorethan1indelxS($$$$);#priority substituation > deletion > indertion
                                #INPUT: $1 query sequence
                                #       $2 reference sequence
                                #       $3 threshold of mismatch numbers 
                                #       $4 'l' or 'r'
                                #OUTPUT: if 2 strings have more than $3 mismatches or more than 1 indel, OUTPUT$1=0; otherwise, return OUTPUT$1=1
                                #        if input$4 = 'l', location of the left character in input$1 matching the reference
                                #        if input$4 = 'r', location of the right character in input$1 matching the reference

sub matchScore($$$); #number of matching bp - number of mismatch - 2*number of indel. If there is a mismatch at no more than $3bp from 3' end, Score - $3 - 1
                     #INPUT:$1: sequence
                     #      $2: sam column 5 to express how many insertions
                     #      $3: if a mismatch occurs no more than $3bp from 3' end, dicard these bps and the mismatch one
sub mismatch1($%);
sub nomorethan1indel($$); #INPUT: 2 strings of identical length
sub nomorethan1indelxS($$$); #INPUT: 2 strings of identical length
                             #       $3 threshold of mismatch numbers 
                             #OUTPUT: if 2 strings have more than $3 mismatches or more than 1 indel, return 0; otherwise, return 1
sub nomorethan1internalSD($$);
sub nomorethan1S($$); #INPUT: 2 strings of identical length
                      #OUTPUT: if 2 strings have no more than 1 mismatch, 1; otherwise, 0
sub nomorethanxS($$$);
sub oneindelmismatch($$$$);
sub readFAfile($); #INPUT: name of a .txt or .fa file
                   #OUTPUT: a string of [ACGTN]
sub revcomDNA($); #reverse complementary DNA
sub revcomFASTQ($$); #reverse complement DNA and QC in FASTQ file $1 and save to file $2
sub rmHybrid($%); #remove hybrids. A hybrid is a pair of barcodes whose read number is LOWER than those of each one pairing with other barcodes
                  #in case of A10B, C1B, C2D, discard C1B. But in case of A10B, C1B, C1D, do not discard anything
                  #If a hybrid exists as a barcode pair in pre-pool $1, print warning
sub rmHybrid1(%); #remove hybrids. A hybrid is a pair of barcodes whose read number is LOWER than those of each one pairing with other barcodes
                  #in case of A10B, C1B, C2D, discard C1B. But in case of A10B, C1B, C1D, do not discard anything
sub rmN($%); #For non 1:1 barcode-pair, only keep the most enriched barcode-pair
             #If a non-most-enriched barcode-pair exists as a barcode pair in pre-pool $1, print warning
sub rmLend($$$$$); #INPUT $1 query sequence
                   #      $2 reference sequence
#If length($1) == length($2), it can't detect insertion in $2
                   #      $3 distances of left end of the fragment from left end of $1 to align to $2
                   #      $4 length of $1 fragment to substr() to align to $2
                   #      $5 maximal mismatch allowed, 1 by default
                   #OUTPUT location of right end of INPUT$2 aligned to INPUT$1, 0 if there are > INPUT$5 mismatches or >1bp INDEL
sub rmLendBridge($$$$$$);
sub rmRend($$$$$); #INPUT $1 query sequence
                   #      $2 reference sequence
#If length($1) == length($2), it can't detect insertion in $2
                   #      $3 distances of left end of the fragment from right end of $1 to align to $2
                   #      $4 length of $1 fragment to substr() to align to $2
                   #      $5 maximal mismatch allowed, 1 by default
                   #OUTPUT location of left end of INPUT$2 aligned to INPUT$1, 0 if there are > INPUT$5 mismatches or >1bp INDEL
sub rmRendBridge($$$$$$);
sub rmSite($$$$$$); #OUTPUT: remove restriction-site-like $2. restriction-site-like: no more than 1bp mismatch from $2
sub SQmean($); #INPUT: a string illumina sequencing report
               #OUTPUT: mean of Quality Score of all characters
sub SQtrimEnd3($$$);
sub strDistance($$);

sub betterMatch($$$) {
    my ($match1, $match2, $threshold) = @_;
    my ($mq, $mq2) = ();
    $match1 eq $match2 && return 0;

    if($match1 =~ /^(\d+)$/) {
	$mq = $1;
	$match2 =~ /^(\d+)$/ || die;

	if($mq - $1 >= $threshold) {
	    return 1;
	} elsif($1 - $mq >= $threshold) {
	    return -1;
	} else {
	    return 0;
	}
    } elsif($match1 =~ /^(\d+)_(\d+)$/) {
	$mq = $1;
	$mq2 = $2;
	$match2 =~ /^(\d+)_(\d+)$/ || die;

	if($mq+$mq2-$1-$2 >= $threshold && $mq >= $1 && $mq2 >= $2) {
	    return 1;
	} elsif($mq+$mq2-$1-$2 <= $threshold && $mq <= $1 && $mq2 <= $2) {
	    return -1;
	} else {
	    return 0;
	}
    } else {
	die;
    }
}

#sub betterMatch($$$) {
#    my ($match1, $match2, $threshold) = @_;
#    my ($mq, $mq2) = ();
#    $match1 eq $match2 && return 0;

#    if($match1 =~ /^([ACGT\^\d]+)_(\d+)$/) {
#	$mq = $2;
#	$match2 =~ /^([ACGT\^\d]+)_(\d+)$/ || die;

#	if($mq - $2 >= $threshold) {
#	    return 1;
#	} elsif($2 - $mq >= $threshold) {
#	    return -1;
#	} else {
#	    return 0;
#	}
#    } elsif($match1 =~ /^[ACGT\^\d]+_(\d+)\*[ACGT\^\d]+_(\d+)$/) {
#	$mq = $1;
#	$mq2 = $2;
#	$match1 =~ /^[ACGT\^\d]+_(\d+)\*[ACGT\^\d]+_(\d+)$/ || die;

#	if($mq+$mq2-$1-$2 >= $threshold && $mq >= $1 && $mq2 >= $2) {
#	    return 1;
#	} elsif($mq+$mq2-$1-$2 <= $threshold && $mq <= $1 && $mq2 <= $2) {
#	    return -1;
#	} else {
#	    return 0;
#	}
#    } else {
#	die;
#    }
#}

sub bowtieMatchLeft($$) {
    my($sam, $maxS) = @_;
    my($t) = (0);

#    $sam =~ /\d\d[ID]/ && die;
    $sam =~ s/\d[ID]//g;
    if($sam =~ s/^(\d+)S//) {
	$1 > $maxS && return 0;
    }

    while($sam =~ s/^(\d+)M//) {
	$t += $1;
    }

    return $t;
}

sub bowtieMatchRight($$) {
    my($sam, $maxS) = @_;
    my($t) = (0);

#    $sam =~ /\d\d[ID]/ && print "$sam\n";
    $sam =~ s/\d[ID]//g;
    if($sam =~ s/(\d+)S$//) {
	$1 > $maxS && return 0;
    }

    while($sam =~ s/(\d+)M$//) {
	$t += $1;
    }

    return $t;
}

sub checkSAM($$$) {
    my ($sam, $fragLen, $line) = @_;
    my @lines = split "\t", $line;

    $lines[1] != $sam && die "$lines[1]*$sam $line*";
#    $sam ne '321' && $lines[5] =~ /\d\d[SID]/i && die "$_$lines[5]\n", $FILE."_$sam.sam";

	if($sam =~ /^67|131|323$/) {
	    /CP\s*$/ || die "$sam";
	    $lines[6] eq '=' || die "$sam*$lines[6]";
	    $lines[8] > $fragLen && die;
	} elsif($sam =~ /^65|129|321$/) {
	    /CP\s*$/ && die "$sam";
	    /UP\s*$/ && next;
	    /DP\s*$/ || die "$sam";
	    $lines[6] eq '=' || next;
	    $lines[8] > $fragLen && next;
	} elsif($sam =~ /^0|64|256$/) {
	    $line =~ /UU\s*$/ || die "wrong end $line\n";
	} else {
	    die;
	}

    $lines[0] =~ /^([^\s]+)(BARCODE)([ACGTN]+)/ || die;
    return ($1.$2.$3, @lines);
}

sub checkSite($$$$$$) {
    my($read, $SITE, $BRIDGELEN, $BARCODELEN, $BARCODELENMAX, $LorR) = @_;
    $SITE || return ();
    my $SITELEN = length($SITE);

    if($LorR eq 'l') {
	if($read =~ s/^([ACGTN]{1,$BARCODELENMAX}).{$BRIDGELEN}($SITE)/$2/) {
#	    length($1) > 20 && die "$1\t$BRIDGELEN\t$SITE";
	    return($read, $1);
	} elsif($read =~ s/^([ACGTN]{$BARCODELEN}).{$BRIDGELEN}([ACGTN]{$SITELEN})/$2/) {
	    if(&nomorethan1S($2, $SITE)) {
		return($read, $1);
	    } else {
		return('');
	    }
	} else {
	    die "$read\n$SITE, $BRIDGELEN, $BARCODELENMAX";
	}
    } elsif($LorR eq 'r') {
	if($read =~ s/($SITE).{$BRIDGELEN}([ACGTN]{1,$BARCODELENMAX})$/$1/) {
	    return($read, $2);
	} elsif($read =~ s/([ACGTN]{$SITELEN}).{$BRIDGELEN}([ACGTN]{$BARCODELEN})$/$1/) {
	    if(&nomorethan1S($1, $SITE)) {
		return($read, $2);
	    } else {
		return();
	    }
	} else {
	    die "$read\n$SITE, $BRIDGELEN, $BARCODELENMAX";
	}
    } else {
	die;
    }
}

sub convertfqTOfa($$) {
    my($file, $output) = @_;
    $output || ($output = $file);
    my(@lines) = ();

    open(I, $file.".fq") || open(I, $file.".fastq") || die;
    open(O, ">$output.fa") || die;
    while(<I>) {
	push @lines,$_;

	if(@lines == 4){
	    print O ">$lines[0]$lines[1]";
	    @lines = ();
	} elsif(@lines > 4) {
	    die;
	}
    }
}

sub convertHashTOHashOfArray(%) {
    my %h = @_;
    my (%output, $k, $v) = ();

    while(($k, $v) = each(%h)) {
	exists($output{$k}) ? die : ($output{$k}->[0] = $v);
    }

    return %output;
}

sub convertTOfa($) {
    my($file) = shift;
    my($firstLine) = (0);

   open(I, $file.".txt") || die;
    open(O, ">$file.fa") || die;
    while(<I>) {
	if(/^>/) {
	    $firstLine && die;
	    print O;
	    $firstLine++;
	} elsif(s/^\s*\d*\s*([ACGTN\s]+)$/$1/i) {
	    s/\s//g;

	    unless($firstLine) {
		$firstLine++;
		print O '>'.$file."\n";
	    }

	    print O uc;
	} else {
	    /^\s*$/ || die;
	}
    }
}

sub differby1D($$) {
    my ($query, $reference) = @_;
    my(@m, @mask, $m, $mask) = ();
    my $len = length($query);

    if(length($reference) - length($query) == 1) {
	$mask = $query ^ $reference;
	$mask =~ /^([\0]+)/ ? ($m = length($1)) : ($m = 0);
	$mask = $query ^ substr($reference, 1);

	if($mask =~ /([\0]+)$/) {
	    $m + length($1) < $len || return 1;
	} else {
	    $m < $len || return 1;
	}
    } else {
	die;
    }

    return 0;
}

sub diffnomorethan($$$) {
    my ($a, $b, $threshold) = @_;

    my $cmp = $a ^ $b; # bitwise exclusive or (XOR) of the two strings; the resulting string will have bytes of 0 at all locations where the two letters are the same

    my $mis = 0;

    while ($cmp =~ /[^\0]/g) {   #match non-zero byte only
	$mis++;
	$mis > $threshold && return 0;
    }

    return 1;
}

sub localDistance($$$) {
    my ($str1, $str2, $indel) = @_;
    $indel || ($indel = 0);
    my ($i, $mini, $window, @distance) = ();

    length($str1) > length($str2) && (($str1, $str2) = ($str2, $str1));
    $window = length($str1) + $indel;
    if(length($str2) > $window) {
	for($i = 0; $i < length($str2) - $window; $i++) {
	    $distance[$i] = &strDistance($str1, substr($str2, $i, $window));
	}
	$mini = min(@distance);

	for($i = length($str2) - $window - 1; $i >= 0; $i--) {
	    $distance[$i] > $mini || return ($i, $mini);
	}
    } else {
	return (0, &strDistance($str1, $str2));
    }
}

sub localDistanceR($$$) {

    my ($str1, $str2, $indel) = @_;
    $indel || ($indel = 0);
    my ($i, $mini, $window, @distance) = ();

    length($str1) > length($str2) && (($str1, $str2) = ($str2, $str1));

    $window = length($str1) + $indel;
    if(length($str2) > $window) {
	for($i = 0; $i < length($str2) - $window; $i++) {
	    $distance[$i] = &strDistance($str1, substr($str2, $i, $window));
#print "$i\t$distance[$i]\t", substr($str2, $i, $window), "\n";
	}
	$mini = min(@distance);

	for($i = 0;$i<=length($str2) - $window - 1;  $i++) {
	    $distance[$i] > $mini || return ($i, $mini);
	}
    } else {
	return (0, &strDistance($str1, $str2));
    }
}

sub locationRead($$$$) {
    my ($FILE, $sam, $HINDIII, $fragLen) = @_;
    my (@lines, %locationRead, %locationRead2) = ();
    my ($line, $location, $match, $read, $t) = ();
    $HINDIII =~ /^[y\d]$/ || die $HINDIII;

    open(SAM, $FILE."_$sam.sam") || die $FILE."_$sam.sam";
    while(<SAM>) {
        $line = $_;
        @lines = split;
        $lines[0] =~ /^([^\s]+)(BARCODE)([ACGTN]+)/ || die;
        $read = $1.$2.$3;
        $lines[1] != $sam && die;
        $sam ne '321' && $lines[5] =~ /\d\d[SID]/i && die "$_$lines[5]\n", $FILE."_$sam.sam";
        /10990\:16360/ && print $FILE."_$sam.sam\t$read\n";
#       ($sam eq '256' && $read =~ /10990\:16360/) &&  print $FILE."_$sam.sam\t$read\n";
#       $read =~ /10990\:16360/ && print $FILE."_$sam.sam\t$read";
        if($sam =~ /^67|323$/) {
            /CP\s*$/ || die "$_$FILE";
            $lines[6] eq '=' || die "$_$sam*$lines[6]";
            $lines[8] > $fragLen && die;
        } elsif($sam =~ /^65|321$/) {
            /CP\s*$/ && die "$_*";
            /UP\s*$/ && next;
            /DP\s*$/ || die "$_$FILE";
            $lines[6] eq '=' || next;
            $lines[8] > $fragLen && next;
        } elsif($sam =~ /^0|256$/) {
            $line =~ /UU\s*$/ || die "wrong end $line\n";
        } else {
            die;
        }

        if($sam =~ /^0|256$/) {
            if($lines[5] =~ /^(\d)S\d(M\d[DI])?\d+M/) {
                $t = $1;
                $HINDIII =~ /^\d$/ || next;
                $t > $HINDIII && next; # <=6bp mismatch in noHindIII insertion site is tolerable
            } else {
                $lines[5] =~ /^\d(M\d[DI])?\d+M/ || die; #at least 10bp perfect match at HindIII side, <10bp indel allowed
            }

            $lines[11] =~ /^AS\:i\:(\d+)$/ || die;
            $t = "$1\t$lines[2]"."_$lines[3]";

            if($sam) {
                $locationRead{$read}->[@{$locationRead{$read}}] = $t;
            } else {
                exists($locationRead{$read}) ? die : ($locationRead{$read}->[0] = $t);
            }
        } elsif($sam =~ /^67|323|65|321$/) {
            $lines[8] =~ /^\-(\d+)$/ && next;
            $lines[8] < 20 && die "$read\t", $FILE."_$sam*$lines[8]";

            if($lines[5] =~ /^(\d)S\d(M\d[DI])?\d+M/) {
                $t = $1;
                $HINDIII =~ /^\d$/ || next;
                $t > $HINDIII && next; # <=6bp mismatch in noHindIII insertion site is tolerable
            } else {
                $lines[5] =~ /^\d(M\d[DI])?\d+M/ || die; #at least 10bp perfect match at HindIII side, <10bp indel allowed
            }

            $lines[11] =~ /^AS\:i\:(\d+)$/ || die;
            $t = $1;
            $lines[$#lines-1] =~ /^YS\:i\:(\d+)$/ || die;
            $t .= "_$1\t$lines[2]"."_$lines[3]"."_$lines[8]";

            if($sam =~ /^67|65$/) {
                exists($locationRead{$read}) ? die : ($locationRead{$read}->[0] = $t);
            } elsif($sam =~ /^323|321$/) {
                $locationRead{$read}->[@{$locationRead{$read}}] = $t;
            } else {
                die;
            }
        } else {
            die;
        }
    }

    return %locationRead;
}

sub locationReadMulti($$$$) {
    my ($FILE, $sam, $HINDIII, $fragLen) = @_;
    my $MAXMISMATCH = 2;
    my (@lines, %locationRead) = ();
    my ($align, $line, $location, $match, $md, $mismatch, $multiMismatch, $read, $t) = ();
    $HINDIII =~ /^\d+$/ || die $HINDIII;
    if($sam eq '0') {
	%locationRead = &locationReadMulti($FILE, $sam+256, $HINDIII, $fragLen);
	open(SAM, $FILE."_$sam.sam");
	while(<SAM>) {
	    ($read, @lines) = &checkSAM($sam, $fragLen, $_);

	    if($lines[5] =~ /\d\d[SID]/) { #less than 10bp mismatch patch is allowed
		exists($locationRead{$read}) && delete($locationRead{$read});
	    } else {
		$t = 0;
		$lines[5] =~ s/^(\d+)S// && ($t = $1);
		$lines[5] =~ s/(\d+)S$// && ($t += $1);
		$multiMismatch = 0;

		if(exists($locationRead{$read})) {
		    if($t > $HINDIII) {
			$multiMismatch = 1; # <=6bp mismatch in noHindIII insertion site is tolerable
		    } else {
			$lines[$#lines - 1] =~ s/^MD\:Z\:// || die;
			$mismatch = 0;
			while($lines[$#lines - 1] =~ s/[\^ACGTN]//) {
			    ++$mismatch;
			}
			$lines[$#lines - 1] =~ /^\d+$/ || die "$read,$lines[$#lines - 1]\n" ;
			$mismatch > $MAXMISMATCH && ($multiMismatch = 1); #too many mismatches even in the best hit of multi-hits
		    }
		}

		if($multiMismatch) {
		    delete($locationRead{$read});
		} else {
		    $t = 0;
		    while($lines[5] =~ s/(\d+)[MSD]//) {
			$t += $1;
		    }
		    while($lines[5] =~ s/(\d+)I//) {
			$t -= $1;
		    }
		    $t < 10 && die $read;
		    $lines[5] && die $read;
		    $lines[11] =~ /^AS\:i\:[-]?(\d+)$/ || die;
		    $locationRead{$read}->[@{$locationRead{$read}}] = "$1\t$lines[2]"."_$lines[3]"."_$t";
		}
	    }
	}
	close SAM;
    } elsif($sam eq '256') {
	%locationRead = ();
	open(SAM, $FILE."_$sam.sam") || return ();
	while(<SAM>) {
	    ($read, @lines) = &checkSAM($sam, $fragLen, $_);
	    $lines[5] =~ s/^(\d+)S//;
	    $lines[5] =~ s/(\d+)S$//;
	    $t = 0;
		while($lines[5] =~ s/(\d+)[MSD]//) {
		    $t += $1;
		}
		while($lines[5] =~ s/(\d+)I//) {
		    $t -= $1;
		}
	    $t < 10 && die $read;
	    $lines[5] && die $read;
	    $lines[11] =~ /^AS\:i\:[-]?(\d+)$/ || die;
	    $locationRead{$read}->[@{$locationRead{$read}}] = "$1\t$lines[2]"."_$lines[3]"."_$t";
	}
	close SAM;
    } elsif($sam eq '67') {
	%locationRead = &locationReadMulti($FILE, $sam+256, $HINDIII, $fragLen);
	open(SAM, $FILE."_$sam.sam");
	while(<SAM>) {
	    ($read, @lines) = &checkSAM($sam, $fragLen, $_);
	    $lines[6] eq '=' || die;
	    if($lines[5] =~ /\d\d[SID]/) { #less than 10bp mismatch patch is allowed
		exists($locationRead{$read}) && delete($locationRead{$read});
	    } elsif($lines[8] < 10) {
		exists($locationRead{$read}) && delete($locationRead{$read});
	    } else {
		$multiMismatch = 0;

		if(exists($locationRead{$read})) {
		    $t = 0;
		    $lines[5] =~ s/^(\d+)S// && ($t = $1);
		    $lines[5] =~ s/(\d+)S$// && ($t += $1);

		    if($t > $HINDIII) {
			$multiMismatch = 1; # <=6bp mismatch in noHindIII insertion site is tolerable
		    } else {
			$lines[$#lines - 2] =~ s/^MD\:Z\:// || die $lines[$#lines - 2];
			$mismatch = 0;
			while($lines[$#lines - 2] =~ s/[\^ACGTN]//) {
			    ++$mismatch;
			}
			$lines[$#lines - 2] =~ /^\d+$/ || die "$_*";
			$mismatch > $MAXMISMATCH && ($multiMismatch = 1); #too many mismatches even in the best hit of multi-hits
		    }
		}

		if($multiMismatch) {
		    delete($locationRead{$read});
		} else {
		    $lines[11] =~ /^AS\:i\:(\d+)$/ || die;
		    $t = $1;
		    $lines[$#lines-1] =~ /^YS\:i\:(\d+)$/ || die;
		    $locationRead{$read}->[@{$locationRead{$read}}] = $t."_$1\t$lines[2]"."_$lines[3]"."_$lines[8]";
		}
	    }
	}
	close SAM;
	    open(SAM, $FILE."_".($sam+64).".sam");
	    while(<SAM>) {
		($read, @lines) = &checkSAM($sam+64, $fragLen, $_);
		exists($locationRead{$read}) || next;
		$lines[8] < 0 || die "$read\t$_$lines[8]";

		if($lines[5] =~ /\d\d[SID]/) { #less than 10bp mismatch patch is allowed
		    delete($locationRead{$read});
		} else {
		    $t = 0;
		    $lines[5] =~ s/^(\d+)S// && ($t = $1);
		    $lines[5] =~ s/(\d+)S$// && ($t += $1);
		    if($t > $HINDIII) {
			delete($locationRead{$read}); # <=6bp mismatch in noHindIII insertion site is tolerable
		    } else {
			if($lines[$#lines - 1] =~ s/^MD\:Z\://) {
			    $md = $lines[$#lines - 1];
			} elsif($lines[$#lines - 2] =~ s/^MD\:Z\://) {
			    $md = $lines[$#lines - 2];
			} else {
			    die "$_$lines[$#lines - 3]";
			}

			$mismatch = 0;
			while($md =~ s/[\^ACGTN]//) {
			    ++$mismatch;
			}
			$md =~ /^\d+$/ || die $read;
			$mismatch > $MAXMISMATCH && delete($locationRead{$read}); #too many mismatches even in the best hit of multi-hits
		    }
		}
	    }
		close SAM;
    } elsif($sam eq '65') {
	%locationRead = &locationReadMulti($FILE, $sam+256, $HINDIII, $fragLen);
	open(SAM, $FILE."_$sam.sam");
	while(<SAM>) {
	    ($read, @lines) = &checkSAM($sam, $fragLen, $_);

	    if($lines[5] =~ /\d\d[SID]/) { #less than 10bp mismatch patch is allowed
		exists($locationRead{$read}) && delete($locationRead{$read});
	    } elsif($lines[6] ne '=') {
		exists($locationRead{$read}) && delete($locationRead{$read});
	    } elsif($lines[8] < 10 || $lines[8] > $fragLen) {
		exists($locationRead{$read}) && delete($locationRead{$read});
	    } else {
		$multiMismatch = 0;

		if(exists($locationRead{$read})) {
		    $t = 0;
		    $lines[5] =~ s/^(\d+)S// && ($t = $1);
		    $lines[5] =~ s/(\d+)S$// && ($t += $1);

		    if($t > $HINDIII) {
			$multiMismatch = 1; # <=6bp mismatch in noHindIII insertion site is tolerable
		    } else {
			$lines[$#lines - 2] =~ s/^MD\:Z\:// || die $lines[$#lines - 2];
			$mismatch = 0;
			while($lines[$#lines - 2] =~ s/[\^ACGTN]//) {
			    ++$mismatch;
			}
			$lines[$#lines - 2] =~ /^\d+$/ || die "$_*";
			$mismatch > $MAXMISMATCH && ($multiMismatch = 1); #too many mismatches even in the best hit of multi-hits
		    }
		}

		if($multiMismatch) {
		    delete($locationRead{$read});
		} else {
		    $lines[11] =~ /^AS\:i\:(\d+)$/ || die;
		    $t = $1;
		    $lines[$#lines-1] =~ /^YS\:i\:(\d+)$/ || die;
		    $locationRead{$read}->[@{$locationRead{$read}}] = $t."_$1\t$lines[2]"."_$lines[3]"."_$lines[8]";
		}
	    }
	}
	close SAM;
	    open(SAM, $FILE."_".($sam+64).".sam");
	    while(<SAM>) {
		($read, @lines) = &checkSAM($sam+64, $fragLen, $_);
		exists($locationRead{$read}) || next;
		$lines[8] < 0 || die "$read\t$_$lines[8]";

		if($lines[5] =~ /\d\d[SID]/) { #less than 10bp mismatch patch is allowed
		    delete($locationRead{$read});
		} else {
		    $t = 0;
		    $lines[5] =~ s/^(\d+)S// && ($t = $1);
		    $lines[5] =~ s/(\d+)S$// && ($t += $1);

		    if($t > $HINDIII) {
			delete($locationRead{$read}); # <=6bp mismatch in noHindIII insertion site is tolerable
		    } else {
			if($lines[$#lines - 1] =~ s/^MD\:Z\://) {
			    $md = $lines[$#lines - 1];
			} elsif($lines[$#lines - 2] =~ s/^MD\:Z\://) {
			    $md = $lines[$#lines - 2];
			} else {
			    die "$_$lines[$#lines - 3]";
			}

			$mismatch = 0;
			while($md =~ s/[\^ACGTN]//) {
			    ++$mismatch;
			}
			$md =~ /^\d+$/ || die $read;
			$mismatch > $MAXMISMATCH && delete($locationRead{$read}); #too many mismatches even in the best hit of multi-hits
		    }
		}
	    }
		close SAM;
    } elsif($sam eq '323') {
	%locationRead = ();
	open(SAM, $FILE."_$sam.sam") || return ();
	while(<SAM>) {
	    ($read, @lines) = &checkSAM($sam, $fragLen, $_);
	    $lines[8] || die "$read\t$FILE"."_$sam.sam";
	    $lines[11] =~ /^AS\:i\:(\d+)$/ || die;
	    $t = $1;
	    $lines[$#lines-1] =~ /^YS\:i\:(\d+)$/ || die;
	    $locationRead{$read}->[@{$locationRead{$read}}] = $t."_$1\t$lines[2]"."_$lines[3]"."_$lines[8]";
	}
	close SAM;
    } elsif($sam eq '321') {
	%locationRead = ();
	open(SAM, $FILE."_$sam.sam") || return ();
	while(<SAM>) {
	    ($read, @lines) = &checkSAM($sam, $fragLen, $_);
	    if($lines[6] eq '=') {
		if($lines[8] < $fragLen) {
		    $lines[11] =~ /^AS\:i\:(\d+)$/ || die;
		    $t = $1;
		    $lines[$#lines-1] =~ /^YS\:i\:(\d+)$/ || die;
		    $locationRead{$read}->[@{$locationRead{$read}}] = $t."_$1\t$lines[2]"."_$lines[3]"."_$lines[8]";
		}
	    } else {
		$lines[8] && die;
	    }
	}
	close SAM;
    } else {
	die;
    }

    return %locationRead;
}

sub locationReadNoRestriction($$$) {
    my ($FILE, $sam, $fragLen) = @_;
    my (@lines, %locationRead, %locationRead2) = ();
    my ($align, $line, $location, $match, $read, $t) = ();

    open(SAM, $FILE."_$sam.sam") || die $FILE."_$sam.sam";
    while(<SAM>) {
	$line = $_;
	@lines = split;
	$lines[0] =~ /^([^\s]+)(BARCODE)([ACGTN]+)/ || die;
	$read = $1.$2.$3;
	$lines[1] != $sam && die;
	$sam ne '321' && $lines[5] =~ /\d\d[SID]/i && die "$_$lines[5]\n", $FILE."_$sam.sam";

	if($sam =~ /^67|323$/) {
	    /CP\s*$/ || die "$_$FILE";
	    $lines[6] eq '=' || die "$_$sam*$lines[6]";
	    $lines[8] > $fragLen && die;
	} elsif($sam =~ /^65|321$/) {
	    /CP\s*$/ && die "$_*";
	    /UP\s*$/ && next;
	    /DP\s*$/ || die "$_$FILE";
	    $lines[6] eq '=' || next;
	    $lines[8] > $fragLen && next;
	} elsif($sam =~ /^0|256$/) {
	    $line =~ /UU\s*$/ || die "wrong end $line\n";
	} else {
	    die;
	}

	$lines[5] =~ /\d\d[ISD]/ && next; # less than 10bp mismatch string is allowed
	$align = $lines[5];
	$align =~ s/^(\d)S//;
	$t = $1;
	$align =~ s/(\d)S$//;
	$t += $1;
	$t > 7 && next;

	if($sam =~ /^0|256$/) {
	    $t = 0;
	    while($align =~ s/(\d+)[MSD]//g) {
		$t += $1;
	    }
	    while($align =~ s/(\d+)I//g) {
		$t -= $1;
	    }
	    $align && die $lines[5];
	    $lines[11] =~ /^AS\:i\:(\d+)$/ || die;
	    $t = "$1\t$lines[2]"."_$lines[3]"."_$t";

	    if($sam) {
		$locationRead{$read}->[@{$locationRead{$read}}] = $t;
	    } else {
		exists($locationRead{$read}) ? die : ($locationRead{$read}->[0] = $t);
	    }
	} elsif($sam =~ /^67|323|65|321$/) {
	    $lines[8] =~ /^\-(\d+)$/ && next;
	    $lines[8] < 20 && die "$read\t", $FILE."_$sam*$lines[8]";
	    $lines[11] =~ /^AS\:i\:(\d+)$/ || die;
	    $t = $1;
	    $lines[$#lines-1] =~ /^YS\:i\:(\d+)$/ || die;
	    $t .= "_$1\t$lines[2]"."_$lines[3]"."_$lines[8]";

	    if($sam =~ /^67|65$/) {
		exists($locationRead{$read}) ? die : ($locationRead{$read}->[0] = $t);
	    } elsif($sam =~ /^323|321$/) {
		$locationRead{$read}->[@{$locationRead{$read}}] = $t;
	    } else {
		die;
	    }
	} else {
	    die;
	}
    }

    return %locationRead;
}

sub locnomorethan1indelxS($$$$) {
    my ($query, $reference, $threshold, $lORr) = @_;
    $lORr || ($lORr = 'l');
    my ($lenq, $lenr) = (length($query), length($reference));
    my ($lenD, $lenI) = ($lenr - 1, $lenr + 1);
    $lenq < $lenr && die "$query $reference $lenq $lenr";
    my ($i, $l) = ();

    if($lORr eq 'r') {
	for($i = $lenq - $lenr; $i >= 0; $i--) {
	    diffnomorethan(substr($query, $i, $lenr), $reference, $threshold) && (return $i + $lenr);
	}
	for($i = $lenq - $lenD; $i >= 0; $i--) {
	    differby1D(substr($query, $i, $lenD), $reference) && (return $i + $lenD);
	}
	for($i = $lenq - $lenI; $i >= 0; $i--) {
	    differby1D($reference, substr($query, $i, $lenI)) && (return $i + $lenI);
	}
    } elsif($lORr eq 'l') {
	for($i = 0; $i <= $lenq - $lenr; $i++) {
	    diffnomorethan(substr($query, $i, $lenr), $reference, $threshold) && (return $i + 1);
	}
	for($i = 0; $i <= $lenq - $lenD; $i++) {
	    differby1D(substr($query, $i, $lenD), $reference) && (return $i + 1);
	}
	for($i = 0; $i <= $lenq - $lenI; $i++) {
	    differby1D($reference, substr($query, $i, $lenI)) && (return $i + 1);
	}
    } else {
	die;
    }

    return 0;
}

sub matchScore($$$) {
    my ($str, $sam5, $len) = @_;
    my $matchbp = 0;

    $str =~ /^\^/ && die $str;
    $str =~ /^\^[ACGT]*$/ && die $str;
    $str =~ /^\d/ || die $str;
    $str =~ /\d$/ || die $str;
    $str =~ s/[ACGT][0..4]$//;
    $sam5 =~ /^[\dMSID]+$/ || die;

    while($str =~ /(\d+)/g) {
	$matchbp += $1;
    }

    while($str =~ /[\^ACGT]/g) {
	$matchbp--;
    }

    while($sam5 =~ /I/g) {
	$matchbp -= 2;
    }

    return $matchbp;
}

sub mismatch1($%) {
    my ($endLen, %barcode) = @_;
    my ($barcode, $barcodeLeft, $barcodeRight, $i, $j, $num, $t, @barcode, %convert, %barcodeEnd, %t) = ();

    if($endLen > 0) {
	while(($barcode, $num) = each(%barcode)) {
	    $barcode =~ /^([ACGTN]{$endLen})([ACGTN]+)$/i || die "$endLen*$barcode $num";
	    exists($barcodeEnd{$1}->{$2}) ? die : ($barcodeEnd{$1}->{$2} = $num);
	}

	while(($barcodeLeft, $t) = each(%barcodeEnd)) {
	    %t = &mismatch1(0, %{$t});

	    while(($barcodeRight, $barcode) = each(%t)) {
		exists($convert{$barcodeLeft.$barcodeRight}) ? die : ($convert{$barcodeLeft.$barcodeRight} = $barcodeLeft.$barcode);
		exists($barcode{$barcodeLeft.$barcodeRight}) || die;
		exists($barcode{$barcodeLeft.$barcode}) ? ($barcode{$barcodeLeft.$barcode} += $barcode{$barcodeLeft.$barcodeRight}) : die;
		delete($barcode{$barcodeLeft.$barcodeRight});
	    }
	}
	%barcodeEnd = ();

	while(($barcode, $num) = each(%barcode)) {
	    $barcode =~ /^([ACGTN]+)([ACGTN]{$endLen})$/i || die;
	    exists($barcodeEnd{$2}->{$1}) ? die : ($barcodeEnd{$2}->{$1} = $num);
	}

	while(($barcodeRight, $t) = each(%barcodeEnd)) {
	    %t = &mismatch1(0, %{$t});

	    while(($barcodeLeft, $barcode) = each(%t)) {
		exists($convert{$barcodeLeft.$barcodeRight}) ? die : ($convert{$barcodeLeft.$barcodeRight} = $barcode.$barcodeRight);
		exists($barcode{$barcodeLeft.$barcodeRight}) || die "$barcodeLeft\t$barcodeRight";
		exists($barcode{$barcode.$barcodeRight}) || die;
	    }
	}
 
	while(($t, $barcode) = each(%convert)) {
	    if(length($t) == length($barcode)) {
		&nomorethan1S($t, $barcode) || delete($convert{$t});
	    } elsif(length($t) - length($barcode) == 1) {
		&differby1D($barcode, $t) || die;
            } elsif(length($barcode) - length($t) == 1) {
                &differby1D($t, $barcode) || die;
	    }
	}
    } elsif($endLen) {
	die;
    } else {
	@barcode = sort {$barcode{$b} <=> $barcode{$a}} keys %barcode;
	%convert = ();

	for($i = 0; $i < @barcode; $i++) {
	    $barcode{$barcode[$i]} > 1 || last;

	    for($j = $#barcode; $j > $i; $j--) {

		if(&nomorethan1indelxS($barcode[$i], $barcode[$j], 1)) {
		    exists($convert{$barcode[$j]}) ? die "$i $j $barcode[$i] $barcode[$j] $convert{$barcode[$j]}": ($convert{$barcode[$j]} = $barcode[$i]);
		    @barcode = grep {$_ ne $barcode[$j]} @barcode;
 		}
            }
        }
    }

    return %convert;
}

sub nomorethan1indel($$) {
    my ($query, $reference) = @_;
    $query eq $reference && (return 1);
    length($query) != length($reference) && die;
    my $len = length($query);
    my($i, $k) = ();

	for($k = 1; $k < $len; $k++) {
	    if(substr($query, 0-$k) ne substr($reference, 0-$k)) {
		substr($query, 1, $len-$k) eq substr($reference, 0, $len-$k) && (return 1);
		last;
	    }
	}
	$k == $len && die "$query, $reference\t$k";
 
	for($k = 1; $k < $len; $k++) {
	    if(substr($query, 0, $k) ne substr($reference, 0, $k)) {
		substr($query, $k-1, $len-$k) eq substr($reference, $k) && (return 1);
		last;
	    }
	}
    $k == $len && die "$query, $reference\t$k";

	for($k = 1; $k < $len; $k++) {
	    if(substr($query, 0, $k) ne substr($reference, 0, $k)) {
		substr($query, $k) eq substr($reference, $k-1, $len-$k) && (return 1);
		last;
	    }
	}
    $k == $len && die "$query, $reference\t$k";

	for($k = 1; $k < $len; $k++) {
	    if(substr($query, 0-$k) ne substr($reference, 0-$k)) {
		substr($query, 0, $len-$k) eq substr($reference, 1, $len-$k) && (return -1);
		last;
	    }
	}
	$k == $len && die "$query, $reference\t$k";

    return 0;
}

sub nomorethan1indelxS($$$) {
    my ($str1, $str2, $threshold) = @_;
    my ($len1, $len2) = (length($str1), length($str2));

    if(abs($len1 - $len2) > 1) {
	return 0;
    } elsif($len1 > $len2) {
	return differby1D($str2, $str1);
    } elsif($len1 < $len2) {
	return differby1D($str1, $str2);
    } elsif($threshold) {
	diffnomorethan($str1, $str2, $threshold) && (return 1);
	return nomorethan1indel($str1, $str2);
    } elsif($str1 eq $str2) {
	return 1;
    } else {
	return 0;
    }
}

sub nomorethan1internalSD($$) {
    my ($str1, $str2) = @_;
    length($str1) != length($str2) && die;
    my $len = length($str1);
    my $ds = 0;
    my($k) = ();

    nomorethan1S($str1, $str2) && (return 1);
    nomorethan1S(substr($str1,0,$len-1), substr($str2,1)) && (return 1);#check whether $str2 shifts backward by 1 but its remain sequence perfectly matches to $str1

    for($k = 1; $k < $len; $k++) { #check whether there is one deletion in $str2
	if(substr($str1, 0, $k) ne substr($str2, 0, $k)) {
	    substr($str1, $k) eq substr($str2, $k-1, $len-$k) && (return 1);
	    last;
	}
    }
    $k == $len && die;

    return 0;
}

sub nomorethan1S($$) {
    my ($str1, $str2) = @_;
    length($str1) != length($str2) && die;
    my $len = length($str1);
    my($k) = ();

    for($k = 1; $k < $len; $k++) { #check whether there is no more than one mismatch
	if(substr($str1, 0, $k) ne substr($str2, 0, $k)) {
	    substr($str1, $k) eq substr($str2, $k) && (return 1);
	    last;
	}
    }
    $k == $len && (return 1);

    return 0;
}

sub nomorethanxS($$$) {
    my ($str1, $str2, $threshold) = @_;
    length($str1) != length($str2) && die;
    $str1 eq $str2 && (return 1);
    $threshold eq '0' && (return 0);
    my($len, $k) = (length($str1)-$threshold+1,);

    for($k = 1; $k < $len; $k++) {
	if(substr($str1, 0, $k) ne substr($str2, 0, $k)) {
#	    if($threshold) {
		return nomorethanxS(substr($str1, $k), substr($str2, $k), $threshold-1);
#	    } else {
#		substr($str1, $k) eq substr($str2, $k) && (return 1);
#	    }
	    last;
	}
    }
    $k == $len && (return 1);

    return 0;
}

sub oneindelmismatch($$$$) {
    my ($str1, $str2, $threshold,$lORr) = @_;
    $lORr || ($lORr = 'l');
    my ($i,$tmp) = ();
    return nomorethan1indelxS($str1,$str2,$threshold) if(length($str1) == length($str2));
    return 0 if(abs(length($str1) - length($str2))>1);
    if(abs(length($str1) - length($str2)) == 1){
        length($str1)>length($str2) || (($str1,$str2)=($str2,$str1));
        #print "$str1    $str2\n";
        for($i=0;$i<length($str1);$i++){
            $tmp = $str1;
            #print $tmp,"\n";
            substr($tmp,$i,1) = undef;
            #print $tmp,"\n";
            return 1 if($tmp eq $str2);
        }
        die if($tmp ne substr($str1,0,length($str1)-1));
        return 0;
    }else{
        die "$str1,$str2\n";
    }
     
}

sub readFAfile($) {
    my($file) = shift;
    my($read) = ('');

    open(I, $file) || die;
    while(<I>) {
	s/[\d\s]//g || die;
	/^[ACGT]+$/i || die $_;
	$read .= uc($_);
    }
    return $read;
}

sub revcomDNA($) {
	my $dna = shift;
	my $revcomp = reverse($dna);
	$revcomp =~ tr/ACTGNRYMKSWHDBVactgnrymkswhdbv/TGACNYRKMSWDHVBtgacnyrkmswdhvb/;
	return $revcomp;
}

sub revcomFASTQ($$) {
    my ($fileIN, $fileOUT) = @_;
    my @lines;

    open(FQ, $fileIN) || die;
    open(TMP, ">$fileOUT") || die;

    @lines = ();
    while(<FQ>) {
        push @lines,$_;

        if(@lines == 4){
                chomp($lines[1]);
                chomp($lines[3]);
                $lines[1] = revcomDNA($lines[1])."\n";
                $lines[3] = reverse($lines[3])."\n";
                print TMP @lines;
            @lines = ();
        } elsif(@lines > 4) {
            die scalar @lines;
        }
    }
    close(TMP);
    close(FQ);
}

sub rmHybrid($%) {
    my ($poolFile, %tl) = @_;
    my($count, $barcodeL, $barcodeR, $br, $maxL, $maxR, $num) = (0);
    my(%doubt, %t, %tr) = ();


    while(($barcodeL, $br) = each(%tl)) {
	while(($barcodeR, $num) = each(%{$br})) {
	    exists($tr{$barcodeR}->{$barcodeL}) ? die : ($tr{$barcodeR}->{$barcodeL} = $num);
	}
    }

    while(($barcodeL, $br) = each(%tl)) {
	$maxL = max(values %{$br});

	while(($barcodeR, $num) = each(%{$br})) {
	    $maxR  = max(values %{$tr{$barcodeR}});

	    if($maxL > $num && $maxR > $num) {
		delete($tr{$barcodeR}->{$barcodeL});
		delete($tl{$barcodeL}->{$barcodeR});
		$count += $num;
		$doubt{$barcodeL.'::'.$barcodeR} = $num;
	    }
	}
    }
    print "$count reads are hybrids and discarded\n";

  #  open(I, $poolFile) || die $poolFile;
#    while(<I>) {
#	/^([ACGTN]+\:[ACGTN]+)/ || die "$_\n";
	#exists($doubt{$1}) && print "Warning: the hybrid in the post-pool exists as a barcode pair in the pre-pool: $_";
 #   }

    return %tl;
}

sub rmHybrid1(%) {
    my (%tl) = @_;
    my($count, $barcodeL, $barcodeR, $br, $maxL, $maxR, $num) = (0);
    my(%t, %tr) = ();


    while(($barcodeL, $br) = each(%tl)) {
	while(($barcodeR, $num) = each(%{$br})) {
	    exists($tr{$barcodeR}->{$barcodeL}) ? die : ($tr{$barcodeR}->{$barcodeL} = $num);
	}
    }

    while(($barcodeL, $br) = each(%tl)) {
	$maxL = max(values %{$br});

	while(($barcodeR, $num) = each(%{$br})) {
	    $maxR  = max(values %{$tr{$barcodeR}});

	    if($maxL > $num && $maxR > $num) {
		delete($tr{$barcodeR}->{$barcodeL});
		delete($tl{$barcodeL}->{$barcodeR});
		$count += $num;
	    }
	}
    }
    print "$count reads are hybrids and discarded\n";

    return %tl;
}

sub rmN($%) {
    my ($poolFile, %tl) = @_;
    my($count, $barcodeL, $barcodeR, $br, $i, $maxL, $maxR, $num, $test, $dead) = (0);
    my(@barcode, %doubt, %t, %tr) = ();

    while(($barcodeL, $br) = each(%tl)) {
	%t = %{$br};
	@barcode = sort {$t{$b} <=> $t{$a}} keys %t;
	$test = 0;
	$dead = 0;

	for($i = 1; $i < @barcode; $i++) {
	    if($t{$barcode[0]} > $t{$barcode[$i]}) {
		$count += $tl{$barcodeL}->{$barcode[$i]};
		$doubt{$barcodeL.'::'.$barcode[$i]} = $tl{$barcodeL}->{$barcode[$i]};

		if($tl{$barcodeL}->{$barcode[$i]} > 1 || $tl{$barcodeL}->{$barcode[0]} < 10) {
#		    print "$barcodeL\t$barcode[$i]\t", $tl{$barcodeL}->{$barcode[$i]}, "\n";
		    $test = 1;
		}
		$dead = 1 if($tl{$barcodeL}->{$barcode[$i]} > 1 && $tl{$barcodeL}->{$barcode[0]} < 10*$tl{$barcodeL}->{$barcode[$i]});
		delete($tl{$barcodeL}->{$barcode[$i]});
	    }
	}
	
#	$test && print "$barcodeL\t$barcode[0]\t", $tl{$barcodeL}->{$barcode[0]}, "	$dead\n\n";
	delete($tl{$barcodeL}) if($dead == 1);
    }
 #   print "$count reads are non-max barcode pairs in 1L:NR and discarded\n";

    while(($barcodeL, $br) = each(%tl)) {
	while(($barcodeR, $num) = each(%{$br})) {
	    exists($tr{$barcodeR}->{$barcodeL}) ? die : ($tr{$barcodeR}->{$barcodeL} = $num);
	}
    }

    $count = 0;
    while(($barcodeL, $br) = each(%tr)) {
	%t = %{$br};
	@barcode = sort {$t{$b} <=> $t{$a}} keys %t;
	$test = 0;
	$dead = 0;
	for($i = 1; $i < @barcode; $i++) {
	    if($t{$barcode[0]} > $t{$barcode[$i]}) {
		$count += $tr{$barcodeL}->{$barcode[$i]};
		$doubt{$barcode[$i].'::'.$barcodeL} = $tr{$barcodeL}->{$barcode[$i]};

		if($tr{$barcodeL}->{$barcode[$i]} > 1 || $tr{$barcodeL}->{$barcode[0]} < 10) {
#		    print "non-max read-pair $barcode[$i]\t$barcodeL\t", $tr{$barcodeL}->{$barcode[$i]}, "\n";
		    $test = 1;
		}
		$dead = 1 if($tr{$barcodeL}->{$barcode[$i]} > 1 && $tr{$barcodeL}->{$barcode[0]}<10*$tr{$barcodeL}->{$barcode[$i]});
		delete($tr{$barcodeL}->{$barcode[$i]});
		exists($tl{$barcode[$i]}) || die;

#		if(scalar(keys %{$tl{$barcode[$i]}}) > 1) {
#		    delete($tl{$barcode[$i]}->{$barcodeL});
#		} elsif(exists($tl{$barcode[$i]}->{$barcodeL})) {
#		    delete($tl{$barcode[$i]});
#		} else {
#		    die;
#		}
	    }
	}

#	$test && print "max read-pair $barcode[0]\t$barcodeL\t", $tr{$barcodeL}->{$barcode[0]}, "	$dead\n\n";
	delete($tr{$barcodeL}) if($dead==1);
#	if(scalar(keys %{$tl{$barcode[0]}}) > 1) {
#	    delete($tl{$barcode[0]}->{$barcodeL});
#	} elsif(exists($tl{$barcode[0]}->{$barcodeL})) {
#	    delete($tl{$barcode[0]});
#	} else {
#	    die;
#	}
    }

    %tl = ();
    while(($barcodeL, $br) = each(%tr)) {
	while(($barcodeR, $num) = each(%{$br})) {
	    exists($tl{$barcodeR}->{$barcodeL}) ? die : ($tl{$barcodeR}->{$barcodeL} = $num);
	}
    }
    
#    print "$count reads are non-max barcode pairs in extra NL:1R and discarded\n";

 #   open(I, $poolFile) || die;
  #  while(<I>) {
#	/^([ACGTN]+\:[ACGTN]+)/ || die;
#	exists($doubt{$1}) && print "Warning: this hybrid-like in the post-pool exists as a barcode-pairin the pre-pool: $_";
 #   }

    return %tl;
}

sub rmLend($$$$$) {
    my ($query, $reference, $location, $len, $mismatch) = @_;
    $mismatch || ($mismatch = 1);
    $location =~ /^\d+$/ || ($location = 0);
    $len < length($reference) && return 0;
    my $locationR = $location + $len - length($reference);
    my $m;

    if($query =~ /^([ACGTN]{$location,$locationR}$reference)/) {
	return length($1);
    } else {
	length($query) - $location < length($reference) && return 0;
	$m = &locnomorethan1indelxS(substr($query, $location, $len), $reference, $mismatch, 'r');
	$m && return $m + $location;
    }

    return 0;
}

sub rmLendBridge($$$$$$) {
    my ($query, $reference, $location, $len, $mismatch, $bridgeLen) = @_;
    $bridgeLen > 10 && die "Bridge is too long to make sense";
    my $t = rmLend($query, $reference, $location, $len, $mismatch);
    $t && return $t;

    if($bridgeLen > 1) {
	$t = rmLendBridge($query, substr($reference,1), $location, $len, $mismatch, $bridgeLen - 1);
	$t ? return($t+1) : return 0;
    } elsif($bridgeLen < 1) {
	return 0;
    } else {
	chop($reference);
	$t = rmLend($query, $reference, $location, $len, $mismatch);
	$t ? return ($t+1) : return 0;
    }
}

sub rmRend($$$$$) {
    my ($query, $reference, $location, $len, $mismatch) = @_;
    my $refLen = length($reference);
    $mismatch || ($mismatch = 1);
    $len < $refLen && return 0;
    $location < $refLen  && return 0;
    my ($ll, $locationR) = ($location - $refLen, $location - $len);
    $locationR < 0 && ($locationR = 0);
    my $m;
    $location =~ /^\d+$/ || die;

    if($query =~ /($reference[ACGTN]{$locationR,$ll})$/) {
	return length($1);
    } else {
	length($query) < $location && return 0;
	$m = &locnomorethan1indelxS(substr($query, 0-$location, $len), $reference, $mismatch, 'l');
	$m && return $location - $m + 1;
    }

    return 0;
}

sub rmRendBridge($$$$$$) {
    my ($query, $reference, $location, $len, $mismatch, $bridgeLen) = @_;
    $bridgeLen > 10 && die "Bridge is too long to make sense";
    my $t = rmRend($query, $reference, $location, $len, $mismatch);
    $t && return $t;

    if($bridgeLen > 1) {
	$t = rmRendBridge($query, substr($reference,1), $location, $len, $mismatch, $bridgeLen - 1);
	$t ? return($t+1) : return 0;
    } elsif($bridgeLen < 1) {
	return 0;
    } else {
	$t = rmRend($query, substr($reference,1), $location, $len, $mismatch);
	$t ? return ($t+1) : return 0;
    }
}

sub rmSite($$$$$$) {
    my($read, $SITE, $BRIDGELEN, $BARCODELEN, $BARCODELENMAX, $LorR) = @_;
    $SITE || return ();
    my $SITELEN = length($SITE);

    if($LorR eq 'l') {
	if($read =~ s/^([ACGTN]{1,$BARCODELENMAX}).{$BRIDGELEN}($SITE)//) {
#	    length($1) > 20 && die "$1\t$BRIDGELEN\t$SITE";
	    return($read, $1);
	} elsif($read =~ s/^([ACGTN]{$BARCODELEN}).{$BRIDGELEN}([ACGTN]{$SITELEN})//) {
	    if(&nomorethan1S($2, $SITE)) {
		return($read, $1);
	    } else {
		return();
	    }
	} else {
	    return();
	}
    } elsif($LorR eq 'r') {
	if($read =~ s/($SITE).{$BRIDGELEN}([ACGTN]{1,$BARCODELENMAX})$//) {
	    return($read, $2);
	} elsif($read =~ s/([ACGTN]{$SITELEN}).{$BRIDGELEN}([ACGTN]{$BARCODELEN})$//) {
	    if(&nomorethan1S($1, $SITE)) {
		return($read, $2);
	    } else {
		return();
	    }
	} else {
	    return();
	}
    } else {
	die;
    }
}

sub SQmean($) {
	my $seq = shift;
	my $seqLen = length($seq); 
	my ($sum) = (0);
	$seqLen || die "empty input";

	while($seq =~ s/(.)//) {
	    $sum += ord($1);
	}

	return $sum/$seqLen - 33;
}

sub SQtrimEnd3($$$) {
    my ($seq, $trimLen, $threshold) = @_;
    $trimLen > 0 || die;
    $threshold || return $seq;
    $threshold < 0 && die; #must be wrong parameter
    my ($seqLen, $i) = (length($seq), );
    $seqLen < $trimLen && return '';

    for($i = $seqLen - 1; $i >= $seqLen - $trimLen; $i--) {
	ord(substr($seq, $i, 1)) < $threshold+33 && return &SQtrimEnd3(substr($seq, 0, $i), $trimLen, $threshold);
    }

    return $seq;
}

sub strDistance($$) {
    my ($str1, $str2) = @_;
    my @ar1 = split //, $str1;
    my @ar2 = split //, $str2;

    my @dist;
    $dist[$_][0] = $_ foreach (0 .. @ar1);
    $dist[0][$_] = $_ foreach (0 .. @ar2);
 
    foreach my $i (1 .. @ar1){
        foreach my $j (1 .. @ar2){
            my $cost = $ar1[$i - 1] eq $ar2[$j - 1] ? 0 : 1;
            $dist[$i][$j] = min(
                        $dist[$i - 1][$j] + 1, 
                        $dist[$i][$j - 1] + 1, 
                        $dist[$i - 1][$j - 1] + $cost );
        }
    }
 
    return $dist[@ar1][@ar2];
}

return 1;
END{};








