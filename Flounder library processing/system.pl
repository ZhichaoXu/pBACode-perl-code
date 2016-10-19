#! usr/bin/perl
use strict;
use warnings;

my $useage = "useage: perl system.pl filename.txt configure_pooling.txt\n";
my ($file_name, $configure) = ('','');
unless ($file_name=$ARGV[0] and $configure=$ARGV[1]){
   print $useage;
   exit;
}

open (FN, "< $file_name") || die "can not open $file_name\n";

while(<FN>){
	chomp $_;
#assign file name
system "perl barcodeRMvectorSitePost_side.pl $configure left $_";
system "perl barcodePostpool_side.pl $configure left $_";

}

close FN;

exit;
