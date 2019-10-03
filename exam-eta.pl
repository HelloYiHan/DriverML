#!/usr/bin/perl
use strict;
use autodie;
use Getopt::Long;
our$file_number;
our$date;
my$options=GetOptions(
'f|file=s' => \$file_number,
'd|date=s' => \$date,
);
open ETA,'<',"${date}_subclass-$file_number-eta.tmp";
#open ETA,'<',"subclass-$file_number.tmpeta.tmp";
while(<ETA>){
	chomp;
	if($_<=0){
		#exit(1);
		print(1);
	}
}
close ETA;
#exit(0)
print(0);
