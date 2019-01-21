#!/usr/bin/perl
use autodie;
use strict;
use Getopt::Long;
our%chara;
our$input_file;
our$date;
our$out_name;
my$optionOK=GetOptions(
	'i|input_file=s' =>\$input_file,
        'd|date=s' =>\$date,
);
$out_name=$date."_gene-characteristic.tmp";

open CHAR,'<',"gene-characteristic-unix.txt";

open PLUS,'>',"$out_name";
my$header=<CHAR>;
while(<CHAR>){
	chomp;
	my@qwe=split(/\t/);
	$chara{$qwe[0]}=1;
	print PLUS "$_\n";
}
close CHAR;
our%notinchara;
open INPUTFILE,'<',"$input_file";
#$header=<INPUTFILE>;
while(<INPUTFILE>){
	chomp;
	my@qwe=split(/\t/);
	if(!$chara{$qwe[1]}){
		if($qwe[2]=~/X|Y/){
			$qwe[2]=23;
		}
		$notinchara{$qwe[1]}[0]=$qwe[2];
		$notinchara{$qwe[1]}[1]=$qwe[3];
	}
}
close INPUTFILE;
our%genome;
open GENOME,'<',"genome-characteristic.txt";
$header=<GENOME>;
while(<GENOME>){
	chomp;
	my@qwe=split(/\t/);
	$genome{$qwe[0]}{$qwe[1]}=$qwe[4];
}
close GENOME;
foreach my$qwe(keys %notinchara){
	if(!$genome{$notinchara{$qwe}[0]}){
		print "chr:$notinchara{$qwe}[0] which is in the input file is not in the genome-characteristic file!\n";
		next;
	}
	foreach my$asd(keys %{$genome{$notinchara{$qwe}[0]}}){
		if($notinchara{$qwe}[1] >= $asd && $notinchara{$qwe}[1] <= $asd+100000){
			$notinchara{$qwe}[2]=$genome{$notinchara{$qwe}[0]}{$asd};
		}
	}
}
foreach my$qwe(keys%notinchara){
        if($notinchara{$qwe}[2]){
		print PLUS "$qwe\t$notinchara{$qwe}[0]\t1\t1\t$notinchara{$qwe}[2]\tNaN\t1\tNaN\tNaN\t1\t1\t1\t1\t1\t1\t1\n";
	}
}
close PLUS;
