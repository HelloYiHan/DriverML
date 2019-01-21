#!/usr/bin/perl
use strict;
use autodie;
use FileHandle;
use Getopt::Long;
our($input_file,$class_number);
our$date;
#our$fix_name;
my $optionOK = GetOptions(
	'i|input_file=s'  =>\$input_file,
	'n|class_number=s' => \$class_number,
        'd|date=s' => \$date,
);
#------------------------------------------------------------
our%gene_class;
open CLASS,'<',"${date}_gene-class.tmp";
my$header=<CLASS>;
our$line=0;
while(<CLASS>){
	$line+=1;
	chomp;
	if(/^\d+\s+(\w+\-?\.?\w*)\s+(\d+)/){
		$gene_class{$1}=$2;
	}
	else{
		print "invalid format of gene class in line $line\n";
	}
}
close CLASS;
my%fh;
foreach my$i(1..$class_number){
	open $fh{$i},">${date}_subclass-$i.tmp";
}
open INPUT,'<',"$input_file";
open TMP,'>',"${date}_genes-which-are-in-input-file-but-not-in-the-chara-file.tmp";
our$n=0;#the number of genes which are in input file but not in the characteristic file------#
while(<INPUT>){
	chomp;
	my@element=split(/\t/);#print "$element[1]\n";
	my$class=$gene_class{$element[1]};
	if(exists($gene_class{$element[1]})){
		$fh{$class}->print("$_\n");
	}
	else{
		print TMP "$element[0]\t$element[1]\t$element[2]\n";
		$n+=1;
	}
}
print("The number of genes which are in the input file but not in the characteristic file is $n\n");
close INPUT;
foreach my$i(1..$class_number){
        close $fh{$i};
}
close TMP;
