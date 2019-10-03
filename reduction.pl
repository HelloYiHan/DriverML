#!/usr/bin/perl -w

use strict;
use Carp;
use Getopt::Long;
use English;
use Pod::Usage;
use File::Spec;
use File::Path;
use File::Copy;
use File::Basename;
use Cwd;

my ($input_mt,$exp_mt,$ref_genome);
my ($out_dir,$simul_t, $prefix, $ratio, $ave_prop_CDS, $each_prop_CDS,$pathway);
$simul_t      = 0;
$ratio        = 0.05;
$ave_prop_CDS = 1.0;
my ($help,$man,$version,$usage);
my $optionOK = GetOptions(
	'i|input_mt=s'	        => \$input_mt,
	'g|exp_mt=s'	        => \$exp_mt,
	'f|ref_genome=s'        => \$ref_genome,
);
pod2usage(-verbose=>2) if($man or $usage or $help or $version);
if(!$input_mt) {
	pod2usage(1);
	croak "You need specify input mutation data file";
}
$input_mt = File::Spec->rel2abs($input_mt);
if(!$exp_mt) {
	pod2usage(1);
	croak "You need specify input the predfined gene mutation table";
}
$exp_mt = File::Spec->rel2abs($exp_mt);

if(!$ref_genome) {
	pod2usage(1);
	croak "You need to specify the reference genome";
}
$ref_genome = File::Spec->rel2abs($ref_genome);

croak "The file you sepcified does not exist" unless (-e $input_mt && -e $exp_mt && -e $ref_genome);

croak "The file:$input_mt does not exist" unless (-e $input_mt);

croak "The file:$exp_mt does not exist" unless (-e $exp_mt);

croak "The file:$ref_genome does not exist" unless (-e $ref_genome);

if($each_prop_CDS){
	$each_prop_CDS = File::Spec->rel2abs($each_prop_CDS);
	croak "The file you sepcified does not exist" unless (-e $each_prop_CDS);
}

my $input_base = fileparse($input_mt);

if(!$prefix) {
	$prefix = $input_base;
}

# check for the presence of R
if( !`which R 2> err.log` ){
	unlink("err.log");
	die "R is not found, which is required\n";
}
unlink("err.log");


##--------------------------------------------------------------------------------------------------------
my %mtype_mbase;
%mtype_mbase = (
	AG  => "AT_GC",
	TC  => "AT_GC",
	AC  => "AT_CG",
	TG  => "AT_CG",
	AT  => "AT_TA",
	TA  => "AT_TA",
	CT  => "CG_TA",
	GA  => "CG_TA",
	CA  => "CG_AT",
	GT  => "CG_AT",
	CG  => "CG_GC",
	GC  => "CG_GC"
);

my @mbase = ("AT_GC", "AT_CG", "AT_TA", "CG_TA", "CG_AT", "CG_GC");
my @pointmt = ("silent", "missense", "nonsense", "splicing");

my (%sample_mtype_pointmt);

my (%chr_seq, @chr_name, $chr);

my (@query, @query_match, $qn, $qc, $qm);
my (@db_gene, $dgn, @db_gene_match, $dgm);
my (%gene, @gene_id, $gn,  @strand,  %sample, @sample_id, $sn);
my (%FSindel,  %nFSindel);
my (@mutation_obs);
my ($flank1bp, $ref_var, $mutation);
my (@pathway, $pid, $pname, $pgn);
my (@exp_table, @obs_table, $etn, $otn, $gmt, $smt);

my ($logfile, $gene_exp_file, $gene_obs_file, $pathway_exp_file, $pathway_obs_file);
my ($out_gene_summary, $out_gene_detail, $out_pathway_summary, $out_pathway_detail);


$gene_exp_file    = $input_mt . "_exp.tmp";
$gene_obs_file    = $input_mt . "_obs.tmp";


my ($row, @rowarray, $flag, $sum);
my ($i, $j, $k, $h, $m, $s, $pn, $pt);


$logfile = "reduction". "\.log.tmp";
open(LOG,">>$logfile") || die "Cannot creat the file $logfile: $!";

open(IN, "$input_mt") || die "Cannot open the file $input_mt: $!";
@query = (); $qn = 0; $flag=0;

while(<IN>){
	chomp;
	@rowarray = split(/\s+/);
	$rowarray[2] =~ s/^chr//i;
	$rowarray[6] =~ tr/A-Z/a-z/; $rowarray[6] =~ s/splic\w+/splicing/i; $rowarray[6] =~ s/nonsynonymous/missense/i; $rowarray[6] =~ s/synonymous/silent/i; $rowarray[6] =~ s/inframeshift/nFS_indel/i; $rowarray[6] =~ s/nonframeshift/nFS_indel/i;$rowarray[6] =~ s/frameshift/FS_indel/i; $rowarray[6] =~ s/fs_indel/FS_indel/i;$rowarray[6] =~ s/nfs_indel/nFS_indel/i;
	if( (($rowarray[6] =~ /silent/ || $rowarray[6] =~ /missense/ || $rowarray[6] =~ /nonsense/ || $rowarray[6] =~ /splicing/) && $rowarray[4]=~ /[ACGT]/ && $rowarray[5]=~ /[ACGT]/) || $rowarray[6] =~ /FS_indel/){
		push @query,[@rowarray];
		$qn++;
	}
	else{
		if($flag==0){ print LOG "\nGenes do not have any of one defined mutations: silent, missense, nonsense, splicing, FS_indel and nFS_indel\:\n"; }
		print LOG "@rowarray\n";
		$flag++;
	}
}
close(IN) || die "Cannot close the file $input_mt: $!";
$qc = scalar(@rowarray);
@query = sort custom_c1 @query;

open(IN, "$exp_mt") || die "Cannot open the file $exp_mt: $!";
$row = <IN>;
@db_gene = (); $dgn = 0;
while(<IN>){
	chomp;
	@rowarray = split(/\s+/);
	push @db_gene,[@rowarray];
	$dgn++;
}
close(IN) || die "Cannot close the file $exp_mt: $!";
@db_gene = sort custom_c0 @db_gene;


@query_match = (); $qm = 0; %sample = (); %gene = ();
$flag = 0; $j = 0;
for($i=0;$i<$qn;$i++){
	$k = 0;
	for(;$j<$dgn;$j++){
		if($query[$i][1] eq $db_gene[$j][0]){
			@rowarray = ();
			for($h=0;$h<$qc;$h++){ $rowarray[$h] = $query[$i][$h]; }
			push @query_match,[@rowarray];
			$qm++;
			$k++;

			$sample{$rowarray[0]}=0;
			$gene{$rowarray[1]}=0;

			last;
		}
		elsif($query[$i][1] lt $db_gene[$j][0]){ last; }
	}
	if($k==0){
		if($flag==0){ print LOG "\nQueried genes are not found in the database $exp_mt\:\n"; }
		print LOG "$query[$i][1]\n";
		$flag++;
	}
}

@gene_id = sort keys (%gene);
$gn = scalar(@gene_id);
@sample_id = sort keys(%sample);
$sn = scalar(@sample_id);

@db_gene_match = (); $dgm = 0;

open(TMP, ">$gene_exp_file") || die "Cannot open the file $gene_exp_file: $!";

@strand = ();
$j = 0;
for($i=0;$i<$gn;$i++){
	for(;$j<$dgn;$j++){
		if($gene_id[$i] eq $db_gene[$j][0]){
			push(@strand,$db_gene[$j][4]);
			@rowarray = ();
			push(@rowarray,$db_gene[$j][0]);
			for($h=5;$h<42;$h++){ push(@rowarray,$db_gene[$j][$h]); }
			push @db_gene_match,[@rowarray];
			$dgm++;

			print TMP "$db_gene[$j][0]";
			for($h=5;$h<42;$h++){ print TMP "\t$db_gene[$j][$h]"; }
			print TMP "\n";

			last;
		}
		elsif($gene_id[$i] lt $db_gene[$j][0]){ last; }
	}
}
close(TMP) || die "Cannot close the file $gene_exp_file: $!";

# read reference sequence

%chr_seq = (); @chr_name = ();

open(IN, "$ref_genome") || die "Cannot open the file $ref_genome: $!";
while(<IN>){
	if(/\>chr(\w+)\s+/ || /\>(\w+)\s+/) { $chr=$1; $chr_seq{$chr} = ""; next;}
	else { $chr_seq{$chr} .= $_; }
}
close(IN) || die "Cannot close the file $ref_genome: $!";

@chr_name = keys(%chr_seq);
foreach $chr (@chr_name) {
	$chr_seq{$chr} =~ s/\s+//g; $chr_seq{$chr} =~ s/\d+//g; $chr_seq{$chr} =~ tr/acgt/ACGT/;
}

# print observed mutation table

open(TMP, ">$gene_obs_file") || die "Cannot open the file $gene_obs_file: $!";

$j = 0;
#	d-henness: commented out mutation_obs b/c it does not seem to be used and inflates to over 70GB durring the course of a run.
#	I was able to reproduce the test results with this change.
#@mutation_obs = ();
for($i=0;$i<$gn;$i++){####$gn: number of matched gene####
	%sample_mtype_pointmt = (); %FSindel = (); %nFSindel = ();
	for($s=0;$s<$sn;$s++){####$sn: number of matched sample####
		$FSindel{$sample_id[$s]} = 0;  $nFSindel{$sample_id[$s]} = 0;
		$sample_mtype_pointmt{$sample_id[$s]}={
			 AT_GC => { silent => 0, missense => 0, nonsense => 0, splicing => 0 },
			 AT_CG => { silent => 0, missense => 0, nonsense => 0, splicing => 0 },
			 AT_TA => { silent => 0, missense => 0, nonsense => 0, splicing => 0 },
			 CG_TA => { silent => 0, missense => 0, nonsense => 0, splicing => 0, silent_CpG => 0, missense_CpG => 0, nonsense_CpG => 0, splicing_CpG => 0},
			 CG_AT => { silent => 0, missense => 0, nonsense => 0, splicing => 0, silent_CpG => 0, missense_CpG => 0, nonsense_CpG => 0, splicing_CpG => 0},
			 CG_GC => { silent => 0, missense => 0, nonsense => 0, splicing => 0, silent_CpG => 0, missense_CpG => 0, nonsense_CpG => 0, splicing_CpG => 0}
		}
	}
	for(;$j<$qm;$j++){
		if($gene_id[$i] eq $query_match[$j][1]){
			 if($query_match[$j][6] =~ /^FS_indel/){ $FSindel{$query_match[$j][0]}++; }
			 elsif($query_match[$j][6] =~ /nFS_indel/){ $nFSindel{$query_match[$j][0]}++; }
			 else {

				 $flank1bp = substr($chr_seq{$query_match[$j][2]}, $query_match[$j][3]-2, 3);
				 if(length($query_match[$j][4])==1){
				 	$ref_var = $query_match[$j][4] . $query_match[$j][5];
				 }
			 	 else{
				 	$ref_var = substr($query_match[$j][4],0,1) . substr($query_match[$j][5],0,1);
				 }
				 if($strand[$i] eq "-"){
					 $ref_var =~ tr/ACGT/TGCA/;
					 $flank1bp = reverse($flank1bp); $flank1bp =~tr/ACGT/TGCA/;
				 }
				 $mutation = $query_match[$j][6];
				 if($flank1bp =~/CG/){ $mutation .= "_CpG"; }

				 $sample_mtype_pointmt{$query_match[$j][0]}{$mtype_mbase{$ref_var}}{$mutation}++;
			 }
		}
		if($gene_id[$i] lt $query_match[$j][1]){ last; }
	}
	@rowarray = ();
	push (@rowarray,$gene_id[$i]);
	print TMP "$gene_id[$i]";
	for($s=0;$s<$sn;$s++){
		for($k=0;$k<4;$k++){
			for($h=0;$h<6;$h++){
				$pn = $sample_mtype_pointmt{$sample_id[$s]}{$mbase[$h]}{$pointmt[$k]};
				print TMP "\t$pn";
				push (@rowarray,$pn);
			}
			$pt = $pointmt[$k] . "_CpG";
			for($h=3;$h<6;$h++){
				$pn = $sample_mtype_pointmt{$sample_id[$s]}{$mbase[$h]}{$pt};
				print TMP "\t$pn";
				push (@rowarray,$pn);
			}
		}
		print TMP "\t$FSindel{$sample_id[$s]}\t$nFSindel{$sample_id[$s]}";
		push (@rowarray,$FSindel{$sample_id[$s]});
		push (@rowarray,$nFSindel{$sample_id[$s]});

	}
	print TMP "\n";
  #	push @mutation_obs,[@rowarray];

}
close(TMP) || die "Cannot close the file $gene_obs_file: $!";

sub custom_c0 {
	$a->[0] cmp $b->[0];
}

sub custom_c1 {
	$a->[1] cmp $b->[1];
}

sub generate_random_string
{
	my ($length_of_randomstring)=shift;

	my @chars=('a'..'z','A'..'Z','0'..'9','_');
	my $random_string;
	foreach (1..$length_of_randomstring)
	{
		$random_string.=$chars[rand @chars];
	}
	return $random_string;
}
