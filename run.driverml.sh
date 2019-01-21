#!/bin/bash

help_file="
\n
\n
\tThis is a free application used for identifying cancer driver genes.\n
\n
\tUSAGE:path/run.driverml.sh -w <path to the program> -i <input file> -f <reference genome> -p <prior file> [options]\n
\n
\tRequired Arguments:\n
\n
\t-w/--pathway\t\t\tThe path to the program.\n
\n
\t-i/--input\t\t\tThe input mutation data file which is in Mutation Annotation Format (MAF). The MAF specifications could\n\t\t\t\t\tbe seen on the manual of this application or NCI Wiki.\n
\n
\t-f/--reference_file\t\tThe human reference genome file for input data.\n
\n
\t-r/--reference_training\t\tThe human reference genome file for training data.\n
\n
\tOptions:\n
\n
\t-g/--mutation_table\t\tThe predefined gene mutation table. Default: GRCh38.\n
\n
\t-y/--tumor_type\t\t\tTraining mutation data. Default: Pan-caner.\n
\n
\t-m/--multicore <int>\t\tSet the number of parallel instances to be run concurrently. Default: 20.\n
\n
\t-o/--output <chr>\t\tThe name of the output file. Default: summary.\n
\n
\t-t/--simulation_time <int>\tThe times in monte carlo simulation. Default: 10000.\n
\n
\t-n/--interpolation_number <int>\tThe max interpolation number used to estimate gene characteristics which are unknown. Default: 100.\n
\n
\t-c/--cluster <int>\t\tThe max cluster number used to estimate background mutation rate. Default: 50.\n
\n
\t-p/--prior\t\t\tSet the prior information. Default: Non-TCGA genes from DriverDB and IntOGen databases.\n
\n
\t-d/--indelratio\t\t\tThe ratio of point to indel mutation in background. Default: 0.05.\n
\n
\tOther:\n
\n
\t-h/--help\t\t\tDisplay the help file.\n
\n
\t-v/--version\t\t\tDisplay version information.\n
"

# initialization of variables
path=
input_file=
mutation_table_file=hg38_refGene.exp
reference_file=
reference_training=
prior=
multi_core=4
monte_carlo_times=10000
indel_ratio=0.05
interpolation_number=100
cluster_number=50
eps=1
tumortype=
output_prefix=summary
date=$(date +%m_%d_%H_%M_%S_%N)
date_pre=pre_${date}
# read the arguments
TEMP=`getopt -o w:i:g:f:p:y:m:o:t:r:d:n:c:e:hv --long pathway:,input:,mutation_table:,reference_genome:,prior:,tumor_type:,multicore:,output:,simulation_time:,reference_training:,indelratio:,interpolation:,cluster:,eps:,help,version -- "$@"`
eval set -- "$TEMP"

# extract arguments into variables.
while true ; do
    case "$1" in
	-w|--pathway)
		path=$2 ; shift 2 ;;
        -i|--input)
		input_file=$2 ; shift 2 ;;
        -g|--mutation_table)
		mutation_table_file=$2 ; shift 2 ;;
        -f|--reference_genome)
		reference_file=$2 ; shift 2 ;;
	-p|--prior)
		prior=$2 ; shift 2 ;;
	-y|--tumor_type)
		tumortype=$2 ; shift 2 ;;
	-m|--mluticore)
		multi_core=$2 ; shift 2 ;;
	-o|--output)
		output_prefix=$2 ; shift 2 ;;
        -t|--simulation_time)
		monte_carlo_times=$2 ; shift 2 ;;
        -r|--reference_training)
		reference_training=$2 ; shift 2 ;; 
	-d|--indelratio)
		indel_ratio=$2 ; shift 2 ;;
	-n|--interpolation)
		interpolation_number=$2 ; shift 2 ;;
	-c|--cluster)
		cluster_number=$2 ; shift 2 ;;
        -e|--eps)
                eps=$2 ; shift 2 ;;
        -h|--help)
		echo -e $help_file; exit 1 ;;
        -v|--version)
		echo '1.0'; exit 0 ;;
	--) shift ; break ;;
        *) echo "Internal error!" ; exit 1 ;;
    esac
done

# warning for required parameters
if [ -z $path ]; then
	echo -e $help_file
	exit 0
elif [ -z $input_file ]; then
	echo -e $help_file
	exit 0
elif [ -z $reference_training ]; then
	echo -e $help_file
        exit 0
elif [ -z $reference_file ]; then
        echo -e $help_file
        exit 0
fi

if [ -z $prior ]; then
	prior=${path}/prior/PAN
else
	prior_t=${path}/$prior
	prior=$prioy_t
fi

if [ -z $tumortype ];then
	training=${path}/training/PAN
else
	training=${path}/$tumortype
fi

#pre-processing
###################################
Hugo_Symbol=$(awk -F "\t" 'NR==1{gsub(/\t/,"\n");print $0}' $input_file | grep "Hugo_Symbol" -n | awk -F ":" '{print $1}')
Chromosome=$(awk -F "\t" 'NR==1{gsub(/\t/,"\n");print $0}' $input_file | grep "Chromosome" -n | awk -F ":" '{print $1}')
Start_Position=$(awk -F "\t" 'NR==1{gsub(/\t/,"\n");print $0}' $input_file | grep "Start_Position" -n | awk -F ":" '{print $1}')
Variant_Classification=$(awk -F "\t" 'NR==1{gsub(/\t/,"\n");print $0}' $input_file | grep "Variant_Classification" -n | awk -F ":" '{print $1}')
Reference_Allele=$(awk -F "\t" 'NR==1{gsub(/\t/,"\n");print $0}' $input_file | grep "Reference_Allele" -n | awk -F ":" '{print $1}')
Tumor_Seq_Allele2=$(awk -F "\t" 'NR==1{gsub(/\t/,"\n");print $0}' $input_file | grep "Tumor_Seq_Allele2" -n | awk -F ":" '{print $1}')
Tumor_Sample_Barcode=$(awk -F "\t" 'NR==1{gsub(/\t/,"\n");print $0}' $input_file | grep "Tumor_Sample_Barcode" -n | awk -F ":" '{print $1}')
Variant_Type=$(awk -F "\t" 'NR==1{gsub(/\t/,"\n");print $0}' $input_file | grep "Variant_Type" -n | awk -F ":" '{print $1}')

if [ -n "$Hugo_Symbol" ] && [ -n "$Chromosome" ] && [ -n "$Start_Position" ] && [ -n "$Variant_Classification" ] && [ -n "$Reference_Allele" ] && [ -n "$Tumor_Seq_Allele2" ] && [ -n "$Tumor_Sample_Barcode" ] && [ -n "$Variant_Type" ];then

awk -F "\t" -v Hugo="$Hugo_Symbol" -v Chromo="$Chromosome" -v Start="$Start_Position" -v Variant="$Variant_Classification" -v Reference="$Reference_Allele" -v Tumor_Seq="$Tumor_Seq_Allele2" -v Tumor_Sample="$Tumor_Sample_Barcode" -v Variant_T="$Variant_Type" '{print $Hugo"\t"$Chromo"\t"$Start"\t"$Variant"\t"$Reference"\t"$Tumor_Seq"\t"$Tumor_Sample"\t"$Variant_T}' $input_file > ${date}_input_1.tmp
else
echo "The format of mutation file is not acceptable. MAF format is required."
exit 0
fi

awk -F "\t" -v row=1 'NR>=row && $4~/Frame_Shift_Del/{print $7"\t"$1"\t"$2"\t"$3"\t"$5"\t"$6"\t""Frame_Shift_Del""\t"$8}NR>=row && $4~/Frame_Shift_Ins/{print $7"\t"$1"\t"$2"\t"$3"\t"$5"\t"$6"\t""Frame_Shift_Ins""\t"$8}NR>=row && $4~/In_Frame_Del/{print $7"\t"$1"\t"$2"\t"$3"\t"$5"\t"$6"\t""In_Frame_Del""\t"$8}NR>=row && $4~/In_Frame_Ins/{print $7"\t"$1"\t"$2"\t"$3"\t"$5"\t"$6"\t""In_Frame_Ins""\t"$8}NR>=row && $4~/Missense_Mutation/{print $7"\t"$1"\t"$2"\t"$3"\t"$5"\t"$6"\t""Missense_Mutation""\t"$8}NR>=row && $4~/Nonsense_Mutation/{print $7"\t"$1"\t"$2"\t"$3"\t"$5"\t"$6"\t""Nonsense_Mutation""\t"$8}NR>=row && $4~/Nonstop_Mutation/{print $7"\t"$1"\t"$2"\t"$3"\t"$5"\t"$6"\t""Nonstop_Mutation""\t"$8}NR>=row && $4~/Silent/{print $7"\t"$1"\t"$2"\t"$3"\t"$5"\t"$6"\t""Silent""\t"$8}NR>=row && $4~/Splice_Site/ && $8~/SNP|DNP|TNP/{print $7"\t"$1"\t"$2"\t"$3"\t"$5"\t"$6"\t""Splice_Site""\t"$8}NR>=row && $4~/Translation_Start_Site/ && $8~/SNP|DNP|TNP/{print $7"\t"$1"\t"$2"\t"$3"\t"$5"\t"$6"\t""Translation_Start_Site""\t"$8}' ${date}_input_1.tmp > ${date}_input_intermediate2_file.tmp

awk -F "\t" '{gsub(/chr/,"",$3);print $0}' ${date}_input_intermediate2_file.tmp > ${date}_input_intermediate_file.tmp

awk -F "\t" -v row=1 'NR>=row && $7~/Frame_Shift_Del|Frame_Shift_Ins/{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t""Fs_indel"}NR>=row && $7~/In_Frame_Del|In_Frame_Ins/{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t""nFs_indel"}NR>=row && $7~/Missense_Mutation/{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t""missense"}NR>=row && $7~/Nonsense_Mutation/{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t""nonsense"}NR>=row && $7~/Silent/{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t""silent"}NR>=row && $7~/Splice_Site/ && $8~/SNP|DNP|TNP/{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t""splicing"}' ${date}_input_intermediate_file.tmp > ${date}_input_intermediate_file_pre.tmp

awk -F "\t" '$3~/^[1-9]$|^[1-2][0-9]$|^X$|^Y$/ && $2!~/^ENSG/ && $2!~/^LOC/{print $0}' ${date}_input_intermediate_file_pre.tmp > ${date}_input_file.tmp

$path/add-genes-to-characteristic.pl -i ${date}_input_file.tmp -d $date

#Hugo_Symbol=$(awk -F "\t" 'NR==1{gsub(/\t/,"\n");print $0}' $training | grep "Hugo_Symbol" -n | awk -F ":" '{print $1}')
#Chromosome=$(awk -F "\t" 'NR==1{gsub(/\t/,"\n");print $0}' $training | grep "Chromosome" -n | awk -F ":" '{print $1}')
#Start_Position=$(awk -F "\t" 'NR==1{gsub(/\t/,"\n");print $0}' $training | grep "Start_Position" -n | awk -F ":" '{print $1}')
#Variant_Classification=$(awk -F "\t" 'NR==1{gsub(/\t/,"\n");print $0}' $training | grep "Variant_Classification" -n | awk -F ":" '{print $1}')
#Reference_Allele=$(awk -F "\t" 'NR==1{gsub(/\t/,"\n");print $0}' $training | grep "Reference_Allele" -n | awk -F ":" '{print $1}')
#Tumor_Seq_Allele2=$(awk -F "\t" 'NR==1{gsub(/\t/,"\n");print $0}' $training | grep "Tumor_Seq_Allele2" -n | awk -F ":" '{print $1}')
#Tumor_Sample_Barcode=$(awk -F "\t" 'NR==1{gsub(/\t/,"\n");print $0}' $training | grep "Tumor_Sample_Barcode" -n | awk -F ":" '{print $1}')
#Variant_Type=$(awk -F "\t" 'NR==1{gsub(/\t/,"\n");print $0}' $training | grep "Variant_Type" -n | awk -F ":" '{print $1}')
#
#if [ -n "$Hugo_Symbol" ] && [ -n "$Chromosome" ] && [ -n "$Start_Position" ] && [ -n "$Variant_Classification" ] && [ -n "$Reference_Allele" ] && [ -n "$Tumor_Seq_Allele2" ] && [ -n "$Tumor_Sample_Barcode" ] && [ -n "$Variant_Type" ];then
#
#awk -F "\t" -v Hugo="$Hugo_Symbol" -v Chromo="$Chromosome" -v Start="$Start_Position" -v Variant="$Variant_Classification" -v Reference="$Reference_Allele" -v Tumor_Seq="$Tumor_Seq_Allele2" -v Tumor_Sample="$Tumor_Sample_Barcode" -v Variant_T="$Variant_Type" '{print $Hugo"\t"$Chromo"\t"$Start"\t"$Variant"\t"$Reference"\t"$Tumor_Seq"\t"$Tumor_Sample"\t"$Variant_T}' $training > ${date_pre}_input_1.tmp
#else
#echo "The format of mutation file is not acceptable. MAF format is required."
#exit 0
#fi
#
#awk -F "\t" -v row=1 'NR>=row && $4~/Frame_Shift_Del/{print $7"\t"$1"\t"$2"\t"$3"\t"$5"\t"$6"\t""Frame_Shift_Del""\t"$8}NR>=row && $4~/Frame_Shift_Ins/{print $7"\t"$1"\t"$2"\t"$3"\t"$5"\t"$6"\t""Frame_Shift_Ins""\t"$8}NR>=row && $4~/In_Frame_Del/{print $7"\t"$1"\t"$2"\t"$3"\t"$5"\t"$6"\t""In_Frame_Del""\t"$8}NR>=row && $4~/In_Frame_Ins/{print $7"\t"$1"\t"$2"\t"$3"\t"$5"\t"$6"\t""In_Frame_Ins""\t"$8}NR>=row && $4~/Missense_Mutation/{print $7"\t"$1"\t"$2"\t"$3"\t"$5"\t"$6"\t""Missense_Mutation""\t"$8}NR>=row && $4~/Nonsense_Mutation/{print $7"\t"$1"\t"$2"\t"$3"\t"$5"\t"$6"\t""Nonsense_Mutation""\t"$8}NR>=row && $4~/Nonstop_Mutation/{print $7"\t"$1"\t"$2"\t"$3"\t"$5"\t"$6"\t""Nonstop_Mutation""\t"$8}NR>=row && $4~/Silent/{print $7"\t"$1"\t"$2"\t"$3"\t"$5"\t"$6"\t""Silent""\t"$8}NR>=row && $4~/Splice_Site/ && $8~/SNP|DNP|TNP/{print $7"\t"$1"\t"$2"\t"$3"\t"$5"\t"$6"\t""Splice_Site""\t"$8}NR>=row && $4~/Translation_Start_Site/ && $8~/SNP|DNP|TNP/{print $7"\t"$1"\t"$2"\t"$3"\t"$5"\t"$6"\t""Translation_Start_Site""\t"$8}' ${date_pre}_input_1.tmp > ${date_pre}_input_intermediate2_file.tmp
#
#awk -F "\t" '{gsub(/chr/,"",$3);print $0}' ${date_pre}_input_intermediate2_file.tmp > ${date_pre}_input_intermediate_file.tmp
#
#awk -F "\t" -v row=1 'NR>=row && $7~/Frame_Shift_Del|Frame_Shift_Ins|Frame_Shift_Indel/{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t""Fs_indel"}NR>=row && $7~/In_Frame_Del|In_Frame_Ins|In_Frame_Indel/{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t""nFs_indel"}NR>=row && $7~/Missense_Mutation/{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t""missense"}NR>=row && $7~/Nonsense_Mutation/{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t""nonsense"}NR>=row && $7~/Silent/{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t""silent"}NR>=row && $7~/Splice_Site/{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t""splicing"}' ${date_pre}_input_intermediate_file.tmp > ${date_pre}_input_intermediate_file_pre.tmp
#
#awk -F "\t" '$3~/^[1-9]$|^[1-2][0-9]$|^X$|^Y$/ && $2!~/^ENSG/ && $2!~/^LOC/{print $0}' ${date_pre}_input_intermediate_file_pre.tmp > ${date_pre}_input_file.tmp
#
#$path/add-genes-to-characteristic.pl -i ${date_pre}_input_file.tmp -d $date_pre

##################################
#awk -F "\t" -v row=1 'NR>=row && $7~/Frame_Shift_Del|Frame_Shift_Ins/{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t""Fs_indel"}NR>=row && $7~/In_Frame_Del|In_Frame_Ins/{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t""nFs_indel"}NR>=row && $7~/Missense_Mutation/{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t""missense"}NR>=row && $7~/Nonsense_Mutation/{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t""nonsense"}NR>=row && $7~/Silent/{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t""silent"}NR>=row && $7~/Splice_Site/ && $8~/SNP|DNP|TNP/{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t""splicing"}' $input_file > ${date}_input_intermediate_file.tmp

#awk -F "\t" '$3~/^[1-9]$|^[1-2][0-9]$|^X$|^Y$/ && $2!~/^ENSG/ && $2!~/^LOC/{print $0}' ${date}_input_intermediate_file.tmp > ${date}_input_file.tmp

#$path/add-genes-to-characteristic.pl -i ${date}_input_file.tmp -d $date

awk -F "\t" -v row=1 'NR>=row && $7~/Frame_Shift_Del|Frame_Shift_Ins|Frame_Shift_Indel/{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t""Fs_indel"}NR>=row && $7~/In_Frame_Del|In_Frame_Ins|In_Frame_Indel/{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t""nFs_indel"}NR>=row && $7~/Missense_Mutation/{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t""missense"}NR>=row && $7~/Nonsense_Mutation/{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t""nonsense"}NR>=row && $7~/Silent/{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t""silent"}NR>=row && $7~/Splice_Site/{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t""splicing"}' $training > ${date_pre}_input_intermediate_file.tmp

awk -F "\t" '$3~/^[1-9]$|^[1-2][0-9]$|^X$|^Y$/ && $2!~/^ENSG/ && $2!~/^LOC/{print $0}' ${date_pre}_input_intermediate_file.tmp > ${date_pre}_input_file.tmp

$path/add-genes-to-characteristic.pl -i ${date_pre}_input_file.tmp -d $date_pre

# find the best parameter
a=1
b=$cluster_number
while [ 1 -gt 0 ]
do
if [ $a -eq $b ];then
cluster_final=$a
break
elif [ $(($b-$a)) -gt 1 ];then
c=$((($a+$b)/2))
elif [ $b -eq $c ];then
cluster_final=$a
break
else
c=$b
fi
$path/cluster.r $interpolation_number $c ${date}
$path/pre-drgap.pl -n $c -i ${date}_input_file.tmp -d $date
for t in $(seq 1 $c)
do
$path/reduction.pl -i ${date}_subclass-$t.tmp -g ${path}/$mutation_table_file -f $reference_file
$path/write_eta.r ${date}_subclass-$t.tmp_exp.tmp ${date}_subclass-$t.tmp_obs.tmp
done
eta_zero=0
for q in $(seq 1 $c)
do
eta_zero=$(($eta_zero+$(perl $path/exam-eta.pl -f $q -d $date)))
done
if [ $eta_zero -gt 0 ]; then
b=$c
else
a=$c
fi
done
echo "Final cluster number is $cluster_final"

#manipulate data format
$path/cluster.r $interpolation_number $cluster_final ${date}
$path/pre-drgap.pl -n $cluster_final -i ${date}_input_file.tmp -d $date

g++ -fopenmp -o ${date}_monte_carlo.out.tmp $path/monte_carlo_sim.cpp

for p in $(seq 1 1)
do

$path/cluster.r $interpolation_number 1 ${date_pre}

$path/pre-drgap.pl -n 1 -i ${date_pre}_input_file.tmp -d $date_pre

$path/reduction.pl -i ${date_pre}_subclass-$p.tmp -g ${path}/hg19_refGene.exp -f $reference_training

done

#statistic test

for p in $(seq 1 $cluster_final)
do
$path/reduction.pl -i ${date}_subclass-$p.tmp -g ${path}/$mutation_table_file -f $reference_file

$path/sta.r $p $multi_core $monte_carlo_times 1 $indel_ratio $prior $date $eps

./${date}_monte_carlo.out.tmp $multi_core $monte_carlo_times ${date}_${p}_eta.tmp ${date}_${p}_N.tmp ${date}_${p}_lrt.tmp ${date}_${p}_para.tmp $date $p

done

# aggregate results

$path/assemble.r $cluster_final 1 $output_prefix $date

# remove temporary files
echo $date
rm ${date}*.tmp
rm ${date_pre}*.tmp
