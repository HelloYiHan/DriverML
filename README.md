# DriverML
DriverML integrates the Rao’s score test and supervised machine learning to identify cancer driver genes. Rigorous and unbiased benchmark analysis and comparisons of DriverML with 20 other existing tools in 31 independent datasets from The Cancer Genome Atlas (TCGA) show that DriverML is robust and powerful among various datasets and outperforms the other tools with a better balance of precision and recall. 

## Access
DriverML is free for non-commerical use only.

## Installation
unzip DriverML-master.zip  
cd DriverML-master  
tar jxvf training.tar.bz2  
chmod +x *.pl *.r *.sh

## Running
Requirements: Perl and R are required on user’s Linux environment.  

Usage: path/run.driverml.sh	-w <AbsolutePath/DriverML-master> –i <AbsolutePath/mutation.txt> -f <AbsolutePath/reference file for input data> -r <AbsolutePath/reference file for training data> [options]  

Required arguments:  

| Alias| Description |
|:---------------|:-----------------------------------------------|
| -w/--path	     |       The absolute path to the DriverML_v1.0.0. |
| -i/--input	   |     The list of tumor mutations to be analyzed. It should be put in the DriverML-Master/ directory.|  
|-f/--reference_genome|	The reference genome for input(-i) data.|
|-r/--reference_training	|The human reference genome file for training data. Deafuat: hg19 reference file|  

Options:  

| Alias| Description |
|:---------------|:-----------------------------------------------|
|-g/--mutation_table|	The predefined gene mutation table which could be found in the package of this application. It could be either hg19_refGene.exp or hg38_refGene.exp according to -i input file. Default: hg38_refGene.exp|
|-y/--tumor_type|	Training mutation data.Default: Pan-cancer|
|-m/--multicore	|Set the number of parallel instances to be run concurrently. Default: 4|
|-t/--simulation_time|	Set the number for Monte Carlo Simulation. Default: 2500|
|-o/--output|	Set the prefix of the output file. Default: summary.|
|-c/--cluster_number|	Set the upper limit of cluster number for computing BMR. Default: 1|
|-p/--prior|	Set the prior information. Default: Non-TCGA genes from DriverDB and IntOGen databases.|
|-n/--interpolation_number| Set the upper limit of interpolation number for making gene clusters. Default: 100|
|-d/--indel_ratio	|The ratio of point mutation to indel in the background. Default: 0.05|
|-h/--help|	Help information.|
|-v/--version|	Software version.|

## Notes
The mutation file needs to be MAF format. Eight columns named Hugo_Symbol, Chromosome, Start_Position, Variant_Classification, Variant_Type, Reference_Allele, Tumor_Seq_Allele2 and Tumor_Sample_Barcode are required. The first row in the mutation file must be the header including column names above (case sensitive).  
Eight required columns are as below:  
* Tumor sample barcode is the sample ID in which this mutation occurs.  
* Hugo symbol requires HUGO symbols of genes.   
* Chromosome that contains the gene.  
* Start position requires the lowest numeric position of the reported variant on the genomic reference sequence (1-based coordinate system).   
* Reference allele is the plus strand reference allele at this position. Include the sequence deleted for a deletion, or "-" for an insertion.  
* Tumor allele2 is tumor sequencing (discovery) allele 2.  
* Variant classification describes the translational effect of a variant allele. It is one of Silent, Missense_Mutation, Nonsense_Mutation, Splice_Site, Frame_Shift_Del, Frame_Shift_Ins, In_Frame_Del, and In_Frame_Ins. Other mutation types will not be analyzed.  
* Variant type is Type of mutation, such as SNP, DNP, TNP, INS, and DEL.  
The mutation file needs to be Detailed information about the MAF format could be found at https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification. An example file could be seen in /yourpath/DriverML_v1.0.0/example_data/UVM.txt.  
* If the tumor allele2 is not mutated, it should be replaced by the tumor allele1.

## Computational results in our manuscript
**ExInAtor is a good software and designed specially for the detection of driver lncRNA genes using ‘whole-genome sequencing’ mutation data, although it is also applicable to the identification of protein coding genes. Since the whole-exome mutation data was used in the DriverML manuscript, the results of ExInAtor in our manuscript do not represent its optimal performance.**

## Output
The output of DriverML is a summary of putative driver genes, including the numbers of each mutation type, the value of the statistic, p-value, and FDR adjusted p-value.    
## Interpretation of the Output
There are 11 columns in the output file. The first column is the gene symbol. The second to eighth columns are numbers of mutations. The ninth column(LRT) is values of the statistics. The tenth and eleventh columns are the P-values and adjusted P-values(Benjamini-Hochberg Procedure). The bigger the statistics, the more likely that gene is a driver. You could refer to the adjusted P-values for possibilities. The negative values of LRT represent they are not likely to be drivers and you could ignore them.
## Example
nohup /AbsolutePath/DriverML-master/run.driverml.sh -w /AbsolutePath/DriverML-master -i example/UVM.txt -f /AbsolutePath/GRCh38.fa -r /AbsolutePath/hg19.fa -m 10 -o UVM-summary.txt > UVM-nohup.out
## Citation
Yi Han, Juze Yang, Xinyi Qian, Wei-Chung Cheng, Shu-Hsuan Liu, Xing Hua, Liyuan Zhou, Yaning Yang, Qingbiao Wu, Pengyuan Liu, Yan Lu, DriverML: a machine learning algorithm for identifying driver genes in cancer sequencing studies, Nucleic Acids Research, Volume 47, Issue 8, 07 May 2019, Page e45, https://doi.org/10.1093/nar/gkz096
## Contact
If you have any questions, please do not hesitate to contact us.
* The email of the  developer is yihan@zju.edu.cn & 250147506@qq.com.
* If you use WeChat, please add me with the ID: hyperwell. 如果你使用微信，请加我微信：hyperwell。
## Last update
Wednesday January 6, 2021
