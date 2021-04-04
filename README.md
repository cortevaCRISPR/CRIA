# CRIA

CRIspr Analysis (CRIA) is a set of bash/python scripts which analyzes NGS reads from amplicon sequencing experiments to identify mutations (SNPs and InDels) in your genome edited samples. The NGS reads are trimmed (AdaterRemoval, seqtk, trimmomatic), mapped to the amplicon sequence/reference (BWA) and are called for variants (samtools mpileup). The different variants observed in the samples with respect to the reference are color coded for visualization and are quantitatively reported.

## 1. System Requirements

CRIA expects the following programs to be installed in your environment.

  * AdapterRemoval
  * seqtk
  * trimmomatic
  * BBMap
  * BWA
  * samtools
  * fastx toolkit

And the below python libraries:

  * os,sys,shutil,math
  * pandas
  * uuid
  * xlsxwriter
 
 ## 2. Initial Setup
 
Cria.sh is the master analysis script which makes internal calls to 2 python scripts mutationColoring.py and visualizationPP.py. Cria.sh needs to be provided with parameters and input files through commandline arguments. The software/tools locations used by Cria are hard coded in Cria.sh and needs to be modified based on your environment. The 2 input files needed by Cria are:
  
    1. Input File - Should be named as Reference_TargetSite_Input.txt. This contains the fastq file names and the type of sample (Control or a Sample Rep)
    2. Information File - Should be named as Reference_TargetSite_Information.txt. This contains the sequences of the reference/amplicon, forward primer and the          target site
 
 Both the input file and information file needs to strictly adhere to the naming guidelines for Cria.sh to work appropriately.
 
 <ins>Input File Guidelines:</ins>
 
There are 2 columns in this file. The first is the name of the sample fastq file (without the .fq) and the second is the label. Label denotes if the sample is the  control sample or the treated sample. The label should follow the naming convention - Reference_TargetSite_**Control**\_OptionalFreeText to denote the control sample and Reference_TargetSite_**Rep**\_OptionalFreeText to denote the treated sample. The optional free text is ignored by the script. The terms "Control" and "Rep" in the 2nd column are keywords and should not be modified. Every unique Reference_TargetSite combination is considered as an experiment. And each experiment should have 1 control sample and 1 or more treated samples. You can analyze multiple experiments using Cria.      
 
 |Sample Fastq Filename | Sample Label|
 |:---------------------:|:-----------:|
 |AM1920033|HC69_M2_Control_S1|
 |AM1920002|HC69_M2_Rep_S1|
 |AM1920018|HC69_M2_Rep_S2|
 |Sample1|B73_Target11_Control_Sample1|
 |Sample2|B73_Target11_Rep_Sample2|
 
 <ins>Information File Guidelines:</ins>
 
The information file is where the amplicon/reference sequence, forward primer sequence and target site sequences are provided in a fasta format. The sequences are  expected to be in UPPER CASE. The fasta header for these sequences need to follow the below conventions.

 * ReferenceName_Reference - eg: HC69_Reference ("Reference" is a keyword) 
 * ReferenceName_TargetSite_fwdPrimer - eg: HC69_M2_fwdPrimer ("fwdPrimer" is a keyword)
 * ReferenceName_TargetSite_TargetSite - eg: HC69_M2_TargetSite ("TargetSite" is a keyword")

Please make sure the reference name and target site name used in the information and input file match with each other. Any discrepancies between them causes issues for Cria.  

<ins>NGS Sample Fastq Reads:<ins>
 
The NGS reads are expected to start with the forward primer (some additional bases may be present before the forward primer which will be removed by the script during execution) followed by the amplicon sequence containing the target site. Only reads completely containing the target site sequence are taken into analysis. Cria expects the NGS reads to contain atleats 5 bp following the target site to call a target site completely contained by the read. 

All your sample NGS fastq files must be placed in the folder "Raw_Data" before you start the execution of Cria. The folder is expected to be in the same folder as Cria.sh and other scripts. The fastq files should have the extension ".fq".  


## 3. Execution

sh Cria.sh Reference_TargetSite_Input.txt ProjectName AnalysisDirectory ExpectedEditSequence NumersOfBasePairsToTrim MutationThreshold% IncludeSNPsFlag ReadsType 

<ins>Example command:</ins>

sh Cria.sh HC69_M2_Input.txt HC69_M2 /path/to/your/cria/folder NNNN 0 0.005 0 SE

The Cria.sh script expects 8 commandline parameters to be passed to it.

1. Reference_TargetSite_Input.txt - The input file with the sample fastq names and labels
2. ProjectName - Any name that makes sense to you
3. AnalysisDirectory - Location of the CRIA scripts and Raw_Data folder. This is where the results will be placed
4. ExpectedEditSequence - "NNNN" if random edit expected (SDN1) or the actual edit sequence if SDN2 or SDN3 (eg: "ATCAGACG")
5. NumberOfBasePairsToTrim - Number of extra basepairs before the start of the forward primer sequence that needs to be removed (eg: 6). If there are none, specify 0
6. MutationThreshold% - Only mutations that passes this Threshold% will be reported by Cria (eg: 5 to denote 5 %, 0.1 to denote 0.1%)
7. IncludeSNPsFlag - Specify 1 to Report SNP's in addition to edits and 0 for no SNP's
8. ReadsType - Specify SE for single end reads and PE for paired end reads. PE reads are merged into a single read by BBmerge before analysis by Cria.  
     
## 4. Results

The final results file is an excel file named Results_ProjectName.xlsx. Multiple sheets containing the analysis results are present in the excel file. 

  1. 0_CRIA-Settings -  Basic execution details of Cria are provided
  2. 1_Exp-Summary - A table with the sample name, experiment name, NHEJ%, HDR% (if applicable) and total mutated # of reads found in that sample is reported
  3. 2_CRIA-Stats - A statistics sheet with information on average/min/max read lengths, # raw reads/reads with start of target site/reads with target site covered are reported. This information could be used for ant troubleshooting purposes. 
  4. 3_All-Mutations - A sheet containing all the mutation data across all samples and experiments are provided. 
  5. Individual Sample Mutation Data - Each sample gets its own sheet showing the smaple information, expected target sequence, what kind of allelic change (insertions, deletions, SNPs) was observed, the actual allelic sequence showing the change (wild type target sequence denoted in green and any change denoted in red), numbers and % of total reads, allelic reads, control reads and mutation % are provided.  
  
A lot of temporary and intermediate files will be generated by Cria. But they all will be removed at the end of the analysis. Some intermediate files are retained in the results folder for troubleshooting purposes. 

## 5. Demo

The demo data is present under the folder DEMO. It should take less than 5 minutes to run the script and get the results. The following command was executed to get the demo data.

sh Cria.sh HC69_M2_Input.txt HC69_M2 /path/to/your/cria/folder NNNN 0 0.005 0 SE

