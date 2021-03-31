# CRIA

CRIspr Analysis (CRIA) is a set of bash/python scripts which analyzes NGS reads from amplicon sequencing experiments to identify mutations (SNPs and InDels) in your genome edited samples. The NGS reads are trimmed (AdaterRemoval, seqtk, trimmomatic), mapped to the amplicon sequence/reference (BWA) and are called for variants. The different variants observed in the samples with respect to the reference are color coded for visualization and are quantitatively reported.

## 1. System Requirements

CRIA expects the following programs to be installed in your environment.

  * AdapterRemoval
  * seqtk
  * trimmomatic
  * BWA
  * samtools
  * fastx toolkit

And the below python libraries:

  * os,sys,shutil,math
  * pandas
  * uuid
  * xlsxwriter
 
 ## 2. Initial Setup
 
Cria.sh is the master analysis script which makes internal calls to 2 python scripts mutationColoring.py and visualizationPP.py. Cria.sh needs to be provided with parameters and input files through commandline arguments. The software locations are hard coded in Cria.sh which needs to be modified based on your environment. The 2 input files needed by Cria are:
  
    1. Input File - Should be named as Reference_TargetSite_Input.txt. This contains the fastq file names and the type of sample (Control or a Sample Rep)
    2. Information File - Should be named as Reference_TargetSite_Information.txt. This contains the sequences of the reference/amplicon, forward primer and the          target site
 
 Both the input file and information file needs to strictly adhere to the namning guidelines for Cria.sh to work appropriately.
 
 Input File Guidelines:
 
There are 2 columns in this file. The first is the name of the sample fastq file (without the .fq) and the second is the label. Label denotes if the sample is the  control sample or the treated sample. The label should follow the naming convention - Reference_TargetSite_Control_OptionalFreeText to denote the control sample and Reference_TargetSite_Rep_OptionalFreeText to denote the treated sample. The optional free text are ignored by the script. The terms are "Control" and "Rep" are keywords and should not be modified. Every unique Reference_TargetSite combination is considered as an experiment. And each experiment should have 1 control sample and 1 or more treated samples. You can analyze multiple experiments using Cria.      
 
 |Sample Fastq Filename | Sample Label|
 |:---------------------:|:-----------:|
 |AM1920033|HC69_M2_Control_S1|
 |AM1920002|HC69_M2_Rep_S1|
 |AM1920018|HC69_M2_Rep_S2|
 |Sample1|B73_Target11_Control_Sample1|
 |Sample2|B73_Target11_Rep_Sample2|
 


 
