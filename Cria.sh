#!/bin/bash

################## Commandline Parameters ############################

input=$1
project=$2
dir=$3
edit=$4
progname=$(basename $0)

trim=$5
threshold=$6
includeSNP=$7
type=${8}

##################### Hard-Coded Paramters #############################

trimmomatic="/home/nj8612/bin/trimmomatic-0.33" 
AdapterRemoval="/home/nj8612/bin/adapterremoval-2.1.7/build/AdapterRemoval"
bwa="/mnt/common/rh6/ngsdev/bwa-0.7.8/bwa"
samtools="/mnt/common/rh6/ngsdev/samtools-0.1.19/samtools"
fastx_collapser="/scale03/fs0/gpfs0/research/astcrispr/miniconda2/envs/sushiseqProd/bin/fastx_collapser"
seqtk="/scale03/fs0/gpfs0/research/astcrispr/miniconda2/envs/CriaProd/bin/seqtk"
bbmap="/scale03/fs0/gpfs0/research/astcrispr/miniconda2/envs/CriaProd/bin"
data_dir=$dir"/Raw_Data"
default_cutoff=0.001
realfold=30
version="2.0"
date1=$(date +"%Y-%m-%d %H:%M:%S")

if (( $(echo "$threshold <= $default_cutoff" |bc -l) )); then
	threshold=0.001
fi

cutoff=$(echo $threshold/100 | bc -l)

cd $dir

############### Initial File Setup #########################################

echo -e "Date:"'\n'"CRIA Analysis Version:"'\n'"Threshold %:"'\n'"Allele-to-Control Fold Change:" > VariablesHeader

sed -i -e 's/\r$//' -e 's/\ /_/g' "$input"

sed 's/_Rep.*//g' "$input" > "$project"_Samples.txt

cut -f2 "$project"_Samples.txt | sort | uniq | awk '!/Control/' > "$project"_Constructs.txt

cut -f1 "$project"_Samples.txt | sort | uniq > "$project"_SampleList.txt


while read -r sample; do

		if [[ ! -f Raw_Data/"$sample".fq ]]; then

			echo "$sample" >> "$project"_gsdList.txt

		else

			:

		fi

done < "$project"_SampleList.txt

rm  "$project"_SampleList.txt

if [[ -f "$project"_gsdList.txt ]]; then

	while read -r sample; do
		if [[ $type == "PE" ]]; then
		
			$bbmap/bbmerge.sh in1="$data_dir"/"$sample"_paired_r1.fq in2="$data_dir"/"$sample"_paired_r2.fq out="$data_dir"/"$sample".fq outu1="$data_dir"/"$sample"_um_r1.fq outu2="$data_dir"/"$sample"_um_r2.fq
			rm "$data_dir"/"$sample"_paired_r1.fq
			rm "$data_dir"/"$sample"_paired_r2.fq
			rm "$data_dir"/"$sample"_single.fq
			rm "$data_dir"/"$sample"_um_r1.fq
			rm "$data_dir"/"$sample"_um_r2.fq
		fi

	done < "$project"_gsdList.txt

	rm  "$project"_gsdList.txt

else

	:

fi



awk '$2 ~ "Control" {print $1}' "$input" > "$project".allcontrols

if [[ -s "$project".allcontrols ]]; then

	while read -r control; do
		
		grep "$control" "$dir/$input" | cut -f2 > "$control".samplename
		echo "$control" > "$control".sampleid
		paste "$control".sampleid "$control".samplename > "$control".completename

		grep -v -e ^@ -e ^+ "$data_dir"/"$control".fq | awk '{print $0, "=", length($0)}' | awk '{print $3}' | \
		  awk '{if(min==""){min=max=$1}; if($1>max) {max=$1}; if($1< min) {min=$1}; total+=$1; count+=1} \
		  END{printf ("%.5f\t%.5f\t%.5f\n",total/count,min,max)}' > "$control".distribution

		awk 'END{print NR / 4}' "$data_dir"/"$control".fq > "$control".TotalReads


######################## AdapterRemoval ####################


		$AdapterRemoval --threads 8 --file1 "$data_dir"/"$control".fq --minlength 24 \
		  --adapter1 CTGTCTCTTATACACATCTCCGAGCCCACGAGAC --output1 "$control"_trimmed1.fq \
		  --discard "$control"_discard.fq --settings "$control".settings


######################## TRIMMING #############################


		if [[ "$trim" != "no" ]]; then
			$seqtk trimfq -b "$trim" "$control"_trimmed1.fq > "$control"_trimmed.fq

			rm  "$control"_trimmed1.fq

		else

			mv "$control"_trimmed1.fq "$control"_trimmed.fq

		fi

		awk 'END{print NR / 4}' "$control"_trimmed.fq > "$control".AdapterTrimReads

		paste "$control".completename "$control".distribution "$control".TotalReads "$control".AdapterTrimReads > "$control".ReadStatstmp

		rm  "$control"_discard.fq "$control".settings
		rm  "$control".distribution
		rm  "$control".TotalReads
		rm  "$control".AdapterTrimReads
		rm  "$control".sampleid
		rm  "$control".samplename
		rm  "$control".completename*

	done < "$project".allcontrols

else

	:

fi

rm  "$project".allcontrols


###################### Analysis ################################

#### "sets" comes from "$project"_Constructs.txt ####

while read -r sets; do

	awk -v sets="$sets" '$2==sets' "$project"_Samples.txt > "$sets".sets

	cut -f1 "$sets".sets > "$sets".samplelist
	echo "$sets" > "$sets".name

	awk -F_ '{print $1}' "$sets".name > "$sets".construct

	awk -v sets="$sets" -F_ '{
		if (NF <= 2) {
			print $2 > sets".ts1";
		} else {
			exit;
		}
	}' "$sets".name

	construct=$(< "$sets".construct)

	grep -m1 -A1 -i ^\>"$construct"_Reference$ "$project"_Information.txt | sed '1d' > "$construct".Reference

	export ref=`grep -m1 -A1 -i ^\>"$construct"_Reference$ "$project"_Information.txt | sed '1d'`

	sed "1i\>${construct}" "$construct".Reference > "$construct".fa

	grep "$sets""_Control" "$project"_Samples.txt | cut -f1 > "$sets".control

	if [[ -s "$sets".control ]]; then
		control=$(< "$sets".control)
		varcontrol=yes
		varfold="$realfold"
	else
		control=no
		varcontrol=no
		varfold=N/A
	fi

	mkdir -p "$project"_"$construct"

	ts1=$(< "$sets".ts1)
	echo "$ts1" > "$ts1".name

	grep -m1 -A1 -i ^\>"$construct"_"$ts1"_TargetSite$ "$project"_Information.txt | sed '1d' > "$ts1".TargetSite
	targetSite=$(< "$ts1".TargetSite)
	grep -m1 -A1 -i ^\>"$construct"_"$ts1"_fwdPrimer$ "$project"_Information.txt | sed '1d' > "$ts1".fwdPrimer
	grep -o -b -f "$ts1".fwdPrimer "$construct".Reference | cut -d':' -f1 | \
	  awk '{print $1 + 1}' > "$ts1".fwdstart
	sed 's/^/^/' "$ts1".fwdPrimer > "$ts1".fwdPrimer1

	r1sitelength=$(awk '{print length($1)}' "$ts1".TargetSite)
	grep -o -b -f "$ts1".TargetSite "$construct".Reference | head -1 | cut -d':' -f1 | \
	  awk -v len="$r1sitelength" -v ref="$construct" -v first="$ts1" \
	  '{print $1 + 1 > first".targetstart"}{print $1 + len > first".targetend"} \
	  {print ($1 + len) + 10 > first".SRendcut"}' OFS='\t'

	r1fwdts1sc=$(< "$ts1".fwdstart)
	r1revts1ec=$(< "$ts1".SRendcut)

	cut -c "$r1fwdts1sc"-"$r1revts1ec" "$construct".Reference > "$ts1".WholeFwdSite
	sed '/^$/d' "$ts1".WholeFwdSite > "$ts1".WholeFwdSite1
	fold -w1 "$ts1".WholeFwdSite1 > "$ts1".WholeFwdSite2

	paste "$ts1".fwdstart "$ts1".targetstart "$ts1".targetend "$ts1".SRendcut > "$ts1".positions

	awk -v first="$ts1" '{print $2 - $1 > first".fwdDistS"}{print $3 - $1 > first".fwdDistE"}' "$ts1".positions

	awk '{print $1 + 5}' "$ts1".fwdDistE > "$ts1".fwdDistEadd

	if [[ "$control" != "no" ]]; then

		grep -B1 -A2 -f "$ts1".fwdPrimer1 "$control"_trimmed.fq | sed '/^--$/ d' > "$control"_"$ts1"_primerstart.fq
		awk 'END{print NR / 4}' "$control"_"$ts1"_primerstart.fq > "$control"_"$ts1".PrimerStartReads
		paste "$control".ReadStatstmp "$control"_"$ts1".PrimerStartReads > "$control".ReadStats1

		awk '{printf("%s\t%s\t%s\t%s\t%s\t%s\t%s%.5f%s\t%s%.5f%s\n",$1,$2,$3,$4,$5,$6,$7" (",($7/$6)*100,"%)",$8" (",($8/$6)*100,"%)")}' \
		  "$control".ReadStats1 > "$control".ControlReads

		sed -r 's/@|\(0\.000%\)|\(-nan%\)//g' "$control".ControlReads > "$control".ControlReads1
		sed -i '1i\\' "$control".ControlReads1

		echo -e "SampleID"'\t'"SampleName"'\t'"TargetSite"'\t'"TargetSiteSequence"'\t'"AlleleSequence"'\t'"AlleleChange"'\t'"AlleleReadSequence" \
		  '\t'"#AlleleReads"'\t'"%AlleleReads"'\t'"#TotalAlleleReads"'\t'"#ControlReads"'\t'"Meets_${realfold}x_threshold"'\t'"Mutation%" > AlleleHeader

		rm  "$control".ControlReads
		rm  "$control".ReadStatstmp
		rm  "$control".ReadStats1
		rm  "$control"_"$ts1".PrimerStartReads

	else
		
		echo -e "SampleID"'\t'"SampleName"'\t'"TargetSite"'\t'"TargetSiteSequence"'\t'"AlleleSequence"'\t'"AlleleChange"'\t'"AlleleReadSequence" \
		  '\t'"#AlleleReads"'\t'"%AlleleReads"'\t'"#TotalAlleleReads"'\t'"#ControlReads"'\t'"Meets_${realfold}x_threshold"'\t'"Mutation%" > AlleleHeader

	fi

	sed "1i\>${ts1}" "$ts1".WholeFwdSite1 > "$ts1".fa
	awk -v ts="$ts1" 'NR==2{print ts "\t" "0" "\t" (length($1)) "\t" ts"_Amplicon"}' \
	  "$ts1".fa > "$ts1"_Amplicon.bed

	$bwa index "$ts1".fa
	$samtools faidx "$ts1".fa
	
	rm  "$ts1".fwdPrimer
	rm  "$ts1".fwdstart
	rm  "$ts1".targetend
	rm  "$ts1".WholeFwdSite
	rm  "$ts1".WholeFwdSite1
	rm  "$ts1".positions
	rm  "$ts1".targetstart
	rm  "$ts1".SRendcut
	rm  "$construct".Reference


#### sample is from "$sets".samplelist ####
	while read -r sample; do

		grep "$sample" "$input" | cut -f2 > "$sample".samplename
		echo "$sample" > "$sample".sampleid
		paste "$sample".sampleid "$sample".samplename > "$sample".completename

		grep -v -e ^@ -e ^+ "$data_dir"/"$sample".fq | awk '{print $0, "=", length($0)}' | awk '{print $3}' | \
		  awk '{if(min==""){min=max=$1}; if($1>max) {max=$1}; if($1< min) {min=$1}; total+=$1; count+=1} \
		  END{printf ("%.5f\t%.5f\t%.5f\n",total/count,min,max)}' > "$sample".distribution

		awk 'END{print NR / 4}' "$data_dir"/"$sample".fq > "$sample".TotalReads

		$AdapterRemoval --threads 8 --file1 "$data_dir"/"$sample".fq --minlength 24 \
		  --adapter1 CTGTCTCTTATACACATCTCCGAGCCCACGAGAC \
		  --output1 "$sample"_trimmed1.fq --discard "$sample"_discard.fq --settings "$sample".settings

		if [[ "$trim" != "no" ]]; then

			$seqtk trimfq -b "$trim" "$sample"_trimmed1.fq > "$sample"_trimmed.fq

			rm  "$sample"_trimmed1.fq

		else

			mv "$sample"_trimmed1.fq "$sample"_trimmed.fq

		fi

		awk 'END{print NR / 4}' "$sample"_trimmed.fq > "$sample".AdapterTrimReads

		paste "$sample".completename "$sample".distribution "$sample".TotalReads "$sample".AdapterTrimReads > "$sample".tmp1

		rm  "$sample"_discard.fq "$sample".settings
		rm  "$sample".sampleid
		rm  "$sample".distribution
		rm  "$sample".TotalReads
		rm  "$sample".AdapterTrimReads
		rm  "$sample".samplename

		grep -B1 -A2 -f "$ts1".fwdPrimer1 "$sample"_trimmed.fq | sed '/^--$/ d' > "$sample"_"$ts1"_primerstart.fq
		awk 'END{print NR / 4}' "$sample"_"$ts1"_primerstart.fq > "$sample"_"$ts1".PrimerStartReads

		ts1start=$(< "$sample"_"$ts1".PrimerStartReads)

		if [[ "$ts1start" -ne 0 ]]; then

			minlen=$(< "$ts1".fwdDistEadd)
			java -jar "$trimmomatic"/trimmomatic-0.33.jar SE -threads 8 -phred64 "$sample"_"$ts1"_primerstart.fq "$sample"_"$ts1"_minlen.fq \
			  MINLEN:"$minlen"
			awk 'END{print NR / 4}' "$sample"_"$ts1"_minlen.fq > "$sample"_"$ts1".MinLenReads

		else
			echo "0" > "$sample"_"$ts1".MinLenReads
		fi

		ts1min=$(< "$sample"_"$ts1".MinLenReads)

		if [[ "$ts1min" -ne 0 ]]; then
		
			croplen=$(< "$ts1".fwdDistEadd)
			java -jar "$trimmomatic"/trimmomatic-0.33.jar SE -threads 8 -phred64 "$sample"_"$ts1"_minlen.fq "$sample"_"$ts1"_crop.fq \
			  CROP:"$croplen"
			awk 'END{print NR / 4}' "$sample"_"$ts1"_crop.fq > "$sample"_"$ts1".FinalTargetReads

		else
			echo "0" > "$sample"_"$ts1".FinalTargetReads
		fi

		ts1crop=$(< "$sample"_"$ts1".FinalTargetReads)
		
		paste "$sample"_"$ts1".PrimerStartReads "$sample"_"$ts1".MinLenReads "$sample"_"$ts1".FinalTargetReads > "$sample"_"$ts1".Reads

		rm  "$sample"_"$ts1"_primerstart.fq
		rm  "$sample"_"$ts1"_minlen.fq


		 if [[ "$ts1crop" -ne 0 && "$ts1start" -ne 0 ]]; then
		
			$fastx_collapser -i "$sample"_"$ts1"_crop.fq -o "$sample"_"$ts1"_collapse.fa -Q33 

			sed ':a;N;/^>/s/\n/\t/g' "$sample"_"$ts1"_collapse.fa > "$sample"_"$ts1"_collapse.tmpfa

			awk -v threshold="$cutoff" '{printf ("%.5f\n", $1 * threshold)}' "$sample"_"$ts1".MinLenReads > "$sample"_"$ts1".cutoff

			sed -i 's/-/\t/' "$sample"_"$ts1"_collapse.tmpfa

			SEcutoff=$(< "$sample"_"$ts1".cutoff)
			awk -v r1co="$SEcutoff" '{if ($2 >= r1co) print $1 "\t" $3 "\t" $2}' \
			  "$sample"_"$ts1"_collapse.tmpfa > "$sample"_"$ts1"_collapse.tmp1fa

			if [[ -f "$sample"_"$ts1"_collapse.tmp1fa ]]; then
			
				sed 's/>//' "$sample"_"$ts1"_collapse.tmp1fa | \
				  awk -v sample="$sample" -v ts1="$ts1" '{i=sprintf("%06d",i+1)}{print > sample"_"i"_"ts1".SEalleles"; close(sample"_"i"_"ts1".SEalleles")}'
				find -maxdepth 1 -type f -name "${sample}_*.SEalleles" | sort | sed -e 's/\.\///' -e 's/\.SEalleles//' > \
				  "$sample"_"$ts1"_SEAlleles.txt

				 while read -r allele; do

					awk '{print $3}' "$allele".SEalleles > "$allele".SEallelecount
					awk '{print $2}' "$allele".SEalleles > "$allele".SEcompare

					sed 's/^/^/' "$allele".SEcompare > "$allele".SEcompareLR
					grep -f "$allele".SEcompareLR "$sample"_trimmed.fq > "$allele"_SE.longreads
					sort "$allele"_SE.longreads | uniq -c | sort -rnk1 | sed -e 's/^[ \t]*//' -e 's/\ /\t/' | cut -f2 | head -n1 > "$allele"_SE.longreads1
					fold -w1 "$allele"_SE.longreads1 > "$allele"_SE.longreads2
				
					SEts1length=$(awk 'END{print NR}' "$allele"_SE.longreads2)
					SEts1WFSlength=$(awk 'END{print NR}' "$ts1".WholeFwdSite2)
				
					if [[ "$SEts1WFSlength" -le "$SEts1length" ]]; then

						awk -v len="$SEts1WFSlength" 'NR <= len' "$allele"_SE.longreads2 > "$allele"_SE.longreads2len
						cp "$ts1".WholeFwdSite2 "$allele".WholeFwdSite2

					else

						awk -v len="$SEts1length" 'NR <= len' "$ts1".WholeFwdSite2 > "$allele".WholeFwdSite2
						cp "$allele"_SE.longreads2 "$allele"_SE.longreads2len

					fi

					if diff "$allele".WholeFwdSite2 "$allele"_SE.longreads2len >/dev/null; then
					  
						echo "Wild Type" > "$allele"_SE.chardesc
						ts1tss=$(< "$ts1".fwdDistS)
						ts1tes=$(awk '{print $1 +1}' "$ts1".fwdDistE)
						awk -v trgtss="$ts1tss" -v trgtes="$ts1tes" '{if (NR == trgtss || NR == trgtes) print $1"\n""|"; else print $1}' \
						  "$allele"_SE.longreads2len > "$allele"_SE.longreads3
						tr -d '\n' < "$allele"_SE.longreads3 > "$allele"_SE.longreads4

						cp "$ts1".TargetSite "$allele"_SE.targetseq

						echo "$allele"_SE.targetseq >> "$sample"_tar_parselist.txt
						
						if [[ "$control" != "no" ]]; then

							sed -e 's/^/^/' -e 's/$/$/' "$allele"_SE.longreads1 > "$allele"_SE.longreads1control
							grep -cf "$allele"_SE.longreads1control "$control"_"$ts1"_primerstart.fq \
							  > "$allele"_SE.controlcount

							allelectrlSEts1=$(< "$allele".SEallelecount)
							crtlfoldSEts1=$(awk -v fold="$realfold" '{print $1 * fold}' "$allele"_SE.controlcount)

							if [[ "$crtlfoldSEts1" -ge "$allelectrlSEts1" ]]; then
								echo "no" > "$allele"_SE.controlfold
							else
								echo "yes" > "$allele"_SE.controlfold
							fi

							paste "$allele"_SE.longreads4 "$allele"_SE.chardesc "$allele"_SE.longreads1 \
							  "$allele".SEallelecount "$sample"_"$ts1".FinalTargetReads "$allele"_SE.controlcount \
							  "$allele"_SE.controlfold > "$allele"_SE.char

							rm  "$allele"_SE.controlcount
							rm  "$allele"_SE.controlfold

						else

							paste "$allele"_SE.longreads4 "$allele"_SE.chardesc "$allele"_SE.longreads1 \
							  "$allele".SEallelecount "$sample"_"$ts1".FinalTargetReads > "$allele"_SE.char

						fi

						echo "$allele"_SE.char >> "$sample"_char_parselist.txt

					else
					
						awk '{print ">" FILENAME "\n" $1}' "$allele"_SE.longreads1 > "$allele"_SE.fa

						$bwa mem "$ts1".fa "$allele"_SE.fa > "$allele"_SE.sam

						$samtools view -S -b "$allele"_SE.sam > "$allele"_SE.bam
						$samtools mpileup -f "$ts1".fa -l "$ts1"_Amplicon.bed "$allele"_SE.bam > "$allele"_SE.mpileup
						$samtools sort "$allele"_SE.bam "$allele"_SE_sort
						mv "$allele"_SE_sort.bam "$allele"_SE.bam
						$samtools index "$allele"_SE.bam

						SEts1tss=$(awk '{print $1 + 1}' "$ts1".fwdDistS)
						SEts1tes=$(awk '{print $1 + 1}' "$ts1".fwdDistE)
						awk -v tss="$SEts1tss" -v tes="$SEts1tes" -F'\t' '{if ($2 == tss || $2 == tes) \
						  print $3"\t"$5"|"; else print $3"\t"$5}' "$allele"_SE.mpileup > "$allele"_SE.chartmp

						awk '{if ($2 == ".") print $1; else print $0}' "$allele"_SE.chartmp > "$allele"_SE.chartmp1

						sed -e 's/\.//' -e 's/[0-9]*//g' -e 's/\t//' "$allele"_SE.chartmp1 > "$allele"_SE.chartmp2

						sed 's/\(.\)/\1 /g;s/ $//' "$allele"_SE.chartmp2 > "$allele"_SE.chartmp3



						awk '{
							if ($2 ~ /[+]/ && $0 !~ /[\|]$/) {
								printf "%s%s\n",$1," "$2;
							} else if ($2 ~ /[-]/ && $0 !~ /[\|]$/) {
								printf "%s\n",$1;
							} else if ($2 ~ /[+]/ && $0 ~ /[\|]$/) {
								printf "%s%s%s\n",$1,$NF," "$2;
							} else if ($2 ~ /[-]/ && $0 ~ /[\|]$/) {
								printf "%s%s\n",$1,$NF;
							} else {
								print $0;
							}
							}' "$allele"_SE.chartmp3 > "$allele"_SE.chartmp4



						SEts1pos=$(awk '{print $1 + 2}' "$ts1".fwdDistE)
						SEts1pos2=$(awk '{print $1 + 1}' "$ts1".fwdDistS)
						awk -v pos="$SEts1pos" -v pos2="$SEts1pos2" '{
							if ($2 ~ /[+-]/ && $0 !~ /[\|]$/ && NR < pos && NR >= pos2) {
								for (i=2;i<=NF;i++)
								printf "%s",$i;
								printf "; ";
							} else if ($2 ~ /[+-]/ && $0 !~ /[\|]$/ && NR < pos2) {
								printf "%s","Indel outside Target";
								printf "; ";
							} else if ($2 ~ /[+-]/ && $0 !~ /[\|]$/ && NR >= pos) {
								printf "%s","Indel outside Target";
								printf "; ";
							} else if ($2 ~ /[+-]/ && $0 ~ /[\|]$/) {
								for (i=2;i<NF;i++)
								printf "%s",$i;
								printf "; ";
							} else if ($2 ~ /[A-Z]/ && NR < pos && NR >= pos2) {
								printf "%s%s"," "$1" to ",$2" ";
								printf "; ";
							} else if ($2 ~ /[A-Z]/ && NR < pos2) {
								printf "%s","SNP outside Target";
								printf "; ";
							} else if ($2 ~ /[A-Z]/ && NR >= pos) {
								printf "%s","SNP outside Target";
								printf "; ";
							} else if ($2 == "$" && NR < pos && NR >= pos2) {
								printf "%s","Indel extends beyond Target";
								printf "; ";
							} else if ($2 == "$" && NR < pos2) {
								printf "%s","Indel before Target";
								printf "; ";
							} else {
								next;
							}
							}' "$allele"_SE.chartmp3 > "$allele"_SE.chardesc

						sed -i "s/^/\'/" "$allele"_SE.chardesc
						sed -i -e '$a\' "$allele"_SE.chardesc
						awk '{gsub(/; $/,""); print}' "$allele"_SE.chardesc > "$allele"_SE.chardesc1



						SEts1pos1=$(awk '{print $1 + 2}' "$ts1".fwdDistS)
						SEts1pos3=$(< "$ts1".fwdDistS)
						awk -v pos1="$SEts1pos1" -v pos3="$SEts1pos3" '{
							if ($2 ~ /[+]/ && $0 !~ /[\|]$/ && NR != pos3) {
								print $1 "\n" $2;
							} else if ($2 ~ /[+]/ && $0 ~ /[\|]$/ && NR < pos1) {
								print $NF;
								for (i=1;i<NF;i++)
								print $i;
							} else if ($2 ~ /[+]/ && $0 ~ /[\|]$/ && NR > pos1) {
								for (i=1;i<=NF;i++)
								printf "%s\n",$i;
							} else if ($2 ~ /[A-Z]/ && $0 !~ /[\|]$/) {
								print $2;
							} else if ($2 ~ /[A-Z]/ && $0 ~ /[\|]$/ && NR < pos1) {
								print $NF;
								print $2;
							} else if ($2 ~ /[A-Z]/ && $0 ~ /[\|]$/ && NR > pos1) {
								print $2;
								print $NF;
							} else if ($2 == "*" && $0 !~ /[\|]$/) {
								print "-";
							} else if ($2 == "*" && $0 ~ /[\|]$/ && NR < pos1) {
								print $NF;
								print "-";
							} else if ($2 == "*" && $0 ~ /[\|]$/ && NR > pos1) {
								print "-";
								print $NF;
							} else if ($2 == "$" && $0 !~ /[\|]$/ && NR < pos1) {
								print $1 "\n" "[";
							} else if ($2 == "$" && $0 ~ /[\|]$/ && NR < pos1) {
								print $NF "\n" $1 "\n" "[";
							} else if ($2 == "$" && $0 ~ /[\|]$/ && NR > pos1) {
								print $1 "\n" $NF "\n" "[";
							} else if ($2 == "|" && NR < pos1) {
								print $2 "\n" $1;
							} else if ($2 == "|" && NR > pos1) {
								print $1 "\n" $2;
							} else {
								print $1;
							}
						}' "$allele"_SE.chartmp4 > "$allele"_SE.chartmp5



						tr -d '\n' < "$allele"_SE.chartmp5 | awk '{print $0}' > "$allele"_SE.chartmp6



						awk 'BEGIN{FS=""} {
							if ($1 ~/[-]/) {
								next;
							} else {
								print $1;
							}
						}' "$allele"_SE.chartmp5 | tr -d '\n' | awk '{print $0}' > "$allele"_SE.targetseqtmp



						grep -Po "(?<=|)[^|]*(?=|)" "$allele"_SE.targetseqtmp | awk 'NR==2' > "$allele"_SE.targetseq

						echo "$allele"_SE.targetseq >> "$sample"_tar_parselist.txt


						if [[ "$control" != "no" ]]; then

							sed -e 's/^/^/' -e 's/$/$/' "$allele"_SE.longreads1 > "$allele"_SE.longreads1control
							grep -cf "$allele"_SE.longreads1control "$control"_"$ts1"_primerstart.fq \
							  > "$allele"_SE.controlcount

							allelectrlSEts1=$(< "$allele".SEallelecount)
							crtlfoldSEts1=$(awk -v fold="$realfold" '{print $1 * fold}' "$allele"_SE.controlcount)

							if [[ "$crtlfoldSEts1" -ge "$allelectrlSEts1" ]]; then
								echo "no" > "$allele"_SE.controlfold
							else
								echo "yes" > "$allele"_SE.controlfold
							fi

							paste "$allele"_SE.chartmp6 "$allele"_SE.chardesc1 "$allele"_SE.longreads1 \
							  "$allele".SEallelecount "$sample"_"$ts1".FinalTargetReads "$allele"_SE.controlcount \
							  "$allele"_SE.controlfold > "$allele"_SE.char

							rm  "$allele"_SE.controlcount
							rm  "$allele"_SE.controlfold

						else

							paste "$allele"_SE.chartmp6 "$allele"_SE.chardesc1 "$allele"_SE.longreads1 \
							  "$allele".SEallelecount "$sample"_"$ts1".FinalTargetReads > "$allele"_SE.char

						fi

						echo "$allele"_SE.char >> "$sample"_char_parselist.txt

						rm  "$allele"_SE.bam*
						rm  "$allele"_SE.fa
						rm  "$allele"_SE.sam
						rm  "$allele"_SE.mpileup
						rm  "$allele"_SE.chartmp*
						rm  "$allele"_SE.targetseqtmp

					fi

					rm  "$allele"_SE.longreads1
					rm  "$allele".SEalleles
					rm  "$allele".SEallelecount
					rm  "$allele".SEcompare
					rm  "$allele".SEcompareLR
					rm  "$allele"_SE.longreads*
					rm  "$allele"_SE.chardesc*
					rm  "$allele".WholeFwdSite2
				
				 done < "$sample"_"$ts1"_SEAlleles.txt

				mapfile -t < "$sample"_char_parselist.txt
				cat "${MAPFILE[@]}" > "$sample"_"$ts1"_SE.char1

				mapfile -t < "$sample"_tar_parselist.txt
				cat "${MAPFILE[@]}" > "$sample"_"$ts1"_SE.tarseq

				if [[ "$control" != "no" ]]; then

					awk -F'\t' '{printf("%s\t%s\t%s\t%s\t%.5f%s\t%s\t%s\n",$1,$2,$3,$4,($4/$5)*100,"%\t"$5,$6,$7)}' \
					  "$sample"_"$ts1"_SE.char1 > "$sample"_"$ts1"_SE.characterization

				else

					awk -F'\t' '{printf("%s\t%s\t%s\t%s\t%.5f%s\n",$1,$2,$3,$4,($4/$5)*100,"%\t"$5)}' \
					  "$sample"_"$ts1"_SE.char1 > "$sample"_"$ts1"_SE.characterization

				fi

				rm  "$sample"_"$ts1"_collapse.tmp1fa
				rm  "$sample"_"$ts1"_SEAlleles.txt
				rm  "$sample"_*_"$ts1"_SE.char
				rm  "$sample"_*_"$ts1"_SE.targetseq
				rm  "$sample"_tar_parselist.txt
				rm  "$sample"_char_parselist.txt
			
			else
			
				echo -e "No reads cover target site\tN/A\tN/A\t0\tN/A" > "$sample"_"$ts1"_SE.char1
				paste "$sample"_"$ts1"_SE.char1 "$sample"_"$ts1".FinalTargetReads > "$sample"_"$ts1"_SE.characterization

			fi

			awk '{print $1, "Read1"}' "$ts1".name > "$ts1".forwardname

			paste "$ts1".forwardname "$sample"_"$ts1"_SE.tarseq "$sample"_"$ts1"_SE.characterization > "$sample"_"$ts1"_SE.compile

			awk '{FS=OFS="\t"}{if(NR==1){print;next}; $1="@"; print}' "$sample"_"$ts1"_SE.compile > "$sample"_"$ts1"_SE.compile1

			rm  "$sample"_"$ts1"_crop.fq
			rm  "$sample"_"$ts1"_collapse.fa
			rm  "$sample"_"$ts1"_collapse.tmpfa
			rm  "$sample"_"$ts1".cutoff
			rm  "$sample"_"$ts1"_SE.compile
			rm  "$sample"_"$ts1"_SE.char1
			rm  "$sample"_"$ts1"_SE.characterization
			rm  "$ts1".forwardname
			rm  "$sample"_"$ts1"_SE.tarseq

		elif [[ "$ts1crop" -eq 0 && "$ts1start" -ne 0  ]]; then
		
			echo -e "N/A"'\t'"No reads cover target site"'\t'"N/A"'\t'"0"'\t'"N/A" > "$sample"_"$ts1"_SE.characterization
			awk '{print $1, "5'\'' to 3'\''"}' "$ts1".name > "$ts1".forwardname
			paste "$ts1".forwardname "$ts1".TargetSite "$sample"_"$ts1"_SE.characterization "$sample"_"$ts1".PrimerStartReads > "$sample"_"$ts1"_SE.compile1

			rm  "$ts1".forwardname
			rm  "$sample"_"$ts1"_SE.characterization

		else
		
			echo -e "N/A"'\t'"No reads begin with primer sequence"'\t'"N/A"'\t'"0"'\t'"N/A"'\t'"0" > "$sample"_"$ts1"_SE.characterization
			awk '{print $1, "5'\'' to 3'\''"}' "$ts1".name > "$ts1".forwardname
			paste "$ts1".forwardname "$ts1".TargetSite "$sample"_"$ts1"_SE.characterization > "$sample"_"$ts1"_SE.compile1

			rm  "$ts1".forwardname
			rm  "$sample"_"$ts1"_SE.characterization

		fi

		paste "$sample".tmp1 "$sample"_"$ts1".Reads > "$sample"_"$sets".CompiledReads

		awk '{printf("%s\t%s\t%s\t%s\t%s\t%s\t%s%.5f%s\t%s%.5f%s\t%s%.5f%s\t%s%.5f%s", \
		  $1,$2,$3,$4,$5,$6,$7" (",($7/$6)*100,"%)",$8" (",($8/$6)*100,"%)",$9" (",($9/$6)*100,"%)",$10" (",($10/$6)*100,"%)")}' \
		  "$sample"_"$sets".CompiledReads > "$sample"_"$sets".CompiledReads1

		sed -r 's/@|\(0\.000%\)|\(-nan%\)//g' "$sample"_"$sets".CompiledReads1 > "$sample"_"$sets".CompiledReads2

		awk 'END{print NR}' "$sample"_"$ts1"_SE.compile1 > "$sample".range

		maxrange=$(< "$sample".range)

		awk -v mxrg="$maxrange" -F'\t' '{
			NR+1
				print $0;
			for (i=NR+1;i<=mxrg;i++)
				print "@\t@";
			}' "$sample".completename > "$sample".completename2

		paste "$sample".completename2 "$sample"_"$ts1"_SE.compile1 > "$sample"_SE_"$sets".completealldatacompile
		sed 's/@//g' "$sample"_SE_"$sets".completealldatacompile > "$sample"_SE_"$sets".completealldatacompile1

		cat AlleleHeader "$sample"_SE_"$sets".completealldatacompile1 > "$sample"_"$sets".AlleleReport
		awk '{gsub(/;$/,""); print}' "$sample"_"$sets".AlleleReport > "$sample"_"$sets".AlleleReport1
		sed 's/\t/,/g' "$sample"_"$sets".AlleleReport1 > "$sample"_"$sets"_AlleleReport.csv

		mv "$sample"_"$sets"_AlleleReport.csv "$project"_"$construct"/


		rm  "$sample"_"$ts1".MinLenReads
		rm  "$sample"_"$ts1".PrimerStartReads
		rm  "$sample"_"$ts1".Reads
		rm  "$sample"_"$ts1"_SE.compile1
		rm  "$sample"_"$ts1".FinalTargetReads
		rm  "$sample"_"$sets".CompiledReads
		rm  "$sample"_"$sets".CompiledReads1
		rm  "$sample".completename*
		rm  "$sample".range
		rm  "$sample".tmp1
		rm  "$sample"_trimmed.fq
		rm  "$sample"_SE_"$sets".completealldatacompile*
		rm  "$sample"_"$sets".AlleleReport*
		
		if [[ "$includeSNP" -eq 1 ]]; then
			python visualizationPP.py "$project"_"$construct"/"$sample"_"$sets"_AlleleReport.csv "$edit" "$ref" "$targetSite" 1 > "$project"_"$construct"/"$sample"_stats.txt
			echo -e "python visualizationPP.py "$dir"/"$project"_"$construct"/"$sample"_"$sets"_AlleleReport.csv "$edit" "$ref" "$targetSite" "1" > "$dir"/"$project"_"$construct"/"$sample"_stats.txt" >> "$dir"/VVcmd.txt

		else
			python visualizationPP.py "$project"_"$construct"/"$sample"_"$sets"_AlleleReport.csv "$edit" "$ref" "$targetSite" 0 > "$project"_"$construct"/"$sample"_stats.txt
			echo -e "python visualizationPP.py "$dir"/"$project"_"$construct"/"$sample"_"$sets"_AlleleReport.csv "$edit" "$ref" "$targetSite" "0" > "$dir"/"$project"_"$construct"/"$sample"_stats.txt" >> "$dir"/VVcmd.txt
		fi
		
	done < "$sets".samplelist

	awk 'FNR==1{print ""}1' *_"$sets".CompiledReads2 > "$project"_"$sets".CompiledReads
	
	if [[ "$control" != "no" ]]; then

		cat "$project"_"$sets".CompiledReads "$control".ControlReads1 > "$project"_"$sets".CompiledReads1
		cat "$project"_"$sets".CompiledReads1 > "$project"_"$sets".StatsReport

		rm  "$control".ControlReads1
		rm  "$project"_"$sets".CompiledReads1
		rm  "$control"_"$ts1"_primerstart.fq
		rm  "$control"_trimmed.fq

	else

		cat "$project"_"$sets".CompiledReads > "$project"_"$sets".StatsReport

	fi

	mv "$project"_"$sets".StatsReport "$project"_"$sets"_"$date"_StatsReport.txt

	echo -e "$date1"'\n'"$version"'\n'"$threshold"'\n'"$varfold" > "$project"_"$sets".variables
	paste VariablesHeader "$project"_"$sets".variables > "$project"_"$sets"_"$date"_Settings.txt

	mv "$project"_"$sets"_"$date"_StatsReport.txt "$project"_"$construct"/
	mv "$project"_"$sets"_"$date"_Settings.txt "$project"_"$construct"/


	rm  "$ts1".fa*
	rm  "$ts1".name
	rm  "$ts1".TargetSite
	rm  "$ts1".fwdPrimer1
	rm  "$ts1".fwdDistE
	rm  "$ts1".fwdDistS
	rm  "$ts1".fwdDistEadd
	rm  "$ts1"_Amplicon.bed
	rm  "$ts1".WholeFwdSite2
	rm  AlleleHeader
	rm  "$construct".fa*
	rm  "$sets".sets
	rm  "$sets".ts*
	rm  "$sets".name
	rm  "$sets".construct
	rm  "$sets".samplelist
	rm  "$sets".control
	rm  "$project"_"$sets".CompiledReads
	rm  *_"$sets".CompiledReads2
	rm  "$project"_"$sets".variables


done < "$project"_Constructs.txt


rm "$project"_Samples.txt
rm "$project"_Constructs.txt
rm VariablesHeader

sh postProcessCria.sh "$project"_"$construct" $dir

echo "Analysis Complete!"


