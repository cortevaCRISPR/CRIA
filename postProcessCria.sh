#!/bin/bash

export project=$1
export dir=$2


find $dir -type f -name "VVcmd.txt" -exec cat {} \; >> $dir/Vcmd.txt
find $dir -type f -name "Filtered*" | xargs cp -t $dir
find $dir -type f -name "*StatsReport.txt" | xargs cp -t $dir
find $dir -type f -name "*_Settings.txt" -exec cat {} + | head -n 4 > $dir/Settings.txt

echo -e "SampleID"'\t'"SampleName"'\t'"AverageReadLength"'\t'"MinimumReadLength"'\t'"MaximumReadLength"'\t'"TotalRawReads"'\t'"TrimmedReads" \ '\t'"#TS_StartSite"'\t'"#TS_CoverSite"'\t'"#TS_Analyzed" > $dir/ReadsHeader
cat $dir/ReadsHeader $dir/*StatsReport.txt > $dir/RunStats.txt

echo -e "Biotracker\tExperiment\tNHEJ%\tHDR%\tMutatedReads" > $dir/AllStats.txt
find $dir -type f -name "*_stats.txt" -exec cat {} \; >> $dir/AllStats.txt

python mutationColoring.py $dir $project

rm $dir/*StatsReport.txt
rm $dir/ReadsHeader
rm $dir/VVcmd.txt
rm $dir/Vcmd.txt
rm $dir/Filtered*
rm $dir/AllStats.txt
rm $dir/RunStats.txt
rm $dir/Settings.txt
#rm $dir/*txt


	