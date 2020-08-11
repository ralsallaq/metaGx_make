#!/bin/bash
module load python/3.5.2
readScratch=$1
analysisDir=$2
assemblyChoice=$3
ecpwyF=$4
que=$5
ncpus=1
tag="bileAcids"
echo scripts/"$tag"_ecpwy_analysis.py
echo "analyzing predicted/annotated genes for existence of $tag pathways"
files=`ls -p $analysisDir/eggnog_map.*/egmap.emapper.annotations.Wcov.csv`
nSamples=`echo $files|wc -w`
echo "Processing $nSamples samples having paired R1 and R2 reads"

function func1 {
local files=$1
local analysisDir=$2
local ecpwyF=$3
local tag=$4

#outF=$analysisDir/allsamples_pwys_missingEnzymes_$tag.csv
outF=$analysisDir/allsamples_pwys_foundIDs_$tag.csv
python scripts/"$tag"_ecpwy_analysis.py -o $outF -i $files -ec $ecpwyF
}

export LSB_JOB_REPORT_MAIL="N"
export -f func1
bsub -q $que  -n $ncpus -R "span[hosts=1] rusage[mem=5000]" -P microbiome -J "bileAcidsPWY" -oo "$readScratch/bileAcidsPWY.log" "func1 '$files' '$analysisDir' '$ecpwyF' '$tag'"
