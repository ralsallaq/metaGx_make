#!/bin/bash
module load python/3.5.2
readScratch=$1
analysisDir=$2
assemblyChoice=$3
ecpwyF=$4
echo $1 $2 $3 $4 $5 $6 $7
if [[ ! $5 -le 1 ]]; then
    threshold=`bc <<< "scale=2; $5/100"`
    echo "threshold is $threshold"
else
    echo "threshold must be a percentage e.g. 10, 30, 75 ,..etc"
    exit 1
fi
que=$6
ncpus=1 #or $7
echo "analyzing predicted/annotated genes for existence of butyrate pathways"
files=`ls -p $analysisDir/eggnog_map.*/egmap.emapper.annotations.Wcov.csv`
nSamples=`echo $files|wc -w`
echo "Processing $nSamples samples having paired R1 and R2 reads"

function func1 {
local files=$1
local analysisDir=$2
local ecpwyF=$3
local threshold=$4

outF=$analysisDir/allsamples_butyrate_featureTab.csv
python ./scripts/butyrate_ecpwy_featureTab.py -o $outF -i $files -ec $ecpwyF -t $threshold
}

export LSB_JOB_REPORT_MAIL="N"
export -f func1
bsub -q $que  -n $ncpus -R "span[hosts=1] rusage[mem=9000]" -P microbiome -J "butyratePWY" -oo "$readScratch/butyrateFT.log" "func1 '$files' '$analysisDir' '$ecpwyF' '$threshold'"
