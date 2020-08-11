#!/bin/bash
module load fastqc/0.11.5
dataDir=$1 
readScratch=$2
que=$3
ncpus=$4

if [ ! -d qualityChecks ]; then
    mkdir qualityChecks
fi

files=`ls -p $readScratch/*.reads.fq.gz`
#select 3 arbitrary file numbers
nSamples=3

function func1 {
local runN=$1
local files=$2
local ncpus=$4
file=`echo $files | cut -f $runN -d " "`
base=`basename $file .fq.gz`
if [ ! -f qualityChecks/"$base"_fastqc.html ]; then
    fastqc $file -o qualityChecks
else
    echo "the file qualityChecks/"$base"_fastqc.html already exist nothing to do"
fi

}
export LSB_JOB_REPORT_MAIL="N"
export -f func1
bsub -q $que  -n 1 -R "span[hosts=1] rusage[mem=3000]" -P microbiome -J "chkQ[1-$nSamples]" -oo "$readScratch/chkQ_%I.log" "func1 \$LSB_JOBINDEX '$files' '$ncpus'"
