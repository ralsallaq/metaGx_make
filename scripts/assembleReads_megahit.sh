#!/bin/bash
activate=/hpcf/apps/conda2/install/4.1.1/bin/activate
echo "activate megahit env to get megahit"
source $activate megahit
readScratch=$1
analysisDir=$2
que=$3
ncpus=$4
files=`ls -p $readScratch/*.temp.fq.gz`
nSamples=`echo $files|wc -w`
echo "Processing $nSamples samples having paired R1 and R2 reads"

function func1 {
local runN=$1
local files=$2
local analysisDir=$3
local ncpus=$4
file=`echo $files | cut -f $runN -d " "`
scratch=`dirname "$file"`
sname=`basename $file .temp.fq.gz`
echo "processing reads from sample $sname in file $file"
echo "assemble the reads into contigs using megahit"
#megahit -t $ncpus --k-min 45 --k-max 151 --k-step 26 --min-count 2 -m 0.5 -r $analysisDir/$sname.merged.fq.gz --12 $analysisDir/$sname.qtrimmed.fq -o $analysisDir/$sname.megahit_out
megahit -t $ncpus --k-min 45 --k-max 151 --k-step 26 --min-count 2 -r $analysisDir/$sname.merged.fq.gz --12 $analysisDir/$sname.qtrimmed.fq -o $analysisDir/$sname.megahit_out

}
export LSB_JOB_REPORT_MAIL="N"
export -f func1
bsub -q $que  -n $ncpus -R "span[hosts=1] rusage[mem=11000]" -P microbiome -J "assemble_spades[1-$nSamples]" -oo "$readScratch/spades.%I.log" "func1 \$LSB_JOBINDEX '$files' '$analysisDir' '$ncpus'"
