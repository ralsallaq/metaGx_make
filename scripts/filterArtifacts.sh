#!/bin/bash
echo "first load java"
module load java/1.8.0_181
activate=/hpcf/apps/conda2/install/4.1.1/bin/activate
echo "activate metaG env to get bbtools"
source $activate metaG
readScratch=$1
analysisDir=$2
que=$3
ncpus=$4
files=`ls -p $analysisDir/*.noAdpt.fq.gz`
nSamples=`echo $files|wc -w`
echo "Processing $nSamples samples having paired R1 and R2 reads"
echo "removing artifacts synthetic molecules and artifacts"
if [ ! -d resources ];then
    echo "missing bbtools resources directory; you need to move it here ... stopping"
    exit 1
fi

function func1 {
local runN=$1
local files=$2
local readScratch=$3
local ncpus=$5
file=`echo $files | cut -f $runN -d " "`
analysisDir=`dirname "$file"`
sname=`basename $file .noAdpt.fq.gz`
echo "processing reads from sample $sname in file $file"
echo "removing artifacts synthetic molecules and artifacts"
if [ ! -f $analysisDir/$sname.filtered.fq.gz ]; then
    bbduk.sh in=$file out=$analysisDir/$sname.filtered.fq.gz ref=artifacts,phix  k=31 hdist=1 ordered cardinality stats=$analysisDir/$sname.filtered.stats.txt 
else
    echo "file $analysisDir/$sname.filtered.fq.gz already exists; nothing to do"
fi
#if [ -f $analysisDir/$sname.filtered.fq.gz ]; then
#    echo "create symbolic link to the filtered file for the next step"
#    unlink $readScratch/$sname.temp.fq.gz; ln -s $analysisDir/$sname.filtered.fq.gz $readScratch/$sname.temp.fq.gz
#fi
}
export LSB_JOB_REPORT_MAIL="N"
export -f func1
bsub -q $que  -n $ncpus -R "span[hosts=1] rusage[mem=10000]" -P microbiome -J "filterArtifacts[1-$nSamples]" -oo "$readScratch/filterArtifacts.%I.log" "func1 \$LSB_JOBINDEX '$files' '$readScratch' '$ncpus'"
