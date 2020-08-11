#!/bin/bash
activate=/hpcf/apps/conda2/install/4.1.1/bin/activate
echo "activate metaG env to get bbtools"
source $activate metaG
readScratch=$1
analysisDir=$2
que=$3
files=`ls -p $analysisDir/*.decontaminated.fq`
nSamples=`echo $files|wc -w`
echo "Processing $nSamples samples having paired R1 and R2 reads"
echo "quality trim reads"
ncpus=8

function func1 {
local runN=$1
local files=$2
local readScratch=$3
local ncpus=$5
file=`echo $files | cut -f $runN -d " "`
analysisDir=`dirname "$file"`
sname=`basename $file .decontaminated.fq`
echo "processing reads from sample $sname in file $file"
#This phase handles overlapping reads,
#and also nonoverlapping reads, if there is sufficient coverage and sufficiently short inter-read gaps
#For very large datasets, "prefilter=1" or "prefilter=2" can be added to conserve memory.
echo "quality trim reads with quality less than 10 and length less than 50"
bbduk.sh in=$file out=$analysisDir/$sname.qtrimmed.fq.gz qtrim=r trimq=10 minlen=50 ordered
echo "create symbolic link in scratch"
if [ -f $analysisDir/$sname.qtrimmed.fq.gz ]; then
    unlink $readScratch/$sname.temp.fq.gz; ln -s $analysisDir/$sname.qtrimmed.fq.gz $readScratch/$sname.temp.fq.gz 
fi

}
export LSB_JOB_REPORT_MAIL="N"
export -f func1
bsub -q $que  -n $ncpus -R "span[hosts=1] rusage[mem=10000]" -P microbiome -J "qtrimmed[1-$nSamples]" -oo "$readScratch/qtrimmed.%I.log" "func1 \$LSB_JOBINDEX '$files' '$readScratch' '$ncpus'"
