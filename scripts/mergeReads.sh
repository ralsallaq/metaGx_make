#!/bin/bash
activate=/hpcf/apps/conda2/install/4.1.1/bin/activate
echo "activate metaG env to get bbtools"
source $activate metaG
readScratch=$1
analysisDir=$2
files=`ls -p $readScratch/*.temp.fq.gz`
nSamples=`echo $files|wc -w`
echo "Processing $nSamples samples having paired R1 and R2 reads"
echo "merging overlapping paired reads and quality trimming the ones that do not overlap"
ncpus=8

function func1 {
local runN=$1
local files=$2
local analysisDir=$3
local ncpus=$5
file=`echo $files | cut -f $runN -d " "`
scratch=`dirname "$file"`
sname=`basename $file .temp.fq.gz`
echo "processing reads from sample $sname in file $file"
echo "First merging paired reads and finding which ones overlap (merged) and those which do not (unmerged)"
#This phase handles overlapping reads,
#and also nonoverlapping reads, if there is sufficient coverage and sufficiently short inter-read gaps
#For very large datasets, "prefilter=1" or "prefilter=2" can be added to conserve memory.
bbmerge-auto.sh in=$file out=$analysisDir/$sname.merged.fq.gz outu=$analysisDir/$sname.unmerged.fq.gz strict k=61 extend2=80 rem ordered ihist=$analysisDir/$sname.ihist_merge.txt

if [ -f $analysisDir/$sname.unmerged.fq.gz ]; then
    echo "second quality trim the unmerged"
    bbduk in=$analysisDir/$sname.unmerged.fq.gz out=$analysisDir/$sname.qtrimmed.fq.gz qtrim=r trimq=10 minlen=50 ordered
fi

}
export LSB_JOB_REPORT_MAIL="N"
export -f func1
bsub -q "standard"  -n $ncpus -R "span[hosts=1] rusage[mem=10000]" -P microbiome -J "filterArtifacts[1-$nSamples]" -oo "$readScratch/filterArtifacts.%I.log" "func1 \$LSB_JOBINDEX '$files' '$analysisDir' '$ncpus'"
