#!/bin/bash
activate=/hpcf/apps/conda2/install/4.1.1/bin/activate
echo "activate metaG env to get quast"
source $activate metaG
readScratch=$1
analysisDir=$2
que=$3
ncpus=$4
shift 
shift
shift
shift
assemblies=$@
echo "evaluating $assemblies assemblies"
files=`ls -p $readScratch/*.reads.fq.gz`
nSamples=`echo $files|wc -w`
echo "Processing $nSamples samples having paired R1 and R2 reads"

function func1 {
local runN=$1
local files=$2
local assemblies=$3
local analysisDir=$4
local ncpus=$5
file=`echo $files | cut -f $runN -d " "`
scratch=`dirname "$file"`
sname=`basename $file .reads.fq.gz`
dircs=()
for ass in $assemblies;do
    if [ $ass = "spades" ]; then
        dircs+=("$analysisDir/$sname.spades_out/scaffolds.fasta")
    fi
    if [ $ass = "megahit" ]; then
        dircs+=("$analysisDir/$sname.megahit_out/contigs.fa")
    fi
done
echo "the following files will be evaluated using quast ${dircs[@]}"
if [ ! -f $analysisDir/$sname.quast/report.html ]; then
    quast.py  -t $ncpus -f --mgm -o $analysisDir/$sname.quast  ${dircs[@]} 
fi
}
export LSB_JOB_REPORT_MAIL="N"
export -f func1
bsub -q $que  -n $ncpus -R "span[hosts=1] rusage[mem=11000]" -P microbiome -J "evaluateAssm[1-$nSamples]" -oo "$readScratch/evaluate.%I.log" "func1 \$LSB_JOBINDEX '$files' '$assemblies' '$analysisDir' '$ncpus'"
