#!/bin/bash
activate=/hpcf/apps/conda2/install/4.1.1/bin/activate
echo "activate metaG env to get diamond"
source $activate metaG
readScratch=$1
analysisDir=$2
assemblyChoice=$3
echo "mapping prokka inferred CDS to uniref90"
files=`ls -p $readScratch/*.temp.fq.gz`
nSamples=`echo $files|wc -w`
echo "Processing $nSamples samples having paired R1 and R2 reads"
que=$4
ncpus=$5

function func1 {
local runN=$1
local files=$2
local assemblyChoice=$3
local analysisDir=$4
local ncpus=$5
uniref90DB=/research/home/margogrp/ralsalla/metaG/databases/uniref90.9.2019.dmnd
file=`echo $files | cut -f $runN -d " "`
scratch=`dirname "$file"`
sname=`basename $file .temp.fq.gz`

echo "map annotated peptides using diamond algorithm to uniref"
echo "diamond is much faster than blasp (blast proteins)"
mkdir $analysisDir/uniref90_map.$sname
echo "the prokka annotation file *.faa encompasses the translated protein coding sequences CDS" 
diamond blastp --query $analysisDir/prokka_Annot.$sname/prokka.faa --db $uniref90DB  --daa $analysisDir/uniref90_map.$sname/unirefmap
}
export LSB_JOB_REPORT_MAIL="N"
export -f func1
bsub -q $que  -n $ncpus -R "span[hosts=1] rusage[mem=10000]" -P microbiome -J "mapAnnot[1-$nSamples]" -oo "$readScratch/uniref.%I.log" "func1 \$LSB_JOBINDEX '$files' '$assemblyChoice' '$analysisDir' '$ncpus'"
