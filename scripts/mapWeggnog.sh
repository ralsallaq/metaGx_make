#!/bin/bash
module load python/2.7.15-rhel7
readScratch=$1
analysisDir=$2
assemblyChoice=$3
echo "determining coverage and taxonomic make up for $assemblyChoice assembly"
files=`ls -p $readScratch/*.reads.fq.gz`
nSamples=`echo $files|wc -w`
echo "Processing $nSamples samples having paired R1 and R2 reads"
ncpus=8

function func1 {
local runN=$1
local files=$2
local assemblyChoice=$3
local analysisDir=$4
local ncpus=$5
eggnogDataDir=/research/home/margogrp/ralsalla/metaG/eggnog-mapper/data/
emapperEX=/research/home/margogrp/ralsalla/metaG/eggnog-mapper/emapper.py
echo "running emapperEX from the directory $emapperEX"
file=`echo $files | cut -f $runN -d " "`
scratch=`dirname "$file"`
sname=`basename $file .reads.fq.gz`

echo "map annotated peptides from sample "$sname" using eggnog mapper with the diamond algorithm"
echo "eggnog assigns orhtolog-groups (OGs=cluster of proteins that have the same function but from different species) to peptides"
if [ ! -f  analysis/eggnog_map.$sname/egmap.emapper.annotations ]; then
mkdir $analysisDir/eggnog_map.$sname
echo "the prokka CDS gene prediction file *.faa encompasses the translated protein coding sequences CDS; the --resume option is for generating output from previous run" 
$emapperEX --resume -i $analysisDir/prokka_Annot.$sname/prokka.faa -m diamond --data_dir $eggnogDataDir --cpu $ncpus -o $analysisDir/eggnog_map.$sname/egmap 
fi

}
export LSB_JOB_REPORT_MAIL="N"
export -f func1
bsub -q "standard"  -n $ncpus -R "span[hosts=1] rusage[mem=12000]" -P microbiome -J "mapAnnot[1-$nSamples]" -oo "$readScratch/eggnog.%I.log" "func1 \$LSB_JOBINDEX '$files' '$assemblyChoice' '$analysisDir' '$ncpus'"
