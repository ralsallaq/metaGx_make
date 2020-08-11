#!/bin/bash
activate=/hpcf/apps/conda2/install/4.1.1/bin/activate
echo "activate metaG env to get bwa"
source $activate metaG
module load samtools/1.9
readScratch=$1
analysisDir=$2
que=$3
ncpus=$4
dbSelect=$5
dbIndexDir=$6
geneFractionThreshold=$7 #integer between 1-100 (use 20 for low depth and 80 for good depth)
echo $1 $2 $3 $4 $5 $6 $7

files=`ls -p $analysisDir/*.qtrimmed.fq.gz`
nSamples=`echo $files|wc -w`
echo "Processing $nSamples samples having paired R1 and R2 reads"

function func1 {
local runN=$1
local files=$2
local readScratch=$3
local dbSelect=$4
local ncpus=$5
local dbIndexDir=$6
local geneFractionThreshold=$7 
file=`echo $files | cut -f $runN -d " "`
analysisDir=`dirname "$file"`
sname=`basename $file .qtrimmed.fq.gz`
echo $runN $readScratch $dbSelect $ncpus $dbIndexDir 
echo $sname
if [ ! -d AMR ]; then
    mkdir AMR
fi
if [[ "$dbSelect" == "drugs" ]]; then
    echo "here drugs"
    ln -s $dbIndexDir/megaresV2_$dbSelect.* ./ 
    amrDB="/research/home/margogrp/ralsalla/megares_v2.00/megares_drugs_database_v2.00.fasta"
    bwa_index_prefix="megaresV2_$dbSelect"
    amr_annotations="/research/home/margogrp/ralsalla/megares_v2.00/megares_"$dbSelect"_annotations_v2.00.csv"
fi
if [[ "$dbSelect" == "modified" ]]; then
    echo "here modified"
    ln -s $dbIndexDir/megaresV2_$dbSelect.* ./ 
    amrDB="/research/home/margogrp/ralsalla/megares_v2.00/megares_modified_database_v2.00.fasta"
    bwa_index_prefix="megaresV2_$dbSelect"
    amr_annotations="/research/home/margogrp/ralsalla/megares_v2.00/megares_"$dbSelect"_annotations_v2.00.csv"
fi
echo "define some parameters"
threshold="$geneFractionThreshold"
rarefactionMin=5
rarefactionMax=100
rarefactionSkip=5
numIterationsForSampling=1


echo "processing reads from sample $sname in file $file"
echo "First align reads to indexed MegaRes DB"
echo "First align reads to indexed MegaRes DB with bwa using maximum exact match (MEM)"
#the option -p is for interleaved paired files
echo "bwa idxbase=$bwa_index_prefix"
if [ ! -f $readScratch/$sname.megaResalgn.sam ]; then
    bwa mem -t $ncpus -o $readScratch/$sname.megaResalgn.sam -p $bwa_index_prefix $file 
else
    echo "the file $readScratch/$sname.megaResalgn.sam already exists; nothing to do"
fi

echo "next will run resistome analysis, default threshold here is 80"
echo "ref=$amrDB; annot_fp=$amr_annotations; sam_fp=$readScratch/$sname.megaResalgn.sam"
resistome \
    -ref_fp $amrDB \
    -annot_fp $amr_annotations \
    -sam_fp $readScratch/$sname.megaResalgn.sam \
    -gene_fp AMR/$sname.gene.tsv \
    -group_fp AMR/$sname.group.tsv \
    -class_fp AMR/$sname.class.tsv \
    -mech_fp AMR/$sname.mechanism.tsv \
    -t $threshold
#wait
echo "next will run rarefaction with min rarefaction level=5 and max rarefaction level=100 to get a distribution"
echo "ref=$amrDB; annot_fp=$amr_annotations; sam_fp=$readScratch/$sname.megaResalgn.sam"
rarefaction \
    -ref_fp $amrDB \
    -sam_fp $readScratch/$sname.megaResalgn.sam \
    -annot_fp $amr_annotations \
    -gene_fp AMR/$sname.gene.rareFD.tsv \
    -group_fp AMR/$sname.group.rareFD.tsv \
    -class_fp AMR/$sname.class.rareFD.tsv \
    -mech_fp AMR/$sname.mech.rareFD.tsv \
    -min $rarefactionMin \
    -max $rarefactionMax \
    -skip $rarefactionSkip \
    -samples $numIterationsForSampling \
    -t $threshold

wait
echo "next will run SNPFinder"
echo "this results into a tsv file per sample with three columns: Gene, Haplotype Pattern, and Occurrence"
snpfinder \
    -amr_fp $amrDB \
    -sampe $readScratch/$sname.megaResalgn.sam \
    -out_fp AMR/$sname.tsv
}
export LSB_JOB_REPORT_MAIL="N"
export -f func1
#func1 1 "$files" "$readScratch" "$dbSelect" "$ncpus" "$dbIndexDir"
bsub -q $que -n $ncpus -R "span[hosts=1] rusage[mem=10000]" -P microbiome -J "runResistome[1-$nSamples]" -oo "$readScratch/resistome_%I.log" "func1 \$LSB_JOBINDEX '$files' '$readScratch' '$dbSelect' '$ncpus' '$dbIndexDir' '$geneFractionThreshold'"
