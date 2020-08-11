#!/bin/bash
activate=/hpcf/apps/conda2/install/4.1.1/bin/activate
echo "activate metaG env to get diamond"
source $activate metaG
module load seqtk/1.2-r101c
readScratch=$1
analysisDir=$2
assemblyChoice=$3
echo "mapping prokka inferred CDS to uniref90"
files=`ls -p $analysisDir/*.qtrimmed.fq.gz`
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
local scratch=$6
uniref90DB=/research/home/margogrp/ralsalla/metaG/databases/uniref90.9.2019.dmnd
file=`echo $files | cut -f $runN -d " "`
#scratch=`dirname "$file"`
sname=`basename $file .qtrimmed.fq.gz`
if [ $assemblyChoice=="spades" ]; then
    assmblyFs=$analysisDir/$sname.spades_out/scaffolds.fasta
elif [ $assemblyChoice=="megahit" ]; then
    assmblyFs=$analysisDir/$sname.megahit_out/contigs.fa
fi
echo "the assembly file is $assmblyFs and it has the following number of contigs"
echo `grep -c "^>" $assmblyFs`

echo "map qtrimmed reads using diamond algorithm to uniref"
echo "diamond is much faster than blasx (match DNA sequence against a protein database fasta)"
if [ ! -d $analysisDir/uniref90_map.$sname ]; then  
    mkdir $analysisDir/uniref90_map.$sname
fi
##blasting qtrimmed files against uniref90
#echo "first convert the fastq file to fasta"
#seqtk seq -a $file > $scratch/$sname.qtrimmed.fa 
#echo "then blastx of qtrimmed reads against uniref90 using diamond"
#diamond blastx --threads $ncpus --query $scratch/$sname.qtrimmed.fa --db $uniref90DB  --daa $analysisDir/uniref90_map.$sname/unirefmapReads

#blasting contigs against uniref90
echo "using blastx of contigs in $assmblyFs against uniref90 using diamond" 
diamond blastx --threads $ncpus --query $assmblyFs --db $uniref90DB  --daa $analysisDir/uniref90_map.$sname/unirefmapContigs
}
export LSB_JOB_REPORT_MAIL="N"
export -f func1
bsub -q $que  -n $ncpus -R "span[hosts=1] rusage[mem=10000]" -P microbiome -J "mapReads[1-$nSamples]" -oo "$readScratch/unirefReads.%I.log" "func1 \$LSB_JOBINDEX '$files' '$assemblyChoice' '$analysisDir' '$ncpus' '$readScratch'"
