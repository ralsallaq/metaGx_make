#!/bin/bash
activate=/hpcf/apps/conda2/install/4.1.1/bin/activate
echo "activate metaG env to get spades/metaspades"
source $activate metaG
module load samtools/1.9
readScratch=$1
analysisDir=$2
que=$3
ncpus=$4
files=`ls -p $analysisDir/*.qtrimmed.fq.gz`
nSamples=`echo $files|wc -w`
echo "Processing $nSamples samples having paired R1 and R2 reads"

function func1 {
local runN=$1
local files=$2
local analysisDir=$3
local ncpus=$4
file=`echo $files | cut -f $runN -d " "`
scratch=`dirname "$file"`
sname=`basename $file .qtrimmed.fq.gz`
echo "processing reads from sample $sname in file $file"
echo "assemble the reads into contigs using metaspades"
#--12 is for file with interlaced forward and reverse paired-end reads, option --merged (not used here) is for merged paired-end reads 
#spades.py --meta -t $ncpus --phred-offset 33 -k 25,55,95,125 -s $analysisDir/$sname.merged.fq.gz --12 $analysisDir/$sname.qtrimmed.fq -o $analysisDir/$sname.spades_out 
#this will run spades assembler with error correction (default)
#spades.py --meta -t $ncpus --phred-offset 33 -k 25,55,95,125 --merged $analysisDir/$sname.merged.fq.gz --12 $analysisDir/$sname.qtrimmed.fq -o $analysisDir/$sname.spades_out 
#this will run spades assembler without error correction
#spades.py --meta --only-assembler -t $ncpus --phred-offset 33 -k 25,55,95,125 --12 $analysisDir/$sname.qtrimmed.fq -o $analysisDir/$sname.spades_out 

#echo "allowing the error correction algorithm and the assembler from spades"
echo "only the assembler from spades"
if [ ! -f $analysisDir/$sname.spades_out/scaffolds.fasta ]; then
spades.py --meta --only-assembler -t $ncpus --phred-offset 33 -k 25,55,95  --12 $file -o $analysisDir/$sname.spades_out 
fi

}
export LSB_JOB_REPORT_MAIL="N"
export -f func1
bsub -q $que  -n $ncpus -R "span[hosts=1] rusage[mem=12000]" -P microbiome -J "assemble_spades[1-$nSamples]" -oo "$readScratch/spades.%I.log" "func1 \$LSB_JOBINDEX '$files' '$analysisDir' '$ncpus'"
