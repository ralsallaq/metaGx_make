#!/bin/bash
activate=/hpcf/apps/conda2/install/4.1.1/bin/activate
echo "activate metaG env to get bwa"
source $activate metaG
module load samtools/1.9
readScratch=$1
analysisDir=$2
que=$3
ncpus=$4
files=`ls -p $analysisDir/*.filtered.fq.gz`
nSamples=`echo $files|wc -w`
bwa_index_prefix="human_genome_bwa"
echo "Processing $nSamples samples having paired R1 and R2 reads"

function func1 {
local runN=$1
local files=$2
local readScratch=$3
local bwa_index_prefix=$4
local ncpus=$5
file=`echo $files | cut -f $runN -d " "`
analysisDir=`dirname "$file"`
sname=`basename $file .filtered.fq.gz`
echo "processing reads from sample $sname in file $file"
echo "First align reads to indexed human genome"
echo "First align reads to indexed human genome with bwa using maximum exact match (MEM)"
#the option -p is for interleaved paired files
echo "bwa idxbase=$bwa_index_prefix"
bwa mem -t $ncpus -o $readScratch/$sname.alignment.sam -p $bwa_index_prefix $file 

echo "Second: convert sam to bam"
#-b :output in the bam format -h include the header in the output
if [ -f $readScratch/$sname.alignment.sam ]; then
    samtools view -bh -o $readScratch/$sname.alignment.bam $readScratch/$sname.alignment.sam 
fi

echo "Third extract unaligned pairs (i.e. non-human sequences)"
#samtools fastq  --threads $ncpus  -f 12 $analysisDir/$sname.alignment.bam  > $readScratch/$sname.decontaminated.fq
if [ -f $readScratch/$sname.alignment.bam ]; then
    # option -f 4  sometimes failed to output pairs
    samtools fastq  --threads $ncpus  -f 12 $readScratch/$sname.alignment.bam  > $analysisDir/$sname.decontaminated.fq
fi

#echo "create symbolic links to files in analysisDir for the next step"
#if [ -f $analysisDir/$sname.decontaminated.fq ]; then
#    unlink $readScratch/$sname.temp.fq; ln -s $analysisDir/$sname.decontaminated.fq $readScratch/$sname.temp.fq 
#fi
}
export LSB_JOB_REPORT_MAIL="N"
export -f func1
bsub -q $que -n $ncpus -R "span[hosts=1] rusage[mem=10000]" -P microbiome -J "decontaminate[1-$nSamples]" -oo "$readScratch/decontaminate%I.log" "func1 \$LSB_JOBINDEX '$files' '$readScratch' '$bwa_index_prefix' '$ncpus'"
