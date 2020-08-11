#!/bin/bash
activate=/hpcf/apps/conda2/install/4.1.1/bin/activate
echo "activate metaG env to get bwa"
source $activate metaG
module load samtools/1.9
readScratch=$1
analysisDir=$2
assemblyChoice=$3
que=$4
ncpus=$5
echo "mapping reads onto $assemblyChoice assembly and get BAM files for sample reads"
files=`ls -p $analysisDir/*.qtrimmed.fq.gz`
nSamples=`echo $files|wc -w`
echo "Processing $nSamples samples having paired R1 and R2 reads"

function func1 {
local runN=$1
local files=$2
local assemblyChoice=$3
local analysisDir=$4
local ncpus=$5
file=`echo $files | cut -f $runN -d " "`
#scratch=`dirname "$file"`
sname=`basename $file .qtrimmed.fq.gz`
echo "Processing sample $sname"
if [[ $assemblyChoice = "spades" ]]; then
    assmblyFs="$analysisDir/$sname.spades_out/scaffolds.fasta"
elif [[ $assemblyChoice = "megahit" ]]; then
    assmblyFs="$analysisDir/$sname.megahit_out/contigs.fa"
fi
assemblyDir=`dirname $assmblyFs`

echo " The following steps were from https://inf-biox121.readthedocs.io/en/2016/Assembly/practicals/03_Mapping_reads_to_an_assembly.html"

echo "cd to the assembly directory and First index the assembly fasta file"
cd $assemblyDir
echo "now in $PWD"
# -p is a prefix we do not needed it as we are inside a directory
#bwa index -a bwtsw -p $sname.bwaIndex $assmblyFs
if [ ! -f scaffolds.fasta.ann ];then
    bwa index -a bwtsw  $assmblyFs
fi
wait
echo "Second align qtrimmed reads onto the indexed assembly using maximum exact match (MEM)"
if [ ! -f bwa/align.sam ]; then
    mkdir bwa
    wait
    cd bwa
    echo "now in $PWD"
    wait
    #the option -p is for interleaved paired files
    #-p here for interleaved fq files
    #the - for both samtools commands indicate that instead of using a file as input, the input comes from a pipe (technically, from ‘standard in’, or ‘STDIN’).
    #bwa mem -t $ncpus ../$assmblyFs -p  $file | samtools -view -buS - | samtools sort - -o map_pe.sorted.bam 
    #here I am following http://merenlab.org/tutorials/assembly-based-metagenomics/
    echo "the assembly fasta is $assmblyFs"
    
    #bwa mem -t $ncpus $assmblyFs -p  $file | samtools -view -F 4 -buS - | samtools sort - -o map_pe.sorted.bam 
    
    bwa mem -t $ncpus $assmblyFs -p  $file -o align.sam 
fi

#wait
#echo "extract the mapped reads from the alignment in bam format"
#samtools view -buS -F 4 -o map_pe.bam align.sam
#wait
#echo "sort the reads in the bam file"
#samtools sort map_pe.bam -o map_pe.sorted.bam
#wait
#echo "Then index the sorted bam file"
#samtools index map_pe.sorted.bam
#
##echo "Fourth extract mapped and unmapped pairs"
####samtools fastq  --threads $ncpus  -f 12 $scratch/$sname.alignment.bam  > $analysisDir/$sname.decontaminated.fq
##samtools fastq  --threads $ncpus  -f 4 map_pe.sorted.bam  > unmapped.fq
##wait
##samtools fastq  --threads $ncpus  -F 4 map_pe.sorted.bam  > mapped.fq

cd ../../
echo "now in $PWD"
echo "done"

}
export LSB_JOB_REPORT_MAIL="N"
export -f func1
bsub -q $que  -n $ncpus -R "span[hosts=1] rusage[mem=10000]" -P microbiome -J "mapReads[1-$nSamples]" -oo "$readScratch/mapReads.%I.log" "func1 \$LSB_JOBINDEX '$files' '$assemblyChoice' '$analysisDir' '$ncpus'"
