#!/bin/bash
activate=/hpcf/apps/conda2/install/4.1.1/bin/activate
echo "activate metaG env to get diamond"
source $activate metaG
module load seqtk/1.2-r101c
module load samtools/1.9
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
file=`basename $(echo $files | cut -f $runN -d " ")`
sname=`basename $file .qtrimmed.fq.gz`
echo "processing $file from sample $sname"
if [ $assemblyChoice=="spades" ]; then
    assmblyFs=$analysisDir/$sname.spades_out/scaffolds.fasta
elif [ $assemblyChoice=="megahit" ]; then
    assmblyFs=$analysisDir/$sname.megahit_out/contigs.fa
fi

echo "the assembly file is $assmblyFs and it has the following number of contigs"
echo `grep -c "^>" $assmblyFs`

echo "prokka predicted from these contigs the following number of CDS genes"
echo `grep -c "^>" $analysisDir/prokka_Annot.$sname/prokka.faa`

echo "map qtrimmed reads using diamond algorithm onto the CDS genes"
echo "diamond is much faster than blasx (match DNA sequence against a protein database fasta)"
if [ ! -d $analysisDir/prokka_Annot.$sname/readsOntoCDS.map ]; then  
    mkdir $analysisDir/prokka_Annot.$sname/readsOntoCDS.map
fi
##blasting qtrimmed files against prokka predicted CDS genes
#echo "first convert the fastq file to fasta" #no need diamond can take fastq
#seqtk seq -a $file > $scratch/$sname.qtrimmed.fa 
echo "build a diamond database for prokka"
echo "cd to prokka"
cd $analysisDir/prokka_Annot.$sname/
echo "now inside $PWD"
if [ ! -f readsOntoCDS.map/cdsMappedReads.bam ]; then
    echo "build diamond database prokka.dmnd out of prokka.faa file"
    diamond makedb --in prokka.faa -d prokka
    wait
    echo "construct a tab-delimited file from prokka fasta for use with samtools to convert sam to bam, this file is called prokka.faa.fai"
    samtools faidx prokka.faa
    wait
    echo "then blastx of qtrimmed reads against prokka.faa file which contains predicted CDS. The aligner is diamond. This mapping is to eventually quantify the CDS genes"
    
    #diamond blastx --threads $ncpus --query $file --db $analysisDir/prokka_Annot.$sname/prokka.dmnd  --daa $analysisDir/prokka_Annot.$sname/readsOntoCDS.map/cdsMappedReads 
    echo "--outfmt 101 will produce a sam file for the alignment"
    
    #diamond blastx --threads $ncpus --query $scratch/$sname.qtrimmed.fa --db prokka.dmnd  --outfmt 101  -o readsOntoCDS.map/cdsMappedReads.sam 
    diamond blastx --threads $ncpus --query ../$file --db prokka.dmnd  --outfmt 101  -o readsOntoCDS.map/cdsMappedReads.sam 
    
    echo "convert sam to bam"
    samtools view -bh -t prokka.faa.fai -o readsOntoCDS.map/cdsMappedReads.bam  readsOntoCDS.map/cdsMappedReads.sam
fi
echo "exiting the directory"
cd ../../ 
}
export LSB_JOB_REPORT_MAIL="N"
export -f func1
bsub -q $que  -n $ncpus -R "span[hosts=1] rusage[mem=10000]" -P microbiome -J "mapReads[1-$nSamples]" -oo "$readScratch/mapReadsOnCDS.%I.log" "func1 \$LSB_JOBINDEX '$files' '$assemblyChoice' '$analysisDir' '$ncpus' '$readScratch'"
