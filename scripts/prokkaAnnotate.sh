#!/bin/bash
activate=/hpcf/apps/conda2/install/4.1.1/bin/activate
echo "activate prokka env to get prokka"
source $activate prokka
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
file=`echo $files | cut -f $runN -d " "`
scratch=`dirname "$file"`
sname=`basename $file .reads.fq.gz`
if [ $assemblyChoice=="spades" ]; then
    assmblyDir=$analysisDir/$sname.spades_out/scaffolds.fasta
elif [ $assemblyChoice=="megahit" ]; then
    assmblyDir=$analysisDir/$sname.megahit_out/contigs.fa
fi

echo "annotating assembly of choice ($assemblyChoice) using prokka"
echo "annotation means: finding the genes and doing functional annotation of those"
echo " from  https://metagenomics-workshop.readthedocs.io/en/latest/annotation/functional_annotation.html  -- PROKKA automates the process of locating open reading frames (ORFs) and RNA regions on contigs, translating ORFs to protein sequences, searching for protein homologs and producing standard output files. For gene finding and translation, PROKKA makes use of the program Prodigal. Homology searching (via BLAST and HMMER) is then performed using the translated protein sequences as queries against a set of public databases (CDD, PFAM, TIGRFAM) as well as custom databases that come with PROKKA." 
echo "Outputs: GFF file = a standardized tab-delimited format for genome annotations. GBK = a Genbank file encompassing more detailed description of nucleoatide sequences and the genes encoded in these"
echo "For description of the columns (following all the lines with ## and before fasta formatted sequences) see http://gmod.org/wiki/GFF3; in order to skip the ## lines (which represent the annotated sequence regions) do for example grep -v "^#" analysis/prokka_Annot.1729701_EW1H-0101_S1/prokka.gff|less"
mkdir $analysisDir/prokka_Annot.$sname
#prokka --norrna --notrna --metagenome --outdir $analysisDir/prokka_Annot.$sname --force --centre metaG --prefix prokka --cpus $ncpus $assmblyDir
#prokka --outdir $analysisDir/prokka_Annot.$sname --force --centre metaG --prefix prokka --cpus $ncpus $assmblyDir
if [ ! -f $analysisDir/prokka_Annot.$sname/prokka.faa ]; then
    prokka --outdir $analysisDir/prokka_Annot.$sname --force --centre metaG --norna --notrna --metagenome --prefix prokka --cpus $ncpus $assmblyDir
else
    echo "The file $analysisDir/prokka_Annot.$sname/prokka.faa already exists, nothing to do"
fi

}
export LSB_JOB_REPORT_MAIL="N"
export -f func1
bsub -q "standard"  -n $ncpus -R "span[hosts=1] rusage[mem=12000]" -P microbiome -J "annotate[1-$nSamples]" -oo "$readScratch/prokka.%I.log" "func1 \$LSB_JOBINDEX '$files' '$assemblyChoice' '$analysisDir' '$ncpus'"
