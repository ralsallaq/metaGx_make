#!/bin/bash
module load samtools/1.9
readScratch=$1
analysisDir=$2
que=$3
ncpus=$4
krakenDBDir=$5 ##/scratch_space/ralsalla/databases/kraken2_db/ (using kraken2)
if [[ "$krakenDBDir" != "/scratch_space/ralsalla/databases/kraken2_db/" ]]; then
    echo "to use kraken2 please select the correct database"
    exit 1
fi

files=`ls -p $analysisDir/*.qtrimmed.fq.gz`
nSamples=`echo $files|wc -w`
echo "Processing $nSamples samples having paired R1 and R2 reads"

function func1 {
local runN=$1
local files=$2
local readScratch=$3
local krakenDBDir=$4
local ncpus=$5
file=`echo $files | cut -f $runN -d " "`
analysisDir=`dirname "$file"`
sname=`basename $file .qtrimmed.fq.gz`

activate=/hpcf/apps/conda2/install/4.1.1/bin/activate
echo "activate metaG env to get bbtools"
source $activate metaG
rm -f $analysisDir/"$sname".qtrimmed.R1.fq.gz $analysisDir/"$sname".qtrimmed.R2.fq.gz
reformat.sh in=$file  out1=$analysisDir/"$sname".qtrimmed.R1.fq.gz out2=$analysisDir/"$sname".qtrimmed.R2.fq.gz
wait
#######kraken1 
#if you want the sequences that were classified or not classified run this
#kraken --preload --threads $ncpus --classified-out $analysisDir/"$sname".qtrimmed.classified --unclassified $analysisDir/"$sname".qtrimmed.unclassified --db $krakenDBDir --fastq-input --gzip-compressed --paired  $analysisDir/"$sname".qtrimmed.R1.fq.gz $analysisDir/"$sname".qtrimmed.R2.fq.gz" > $analysisDir/"$sname".krakenClass
#else for only to get classifications

echo "run kraken classification using kraken2"
########kraken 2
kraken2 --use-names  --threads $ncpus --db $krakenDBDir  --gzip-compressed --paired  $analysisDir/"$sname".qtrimmed.R1.fq.gz $analysisDir/"$sname".qtrimmed.R2.fq.gz --report $analysisDir/"$sname".krakenReport > $analysisDir/"$sname".krakenClass
#if [ ! -f $analysisDir/"$sname".krakenClass ]; then
#    kraken2 --use-names  --threads $ncpus --db $krakenDBDir  --gzip-compressed --paired  $analysisDir/"$sname".qtrimmed.R1.fq.gz $analysisDir/"$sname".qtrimmed.R2.fq.gz --report $analysisDir/"$sname".krakenReport > $analysisDir/"$sname".krakenClass
#fi
rm -f $analysisDir/"$sname".krakenLabels 
if [ ! -f $analysisDir/"$sname".krakenLabels ]; then
    grep '^C' $analysisDir/"$sname".krakenClass | sed 's/^C\t//g' >  $analysisDir/"$sname".krakenLabels
fi
echo "the output has the following columns:
1. 'C'/'U': one letter code indicating that the sequence was either classified or unclassified.
2. The sequence ID, obtained from the FASTA/FASTQ header.
3. The taxonomy ID Kraken used to label the sequence; this is 0 if the sequence is unclassified.
4. The length of the sequence in bp.
5. A space-delimited list indicating the LCA mapping of each k-mer in the sequence. For example, "562:13 561:4 A:31 0:1 562:3" would indicate that:
*the first 13 k-mers mapped to taxonomy ID #562
*the next 4 k-mers mapped to taxonomy ID #561
*the next 31 k-mers contained an ambiguous nucleotide
*the next k-mer was not in the database
*the last 3 k-mers mapped to taxonomy ID #562
"

}
export LSB_JOB_REPORT_MAIL="N"
export -f func1

bsub -q $que  -n $ncpus -R "span[hosts=1] rusage[mem=10000]" -P microbiome -J "classifyKraken[1-$nSamples]" -oo "$readScratch/runKraken_%I.log"  "func1 \$LSB_JOBINDEX '$files' '$readScratch' '$krakenDBDir' '$ncpus'" 