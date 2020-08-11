#!/bin/bash
activate=/hpcf/apps/conda2/install/4.1.1/bin/activate
echo "activate metaG env to get quast"
source $activate metaG
readScratch=$1
analysisDir=$2
assemblyChoice=$3
que=$4
ncpus=$5
echo "determining coverage and taxonomic make up for $assemblyChoice assembly"
files=`ls -p $readScratch/*.temp.fq.gz`
nSamples=`echo $files|wc -w`
echo "Processing $nSamples samples having paired R1 and R2 reads"

function func1 {
local runN=$1
local files=$2
local assemblyChoice=$3
local analysisDir=$4
local ncpus=$5
file=`echo $files | cut -f $runN -d " "`
scratch=`dirname "$file"`
sname=`basename $file .temp.fq.gz`
if [ $assemblyChoice=="spades" ]; then
    assmblyDir=$analysisDir/$sname.spades_out/scaffolds.fasta
elif [ $assemblyChoice=="megahit" ]; then
    assmblyDir=$analysisDir/$sname.megahit_out/contigs.fa
fi

echo "determining taxonomic makeup using sketch tables I am saving the top 50 hits"
sendsketch.sh in=$assmblyDir address=silva records=50 out=$analysisDir/taxonomic_makeup.$sname

wait
echo "Calculate the coverage distribution, and capture reads that did not make it into the assembly"
bbmap.sh in=$analysisDir/$sname.noAdpt.fq.gz ref=$assmblyDir nodisk covhist=$analysisDir/$sname.$assemblyChoice.covhist.txt covstats=$analysisDir/$sname.$assemblyChoice.covstats.txt outm=$analysisDir/$sname.$assemblyChoice.assembled.fq.gz outu=$analysisDir/$sname.$assemblyChoice.unassembled.fq.gz maxindel=130 minid=90 qtrim=10 untrim ambig=all 
}
export LSB_JOB_REPORT_MAIL="N"
export -f func1
bsub -q $que  -n $ncpus -R "span[hosts=1] rusage[mem=11000]" -P microbiome -J "postAssembly[1-$nSamples]" -oo "$readScratch/postAssembly.%I.log" "func1 \$LSB_JOBINDEX '$files' '$assemblyChoice' '$analysisDir' '$ncpus'"
