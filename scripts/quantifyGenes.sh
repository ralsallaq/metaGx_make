#!/bin/bash
module purge
activate=/hpcf/apps/conda2/install/4.1.1/bin/activate
echo "activate htseq env to get htseq"
source $activate htseq
module load samtools/1.9
module load picard/2.9.4
module load bedtools/2.25.0 
module load python/2.7.15-rhel7 #has pandas
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

echo "cd to the directory having the mapped reads"
#cd $analysisDir/prokka_Annot.$sname/readsOntoCDS.map/cdsMappedReads
cd $analysisDir/prokka_Annot.$sname/readsOntoCDS.map/
echo "now in $PWD"
wait
echo "sort the reads in the bam file either by read name or by leftmost alignment coordinate; we will sort by the coordinate which is the default"
if [ ! -f cdsMappedReads.sorted.bam ]; then
    samtools sort -@"$ncpus" cdsMappedReads.bam -o cdsMappedReads.sorted.bam
fi
wait
#echo "Then index the sorted bam file"
#samtools index map_pe.sorted.bam
#wait
echo "then we will use picard tools to remove duplicates in the sorted and indexed bam file"
if [ ! -f cdsMappedReads.markedup.bam ]; then
    #the MAX_FILE_HANDLES_FOR_READ_ENDS_MAP is found by the unix command ulimit -n
    java -jar /hpcf/apps/picard/install/2.9.4/picard.jar MarkDuplicates INPUT=cdsMappedReads.sorted.bam OUTPUT=cdsMappedReads.markedup.bam METRICS_FILE=cdsMappedReads.markedup.metrics AS=TRUE VALIDATION_STRINGENCY=LENIENT MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=60000 REMOVE_DUPLICATES=TRUE
fi
wait
#######BEDTOOLS method ##########
echo "Conver the index file (fai file) for prokka.faa to bed using awk, the bed file has the CDS id followed by the start basepair number set here to zero and the end bp number (length)"
rm prokka.bed
if [ ! -f prokka.bed ]; then
    awk 'BEGIN {FS="\t"};{print $1 FS "0" FS $2}' ../../../$analysisDir/prokka_Annot.$sname/prokka.faa.fai > prokka.bed
fi

echo "use bedtools to quantify the CDS genes two outputs are of interest (number of reads covering the gene =depth) and the percentage of the gene covered at that depth"
rm gene_coverage.hist
if [ ! -f gene_coverage.hist ]; then
    #bedtools coverage -abam cdsMappedReads.markedup.bam -b prokka.bed -hist > gene_coverage.hist
    echo "the following is from https://github.com/arq5x/bedtools-protocols/blob/master/bedtools.md#ap2b-assessing-coverage-in-exome-capture-experiments-"
    echo "it will give the following columns: ['cds','depth','num_bases_atDepth','sizeOfGene','percOfCDSCoveredAtdepth']"
    bedtools genomecov -ibam cdsMappedReads.markedup.bam  > gene_coverage.hist
fi

echo "Longer genes will of-course have more mapped reads than shorter genes, to facilitate comparison across samples we need to account for the fact that the total number of reads sequenced (the sequencing depth) may differ significantly between samples/n there are several ways to normalize abundance values in metagenomes/n here we will use Transcripts Per Million method (TPM)/n we need for this: 1. The average read length of the sample 2.The length of all genes"
#echo "we get the gene length from our produced gtf filei; here we extract only the start, stop and gene name fields from the file, then remove the ‘gene_id’ string, print the gene name first followed by the length of the gene, change the separator to tab and store the results in the .genelengths file"
#if [ ! -f genelengths ]; then
#    grep '^gnl' $analysisDir/prokka_Annot.$sname/prokka.gff > prokka_CDS.gtf
#    cut -f4,5,9 prokka_CDS.gtf | sed 's/gene_id //g' | gawk '{print $3,$2-$1+1}' | tr ' ' '\t' > genelengths
#fi
wait
echo "we get the read counts for each gene using bedtools"
rm -f tpm_coverage
if [ ! -f tpm_coverage ]; then
    bedtools coverage -abam cdsMappedReads.markedup.bam -b prokka.bed -counts > bedToolCovOut
    cat bedToolCovOut |awk '{print $1"\t"$4}' > gene_read.count
    echo "we get the gene length from our produced gtf filei; here we extract only the start, stop and gene name fields from the file, then remove the ‘gene_id’ string, print the gene name first followed by the length of the gene, change the separator to tab and store the results in the .genelengths file"
    grep '^gnl' ../../../$analysisDir/prokka_Annot.$sname/prokka.gff > prokka_CDS.gtf
    cat prokka_CDS.gtf|cut -f4,5,9|sed 's/ID=//g;s/;.*$//g'|awk '{print $3,$2-$1+1}' | tr ' ' '\t' > genelengths
    #echo "we get the gene length from our produced bedToolCovOut" #<==this gives the lengths of CDS in aa space not nucleotide space
    #cat bedToolCovOut |awk '{print $1"\t"$3-$2}' > genelengths 
    
    echo "Now we can calculate TPM using the tpm_table.py script from https://github.com/EnvGen/metagenomics-workshop/blob/master/in-house/tpm_table.py"
    python ../../../scripts/tpm_table.py -n $sname -c gene_read.count -i <(echo -e "$sname\t100") -l genelengths > tpm_coverage
fi

cd ../../../
echo "now in $PWD"

echo "combine annotations of eggnog with coverage info to get coverage info for pathways/enzymes"
grep -v '^#' $analysisDir/eggnog_map.$sname/egmap.emapper.annotations |awk 'BEGIN{FS="\t"}{print $0}'>$analysisDir/eggnog_map.$sname/eggnog_map.$sname.annot.reformat

rm $analysisDir/eggnog_map.$sname/egmap.emapper.annotations.Wcov.csv
if [ ! -f $analysisDir/eggnog_map.$sname/egmap.emapper.annotations.Wcov.csv ]; then
    reformatEgg="$analysisDir/eggnog_map.$sname/eggnog_map.$sname.annot.reformat"
    covF="$analysisDir/prokka_Annot.$sname/readsOntoCDS.map/gene_coverage.hist"
    python -c "import pandas as pd
df1 = pd.read_csv('$reformatEgg', header=None, sep='\t')
#the 8th column has something like this ko:K03148,ko:K21029, these are the KEGG Orthologs see for example the following link for enzyme ec:2.3.1.15 https://www.genome.jp/dbget-bin/www_bget?ec:2.3.1.15 and look under Orthology to see KOs
#the 9th column has something like this ko00730,ko01100,ko04122,map00730,map01100,map0 these found in the link https://www.genome.jp/dbget-bin/www_bget?ec:2.3.1.15 for the enzyme with ec 2.3.1.15 under the Pathway so map00561 for example is ec00561 pathway
#
df1.columns=['cds','seedOrth','evalue','score','taxG','proteinName','geneOntology','ec','kEGG_ko','kEGG_PWY','kEGG_module','kEGG_rxn','kEGG_rclass','bRITE','kEGG_TC','CAZy','biGG_rxn','tx_scope','OGs','bestOG','cOG funcat','description']
df2 = pd.read_csv('$covF', header=None, sep='\t')
#df2.columns=['cds','indx_1stbp_protein','indx_lastbp_protein','depth','num_bases_atDepth','sizeOgGene','percOfCDSCoveredAtdepth'] 
df2.columns=['cds','depth','num_bases_atDepth','sizeOfGene','percOfCDSCoveredAtdepth']
print(df1.head(),'\n',df2.head())
merged=df2.merge(df1, on='cds' , how='left')
print(merged.head())
merged.to_csv('$analysisDir/eggnog_map.$sname/egmap.emapper.annotations.Wcov.csv')
            "
else
    echo "file $analysisDir/eggnog_map.$sname/egmap.emapper.annotations.Wcov.csv already exists, nothing to do"
fi


echo "done"

}
export LSB_JOB_REPORT_MAIL="N"
export -f func1
bsub -q $que  -n $ncpus -R "span[hosts=1] rusage[mem=10000]" -P microbiome -J "quanGenes[1-$nSamples]" -oo "$readScratch/quantGenes.%I.log" "func1 \$LSB_JOBINDEX '$files' '$assemblyChoice' '$analysisDir' '$ncpus'"
