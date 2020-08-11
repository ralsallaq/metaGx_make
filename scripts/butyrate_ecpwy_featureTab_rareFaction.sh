#!/bin/bash
module load samtools/1.10
module load python/3.5.2
module load bedtools/2.25.0 
readScratch=$1
analysisDir=$2
assemblyChoice=$3
ecpwyF=$4
echo $1 $2 $3 $4 $5 $6 $7
if [[ ! $5 -le 1 ]]; then
    threshold=`bc <<< "scale=2; $5/100"`
    echo "threshold is $threshold"
else
    echo "threshold must be a percentage e.g. 10, 30, 75 ,..etc"
    exit 1
fi
que=$6
ncpus=1 #or $7
echo "rarefaction analysis of predicted/annotated genes for existence of butyrate pathways"
bamfiles=`find $analysisDir/prokka_Annot.*/* -name cdsMappedReads.markedup.bam`
snames=`echo $bamfiles| tr " " "\n"|sed 's/.*prokka_Annot\.//g;s/\/readsOntoCDS.*$//g'`
nSamples=`echo $snames|wc -w`
echo "Processing $nSamples samples having paired R1 and R2 reads"

function func1 {
local runN=$1
local snames=$2
local analysisDir=$3
local ecpwyF=$4
local threshold=$5

sname=`echo $snames | cut -f $runN -d " "`
echo "processing $sname"
#alignment with all reads aligned
bamfile=`find "$analysisDir"/prokka_Annot."$sname"/* -name cdsMappedReads.markedup.bam`
reformatEgg=`find "$analysisDir"/eggnog_map."$sname"/* -name eggnog_map."$sname".annot.reformat`

echo $bamfile $reformatEgg
tempRAM=$(mktemp -d $readScratch/$sname.XXXXXXXXX)

min=5
max=95
skip=5
if [[ $max -eq 100 ]]; then
    max=95
    append=99.99999
else
    append=""
fi
seq_frac=`seq $min $skip $max`
subsamplingFrac=`echo "$seq_frac" "$append"` 
outF="$analysisDir"/"$sname"_butyrateFT_rareFaction.tsv
rm -f $outF
for f in $subsamplingFrac; do
    frac=`bc <<< "scale=2; $f/100"`
    samtools view -h -s "1$frac" -b "$bamfile" > $tempRAM/sub_"$f".bam 
    bedtools genomecov -ibam $tempRAM/sub_"$f".bam  > $tempRAM/gene_coverage_sub_"$f".hist
    covF=""$tempRAM"/gene_coverage_sub_"$f".hist"
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
merged.to_csv('$tempRAM/egmap.emapper.annotations_sub_$f.Wcov.csv')
            "

    python ./scripts/butyrate_ecpwy_featureTab.py -o "$tempRAM"/sub_"$f"_featureTab.csv  -i "$tempRAM"/egmap.emapper.annotations_sub_"$f".Wcov.csv -ec "$ecpwyF" -t "$threshold"

    distinctFeatures=`echo $(cat "$tempRAM"/sub_"$f"_featureTab.csv|wc -l)-1|bc`

    echo "$f $distinctFeatures"|column -t >>$outF 
done

}

export LSB_JOB_REPORT_MAIL="N"
export -f func1
bsub -q $que  -n $ncpus -R "span[hosts=1] rusage[mem=6000]" -P microbiome -J "butyrateRareF[1-$nSamples]" -oo "$readScratch/butyrateRareF.%I.log" "func1 \$LSB_JOBINDEX '$snames' '$analysisDir' '$ecpwyF' '$threshold'"
