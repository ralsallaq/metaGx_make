#!/bin/bash
activate=/hpcf/apps/conda2/install/4.1.1/bin/activate
echo "activate taxtastic env to get pandas"
source $activate taxtastic
if [ ! -f sampleInfo.csv ]; then
    echo "please prepare sampleInfo.csv file with at least two columns: R1 and R2 encompassing the full paths for reads 1 and 2"
    exit 1
else
    echo "a csv file sampleInfo.csv is written for the read paths using get_sampleInfo.sh script"
    nSamples=`echo $(cat sampleInfo.csv|wc -l)-1|bc`
    echo "Processing $nSamples samples having paired R1 and R2 reads"
fi
echo "Gathering files from sampleInfo.csv file requires python" 
R1F=`echo $(python -c "import pandas as pd; df=pd.read_csv('sampleInfo.csv', index_col=0); print(df['R1'].values)"|sed 's/\[//g;s/\]//g')`
R2F=`echo $(python -c "import pandas as pd; df=pd.read_csv('sampleInfo.csv', index_col=0); print(df['R2'].values)"|sed 's/\[//g;s/\]//g')`
snameS=`echo $(python -c "import pandas as pd; df=pd.read_csv('sampleInfo.csv', index_col=0); print(df['sname'].values)"|sed 's/\[//g;s/\]//g')`

echo $snameS

activate=/hpcf/apps/conda2/install/4.1.1/bin/activate
echo "activate metaG env to get bbtools"
source $activate metaG
#R1F=`cat sampleInfo.csv | awk '{FS=","}{print $2}'`
#R2F=`cat sampleInfo.csv | awk '{FS=","}{print $3}'`
dataDir=$1 
readScratch=$2
que=$3
ncpus=$4

function func1 {
local runN=$1
local R1F=$2
local R2F=$3
local snameS=$4
local readScratch=$5
local ncpus=$6

R1=`echo $R1F | cut -f $runN -d " "|sed "s/'//g"`
R2=`echo $R2F | cut -f $runN -d " "|sed "s/'//g"`
sname=`echo $snameS | cut -f $runN -d " "|sed "s/'//g"`

echo "processing reads $R1 and $R2 from sample $sname"
echo "verify pairing of the reads and write an interleaved file"

if [ ! -f $readScratch/$sname.reads.fq.gz ]; then
    reformat.sh in1=$R1 in2=$R2 vpair out=$readScratch/$sname.reads.fq.gz  
else
    echo "file $readScratch/$sname.reads.fq.gz already exist, nothing to do"
fi

#echo "create symbolic links to these files"
#rm $readScratch/$sname.temp.fq.gz; ln -s $readScratch/$sname.reads.fq.gz $readScratch/$sname.temp.fq.gz 
}
export LSB_JOB_REPORT_MAIL="N"
export -f func1
#func1 3 "$R1F" "$R2F" "$snameS" "$readScratch" "$ncpus"
bsub -q $que  -n 1 -R "span[hosts=1] rusage[mem=12000]" -P microbiome -J "verify[1-$nSamples]" -oo "$readScratch/verifyPairing_%I.log" "func1 \$LSB_JOBINDEX '$R1F' '$R2F' '$snameS' '$readScratch' '$ncpus'"
