#!/bin/bash
#!/bin/bash
activate=/hpcf/apps/conda2/install/4.1.1/bin/activate
echo "activate cutadapt env to get cutadapt"
source $activate cutadapt
readScratch=$1
analysisDir=$2
que=$3
ncpus=$4
files=`ls -p $analysisDir/*.decontaminated.fq`
nSamples=`echo $files|wc -w`
echo "Processing $nSamples samples having paired R1 and R2 reads"
echo "quality trim reads"

function func1 {
local runN=$1
local files=$2
local readScratch=$3
local ncpus=$4
file=`echo $files | cut -f $runN -d " "`
analysisDir=`dirname "$file"`
sname=`basename $file .decontaminated.fq`
echo "processing reads from sample $sname in file $file"
echo "quality trim reads with quality less than 10 and length less than 50"
#echo "remove adapters and discard any reads with quality below 10 or with length <50 base-pairs (the pair discarded if *ANY* of them is too short or its quality is bad)"
#cutadapt --interleaved -j $ncpus -q 10 -m 50 -pair-filter=any -A XXX -o $analysisDir/$sname.qtrimmed.fq.gz $file 

echo "remove adapters and discard any pairs with quality below 10 or with length <50 base-pairs (the pair discarded if *ANY* of them is too short or its quality is bad)"
echo "the options -pair-filter=any  and -pair-filter=both will return only R1 reads so do not use them"
if [ ! -f $analysisDir/$sname.qtrimmed.fq.gz ]; then
    cutadapt --interleaved -j $ncpus -q 10 -m 50 -A XXX -o $analysisDir/$sname.qtrimmed.fq.gz $file 
else
    echo "file $analysisDir/$sname.qtrimmed.fq.gz already exists nothing to do"
fi

#echo "create symbolic link in scratch"
#if [ -f $analysisDir/$sname.qtrimmed.fq.gz ]; then
#    unlink $readScratch/$sname.temp.fq.gz; ln -s $analysisDir/$sname.qtrimmed.fq.gz $readScratch/$sname.temp.fq.gz 
#fi

}
export LSB_JOB_REPORT_MAIL="N"
export -f func1
bsub -q $que  -n $ncpus -R "span[hosts=1] rusage[mem=10000]" -P microbiome -J "qtrimmed[1-$nSamples]" -oo "$readScratch/qtrimmed.%I.log" "func1 \$LSB_JOBINDEX '$files' '$readScratch' '$ncpus'"
