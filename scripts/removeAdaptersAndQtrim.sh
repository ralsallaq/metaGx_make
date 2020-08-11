#!/bin/bash
activate=/hpcf/apps/conda2/install/4.1.1/bin/activate
echo "activate cutadapt env to get cutadapt"
source $activate cutadapt
adaptersF=$1
readScratch=$2
analysisDir=$3
que=$4
if [ ! -d $analysisDir ]; then
    mkdir analysis
fi
if [ ! -f $adapters ];then
    echo "the file adapters.fasta is required but missing, this file is produced by inspecting\n
          output from fastqc applied to a number of input files ... stopping"
    exit 1
fi
files=`ls -p $readScratch/*.temp.fq.gz`
nSamples=`echo $files|wc -w`
echo "Processing $nSamples samples having paired R1 and R2 reads"
ncpus=4

function func1 {
local runN=$1
local files=$2
local adaptersF=$3
local analysisDir=$4
local ncpus=$5
file=`echo $files | cut -f $runN -d " "`
scratch=`dirname "$file"`
sname=`basename $file .temp.fq.gz`
echo "processing reads from sample $sname in file $file"
echo "using cutadapt to remove adapters note that the number of reads will be the same except that some will be shorter due to removing partial sequences from them"
#cutadapt --interleaved -j $ncpus -a $FAdapter -A $RAdapter -o $analysisDir/$sname.noAdpt.fq.gz $file 
#remove adapters and discard any reads with quality below 10 or with length <50 base-pairs (the pair discarded if *ANY* of them is too short or its quality is bad)
cutadapt --interleaved -j $ncpus -q 10 -m 50 -pair-filter=any -a file:$adaptersF -A file:$adaptersF -o $analysisDir/$sname.noAdptQtrim.fq.gz $file 
echo "create symbolic links to these files"
if [ -f $analysisDir/$sname.noAdptQtrim.fq.gz ]; then
unlink $scratch/$sname.temp.fq.gz; ln -s $analysisDir/$sname.noAdptQtrim.fq.gz $readScratch/$sname.temp.fq.gz 
fi
}
export LSB_JOB_REPORT_MAIL="N"
export -f func1
bsub -q $que  -n $ncpus -R "span[hosts=1] rusage[mem=4000]" -P microbiome -J "remAdapt[1-$nSamples]" -oo "$readScratch/remAdapt_%I.log" "func1 \$LSB_JOBINDEX '$files' '$adaptersF' '$analysisDir' '$ncpus'"
