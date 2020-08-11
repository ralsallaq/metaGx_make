#!/bin/bash
readScratch=$1
analysisDir=$2
key=$3
que=$4
ncpus=1
inputFs=`ls -p $readScratch/*.reads.fq.gz`

if [ ! $key == "all" ] && [ ! $key == "decontaminated" ]; then
    files=`ls -p $analysisDir/*.$key.fq.gz`
elif [ $key == "decontaminated" ]; then
    files=`ls -p $analysisDir/*.$key.fq`
else
    noAdptFs=`ls -p $analysisDir/*.noAdpt.fq.gz`
    filtFs=`ls -p $analysisDir/*.filtered.fq.gz`
    decontFs=`ls -p $analysisDir/*.decontaminated.fq`
    #errcFs=`ls -p $analysisDir/*.ecct.fq.gz` 
    qtrimFs=`ls -p $analysisDir/*.qtrimmed.fq.gz`
    #unmergedFs=`ls -p $analysisDir/*.unmerged.fq.gz`
    files=$noAdptFs
fi

nSamples=`echo $files|wc -w`

echo "Processing $nSamples samples having paired R1 and R2 reads"

function func1 {
local files=$1
local key=$2
local analysisDir=$3
local ncpus=$4
local readScratch=$5
nfiles=`echo $files|wc -w`
if [ ! $key == "all" ] && [ ! $key == "decontaminated" ]; then
echo "#sample   #input # $key" > $analysisDir/trackReads.out
    for f in $files; do
        sname=`basename $f .$key.fq.gz`
        inputN=`echo $(zcat $readScratch/$sname.reads.fq.gz|wc -l)/4|bc`
        keyN=`echo $(zcat $analysisDir/$sname.$key.fq.gz|wc -l)/4|bc` 
        echo "$sname    $inputN  $keyN" >> $analysisDir/trackReads.out
    done

elif [ $key == "decontaminated" ]; then
echo "#sample   #input # $key" > $analysisDir/trackReads.out
    for f in $files; do
        sname=`basename $f .$key.fq`
        echo "$f $sname"
        inputN=`echo $(zcat $readScratch/$sname.reads.fq.gz|wc -l)/4|bc`
        noAdptN=`echo $(zcat $analysisDir/$sname.noAdpt.fq.gz|wc -l)/4|bc`
        keyN=`echo $(zcat $f|wc -l)/4|bc`
        echo "$sname  $inputN  $noAdptN  $keyN" >> $analysisDir/trackReads.out
    done
else
echo "#sample   #input  #noAdpt #filtered   #decontaminated #qtrimmed" | column -t > $analysisDir/trackReads.out
    for f in $files; do
        echo "$f $sname"
        sname=`basename $f .noAdpt.fq.gz`
        inputN=`echo $(zcat $readScratch/$sname.reads.fq.gz|wc -l)/4|bc`
        noAdptN=`echo $(zcat $analysisDir/$sname.noAdpt.fq.gz|wc -l)/4|bc`
        filtN=`echo $(zcat $analysisDir/$sname.filtered.fq.gz|wc -l)/4|bc`
        decontN=`echo $(cat $analysisDir/$sname.decontaminated.fq|wc -l)/4|bc`
        qtrimN=`echo $(zcat $analysisDir/$sname.qtrimmed.fq.gz|wc -l)/4|bc`
        echo "$sname    $noAdptN  $noAdptN    $filtN  $decontN    $qtrimN" | column -t >> $analysisDir/trackReads.out
    done
fi

}
export LSB_JOB_REPORT_MAIL="N"
export -f func1
bsub -q $que  -n $ncpus -R "span[hosts=1] rusage[mem=10000]" -P microbiome -J "trackReads" -oo "$analysisDir/trackReads.log" "func1 '$files' '$key' '$analysisDir' '$ncpus' '$readScratch'"
