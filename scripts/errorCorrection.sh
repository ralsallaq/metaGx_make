#!/bin/bash
activate=/hpcf/apps/conda2/install/4.1.1/bin/activate
echo "activate metaG env to get bbtools"
source $activate metaG
readScratch=$1
analysisDir=$2
files=`ls -p $analysisDir/*.qtrimmed.fq.gz`
nSamples=`echo $files|wc -w`
echo "Processing $nSamples samples having paired R1 and R2 reads"
echo "correcting errors using bbtools"
ncpus=8
if [ ! -d resources ];then
    echo "missing bbtools resources directory; you need to move it here ... stopping"
    exit 1
fi

function func1 {
local runN=$1
local files=$2
local analysisDir=$3
local ncpus=$4
file=`echo $files | cut -f $runN -d " "`
sname=`basename $file .qtrimmed.fq.gz`
echo "processing reads from sample $sname in file $file"
echo "phase1 error correct paired reads via overlap before clumping"
bbmerge.sh in=$file out=$analysisDir/$sname.ecco.fq.gz mix vstrict ordered ihist=$analysisDir/$sname.ihist_merge1.txt
wait
echo "phase 2 error correct reads by multiple passes"
if [ -f $scratch/$sname.ecco.fq.gz ]; then
  clumpify.sh in=$analysisDir/$sname.ecco.fq.gz out=$analysisDir/$sname.eccc.fq.gz ecc passes=4 reorder  
fi
echo "phase 3 error correct using tadpole using kmer counts, it is recommended to use k=1/3 of read length in error-correction mode (k=32 for 100 bps reads)"
#Low-depth reads can be discarded here with the "tossjunk", "tossdepth", or "tossuncorrectable" flags.
#tossjunk flag removes reads that cannot be used for assembly
#tossdepth=4 removes reads containing kmers at or below this depth
#mode=correct to correct errors only (do not assemble)
#For very large datasets, "prefilter=1" or "prefilter=2" can be added to conserve memory.
if [ -f $scratch/$sname.eccc.fq.gz ]; then
    #tadpole.sh in=$analysisDir/$sname.eccc.fq.gz out=$analysisDir/$sname.ecct.fq.gz ecc k=32 ordered mode=correct 
    #tadpole.sh in=$analysisDir/$sname.eccc.fq.gz out=$analysisDir/$sname.ecct.fq.gz ecc k=32 ordered mode=correct tossjunk
    tadpole.sh in=$analysisDir/$sname.eccc.fq.gz out=$analysisDir/$sname.ecct.fq.gz ecc k=32 ordered mode=correct tossdepth=4
fi

##phase 4 normalize : tadpole is superior to bbnorm so may be I will not use it
}
export LSB_JOB_REPORT_MAIL="N"
export -f func1
bsub -q "standard"  -n $ncpus -R "span[hosts=1] rusage[mem=10000]" -P microbiome -J "eCorrect[1-$nSamples]" -oo "$readScratch/eCorrect.%I.log" "func1 \$LSB_JOBINDEX '$files' '$analysisDir' '$ncpus'"
