module load python/2.7.15-rhel7
readScratch=$1
analysisDir=$2
outputF=$3
bsub -q standard  -n 1 -R "span[hosts=1] rusage[mem=35000]" -P microbiome -J "krakenSummary" -oo "$readScratch/getRichnessInsamples.log" "python scripts/getReadsClassificationsInsamples.py '$analysisDir' '$outputF'"
