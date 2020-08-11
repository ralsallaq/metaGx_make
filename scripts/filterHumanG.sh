module load python/3.5.2
#$1 is the full path of the unfiltered human genome ftped from NIH
#$2 is the criteria for filtering "EBV" or "MT" or "both"
#$3 is the full path of the filtered genome
#filterHumanG.py reads human genome in gz fasta format and filters away records with specific ids
bsub -q "standard" -n 1 -R "span[hosts=1] rusage[mem=10000]" -P microbiome -J "filterHumanG" -oo "filterHumanG.log" "python filterHumanG.py -i $1 -o $3 -fltr $2" 
