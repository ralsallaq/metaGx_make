module load bwa/0.7.16a

function func1 {
    #ftpIndex=$1
    filtered=$1
    bwaIndex=$2
    idxbase=$3 #prefix of index file
    

    bwa index -a bwtsw -p $idxbase $filtered
    wait
    #create an archive for the bwa index for the filtered genome 
    tar czvf $bwaIndex human_genome_bwa.* 
}
export LSB_JOB_REPORT_MAIL="N"
export -f func1

bsub -q "standard" -n 1 -R "span[hosts=1] rusage[mem=12000]"  -P microbiome -J "bwaIndex" -oo "bwaIndex.log" "func1 $1 $2 $3"
