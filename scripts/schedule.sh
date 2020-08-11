#!/bin/bash
#bash script to schedule the pipeline on the cluster using gnu-make and lsf (IBM)

usage()
{
        echo "usage:./scripts/schedule.sh -h "human_genome_bwa" -a adaptersFasta -b bbtoolsResources -d dataDir -s scratchDir -w analysisDir -q HPC-que -t ncpus"

}


#input idxbase="human_genome_bwa"
#input adaptersF=/research/home/margogrp/ralsalla/plate_164606_wgs-1/adapters.fasta
#input bbresources=/research/home/margogrp/ralsalla/plate_164606_wgs-1/resources
#input dataDir=/research/dept/hart/PI_data_distribution/margogrp/GSF/margogrp_164606_wgs-1/
#input readScratch=/scratch_space/ralsalla/reads_164606/
#input analysisDir=/research/home/margogrp/ralsalla/plate_164606_wgs-1/analysis/
#input #que can be standard or normal or short
#input que="standard"
#input ncpus=4

while getopts ":h: :a: :b: :d: :s: :w: :q: :t:" opt; do
    case $opt in
        h)
            echo "The human genome idxbase is set to $OPTARG"
            idxbase=$OPTARG
            ;;
        a)
            echo "The adapters fasta file is set to $OPTARG"
            adaptersF=$OPTARG
            ;;
        b)
            echo "The bbtools resource directory is set to $OPTARG"
            bbresources=$OPTARG
            ;;
        d)
            echo "The directory where the data is located is set to $OPTARG"
            dataDir=$OPTARG
            ;;
        s)
            echo "The scratch directory for log files and etc is set to $OPTARG"
            readScratch=$OPTARG
            ;;
        w)
            echo "The analysis/work directory is set to $OPTARG"
            analysisDir=$OPTARG
            ;;
        q)
            echo "The HPC que for running the analyses is set to $OPTARG"
            que=$OPTARG
            ;;
        t)
            echo "The number of cores/threads per node is set to $OPTARG"
            ncpus=$OPTARG
            ;;
        \?)
            echo "Invalid option: -$OPTARG"
            exit 1
            ;;
        :)
            echo "Option -$OPTARG reguires an argument"
            exit 1
            ;;
    esac
done

#test if options are passed
if [ -z "$idxbase" ] || [ -z "$adaptersF" ] || [ -z "$bbresources" ] || [ -z "$dataDir" ] || [ -z "$readScratch" ] || [ -z "$analysisDir" ] || [ -z "$que" ] || [ -z "$ncpus" ]; then 
    usage
    exit 1
fi


function checkJoboutput {
local output=$1
echo "This job output is $output"
check=`echo $output | grep 'make\|Job' |awk '{print $NF}'`
echo check is $check
{ # previous command failed try this next: check if make returns up to date
    #check=`echo $output | awk '{print $NF}'`
    #echo "check is $check"
    if [[ $check == "date." ]]; then
        echo "The command returned up to date status"
        continue
    elif [[ -z $check ]]; then
        echo "The command failed"
        exit 1
    else
        false
    fi

} || \
{ #try to see if it is a scheduled job array on the cluster
    if [[ -z $check ]]; then
        echo "The command failed"
        exit 1
    fi

    jj=`echo $output | grep 'Job'|tail -1|cut -f 2 -d "<"|cut -f 1 -d ">"`
    #check if it is a job array or not
    lines=`bjobs -r $jj|wc -l`
    if [[ "$lines" -gt 2 ]]; then
        echo "this seems an array job with $((lines-1)) jobs"

        vars=(`bjobs -A $jj |head -1`)
        vals=(`bjobs -A $jj |tail -1`)
        echo "$vars $vals"
        #assign vals to vars
        for ((i=0;i<${#vars[@]};i++)); do 
            eval $(echo ${vars[$i]}=${vals[$i]})
        done
    
        if [ $PEND -gt 0 ]; then
            status="PEND"
        elif [ $RUN -gt 0 ]; then
            status="RUN"
        elif [ $DONE = $NJOBS ]; then
            status="DONE"
        elif [ $EXIT = $NJOBS ]; then
            status="EXIT"
        fi
    
        echo "will wait until job array $jj finishes"
        while [ "$status" == "PEND" ] || [ "$status" == "RUN" ]; do
                #echo "waiting job array is $status"
                sleep 1
            vars=(`bjobs -A $jj |head -1`)
            vals=(`bjobs -A $jj |tail -1`)
            #assign vals to vars
            for ((i=0;i<${#vars[@]};i++)); do 
                eval $(echo ${vars[$i]}=${vals[$i]})
            done

            if [ $PEND -gt 0 ]; then
                status="PEND"
            elif [ $RUN -gt 0 ]; then
                status="RUN"
            elif [ $DONE = $NJOBS ]; then
                status="DONE"
            elif [ $EXIT = $NJOBS ]; then
                status="EXIT"
            else
                status="Partial EXIT"
            fi
        done
        #sleep 130
        echo "job array $jj is finished with status $status"
        if [[ "$status" = "EXIT" ]]; then
            echo "will try to exit schedule"
            exit 1
        elif [[ "$status" = "Partial EXIT" ]]; then
            echo "It seems some of the jobs in the array have EXIT status"
        fi
    else
        echo "this seems a regular job with $((lines-1)) jobs"
        false
    fi
} || \
{ #try to see if it is a scheduled job on the cluster
    if [[ -z $check ]]; then
        echo "The command failed"
        exit 1
    fi
    jj=`echo $output | grep 'Job'|tail -1|cut -f 2 -d "<"|cut -f 1 -d ">"`

    status=`bjobs -r "$jj"|awk '{print $3}'|head -150|tail -1`
    echo "will wait until job $jj finishes"
    while [ "$status" == "PEND" ] || [ "$status" == "RUN" ]; do
            #echo "waiting job is $status"
            sleep 1
            status=`bjobs -r "$jj"|awk '{print $3}'|tail -1`
    done
    sleep 130
    echo "singular job $jj is finished with status $status"
    if [ "$status" == "EXIT" ]; then
        echo "will try to exit schedule"
        exit 1
    fi
} || \
{ #Now check if the command failed
#check=`echo $output | awk '{print $NF}'`
#echo "check is $check"
#check if $check is empty
if [[ -z $check ]]; then
    echo "The command failed"
    exit 1
fi
} || \
{ #catch seems either is is up to date or it is not running for some reason go to next step
echo "not sure what the command returned, not failed, not up to date and no job number ...check"
exit 1
}
} #end of function

function pipeline {

    local idxbase=$1
    local adaptersF=$2
    local bbresources=$3
    local dataDir=$4
    local readScratch=$5
    local analysisDir=$6
    local que=$7
    local ncpus=$8

    echo "the pipeline is started"
    echo `date`
    echo " the command used:"

    echo "./scripts/schedule.sh -h $idxbase  -a $adaptersF  -b $bbresources  -d $dataDir   -s $readScratch  -w  $analysisDir  -q $que  -t $ncpus"
    echo "==============================================================="

    
    make prepAndlinkHumanG readScratch="$readScratch" analysisDir="$analysisDir"
    sleep 50
    out=`make linkreads dataDir="$dataDir" readScratch="$readScratch" que="$que" ncpus="$ncpus"`
    checkJoboutput "$out"
    out=`make checkquality dataDir="$dataDir" readScratch="$readScratch" que="$que" ncpus="$ncpus"`
    checkJoboutput "$out"
    
    out=`make removeAdapters adaptersF="$adaptersF" readScratch="$readScratch" analysisDir="$analysisDir" que="$que" ncpus="$ncpus"`
    checkJoboutput "$out"
    
    out=`make filterArtifacts readScratch="$readScratch" analysisDir="$analysisDir" que="$que" ncpus="$ncpus"`
    checkJoboutput "$out"
    
    out=`make decontByMap readScratch="$readScratch" analysisDir="$analysisDir" que="$que" ncpus="$ncpus"`
    checkJoboutput "$out"
    
    out=`make qualityTrim readScratch="$readScratch" analysisDir="$analysisDir" que="$que" ncpus="$ncpus"`
    checkJoboutput "$out"
    
    out=`make trackReads readScratch="$readScratch" analysisDir="$analysisDir" que="$que"`
    checkJoboutput "$out"
    
    out=`make assembleTocontigsWspades readScratch="$readScratch" analysisDir="$analysisDir" que="$que" ncpus="$ncpus"`
    checkJoboutput "$out"

    out=`make evaluateAssemblies readScratch="$readScratch" analysisDir="$analysisDir" que="$que" ncpus="$ncpus"`
    checkJoboutput "$out"
    
    #out=`make taxonomicMakeup readScratch="$readScratch" analysisDir="$analysisDir" que="$que" ncpus="$ncpus"`
    #checkJoboutput "$out"
    
    out=`make annotate readScratch="$readScratch" analysisDir="$analysisDir" que="$que" ncpus="$ncpus"`
    checkJoboutput "$out"

    out=`make mapAnnotationsToNOG readScratch="$readScratch" analysisDir="$analysisDir" que="$que" ncpus="$ncpus"`
    checkJoboutput "$out"

    out=`make mapReadsToAnnotations readScratch="$readScratch" analysisDir="$analysisDir"  que="$que"  ncpus="$ncpus"`
    checkJoboutput "$out"

    out=`make quanitfyGenes readScratch="$readScratch" analysisDir="$analysisDir"  que="$que"  ncpus="$ncpus"`
    checkJoboutput "$out"
    
    out=`make butyratePWYs readScratch="$readScratch" analysisDir="$analysisDir"  que="$que"  ncpus="$ncpus"`
    checkJoboutput "$out"

    out=`make bileAcidsPWYs readScratch="$readScratch" analysisDir="$analysisDir"  que="$que"  ncpus="$ncpus"`
    checkJoboutput "$out"

    #out=`make runResistomeAnalysis`
    #checkJoboutput "$out"

    out=`make classifyWKraken readScratch="$readScratch" analysisDir="$analysisDir" que="$que" ncpus="$ncpus"`

}

export LSB_JOB_REPORT_MAIL="N"
export -f checkJoboutput
export -f pipeline

bsub -q "standard" -n 1 -P microbiome -J "schedule_metaG" -oo "schedule.log" "pipeline '$idxbase' '$adaptersF' '$bbresources' '$dataDir' '$readScratch' '$analysisDir' '$que' '$ncpus'"

