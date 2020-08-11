genomeSource=$1

genomeIndexSource=$2

genomeDistFile=$3

bwaIndexDistFile=$4

wget $genomeSource -O $genomeDistFile

wget $genomeIndexSource -O $bwaIndexDistFile  
