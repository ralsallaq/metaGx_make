module load python/2.7.15-rhel7
resistomeTags="gene class mechanism group"
for tag in $resistomeTags; do
    files=`ls -p AMR/*$tag.tsv`
    echo "number of files with tag $tag is" $(echo $files|wc -w)
    python ./scripts/combineResistomeFiles.py -i $files -o AMR/allsamples.$tag.csv
done


#rarefactTags="gene.rareFD class.rareFD mech.rareFD group.rareFD"
#for tag in $rarefactTags; do
#    files=`ls -p AMR/*$tag.tsv`
#    echo "number of files with tag $tag is" $(echo $files|wc -w)
#    python ./scripts/combineResistomeFiles.py -i $files -o AMR/allsamples.$tag.csv
#done

files=`ls -p AMR/*_L???.tsv`
echo "number of output files from SNP finder is " $(echo $files|wc -w)
python ./scripts/combineResistomeFiles.py -i $files -o AMR/allsamples.SNP.csv
