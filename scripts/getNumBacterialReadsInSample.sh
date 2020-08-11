krakenTag="krakenLabels"
files=`ls -p analysis/*.krakenLabels`
nFs=`echo $files|wc -w`
echo "the number of samples = $nFs"
rm -f analysis/krakenClassReadsInSamples.csv
echo "Sample,reads classified">>analysis/krakenClassReadsInSamples.csv
for f in $files;do
    sname=`basename $f|sed -e 's/.krakenLabels//g'`
    nReadsClassified=`cat $f|wc -l`
    echo `echo $sname,$nReadsClassified`>>analysis/krakenClassReadsInSamples.csv
done

