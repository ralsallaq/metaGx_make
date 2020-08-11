files1=`ls -p analysis/*.noAdpt.fq.gz`
files2=`ls -p analysis/*.filtered.fq.gz`
files=`echo "$files1" "$files2"`
echo $files
rm -f $files
