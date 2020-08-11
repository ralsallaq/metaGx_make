echo "samples = "`echo $(cat analysis/trackReads.out|wc -l)-1|bc`
echo "samples with low sequencing read number (<1e6) = "
cat analysis/trackReads.out |awk '{if($2<1e6) print $1"\t"$2}'|column -t
echo "assembled = "$(find analysis/*.spades_out/ -name scaffolds.fasta|grep 'out/scaffolds.fasta'|wc -l)
echo "prokka predicted = "$(find analysis/prokka_Annot.*/* -name prokka.faa|wc -l)
echo "eggeNOG annotated = "$(find analysis/eggnog_map.*/* -name egmap.emapper.annotations|wc -l)
