SHELL:=/bin/bash
#genomeSource:=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz
#genomeIndexSource:=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.bwa_index.tar.gz
#HumanGenomeFilter="both" #can be "EBV" or "MT" or "both" 
#idxbase="human_genome_bwa"
#adaptersF=/research/home/margogrp/ralsalla/plate_164606_wgs-1/adapters.fasta
#bbresources=/research/home/margogrp/ralsalla/plate_164606_wgs-1/resources
#dataDir=/research/dept/hart/PI_data_distribution/margogrp/GSF/margogrp_164606_wgs-1/
#readScratch=/scratch_space/ralsalla/reads_164606/
#analysisDir=/research/home/margogrp/ralsalla/plate_164606_wgs-1/analysis/
##que can be standard or normal or short
#que="standard"
#ncpus=4
##define variables
##createDirs:
##	mkdir working && mkdir working/bwa_index 	
##./working/genome.fna.gz working/bwa_index/index.tgz: $(genomeSource) $(genomeIndexSource) 
##	./scripts/getHGenomeAndIndex.sh $^ $@
#
##HumanGenomeFile:=working/genome.fna.gz
##	
##./working/genome.filtered.fna: working/genome.fna.gz  
##	./scripts/filterHumanG.sh $^ $(HumanGenomeFilter) $@   
##HumanGenomeFile:=./working/genome.filtered.fna.gz
##./working/bwa_index/bwa_index.tgz $(idxbase): $(HumanGenomeFile)
##	./scripts/bwa_indexHumanG.sh $^ $@

prepAndlinkHumanG:
	mkdir humanG && ln -s /research/home/margogrp/ralsalla/metaG/pipeline/working/* humanG/ ; ln -s /research/home/margogrp/ralsalla/metaG/pipeline/human_genome_bwa.* ./ ; mkdir $(readScratch) ; mkdir $(analysisDir) 
linkreads:
	./scripts/verifyPair.sh $(dataDir) $(readScratch) $(que) $(ncpus) 
checkquality:
	./scripts/checkquality.sh $(dataDir) $(readScratch) $(que) $(ncpus)
removeAdapters:
	./scripts/removeAdapters.sh $(adaptersF) $(readScratch) $(analysisDir) $(que) $(ncpus)
filterArtifacts:
	./scripts/filterArtifacts.sh $(readScratch) $(analysisDir) $(que) $(ncpus)
decontByMap:
	./scripts/decontaminate.sh $(readScratch) $(analysisDir) $(que) $(ncpus)
qualityTrim:
	./scripts/qtrimReads.sh $(readScratch) $(analysisDir) $(que) $(ncpus)
trackReads: 
	scripts/trackReads.sh $(readScratch)  $(analysisDir) all $(que) 
assembleTocontigsWspades:
	./scripts/assembleReads_spades.sh $(readScratch) $(analysisDir) $(que) $(ncpus)
assembleTocontigsWmegahit:
	./scripts/assembleReads_megahit.sh $(readScratch) $(analysisDir) $(que) $(ncpus)
assembleTocontigsWtadpole:
	./scripts/assembleReads_tadwrapper.sh $(readScratch) $(analysisDir) $(que) $(ncpus)
#to do:
evaluateAssemblies:
	./scripts/evaluateAssemblies.sh $(readScratch) $(analysisDir) $(que) $(ncpus) "spades"
#pick from the evaluation one assembly for further analysis
assemblyChoice="spades"
taxonomicMakeup:
	./scripts/postAssembly.sh $(readScratch) $(analysisDir) $(assemblyChoice) $(que) $(ncpus)

#annotate the contigs using prokka
annotate:
	./scripts/prokkaAnnotate.sh $(readScratch) $(analysisDir) $(assemblyChoice) $(que) $(ncpus)

#map annotated peptides using eggnog to NOG db
mapAnnotationsToNOG:
	./scripts/mapWeggnog.sh $(readScratch) $(analysisDir) $(assemblyChoice) $(que) $(ncpus)

mapReadsToAnnotations:
	./scripts/mapReadsToAnnotations.sh $(readScratch) $(analysisDir) $(assemblyChoice) $(que) $(ncpus)

#quantify CDS genes to get genes counts (abundances or number of reads mapped): 
quanitfyGenes: 
	./scripts/quantifyGenes.sh $(readScratch) $(analysisDir) $(assemblyChoice) $(que) $(ncpus)
#Get butyrate pathways in the samples
butyratePWYs:
	./scripts/butyrate_ecpwy_analysis.sh $(readScratch) $(analysisDir) $(assemblyChoice) /research/home/margogrp/ralsalla/butyrate_EC_pathways.csv $(que) $(ncpus)

#Get bileAcids pathway(s) in the samples
bileAcidsPWYs:
	./scripts/bileAcids_ecpwy_analysis.sh $(readScratch) $(analysisDir) $(assemblyChoice) /research/home/margogrp/ralsalla/bileSalt_EC_pathways.csv $(que) $(ncpus) 

bileFT:
	./scripts/bileAcids_ecpwy_featureTab.sh $(readScratch) $(analysisDir) $(assemblyChoice) /research/home/margogrp/ralsalla/bileSalt_EC_pathways.csv 75 $(que) $(ncpus)

#Get butyrate features in samples (feature table)
butyrateFT:
	./scripts/butyrate_ecpwy_featureTab.sh $(readScratch) $(analysisDir) $(assemblyChoice) /research/home/margogrp/ralsalla/butyrate_EC_pathways.csv 75 $(que) $(ncpus)

#rarefy for butyrate features; how many distinct features per percentage of reads
butyrateFTRarefaction:
	./scripts/butyrate_ecpwy_featureTab_rareFaction.sh $(readScratch) $(analysisDir) $(assemblyChoice) /research/home/margogrp/ralsalla/butyrate_EC_pathways.csv 75 $(que) $(ncpus)

#map qtrimmed reads to MegaRes drug DB or to modified DB to identify AMR genes (ARG)
runResistomeAnalysis:
	./scripts/alignReadsToMegaRes.sh $(readScratch) $(analysisDir) $(que) $(ncpus) modified /research/home/margogrp/ralsalla/databaseFiles/ 70

#get taxonomy classification using Lowest Common Ancestor with Kraken; to get kraken database go to /scratch_space/ralsalla/databases/
classifyWKraken:
	./scripts/runKraken2.sh $(readScratch) $(analysisDir) $(que) $(ncpus) /scratch_space/ralsalla/databases/kraken2_db/ 

#produce a summary file that connects bacterial richness and number of bacteria archaea, eukaryota and virus reads to number of reads
summaryKraken:
	./scripts/getReadsClassificationsInsamples.sh $(readScratch) $(analysisDir) krakenClassReadsSummaryInSamples.csv

mapAnnotationsToUniRef:
#map annotated peptides using diamond to uniref90 db
	./scripts/mapProkkaFaaToUniref.sh $(readScratch) $(analysisDir) $(assemblyChoice) $(que) $(ncpus)
mapReadsToUniRef:
	./scripts/mapReadsToUniref.sh $(readScratch) $(analysisDir) $(assemblyChoice) $(que) $(ncpus)

summaryForPipeline:
	./scripts/testFinished.sh > summaryForPipeline.out
