#!/usr/bin/env nextflow

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
This is the main document. To modify CPU usage and parameters fine tuning please go to nextflow.config
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


// Enable DSL2
nextflow.enable.dsl=2


// Calling modules

include { GET_DATA } from './modules/get_data.nf'
include { PARSE_GENBANK } from './modules/parse_genbank.nf'
include { REMOVE_REDUNDANCY } from './modules/remove_redundancy.nf'
include { GENE_FASTA_DATABASE } from './modules/gene_fasta_database.nf'
include { GENE_CLUSTERING } from './modules/gene_clustering.nf'



/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Print pipeline metadata: Version and Help
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


def print_help() {
	println "\n\033[1;31mSYNOPSIS\033[0m"

	println "\n\033[1;33mUSAGE\033[0m"
	println "\nnextflow run aBPA.nf --data <PATH> --output <PATH> --tax_id <INT> --config <PATH/FILE> [..OPTIONS..]"
	
	println "\n\033[1;33mMANDATORY\033[0m"
	println "  --data <PATH>		Set data file PATH"
	println "  --output <PATH>		Set output directory PATH"
	println "  --tax_id <INT>		Set taxonomical ID value <INT>"
	println "  --config <PATH>		Set config file PATH"
	
	println "\n\033[1;33mOPTIONS\033[0m"
	println "  --threads <INT>		Set number of threads (default: 10)"
	println "  --completeness <INT/FLOAT>	Set gene completeness/breadth of coverage threshold (default: 50)"
	println "  --coverage <INT/FLOAT>	Set mean depth of coverage threshold (default: 0.5)"
	println "  --genomes <INT>		Set number of genomes to download (default: 100)"
	println "  --clustering <INT/FLOAT>	Set clustering threshold (default 0.95)"
	println "  --core-threshold <FLOAT>	Set core genome threshold (default: 0.01)"
	println "  --clean-mode <STRING>	Set pangenome mode (default: strict)"
	println "  --help			Print help page and exit"
	
	println "\n\033[1;31mDESCRIPTION\033[0m"
	println "\n\033[1;33m--data <PATH>\033[0m"
	println "Please specify the full PATH of your data. Example: /home/user/mydata/data"
	
	println "\n\033[1;33m--output <PATH>\033[0m"
	println "Please specify the full PATH of your output folder. You need to make the folder first before running the program."

	println "\n\033[1;33m--tax_id <INT>\033[0m"
	println "Please specify the taxonomical ID for your bacteria. It should be a discrete and unique number."
	
	println "\n\033[1;33m--config <PATH>\033[0m"
	println "\nPlease set file PATH of your config.tab file. Example: /home/user/me/aBPA/config/config.tab"
	println "config.tab file should contain 3 fields separated by tab. First field should have the sample name, second field softclipping value <INT> and third field group ID."
	println "\nExample: \n	SAMPLE1	5	NONUDG\n	SAMPLE2	2	UDG"

	println "\n\033[1;33m--threads <INT>\033[0m"
	println "Set amount of threads to be used globally."

	println "\n\033[1;33m--completeness <INT>\033[0m"
	println "Set gene breadth of coverage threshold as percentage. Genes that have a value less than <INT> will be considered absent."

	println "\n\033[1;33m--coverage <INT/FLOAT>\033[0m"
	println "Set gene normalized coverage threshold. Currently aBPA is using the simplest statistic for normalization: (Gene mean depth/Global mean depth)."

        println "\n\033[1;33m--genomes <INT/FLOAT>\033[0m"
	println "Set amount of FASTA/GENBANK files to be downloaded. Bear in mind disk space."

        println "\n\033[1;33m--clustering <INT/FLOAT>\033[0m"
	println "Set clustering threshold <INT/FLOAT> for FASTA database.\nA value of 0.9 means any group of sequences with identity values equal or bigger than 0.9 will be clustered together and a consensus representative sequence will be produced."

        println "\n\033[1;33m--core-threshold <INT/FLOAT>\033[0m"
	println "Set threshold for core genome building. Similarly as clustering flag but during pangenome step."

        println "\n\033[1;33m--clean-mode <INT/FLOAT>\033[0m"
	println "Set behaviour of pangenome building. Visit Panaroo documentation to know more about this.\n\n"

    exit 0
}

def version() {
	println "aBPA version 0.2"
	exit 0
}

if (params.help) {
    print_help()
}

if (params.version) {
	version()
}

// Main workflow

workflow {
        if (!params.trusted_data) {
                GET_DATA(params.genomes, tax_id, get_data_parallel)
                fastaFiles = GET_DATA.out.fasta_files
                gffFiles = GET_DATA.out.gbk_files

        } else {
                fastaFiles = Channel.of(files("${params.trusted_data}/*fasta"))
                gffFiles = Channel.of(files("${params.trusted_data}/*gb"))

        }

	PARSE_GENBANK(gffFiles, fastaFiles)

	REMOVE_REDUNDANCY(FASTA_DATABASE.out.validFiles.map { fasta, gb -> fasta , gb}, params.remove_redundancy_parallel)

	GENE_FASTA_DATABASE(REMOVE_REDUNDANCY.out.nonRedundant_files.map { fasta, gb -> gb})

	GENE_CLUSTERING(GENE_FASTA_DATABASE.out.fastaDatabase, params.gene_identity_clustering, params.cd_hit_threads)

	ANNOTATE(GENE_CLUSTERING.out.clusteredDatabase, params.prokka_annotate_threads, REMOVE_REDUNDANCY.out.nonRedundant_files.map { fasta, gb -> fasta, gb},
		params.prokka_annotate_parallel)


	makePangenome(prokkaMakeAnnotations.out.prokkaGFF, pangenomeMode, pangenomeThreshold, threadsGlobal)
	formattingPangenome(makePangenome.out.panSequence)
	blastMe(formattingPangenome.out.panGenomeReference)
        outgroupEntrez(outTax)
        makeReads(outgroupEntrez.out.outgroupFasta)
        outgroupAlignmentFAndiltering(makeReads.out.outgroupReads, formattingPangenome.out.panGenomeReference, threadsGlobal)
        makeOutgroupConsensus(outgroupAlignmentFAndiltering.out.outgroupFastaPostAlignment, formattingPangenome.out.panGenomeReference)
	alignment(reads, formattingPangenome.out.panGenomeReference, threadsGlobal, configFile, missingProb, seedAlignment, gapFraction, minReadLength, maxReadLength)
	alignmentSummary(configFile, alignment.out.postAlignedBams)
	normalizationFunction(alignmentSummary.out.refLenght, alignmentSummary.out.rawCoverage)
	updateNormalization(normalizationFunction.out.geneNormalizedSummary, alignmentSummary.out.completenessSummary)


	if (params.genotyper == "gatk") {
		gatkConsensus(formattingPangenome.out.panGenomeReference, alignmentSummary.out.postAlignmentFiles, formattingPangenome.out.panGenomeReferenceDictionary, formattingPangenome.out.panGenomeReferenceIndex)
		extractedSequencesFasta = gatkConsensus.out.gatkConsensusSequences
		vcfFile = gatkConsensus.out.gatkGenotypes

	} else if (params.genotyper == "bcftools") {
		bcftoolsConsensus(formattingPangenome.out.panGenomeReference, alignmentSummary.out.postAlignmentFiles)
		extractedSequencesFasta = bcftoolsConsensus.out.consensusSequences

	} else {
		error "Invalid option for --genotyper. Please choose 'gatk' or 'bcftools'."
	}

	plotCoveragevsCompleteness(updateNormalization.out.geneNormalizedUpdated, geneCompleteness, normalizedCoverageDown)
        applyCoverageBounds(updateNormalization.out.geneNormalizedUpdated, normalizedCoverageDown, normalizedCoverageUp, geneCompleteness)
	makeMatrix(makePangenome.out.initialMatrix , normalizationFunction.out.globalMeanCoverage, applyCoverageBounds.out.geneNormalizedUpdatedFiltered, geneCompleteness, normalizedCoverageDown, normalizedCoverageUp)
	buildHeatmap(makeMatrix.out.finalCsv, makeMatrix.out.INDEX ,makeMatrix.out.matrix, makeMatrix.out.sampleNames)
	plotCoveragevsCompletenessOnFiltered(applyCoverageBounds.out.geneNormalizedUpdatedFiltered, geneCompleteness,normalizedCoverageDown)
	filterGeneAlignments(makePangenome.out.alignedGenesSeqs, extractedSequencesFasta, fastaDatabase.out.validFasta, downloadGenomes, makeOutgroupConsensus.out.extractedSequencesOutgroupFasta, buildHeatmap.out.blackListed)
	pMauve(fastaDatabase.out.validFasta)
	makeMSA(filterGeneAlignments.out.genesAlnSeq, buildHeatmap.out.maskedMatrixGenesNoUbiquitous, buildHeatmap.out.maskedMatrixGenesOnlyAncient, buildHeatmap.out.maskedMatrixGenesUbiquitous, buildHeatmap.out.genesAbovePercentSeries, filterGeneAlignments.out.sampleNames)
	treeThreshold(makeMSA.out.genesAbovePercentMSA)
	treeUbiquitous(makeMSA.out.maskedMatrixGenesUbiquitousMSA)
	treeNoUbiquitous(makeMSA.out.maskedMatrixGenesNoUbiquitousMSA)
	treeAncient(makeMSA.out.maskedMatrixGenesOnlyAncientMSA)
	xmfaToFasta(pMauve.out.pMauveCoreGenome)
	filterMauveFasta(xmfaToFasta.out.pMauveFastaMSA)
	startingTree(filterMauveFasta.out.concatenatedSeqtkMauveFastaMSA)
	findRecombinationSpots(filterMauveFasta.out.concatenatedSeqtkMauveFastaMSA, startingTree.out.startingTreeMauveFasta, startingTree.out.kappa)
	mapRecombinantsToGenes(findRecombinationSpots.out.recombinationMap, filterMauveFasta.out.concatenatedSeqtkMauveFastaMSA, blastMe.out.panGenomeReferenceDB, prokkaMakeAnnotations.out.prokkaGFF)
	getResults(
	resultsDir, fastaDatabase.out.validFasta , fastaDatabase.out.validGff , fastaDatabase.out.fastaDatabaseLogFile , fastaDatabase.out.theFastaDatabase, 
	clustering.out.clusteredDatabase, clustering.out.clusteringLog, prokkaMakeAnnotations.out.prokkaGFF, prokkaMakeAnnotations.out.prokkaLogfile,  makePangenome.out.panarooLog,
	filterGeneAlignments.out.genesAlnSeq, formattingPangenome.out.panGenomeReference, updateNormalization.out.geneNormalizedUpdated, normalizationFunction.out.globalMeanCoverage,
	alignmentSummary.out.postAlignmentFiles, alignmentSummary.out.refLenght, alignmentSummary.out.rawCoverage, alignmentSummary.out.completenessSummary, buildHeatmap.out.finalMatrix,
	buildHeatmap.out.presenceAbsence, buildHeatmap.out.maskedMatrixGenesOnlyAncient, buildHeatmap.out.maskedMatrixGenesUbiquitous, buildHeatmap.out.maskedMatrixGenesNoUbiquitous,
	buildHeatmap.out.genesAbovePercentSeries, treeThreshold.out.genesAbovePercentMSAIqtree ,treeThreshold.out.genesAbovePercentMSALog , treeThreshold.out.genesAbovePercentMSATreefile,
	treeUbiquitous.out.maskedMatrixGenesUbiquitousMSAIqtree, treeUbiquitous.out.maskedMatrixGenesUbiquitousMSALog, treeUbiquitous.out.maskedMatrixGenesUbiquitousMSATreefile, 
	treeNoUbiquitous.out.maskedMatrixGenesNoUbiquitousMSAIqtree , treeNoUbiquitous.out.maskedMatrixGenesNoUbiquitousMSALog , treeNoUbiquitous.out.maskedMatrixGenesNoUbiquitousMSATreefile,
	treeAncient.out.maskedMatrixGenesOnlyAncientMSAIqtree , treeAncient.out.maskedMatrixGenesOnlyAncientMSALog , treeAncient.out.maskedMatrixGenesOnlyAncientMSATreefile
	)
}





process makePangenome {
	conda "${projectDir}/envs/panaroo.yaml"
	
	input:
	path prokkaGFF, stageAs: 'filteredGFF/*'
	val pangenomeMode
	val pangenomeThreshold
	val threadsGlobal

	output:
	path 'pan_genome_reference.fa' , emit: panSequence
	path 'gene_presence_absence.Rtab' , emit: initialMatrix
	path 'aligned_gene_sequences/*' , emit: alignedGenesSeqs
	path 'makePangenome.log', emit: panarooLog

	script:
	"""
	panaroo -i filteredGFF/*.gff -o ./ --clean-mode $pangenomeMode -a core --core_threshold $pangenomeThreshold -t $threadsGlobal

	cat .command.out >> makePangenome.log
	"""
}


process  formattingPangenome {
	conda "${projectDir}/envs/formattingPangenome.yaml"

	input:
	path panGenomeReference, stageAs: 'pan_genome_reference.fa'

	output:
	path 'panGenomeReference.fasta', emit: panGenomeReference
	path 'panGenomeReference.dict', emit: panGenomeReferenceDictionary
	path 'panGenomeReference.fasta.fai', emit: panGenomeReferenceIndex

	script:
	"""
	seqtk seq $panGenomeReference > panGenomeReference.fasta
	picard CreateSequenceDictionary -R panGenomeReference.fasta
	samtools faidx panGenomeReference.fasta
	"""
}


process alignment {
	conda "${projectDir}/envs/alignment.yaml"

	input:
	path reads
	path panRef, stageAs: 'panGenomeReference.fasta'
	val threadsGlobal
	path configFile
	val missingProb
	val seedAlignment
	val gapFraction
	val minReadLength
	val maxReadLength

	output:
	path '*_DMC_P.bam', emit: postAlignedBams
	path '*_final.fastq', emit: postAlignedReads
	

	script:
	"""
	for sample in $reads/*; do
    		bwa index $panRef
    		name=\$(basename "\$sample")
    		softClip=\$(grep "\$name" $configFile | awk '{print \$2}')
    
    		# Making read groups
    		rg_id="\${name%.fastq*}"  # sample name as id
    		rg_sm="\${name%.fastq*}" # sample name again
    		rg_pl="illumina"        # I dont think this is very important for this pipeline so its going to be just illumina because why not
    		rg_lb="lib1"            # group id
    		rg_pu="unit1"           # not sure what Ill put here

    		bwa aln -l $seedAlignment -n $missingProb -o $gapFraction -t $threadsGlobal $panRef "\$sample" > "\${name%.fastq*}.sai"
    		bwa samse -r "@RG\\tID:\$rg_id\\tSM:\$rg_sm\\tPL:\$rg_pl\\tLB:\$rg_lb\\tPU:\$rg_pu" \
		$panRef "\${name%.fastq*}.sai" "\$sample" > "\${name%.fastq*}.sam"
    
    		samtools view -bS "\${name%.fastq*}.sam" > "\${name%.fastq*}.bam"
    		samtools quickcheck "\${name%.fastq*}.bam"
    		samtools sort -o "\${name%.fastq*}_sorted.bam" -O bam -@ $threadsGlobal "\${name%.fastq*}.bam"
    		samtools index "\${name%.fastq*}_sorted.bam"
    		rm "\${name%.fastq*}.bam"
    		samtools view -b -@ $threadsGlobal -F 4 "\${name%.fastq*}_sorted.bam" > "\${name%.fastq*}_sorted_mappedreads.bam"
    		samtools index "\${name%.fastq*}_sorted_mappedreads.bam"
    		bam trimBam "\${name%.fastq*}_sorted_mappedreads.bam" "\${name%.fastq*}_softclipped.bam" -L "\$softClip" -R "\$softClip" --clip
    		samtools view -q 30 -o "\${name%.fastq*}_qc.bam" "\${name%.fastq*}_softclipped.bam"
    		samtools view -e 'length(seq)>$minReadLength && length(seq)<$maxReadLength' -O BAM -o "\${name%.fastq*}_lg.bam" "\${name%.fastq*}_qc.bam"
    		samtools sort -o "\${name%.fastq*}_DMC_P.bam" -O bam -@ $threadsGlobal "\${name%.fastq*}_lg.bam"
    		samtools coverage "\${name%.fastq*}_DMC_P.bam" > "\${name}_genomicsMetrics.txt"
    		samtools fastq -@ $threadsGlobal "\${name%.fastq*}_DMC_P.bam" > "\${name%.fastq*}_final.fastq"
	done

	rm *sam *sai *_lg.bam *_qc.bam *_sorted_mappedreads.bam*
	cat .command.out >> alignment.log
	"""
}


process alignmentSummary {
	conda "${projectDir}/envs/alignment.yaml"
	
	input:
	path configFile
	path bamfiles, stageAs: 'bam/*'
	

	output:
        path 'postPangenomeAlignment*bam' , emit: postAlignmentFiles
	path 'completenessSummary.tab', emit: completenessSummary
	path '*_refLength.txt', emit: refLenght
	path '*_rawCoverage.txt' , emit: rawCoverage


	script:
	"""
	#!/bin/bash
	awk '{print \$NF}' $configFile | uniq > groups.txt

	while read -r groupID; do
		groupName=\$(echo "\$groupID")
		grep -w "\$groupID" $configFile | awk '{print \$1}' > "\$groupName"ID
	done < groups.txt

        for groupFile in *ID; do
                IDs=\$(basename "\${groupFile%ID}")

		echo "Group ID: \$IDs"

                bamFiles=()

                while read -r sampleName; do
                        bamFiles+=(./bam/"\${sampleName%.fastq*}_DMC_P.bam")
                        echo "Adding BAM file: ./bam/\${sampleName%.fastq*}_DMC_P.bam"
		done < "\$groupFile"

                if [ \${#bamFiles[@]} -eq 1 ]; then
                        cp "\${bamFiles[0]}" postPangenomeAlignment_"\${IDs}".bam
                elif [ \${#bamFiles[@]} -gt 1 ]; then
                        samtools merge postPangenomeAlignment_mergedGroup"\${IDs}".bam "\${bamFiles[@]}"
                fi
        done

	#Getting some stats AND Re group the data so genotyping is simpler
        for i in postPangenomeAlignment*bam; do
                samplename=\$(basename ./bam/"\${i%.bam}")
                samtools index "\$i"
                samtools depth -a "\$i" > "\${samplename}_rawCoverage.txt"
                samtools idxstats "\$i" | awk '{sum += \$2} END {print sum}' > "\${samplename}_refLength.txt"
                samtools coverage "\$i" | awk -v samplename="\$samplename" 'NR>1 {print samplename, \$1, \$6}' | sed -e 's/~/_/g' | sed -e 's/ /\t/g' | sort -k 1 -t \$'\t' >> completenessSummary.tab
		mv "\$i" ./"\${i%.bam}_TMP.bam"
		mv "\$i".bai ./"\${i%.bam.bai}_TMP.bam.bai"
		picard AddOrReplaceReadGroups I="\${i%.bam}_TMP.bam" O="\${samplename}.bam" RGLB="\${samplename}" RGSM="\${samplename}" RGPU=Illumina RGPL=ILLUMINA RGID="\${samplename}" RGDS="\${samplename}"
		samtools index "\${samplename}.bam"
	done

	rm *TMP.bam*

	cat .command.out >> alignmentSummary.log
	"""
}

process normalizationFunction {

	input:
	path refLength, stageAs: 'refLength/*'
	path rawCoverage, stageAs: 'rawCoverage/*'


	output:
	path 'geneNormalizedSummary.txt', emit: geneNormalizedSummary
	path 'globalMeanCoverage.txt' , emit: globalMeanCoverage

	script:
	"""
	#!/bin/bash

	echo -e "sampleID\tgene\tnormalizedGeneSimple\tnormalizedGeneScaled\tnormalizedGenomeSimple\tnormalizedGenomeScaled" > geneNormalizedSummary.txt
	echo -e "sampleID \t sampleCoverage \t refCount \t globalMean"  > globalMeanCoverage.txt

	for i in rawCoverage/*_rawCoverage.txt; do
		name=\$(basename "\${i%_rawCoverage.txt}")

		#Compute global mean coverage
		globalMean=\$(awk -v name="\$name" '{sum += \$3; count++} END {if (count > 0) print sum / count; else print "Something went wrong, check log file"}' "\$i")
		finalCount=\$(awk -v name="\$name" '{fcount++} END {print fcount}' "\$i")
		refCount=\$(cat refLength/"\${name}_refLength.txt")
		echo -e "\$name\t\$finalCount\t\$refCount\t\$globalMean" >> globalMeanCoverage.txt

		#Normalize coverage per gene
		awk -v globalMean="\$globalMean" -v name="\$name" -v sampleCoverage="\$finalCount" -v refCount="\$refCount" '
			{
			if (\$2 > geneLength[\$1]) {
				geneLength[\$1] = \$2
			}
			sumgene[\$1] += \$3
			countgene[\$1]++
			}
		END {
			for (gene in geneLength) {
				geneMean = sumgene[gene] / countgene[gene]
				normalizedGeneScaled = (geneMean / globalMean) * geneLength[gene]
				normalizedGeneSimple = (geneMean / globalMean)
				normalizedGenomeSimple = (geneMean / globalMean) * (geneLength[gene] / sampleCoverage)
				normalizedGenomeScaled = (geneMean / globalMean) * (geneLength[gene] / refCount)
				print name"\t"gene"\t"normalizedGeneSimple"\t"normalizedGeneScaled"\t"normalizedGenomeSimple"\t"normalizedGenomeScaled
			}
		}
		' "\$i" >> geneNormalizedSummary.txt
	done
	"""
}

process updateNormalization {

	input:
	path normalized, stageAs: 'normalized/*'
	path completeness, stageAs: 'completeness/*'


	output:
	path 'geneNormalizedUpdated.tab', emit: geneNormalizedUpdated


	script:
	"""
	#!/bin/bash
	echo -e "sampleID\tgene\tnormalizedGeneSimple\tnormalizedGeneScaled\tnormalizedGenomeSimple\tnormalizedGenomeScaled\tgeneCompleteness" > geneNormalizedUpdated.tab

	sed -i -e 's/~/_/g' normalized/geneNormalizedSummary.txt

	awk 'NR>1{print \$1"XYZ"\$2, \$3, \$4, \$5, \$6}' normalized/geneNormalizedSummary.txt > TMP1

	awk '{print \$1"XYZ"\$2, \$3}' completeness/completenessSummary.tab > TMP2

	while read -r ID completeness;do

		if grep -wq "\${ID}" TMP1; then
			oldLine=\$(grep -w "\${ID}" TMP1)
			specificCompleteness=\$(grep -w "\${ID}" TMP2 | awk '{print \$NF}')
			echo -e "\${oldLine}\t\${specificCompleteness}" >> geneNormalizedUpdated.tab
		fi

	done < TMP2

	sed -i -e 's/XYZ/\t/g' geneNormalizedUpdated.tab
	sed -i -e 's/ /\t/g' geneNormalizedUpdated.tab

	rm TMP1 TMP2 
	"""
}

process applyCoverageBounds {
	
	input:
	path geneNormalizedUpdated, stageAs: 'geneNormalizedUpdated.tab'
	val normalizedCoverageDown
	val normalizedCoverageUp
	val completenessBound

	output:
	path 'geneNormalizedUpdatedFiltered.tab', emit: geneNormalizedUpdatedFiltered

	script:
	"""
	awk 'NR==1{print \$0}' $geneNormalizedUpdated > header
	awk -v UpBound=$normalizedCoverageUp '\$3 < UpBound {print \$0}' $geneNormalizedUpdated > TMP1
	awk -v DownBound=$normalizedCoverageDown '\$3 > DownBound {print \$0}' TMP1 > TMP2
	awk -v completenessBound=$completenessBound '\$NF> completenessBound {print \$0}' TMP2 > TMP3
	cat header TMP3 > geneNormalizedUpdatedFiltered.tab
	"""
}

process plotCoveragevsCompleteness {
	conda "${projectDir}/envs/plot.yaml"
	
	input:
	path geneNormalizedUpdated, stageAs: 'gNS/'
	val completeness
	val coverage

	output:
	path 'plotCoverage_vs_Completeness.png', emit: plotCoverage_vs_Completeness
	
	script:
	"""
	plot_cvg_vs_completeness.py gNS/geneNormalizedUpdated.tab $completeness $coverage
	"""
}





process plotCoveragevsCompletenessOnFiltered {
	conda "${projectDir}/envs/plot.yaml"
	
	input:
	path geneNormalizedUpdatedFiltered, stageAs: 'gNS/*'
	val completeness
	val coverage

	output:
	path 'plotCoverageVsCompletenessFiltered.png', emit: plotCoverageVsCompletenessFiltered
	
	script:
	"""
	plot_cvg_vs_completeness.py gNS/geneNormalizedUpdatedFiltered.tab $completeness $coverage
	mv plotCoverage_vs_Completeness.png ./plotCoverageVsCompletenessFiltered.png
	"""
}

process makeMatrix {
	conda "${projectDir}/envs/plot.yaml"

	input:
	path pangenomeRtab, stageAs: 'pangenome/*'
	path gMC, stageAs: 'gMC/*'
	path normalized, stageAs: 'normalized/*'
	val gCompleteness
	val lowerBound
	val upperBound


	output:
	path 'matrix.tab', emit: matrix
	path '*_final.csv', emit: finalCsv
	path 'sample_names', emit: sampleNames
	path 'INDEX', emit: INDEX

	script:
	"""
	awk 'NR==1{print \$0}' pangenome/gene_presence_absence.Rtab > matrix.tab
	awk 'NR>1 {print \$0}' pangenome/gene_presence_absence.Rtab | sort -k 1 -t \$'\t' >> matrix.tab
	awk 'NR>1 {print \$1}' pangenome/gene_presence_absence.Rtab | sort -k 1 -t \$'\t' > INDEX

	awk 'NR>1 {print \$1}' gMC/globalMeanCoverage.txt > sample_names

	while read -r name; do

		echo -e "Gene\tnormalizedCoverage\tcompleteness" > "\${name}"_index.tmp
		grep -w "\$name" normalized/geneNormalizedUpdatedFiltered.tab | awk '{print \$2, \$3, \$NF}' >> "\${name}"_index.tmp

	done < sample_names


	for i in *_index.tmp; do

		sed -i -e 's/ /\t/g' "\$i"
		lambda.py "\$i" $gCompleteness $lowerBound $upperBound

	done
	"""
}


process buildHeatmap {
	conda "${projectDir}/envs/heatmap.yaml"

	input:
	path fCSV, stageAs: 'fCSV/*'
	path INDEX, stageAs: 'INDEX/*'
	path matrix, stageAs: 'matrix/*'
	path names, stageAs: 'names/*'

	output:
	path 'final_matrix.tab', emit: finalMatrix
	path 'presenceAbsence*.png', emit: presenceAbsence
	path 'maskedMatrixGenesOnlyAncient.txt', emit: maskedMatrixGenesOnlyAncient
	path 'maskedMatrixGenesUbiquitous.txt', emit: maskedMatrixGenesUbiquitous
	path 'maskedMatrixGenesNoUbiquitous.txt', emit: maskedMatrixGenesNoUbiquitous
	path 'sampleOrdernoUbiquitous.txt', emit: sampleOrdernoUbiquitous 
	path 'sampleOrderonlyAncient.txt', emit: sampleOrderonlyAncient
	path 'genesAbovePercentSeries.txt', emit: genesAbovePercentSeries
	path 'blackListedQualityChecked.txt', emit: blackListed

	script:
	"""
	for i in fCSV/*_final.csv; do

		name=\$(basename "\$i")
		sed  -e 's/,/\t/g' "\$i" | awk 'NR>1{print \$0}' > "\${name%_index.tmp_final.csv}"_INDEX.Z

	done


	for i in *_INDEX.Z; do
		name=\$(basename "\$i")
		#Create the FINAL_INDEX file
		echo "\${name%_INDEX.Z}" > "\${name}"_FINAL_INDEX
    
		#Process the INDEX file
		while read -r gene; do
			toprint=\$(echo "\$gene 0")
			if grep -wq "\$gene" "\$i"; then
				grep -w "\$gene" "\$i" >> "\${name}"_FINAL_INDEX
			else
				echo "\$toprint" >> "\${name}"_FINAL_INDEX
			fi
		done < INDEX/INDEX
    
		#Extract the last column
		awk '{print \$NF}' "\${name}"_FINAL_INDEX > "\${name}"_FINALCOLUMN

	done
        # I need to add a quality control step right here. User samples can be false positives sometimes or just super low quality and have 0 genes after filtering
        # Then, black list unwanted samples and exclude them from the final_matrix.tab document.
        for checkSample in *_FINALCOLUMN; do
                sampleName=\$(basename "\${checkSample%_INDEX.Z_FINALCOLUMN}")
                counts=\$(grep -c "0" "\$checkSample")
                totalLines=\$(wc -l "\$checkSample" | awk '{print \$1 - 1}')
                proportion=\$(awk -v absence="\$counts" -v record="\$totalLines" 'BEGIN { print (absence / record ) }')         

                if  (( \$(awk -v p="\$proportion" 'BEGIN { print (p > 0.95) }' ) )); then
                        echo "\$sampleName" >> blackListedQualityChecked.txt
                fi
        done
	
	# Do this only if blacklisted
	if [[ -s blackListedQualityChecked.txt ]]; then
		while read -r removeMe; do
			mv "\${removeMe}_INDEX.Z_FINALCOLUMN" "\${removeMe}LowQualitySample"    
			grep -v "\${removeMe}" names/sample_names > names/sample_names.tmp
			mv names/sample_names.tmp names/sample_names
		done < blackListedQualityChecked.txt    
	else
		touch blackListedQualityChecked.txt
	
	fi

	paste matrix/matrix.tab *_FINALCOLUMN > final_matrix.tab
	tr '\n' ' ' < names/sample_names > names_heatmap

	heatmap.py final_matrix.tab names_heatmap
	
	"""
}


process bcftoolsConsensus {
	conda "${projectDir}/envs/consensus.yaml"

	input:
	path panGenomeRef, stageAs: 'panGenomeRef.fasta'
	path bamFiles, stageAs: 'BAM/*'

	output:
	path 'extractedSequences*.fasta', emit: consensusSequences

	script:
	"""
	for b in BAM/*; do
		basename=\$(basename "\$b")
		bcftools mpileup -f $panGenomeRef "\$b" | bcftools call -c | vcfutils.pl vcf2fq > extractedSequences"\${basename%.bam}".fq
		seqtk seq -a extractedSequences"\${basename%.bam}".fq > extractedSequences"\${basename%.bam}".fasta
	done
	"""
}


process gatkConsensus {
	conda "${projectDir}/envs/gatk.yaml"

	input:
	path panGenomeRef, stageAs: 'panGenomeRef.fasta'
	path bamFiles, stageAs: 'BAM/*'
	path panGenomeRefDictionary, stageAs: 'panGenomeRef.dict'
	path panGenomeReferenceIndex, stageAs: 'panGenomeRef.fasta.fai'
	output:
	path 'extractedSequences*.fasta', emit: gatkConsensusSequences
	path '*GenotypedNormalized.vcf.gz', emit: gatkGenotypes

	script:
	"""

	for b in BAM/*; do
		basename=\$(basename "\$b")
		gatk3 -T UnifiedGenotyper --min_base_quality_score 30 \
			--genotype_likelihoods_model BOTH --annotateNDA \
			--genotyping_mode DISCOVERY --output_mode EMIT_ALL_SITES \
			-I "\$b" -R panGenomeRef.fasta \
			-o "\${basename%.bam}"Genotyped.vcf 

		bcftools norm -f panGenomeRef.fasta "\${basename%.bam}"Genotyped.vcf > "\${basename%.bam}"GenotypedNormalized.vcf
		bgzip -i "\${basename%.bam}"GenotypedNormalized.vcf
		bcftools index "\${basename%.bam}"GenotypedNormalized.vcf.gz
		bcftools consensus -a N -M N -f panGenomeRef.fasta "\${basename%.bam}"GenotypedNormalized.vcf.gz -o "\${basename%.bam}"GenotypedNormalizedConsensus.fasta
		seqtk seq "\${basename%.bam}"GenotypedNormalizedConsensus.fasta > extractedSequences"\${basename%.bam}".fasta
	done
	"""
}







/*
 * First seqtk seq every gene.aln file, and then fix FASTA header. Add samples gene sequences to the right alignments. check sequence lenght and fill with Ns if lenght smaller than modern strains.
 * Finally check number of FASTA headers and if not == to \$genomes then add the missing and fill sequence with Ns.
 */


process filterGeneAlignments {
	conda "${projectDir}/envs/filterGeneAlignments.yaml"

	input:
	path genesAln, stageAs: 'panaroo_genes/*'
        path extractedSequencesFasta, stageAs: 'user_genes/*'
	path fFiles, stageAs: 'modern_data/*'
	val genomes
	path outgroupSeq, stageAs: 'outgroup'
	path blackListed, stageAs: 'blackListed.txt'

	output:
	path 'filteredGenes/*_Filtered.fasta', emit: genesAlnSeq
	path 'sampleNames.txt', emit: sampleNames
	path 'specialCases/*fasta', emit: specialCases

	script:
	"""
	#!/bin/bash

	shopt -s nullglob  #Prevent literal interpretation of wildcards if there are no matchings
	#sanity check before starting the process

	if ! find panaroo_genes/ -name "*.aln.fas" | grep -q .; then
    		echo -e "No files with .aln.fas extension were found in panaroo_genes/. Check previous process. Stopping the pipeline."
    		exit 1
	else
		echo -e "Fasta files with expected extension .aln.fas were found. Proceeding with the process."
	fi

	if [[ -e outgroup ]]; then
	        sed -i -e 's/~/_/g' outgroup
	else
        	echo -e "outgroup file was not found. Stopping the pipeline."
        	exit 1
	fi

	if ! find modern_data/ -name "*fna" | grep -q .; then
    		echo -e "No modern sequences with .fna extension were found in modern_data/. Check previous process. Stopping the pipeline."
    		exit 1
	else
        	echo -e "Fasta files with expected extension .fna were found. Proceeding with the process."
	fi


	# modern_samples_list() will create a text file with modern genomes names
	modern_samples_list() {
	fasta=\$1
        	name=\$(basename "\${fasta%.fna}")
        	echo "\${name}" > "\${name}"_modern_TMP
	}
	export -f modern_samples_list
	find modern_data/ -name "*.fna" | parallel -j 30 modern_samples_list

	cat *_modern_TMP >> modernSampleNames.txt
	rm *_modern_TMP


	# Adding the outgroup to this as it is modern too
	echo outgroup >> modernSampleNames.txt
	echo -e "Standardise fasta suffixes"

	panaroo_fasta_suffix() {
	fasta_file=\$1
	filename=\$(basename "\${fasta_file%.aln.fas}")

		mv "\${fasta_file}" panaroo_genes/"\${filename}.fasta"
	}
	export -f panaroo_fasta_suffix
	find panaroo_genes/ -name "*.aln.fas" | parallel -j 10 panaroo_fasta_suffix

	echo -e "Done\n"

	echo -e "Fixing FASTA headers and sequences formatting with seqtk in existing gene alignments"

	mkdir -p panaroo_parsed
	parsing_panaroo_msa() {
	fasta_file=\$1
                name=\$(basename "\${fasta_file%.fasta}" | sed -e 's/~/_/g') # Also replace  ~ characters with _, with double % to remove two dots
                seqtk seq "\${fasta_file}" | awk '/^>/ {sub(/;.*/, "", \$0)} {print}' > panaroo_parsed/"\${name}_parsing_panaroo.fasta"
	}
	export -f parsing_panaroo_msa
	find panaroo_genes/ -name "*.fasta" | parallel -j 30 parsing_panaroo_msa

	echo -e "Done\n"



	mkdir -p blacklisted
	echo -e "Checking if there are low quality samples"
	if [[ -s blackListed.txt ]]; then
	        while read -r removeMe; do
        	        mv user_genes/*"\${removeMe}.fasta" blacklisted/"\${removeMe}"
                	echo "\${removeMe} has been removed from analysis due to low quality."
        	done < blackListed.txt
	else
	        echo -e "Every sample passed quality checks."
	fi

	#index_and_formatting() send user samplenames to userSampleNames.txt and replaces ~ from gene names with _
	index_and_formatting() {
	sample=\$1
        	sampleName=\$(basename "\${sample%.fasta}")
        	echo "\${sampleName}" >> userSampleNames.txt
		sed -i -e 's/~/_/g' "\${sample}"
	}
	export -f index_and_formatting
	find user_genes/ -name "*.fasta" | parallel -j 30 index_and_formatting

	echo -e "Done\n"


	# Make a file with every sample combined
	cat modernSampleNames.txt userSampleNames.txt > sampleNames.txt


	echo -e "Adding user sample genes sequences to each particular gene MSA and replace gene name with sample name"

	#add_user_sample_sequences() adds user sample gene sequences into each panaroo msa and replaces gene name with user sample name
	add_user_sample_sequences() {
	fasta_file=\$1
	        if [[ -e "\$fasta_file" ]]; then
                	name=\$(basename "\${fasta_file%_parsing_panaroo.fasta}")

                	while read -r sampleName; do
                                grep -w -A 1 "\$name" "user_genes/\${sampleName}.fasta" | awk -v newHeader="\$sampleName" '/^>/ {sub(/^>.*/, ">" newHeader, \$0)} {print}' >> "\${fasta_file}"
                	done < userSampleNames.txt
        	else
                	echo -e "There are no files with .fasta extension in panaroo_parsed/ folder. Stopping the process.\n"
			exit 1
        	fi
	}

	export -f add_user_sample_sequences
	find panaroo_parsed/ -name "*.fasta" | parallel -j 30 add_user_sample_sequences
	echo -e "Done\n"


	echo -e "Adding outgroup gene sequences into each gene MSA file"
	#similarly as add_user_sample_sequences(), add_outgroup_sequences() will add outgroup aligned gene sequences into each gene panaroo msa file
	add_outgroup_sequences() {
	fasta_file=\$1

                name=\$(basename "\${fasta_file%_parsing_panaroo.fasta}")
                grep -w -A 1 "\$name" outgroup | awk -v outgroup="outgroup" '/^>/ {sub(/^>.*/, ">" outgroup, \$0)} {print}' >> "\${fasta_file}"
	}
	export -f add_outgroup_sequences
	find panaroo_parsed/ -name "*_parsing_panaroo.fasta" | parallel -j 30 add_outgroup_sequences

	echo -e "Done\n"


	echo -e "Padding incomplete sequences to the right"

	panaroo_parsed_padding() {
	MSA=\$1
        	name=\$(basename "\${MSA}")
        	seqLength=\$(awk 'NR==2 {print length}' "\$MSA")

        	awk -v seq_length="\$seqLength" '{
	                if (\$0 ~ /^>/) {
                        	print
                	} else {
                	while ( length(\$0) < seq_length) {
	                        \$0 = \$0 "n"
                	}
                	print
                	}
        	}' "\${MSA}" > tmp_"\${name}" && mv tmp_"\${name}" "\${MSA}"
	}
	export -f panaroo_parsed_padding
	find panaroo_parsed/ -name "*parsing_panaroo.fasta" | parallel -j 30 panaroo_parsed_padding

	echo -e "Done\n"



	echo -e "Checking if there are Panaroo headers artifacts"

	check_panaroo_headers_artifacts() {
	file_msa=\$1

	        while read -r sampleName; do
        	        header_name=">\$sampleName"
			matched_header=\$(grep "\$sampleName" "\$file_msa")

                	if [[ -n "\$matched_header" ]]; then

                        	while IFS= read -r matched; do

                                	if [[ "\$matched" != "\$header_name" ]]; then
                                        	sed -i -e "s/\${matched}/\${header_name}/g" "\$file_msa"
                                        	echo -e "\$matched name was found in \$file_msa but it should have been \$header_name instead. Fixed"
                                	fi
                        	done <<< "\$matched_header"
                	fi

        	done < modernSampleNames.txt
	}
	export -f check_panaroo_headers_artifacts
	find panaroo_parsed/ -name "*_parsing_panaroo.fasta" | parallel -j 30 check_panaroo_headers_artifacts
	echo -e "Done\n"



	#modern_append_gaps() will check if there are modern strains missing in panaroo_parsed alignments. If yes, append the sample and padding with -
	echo -e "Padding missing samples. Gaps for modern strains and n for ancient samples"

	modern_append_gaps() {
	fasta_file=\$1

		seq_length=\$(awk 'NR==2 {print length}' "\$fasta_file")

                while read -r strain; do
                        if ! grep -wq "\$strain" "\$fasta_file"; then
                                echo ">\$strain" >> "\$fasta_file"
                                gaps_length=\$(printf '%*s' "\$((seq_length))" | tr ' ' '-')
                                echo "\$gaps_length" >> "\$fasta_file"
                        fi
                done < modernSampleNames.txt
	}
	export -f modern_append_gaps
	find panaroo_parsed/ -name "_parsing_panaroo.fasta" | parallel -j 30 modern_append_gaps



	#ancient_append_missingness() will check if there are ancient strains missing in panaroo_parsed alignments. If yes, append the sample and fill it with n as we don't know the ancient status for that gene.
	ancient_append_missingness() {
	fasta_file=\$1

	        seq_length=\$(awk 'NR==2 {print length}' "\$fasta_file")

		while read -r strain; do
                	if ! grep -wq "\$strain" "\$fasta_file"; then
                        	echo ">\$strain" >> "\$fasta_file"
                                fakeSeq=\$(printf '%*s' "\$((seq_length))" | tr ' ' 'n')
                                echo "\$fakeSeq" >> "\$fasta_file"
                        fi
                done < sampleNames.txt
	}

	export -f ancient_append_missingness
	find panaroo_parsed/ -name "*_parsing_panaroo.fasta" | parallel -j 30 ancient_append_missingness

	echo -e "Done \n"


	# At this point there will be many msa that are finished: number of headers == number of samples AND headers names == all samples names. I need to find them based on this logic and send them
	# to the final folder. There will be some msa with header duplication problems and they need to be send to special_cases.

	mkdir -p sanitised_msa

	sanitised_msa() {
	msa_file=\$1

		name=\$(basename "\${msa_file%_parsing_panaroo.fasta}.fasta")
	        sample_count=\$(grep -c '^>' "\$msa_file")
	        total_sample_count=\$(wc -l < sampleNames.txt)

		if [[ "\$sample_count" -ne "\$total_sample_count" ]]; then  #if they are equal then keep going. if not skip current file.
			return
		fi

		missing_samples=0
		while read -r strain; do
			if ! grep -wq "\$strain" "\$msa_file"; then
				missing_samples=1
				break
			fi
		done < sampleNames.txt

		if [[ "\$missing_samples" -eq 0 ]]; then
			mv "\$msa_file" "sanitised_msa/\${name}"
		fi

	}
	export -f sanitised_msa
	find panaroo_parsed/ -name "*_parsing_panaroo.fasta" | parallel -j 30 sanitised_msa

	# Dealing with fragmented Panaroo gene alignments multi-entries
	mkdir -p special_cases
	
	echo -e "Checking if there are leftovers after sanitising"

	unsanitised=\$(ls panaroo_parsed/*fasta | wc -l )
	if (( unsanitised != 0 )); then
		echo -e "There are \$unsanitised samples that have problems in panaroo_parsed/ directory. Inspecting them. ."
	else
		echo -e "It seems every gene msa was succesfully cleaned. Moving on."
	fi

	fix_duplicated_entries() {
	fasta_file=\$1

               	if [[ -e "\$fasta_file" ]] ; then
                       	name=\$(basename "\${fasta_file%_parsing_panaroo.fasta}")
               		# Identifying repeated headers and save them to .dup
                       	awk '/^>/ {count[\$0]++} END {for (header in count) if (count[header] > 1) print substr(header, 2)}' "\$fasta_file" > special_cases/"\${name}.dup"
	
               		# Extracting sequences for each repeated entry into temporary files
                	while read -r entry; do
                        	awk -v sampleName="\$entry" '
                        	\$0 ~ "^>" sampleName {print_header=1; next}
                        	/^>/ {print_header=0}
                        	print_header {print}' "\$fasta_file" > special_cases/"\${name}_duplicated_\${entry}"
                       done < special_cases/"\${name}.dup"
               # Remove repeated entries and their sequences from the original FASTA file
                       awk -v dpdFile=special_cases/"\${name}.dup" '
                       BEGIN {
                       while (getline < dpdFile) {
                               repeated[\$0] = 1
                               }
                       }
                       /^>/ { header = substr(\$0, 2)
                       if (repeated[header]) {
                               skip = 1
                       } else {
                               skip = 0}
                       }
                       !skip' "\$fasta_file" > special_cases/"\${name}.fasta"

                       for indexSeqs in special_cases/"\${name}_duplicated_"*; do
                               geneName=\$(basename "\${indexSeqs%_duplicated_*}")
                               sampleName=\$(basename "\${indexSeqs##*_duplicated_}")

                               # Read the longest line based on letter count
                               sequence=\$(awk '
                               {gsub(/[^a-zA-Z]/, "", \$0); len=length(\$0)}
                               len > max_length {max_length=len; longest=\$0}
                               END {print longest}' "\$indexSeqs")

                       # Finally add the selected sequence back to the cleaned original FASTA file with the header as well
                               echo ">\${sampleName}" >> special_cases/"\${name}.fasta"
                               echo "\$sequence" >> special_cases/"\${name}.fasta"
                       done

               else
                       echo -e "No files found with fasta extension in special_cases folder"    
       		# Cleaning temporary files
               fi
       }

	export -f fix_duplicated_entries

	if (( unsanitised != 0 )); then
		find panaroo_parsed/ -name "*_parsing_panaroo.fasta" | parallel -j 30 fix_duplicated_entries
	fi

	#finally checking that every seq has the same length

	checking_alignment_lengths() {
	msa_file=\$1
	name=\$(basename "\${msa_file}")

               seq_length=\$(awk 'NR%2 == 0 && length > max { max = length } END { print max }' "\$cleanedFasta")

               awk -v numCols="\$seq_length" '{
                       if (\$0 ~ /^>/) {
                               print
                       } else {
                       while ( length(\$0) < numCols) {
                               \$0 = \$0 "N"
                       }
                       print
                       }
               }' "\$msa_file" > special_cases/tmp_"\${name}" && mv special_cases/tmp_"\${name}" "\${msa_file}"
       }

	export -f checking_alignment_lengths
	find special_cases/ -name "*.fasta" | parallel -j 30 checking_alignment_lengths
	echo -e "Done"


	#now check if these cleaned fasta passed the sanitised_msa() test.
	sanitised_msa_special_cases() {
	msa_file=\$1

        	name=\$(basename "\${msa_file}")
        	sample_count=\$(grep -c '^>' "\$msa_file")
        	total_sample_count=\$(wc -l < sampleNames.txt)

        	if [[ "\$sample_count" -ne "\$total_sample_count" ]]; then  #if they are equal then keep going. if not skip current file.
	                return
	        fi

	        missing_samples=0
	        while read -r strain; do
                	if ! grep -wq "\$strain" "\$msa_file"; then
                        	missing_samples=1
                        	break
                	fi
        	done < sampleNames.txt

        	if [[ "\$missing_samples" -eq 0 ]]; then
	                mv "\$msa_file" "sanitised_msa/\${name}"
	        fi

	}
	export -f sanitised_msa_special_cases
	find special_cases/ -name "*.fasta" | parallel -j 30 sanitised_msa_special_cases

	# final steps: sort everything.

	# Make INDEX file so we can sort every file in the same order before building MSA based on first file.
	firstFile=\$(ls -1 sanitised_msa/ | awk 'NR==1 {print \$0}') && awk '/^>/ {print \$0}' sanitised_msa/"\$firstFile" > sorting_index

	mkdir -p sorted

	echo -e "Sorting MSAs"

	sortingMSA() {
	fasta=\$1

	        name=\$(basename "\${fasta%.fasta}")

	        while read -r header; do
                	awk -v headerName="\$header" '
	                        \$0 ~ headerName {
                        	print \$0     # Print the matched header
                        	getline      # Get the sequence line after the header and store it in line
                        	print \$0     # Print the sequene
                        	}
                	' "\$fasta" >> sorted/"\${name}.fasta"
        	done < sorting_index
	}

	export -f sortingMSA
	find sanitised_msa/ -name "*.fasta" | parallel -j 30 sortingMSA

	echo -e "Done"

	cat .command.out >> filterGeneAlignments.log
	"""
}

/*   Add individual gene sequences to these particular gene.aln files
 *   To do this we need: If there is one or more samples missing in any particular gene.aln file, then we get the gene lenght and we add Ns.
 *   At the end , we will have every gene.aln file filled with every sample. So we can now concatenate every gene.aln file based on INDEX.
 */



process makeMSA {

	input:
	path genesAlnSeq, stageAs: 'genes/*'
	path maskedMatrixGenesNoUbiquitous, stageAs: 'maskedMatrixGenesNoUbiquitous.txt'
	path maskedMatrixGenesOnlyAncient, stageAs: 'maskedMatrixGenesOnlyAncient.txt'
	path maskedMatrixGenesUbiquitous, stageAs: 'maskedMatrixGenesUbiquitous.txt'
	path genesAbovePercentSeries, stageAs: 'genesAbovePercentSeries.txt'
	path sampleNames, stageAs: 'sampleNames.txt'

	output:
	path 'genesAbovePercentMSA.fasta', emit: genesAbovePercentMSA
	path 'maskedMatrixGenesNoUbiquitousMSA.fasta', emit: maskedMatrixGenesNoUbiquitousMSA
	path 'maskedMatrixGenesOnlyAncientMSA.fasta', emit: maskedMatrixGenesOnlyAncientMSA
	path 'maskedMatrixGenesUbiquitousMSA.fasta', emit: maskedMatrixGenesUbiquitousMSA
	

	script:
	"""
	#!/bin/bash

	touch genesAbovePercentMSA.fasta
	touch maskedMatrixGenesUbiquitousMSA.fasta
	touch maskedMatrixGenesOnlyAncientMSA.fasta
	touch maskedMatrixGenesNoUbiquitousMSA.fasta

	sed -i -e 's/~/_/g' genesAbovePercentSeries.txt
	sed -i -e 's/~/_/g' maskedMatrixGenesNoUbiquitous.txt
	sed -i -e 's/~/_/g' maskedMatrixGenesOnlyAncient.txt
	sed -i -e 's/~/_/g' maskedMatrixGenesUbiquitous.txt



	while read -r gene; do
		file="genes/\${gene}_Filtered.fasta"
		if [[ -f "\$file" ]] ; then
			paste genesAbovePercentMSA.fasta "\$file" > TMP; mv TMP genesAbovePercentMSA.fasta
		else
			echo "Warning: File \$file not found, skipping..."
		fi
	done < genesAbovePercentSeries.txt
	sed -i -e 's/\t//g' genesAbovePercentMSA.fasta
	awk -F'>' '/^>/ {print ">" \$2} !/^>/' genesAbovePercentMSA.fasta > TMPg; mv TMPg genesAbovePercentMSA.fasta



	while read -r gene; do
		file="genes/\${gene}_Filtered.fasta"
		if [[ -f "\$file" ]] ; then
			paste maskedMatrixGenesNoUbiquitousMSA.fasta "\$file" > TMP; mv TMP maskedMatrixGenesNoUbiquitousMSA.fasta
		else
			echo "Warning: File \$file not found, skipping..."
		fi
	done < maskedMatrixGenesNoUbiquitous.txt
        sed -i -e 's/\t//g' maskedMatrixGenesNoUbiquitousMSA.fasta
        awk -F'>' '/^>/ {print ">" \$2} !/^>/' maskedMatrixGenesNoUbiquitousMSA.fasta > TMPnU; mv TMPnU maskedMatrixGenesNoUbiquitousMSA.fasta



	while read -r gene; do
		file="genes/\${gene}_Filtered.fasta"
		if [[ -f "\$file" ]] ; then
			paste maskedMatrixGenesOnlyAncientMSA.fasta "\$file" > TMP; mv TMP maskedMatrixGenesOnlyAncientMSA.fasta
		else
			echo "Warning: File \$file not found, skipping..."
		fi
	done < maskedMatrixGenesOnlyAncient.txt
        sed -i -e 's/\t//g' maskedMatrixGenesOnlyAncientMSA.fasta
        awk -F'>' '/^>/ {print ">" \$2} !/^>/' maskedMatrixGenesOnlyAncientMSA.fasta > TMPo2; mv TMPo2 maskedMatrixGenesOnlyAncientMSA.fasta

	while read -r gene; do
		file="genes/\${gene}_Filtered.fasta"
		if [[ -f "\$file" ]] ; then
			paste maskedMatrixGenesUbiquitousMSA.fasta "\$file" > TMP; mv TMP maskedMatrixGenesUbiquitousMSA.fasta
		else
			echo "Warning: File \$file not found, skipping..."
		fi
	done < maskedMatrixGenesUbiquitous.txt

        sed -i -e 's/\t//g' maskedMatrixGenesUbiquitousMSA.fasta
        awk -F'>' '/^>/ {print ">" \$2} !/^>/' maskedMatrixGenesUbiquitousMSA.fasta > TMPU2; mv TMPU2 maskedMatrixGenesUbiquitousMSA.fasta

	"""
}


process pMauve {
	conda "${projectDir}/envs/pMauve.yaml"

	input:
	path gffFiles, stageAs: '*'

	output:
	path 'pMauveAlignment.xmfa', emit: pMauveAlignment
	path 'pMauveAlignment.bbcols', emit: pMauveBbcols
	path 'pMauveAlignment.backbone', emit: pMauveBackbone
	path 'pMauveAlignmentCoreGenome', emit: pMauveCoreGenome

	script:
	"""
	progressiveMauve  --output=pMauveAlignment *
	mv pMauveAlignment ./pMauveAlignment.xmfa
	stripSubsetLCBs pMauveAlignment.xmfa pMauveAlignment.bbcols pMauveAlignmentCoreGenome 500
	"""
}

process xmfaToFasta {
	conda "${projectDir}/envs/biopython.yaml"

	input:
	path coreGenome, stageAs: 'pMauveAlignmentCoreGenome.xmfa'

	output:
	path 'pMauveFastaMSA.fasta', emit: pMauveFastaMSA

	script:
	"""
	convertXmfaToFasta.py pMauveAlignmentCoreGenome.xmfa
	"""
}

process filterMauveFasta {
	conda "${projectDir}/envs/seqtk.yaml"

	input:
	path mauveFastaMSA, stageAs: 'pMauveFastaMSA.fasta'

	output:
	path 'concatenatedSeqtkMauveFastaMSA.fasta', emit: concatenatedSeqtkMauveFastaMSA

	script:
	"""
	seqtk seq pMauveFastaMSA.fasta > seqtkMauveFastaMSA.fasta
	sed -i -e '/^>/ s/\\/.*//' seqtkMauveFastaMSA.fasta 

	awk '
	{
	    if (/^>/) {
	        header = \$1  
	        gsub(/^>/, "", header)  
	    } else {
	        seq[header] = seq[header] \$0  
	    }
	}
	END {
	    for (id in seq) {
	        print ">" id  
	        print seq[id] 
	    }
	}' seqtkMauveFastaMSA.fasta > concatenatedSeqtkMauveFastaMSA.fasta
	
	"""
}

process startingTree {
	conda "${projectDir}/envs/iqtree.yaml"	

	input:
	path concatenatedSeqtkMauveFastaMSA, stageAs: 'concatenatedSeqtkMauveFastaMSA.fasta'

	output:
	path 'startingTreeMauveFasta.treefile', emit: startingTreeMauveFasta
	path 'startingTreeMauveFasta.iqtree', emit: startingTreeMauveFastaLog
	path 'kappaValue', emit: kappa

	script:
	"""
	iqtree -s concatenatedSeqtkMauveFastaMSA.fasta --prefix startingTreeMauveFasta -T 10 -B 1000 -m MFP
	
	awk -F':' '
		/A-C:/ {ACtransversion=\$2 + 0}
		/A-G:/ {AGtransition=\$2 + 0}
		/A-T:/ {ATtransversion=\$2 + 0}
		/C-T:/ {CTtransition=\$2 + 0}
		/C-G:/ {CGtransversion=\$2 + 0}
		/G-T:/ {GTtransversion=\$2 + 0}
		END {
		transitionRate= ((AGtransition + CTtransition)/2)
		transversionRate= (( ACtransversion + ATtransversion + CGtransversion + GTtransversion ) / 4)
		kappa = (transitionRate / transversionRate  )		

		print kappa }' startingTreeMauveFasta.iqtree > kappaValue

	"""
}


process findRecombinationSpots {
	conda "${projectDir}/envs/clonalframe.yaml" 

	input:
	path MSA, stageAs: 'concatenatedSeqtkMauveFastaMSA.fasta'
	path startingTree, stageAs: 'startingTreeMauveFasta.treefile'
	path kappa, stageAs: 'kappaValue'
	

	output:
	path 'recombinantOutputs.importation_status.txt', emit: recombinationMap
	
	script:
	"""
	kappa=\$(cat kappaValue)
	ClonalFrameML startingTreeMauveFasta.treefile concatenatedSeqtkMauveFastaMSA.fasta recombinantOutputs -kappa "\${kappa}"

	"""
}

process mapRecombinantsToGenes {
	conda "${projectDir}/envs/blast.yaml"

	input:
	path recombinationMap, stageAs: 'recombinantOutputs.importation_status.txt'
	path MSA, stageAs: 'concatenatedSeqtkMauveFastaMSA.fasta'
	path db, stageAs: 'database/*'
	path GFF, stageAs: 'gff/*'

	output:
	stdout

	script:
	"""
	#!/bin/bash

	awk '!/^NODE/ && NR>1 {print \$0}' recombinantOutputs.importation_status.txt > filteredRecombinationMap

	while read -r name beg end; do
	    awk -v name="\$name" -v start="\$beg" -v end="\$end" '
	    /^>/ {
	        # Process headers
	        headerSequence = substr(\$0, 2)  # remove the ">" thing to get the header name
	        isHeader = (headerSequence == name)  # check if the header matches the target name
	
	        if (isHeader) {
	            fullSeq = ""  # Reset fullSeq for the new header
	            currentHeader = name "_seq"  # make unique header for this sequence
	        }
	    }
	    !/^>/ && isHeader {
	        fullSeq = fullSeq \$0   
	    }
	    /^>/ && fullSeq != "" {
	        # store fullSeq for each unique header when a new header starts
	        # remove sequences shorter than 30 bp? To avoid mapping uncertainty
	        seqTesting = substr(fullSeq, start, end - start + 1)
	        if (length(seqTesting) >= 30) {
	            seqArray[currentHeader] = seqTesting
        	}
	        fullSeq = ""
	    }
	    END {
	        # this is for the last entry
	        if (fullSeq != "") {
	            seqTesting = substr(fullSeq, start, end - start + 1)
	            if (length(seqTesting) >= 30) {
	                seqArray[currentHeader] = seqTesting
	            }
	        }
	        
	        # print the sequences
	        for (header in seqArray) {
	            printf(">%s\\n%s\\n", header, seqArray[header])
	        }
	    }
	    ' concatenatedSeqtkMauveFastaMSA.fasta >> "\${name%.fna}"_TMP.fasta
	done < filteredRecombinationMap
	
	for i in *_TMP.fasta; do
		name=\$(basename "\${i%_TMP.fasta}")
		awk '
		    BEGIN { count = 0 }  
		    /^>/ { 
		        count++  
		        \$0 = \$0 "_" count  # this just add a counter ID for each header to make them unique
		        print
		    } 
		    !/^>/ { 
		        print  # print the sequence line only
		    }
		' "\$i"  > "\$name"_newfileWithExtractedHeadersAndSequences.fasta
	done

	rm *TMP.fasta
	
	# Mapping sequences to PanGenomeReference

	for sample in *newfileWithExtractedHeadersAndSequences.fasta; do
		name=\$(basename "\${sample%_newfileWithExtractedHeadersAndSequences.fasta}")
		blastn -query "\$sample" -db database/panGenomeReferenceDB -out blastResults"\$name" -outfmt 6
	done	


	echo -e "qseqid\\tsseqid\\tpident\\tlength\\tmismatch\\tgapopen\\tqstart\\tqend\\tsstart\\tsend\\tevalue\\tbitscore\\tseqLength" > blastSummaryResults.tab

	for sample in *newfileWithExtractedHeadersAndSequences.fasta; do
    		awk '/^>/ {name = \$0 ; getline ; seqLength = length(\$0); print name, seqLength}' "\$sample" >> recombinantsIDplusLengths.txt
	done

	sed -i -e 's/>//g' recombinantsIDplusLengths.txt
	
	cat blastResults* >> TMP1.tab
	sed -i -e 's/~/_/g' TMP1.tab

	while read -r sample seqLength; do
    		grep -w "\$sample" TMP1.tab | awk -v value="\$seqLength" '{print \$0, value}' >> blastSummaryResults.tab
	done < recombinantsIDplusLengths.txt

	sed -i -e 's/ /\\t/g' blastSummaryResults.tab

	rm TMP1.tab
	"""
}



process treeThreshold {
	conda "${projectDir}/envs/iqtree.yaml"
	
	input:
	path genesMSA, stageAs: 'genesAbovePercentMSA.fasta'

	output:
	path 'genesAbovePercentMSA.contree', emit: genesAbovePercentMSAContree
	path 'genesAbovePercentMSA.iqtree', emit: genesAbovePercentMSAIqtree
	path 'genesAbovePercentMSA.log', emit: genesAbovePercentMSALog
	path 'genesAbovePercentMSA.treefile', emit: genesAbovePercentMSATreefile

	script:
	"""
	iqtree -s genesAbovePercentMSA.fasta --prefix genesAbovePercentMSA -T 10 -B 1000 -m MFP 
	"""
}

process treeUbiquitous {
        conda "${projectDir}/envs/iqtree.yaml"
        
        input:
        path ubiquitousMSA, stageAs: 'maskedMatrixGenesUbiquitousMSA.fasta'

        output:
        path 'maskedMatrixGenesUbiquitousMSA.contree', emit: maskedMatrixGenesUbiquitousMSAContree
	path 'maskedMatrixGenesUbiquitousMSA.iqtree', emit: maskedMatrixGenesUbiquitousMSAIqtree
	path 'maskedMatrixGenesUbiquitousMSA.log', emit: maskedMatrixGenesUbiquitousMSALog
	path 'maskedMatrixGenesUbiquitousMSA.treefile', emit: maskedMatrixGenesUbiquitousMSATreefile

        script:
        """
        iqtree -s maskedMatrixGenesUbiquitousMSA.fasta --prefix maskedMatrixGenesUbiquitousMSA -T 10 -B 1000 -m MFP
        """
}

process treeNoUbiquitous {
        conda "${projectDir}/envs/iqtree.yaml"
        
        input:
        path noUbiquitousMSA, stageAs: 'maskedMatrixGenesNoUbiquitousMSA.fasta'

        output:
        path 'maskedMatrixGenesNoUbiquitousMSA.contree', emit: maskedMatrixGenesNoUbiquitousMSAContree
	path 'maskedMatrixGenesNoUbiquitousMSA.iqtree', emit: maskedMatrixGenesNoUbiquitousMSAIqtree
	path 'maskedMatrixGenesNoUbiquitousMSA.log', emit: maskedMatrixGenesNoUbiquitousMSALog
	path 'maskedMatrixGenesNoUbiquitousMSA.treefile', emit: maskedMatrixGenesNoUbiquitousMSATreefile

        script:
        """
        iqtree -s maskedMatrixGenesNoUbiquitousMSA.fasta --prefix maskedMatrixGenesNoUbiquitousMSA -T 10 -B 1000 -m MFP
        """
}


process treeAncient {
        conda "${projectDir}/envs/iqtree.yaml"
        
        input:
        path ancientMSA, stageAs: 'maskedMatrixGenesOnlyAncientMSA.fasta'

        output:
        path 'maskedMatrixGenesOnlyAncientMSA.contree', emit: maskedMatrixGenesOnlyAncientMSAConqtree
	path 'maskedMatrixGenesOnlyAncientMSA.iqtree', emit: maskedMatrixGenesOnlyAncientMSAIqtree
	path 'maskedMatrixGenesOnlyAncientMSA.log', emit: maskedMatrixGenesOnlyAncientMSALog
	path 'maskedMatrixGenesOnlyAncientMSA.treefile', emit: maskedMatrixGenesOnlyAncientMSATreefile

        script:
        """
        iqtree -s maskedMatrixGenesOnlyAncientMSA.fasta --prefix maskedMatrixGenesOnlyAncientMSA -T 10 -B 1000 -m MFP
        """
}

process outgroupEntrez {
	conda "${projectDir}/envs/entrez.yaml"

	input:
	val outID
	
	output:
	path '*fna', emit: outgroupFasta

	script:
	"""
	#!/bin/bash
	counter=0
	esearch -db assembly -query "txid${outID}[Organism] AND (latest[filter] AND (complete genome[filter] OR chromosome level[filter]))" | esummary | xtract -pattern DocumentSummary -element FtpPath_RefSeq | while read url; do
        
		if [ "\$counter" -ge 1 ]; then
			break
		fi
	
		if [ -z "\$url" ]; then
			continue
		fi
        
		fname="\$(basename "\$url")"
		downloadMeAndCheck() {
			wget "\$url/\${fname}_genomic.gbff.gz"
	        	wget "\$url/\${fname}_genomic.fna.gz"

			gunzip -f "\${fname}_genomic.gbff.gz"
			gunzip -f "\${fname}_genomic.fna.gz"


			if [ -f "\${fname}_genomic.gbff.gz" ] || [ -f "\${fname}_genomic.fna.gz" ]; then	
				rm -f "\${fname}_genomic.gbff.gz" "\${fname}_genomic.fna.gz"
				return 1
			fi

			return 0
		}

		while ! downloadMeAndCheck; do
			echo "Files were corrupted. Retrying"
			sleep 3
		done
		# Adding to the counter
	        counter="\$((counter + 1))"     
                fna_count=\$(ls -1 *.fna 2>/dev/null | wc -l)
                gbff_count=\$(ls -1 *.gbff 2>/dev/null | wc -l)

                # check that there is already one file
                if [ "\$fna_count" -lt 1 ] || [ "\$gbff_count" -lt 1 ]; then
                        echo "It seems a file was already downloaded. Moving on"
                        break
                else
                        continue
                fi

	done
	"""
}


process makeReads {
	conda "${projectDir}/envs/art.yaml"

	input:
	path outgroupFasta, stageAs: 'outgroupFasta.fasta'
	
	output:
	path 'outgroupReads.fq', emit: outgroupReads


	script:
	"""
	art_illumina -i outgroupFasta.fasta -l 150 -f 100 -o outgroupReads
	"""
}


process outgroupAlignmentFAndiltering {
	conda "${projectDir}/envs/alignment.yaml"
	
	input:
	path outgroupReads, stageAs: 'outgroupReads.fq'
	path panGenomeRef, stageAs: 'panGenomeReferenceSeq.fasta'
	val threadsGlobal

	output:
	path 'outgroupFastaPostAlignment.bam', emit: outgroupFastaPostAlignment

	script:
	"""
	bwa index panGenomeReferenceSeq.fasta
	bwa mem -B 1 -E 1 panGenomeReferenceSeq.fasta outgroupReads.fq -t $threadsGlobal > outgroupFasta.sam
	samtools view -bS outgroupFasta.sam > outgroupFasta.bam
	samtools quickcheck outgroupFasta.bam
	samtools sort -o outgroupFastaSorted.bam -O bam -@ $threadsGlobal outgroupFasta.bam
	samtools index outgroupFastaSorted.bam
	samtools view -b -@ 10 -F 4 outgroupFastaSorted.bam > outgroupFastaSortedMappedreads.bam
	samtools index outgroupFastaSortedMappedreads.bam
	samtools sort -o outgroupFastaPostAlignment.bam -O bam -@ $threadsGlobal outgroupFastaSortedMappedreads.bam
	"""
}


process makeOutgroupConsensus {
	conda "${projectDir}/envs/consensusOutgroup.yaml"

	input:
	path outgroupFastaPostAlignment, stageAs: 'outgroupFastaPostAlignment.bam'
	path panGenomeRef, stageAs: 'panGenomeRef.fasta'

	output:
	path 'extractedSequencesOutgroup.fasta', emit: extractedSequencesOutgroupFasta
	path 'extractedSequencesOutgroup.fq', emit: extractedSequencesOutgroupFastq

	script:
	"""
	bcftools mpileup -f panGenomeRef.fasta outgroupFastaPostAlignment.bam | bcftools call -c | vcfutils.pl vcf2fq > extractedSequencesOutgroup.fq
	seqtk seq -a extractedSequencesOutgroup.fq > extractedSequencesOutgroup.fasta

	"""
}


process blastMe {
	conda "${projectDir}/envs/blast.yaml"

	input:
	path panSeq, stageAs: 'panGenomeReference.fasta'

	output:
	path 'panGenomeReferenceDB*', emit: panGenomeReferenceDB

	script:
	"""
	makeblastdb -in panGenomeReference.fasta -dbtype nucl -out panGenomeReferenceDB
	"""
}


process getResults {

	input:
	path makeDir
	path checkedFastas, stageAs: 'checkedFasta/*'
	path checkedGffs, stageAs: 'checkedGff/*'
	path fastaDatabaseLog, stageAs: 'fastaDatabase.log'
	path fastaDatabaseSeqs, stageAs: 'clusteredSequences.fasta'
	path clusteredSequences, stageAs: 'clusteredNonRedundantGenes.fasta'
	path clusteringLogFile, stageAs: 'clustering.log'
	path prokkaGff, stageAs: 'prokkaGff/*'
	path prokkaLog, stageAs: 'prokka.log'
	path panarooLog, stageAs: 'makePangenome.log'
	path genesMSA, stageAs: 'geneMSA/*'
	path panrefG, stageAs: 'pangenomeReferenceGenome.fasta'
	path geneNormalizedUpdated, stageAs: 'geneNormalizedUpdated.tab'
	path globalMeanCoverage, stageAs: 'globalMeanCoverage.tab'
	path PostAlignmentFiles, stageAs: '*'
	path RefLenghts, stageAs: '*'
	path RawCoverage, stageAs: '*'
	path CompletenessSummary, stageAs: 'completenessSummary.tab'
	path FinalMatrix, stageAs: 'finalMatrix.tab'
	path presenceAbsenceplots, stageAs: '*'
	path MaskedMatrixGenesOnlyAncient, stageAs: 'maskedMatrixGenesOnlyAncient.txt'
	path MaskedMatrixGenesUbiquitous, stageAs: 'maskedMatrixGenesUbiquitous.txt'
	path MaskedMatrixGenesNoUbiquitous, stageAs: 'maskedMatrixGenesNoUbiquitous.txt'
	path GenesAbovePercentSeries, stageAs: 'genesAbovePercentSeries.txt'
	path GenesAbovePercentMSAIqtree, stageAs: 'genesAbovePercentMSA.iqtree'
	path GenesAbovePercentMSALog, stageAs: 'genesAbovePercentMSA.log'
	path GenesAbovePercentMSATreefile, stageAs: 'genesAbovePercentMSA.treefile'
	path MaskedMatrixGenesUbiquitousMSAIqtree, stageAs: 'maskedMatrixGenesUbiquitousMSA.iqtree'
	path MaskedMatrixGenesUbiquitousMSALog, stageAs: 'maskedMatrixGenesUbiquitousMSA.log'
	path MaskedMatrixGenesUbiquitousMSATreefile, stageAs: 'maskedMatrixGenesUbiquitousMSA.treefile'
	path MaskedMatrixGenesNoUbiquitousMSAIqtree, stageAs: 'maskedMatrixGenesNoUbiquitousMSA.iqtree'
	path MaskedMatrixGenesNoUbiquitousMSALog, stageAs: 'maskedMatrixGenesNoUbiquitousMSA.log'
	path MaskedMatrixGenesNoUbiquitousMSATreefile, stageAs: 'maskedMatrixGenesNoUbiquitousMSA.treefile'
	path MaskedMatrixGenesOnlyAncientMSAIqtree, stageAs: 'maskedMatrixGenesOnlyAncientMSA.iqtree'
	path MaskedMatrixGenesOnlyAncientMSALog, stageAs: 'maskedMatrixGenesOnlyAncientMSA.log'
	path MaskedMatrixGenesOnlyAncientMSATreefile, stageAs: 'maskedMatrixGenesOnlyAncientMSA.treefile'

	output:
	stdout

	script:
	"""
	#!/bin/bash

	mkdir -p "${makeDir}/modernData/"

	mv checkedFasta/* "${makeDir}/modernData/"
	mv checkedGff/* "${makeDir}/modernData/"
	mv fastaDatabase.log "${makeDir}/modernData/"

	mkdir -p "${makeDir}/clusteredSequences"
	mv clusteredSequences.fasta "${makeDir}/clusteredSequences/"
	mv clusteredNonRedundantGenes.fasta "${makeDir}/clusteredSequences/"
	mv clustering.log "${makeDir}/clusteredSequences/"

	mkdir -p "${makeDir}/prokkaResults/"

	mv prokkaGff/* "${makeDir}/prokkaResults/"
	mv prokka.log "${makeDir}/prokkaResults/"

	mkdir -p "${makeDir}/pangenomeFiles"

	mv makePangenome.log "${makeDir}/pangenomeFiles/"
	mv geneMSA/* "${makeDir}/pangenomeFiles/"
	mv pangenomeReferenceGenome.fasta "${makeDir}/pangenomeFiles/"

	mkdir -p "${makeDir}/alignmentResults"

	mv geneNormalizedUpdated.tab "${makeDir}/alignmentResults/"
	mv globalMeanCoverage.tab "${makeDir}/alignmentResults/"
	mv *_refLength.txt "${makeDir}/alignmentResults/"
	mv *Coverage.txt "${makeDir}/alignmentResults/"
	mv *bam "${makeDir}/alignmentResults/"
	mv completenessSummary.tab "${makeDir}/alignmentResults/"

	mkdir -p "${makeDir}/matrixResults"
	mv finalMatrix.tab "${makeDir}/matrixResults"

	mkdir -p "${makeDir}/plotsResults"
	mv *png "${makeDir}/plotsResults"

	mkdir -p "${makeDir}/MSAs"
	mv maskedMatrixGenesUbiquitous.txt "${makeDir}/MSAs/"
	mv maskedMatrixGenesOnlyAncient.txt "${makeDir}/MSAs/"
	mv maskedMatrixGenesNoUbiquitous.txt "${makeDir}/MSAs/"
	mv genesAbovePercentSeries.txt "${makeDir}/MSAs/"
	mv *MSA* "${makeDir}/MSAs/"

	"""
}


