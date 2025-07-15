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
include { ANNOTATE } from './modules/annotate.nf'
include { MAKE_PANGENOME } from './modules/make_pangenome.nf'
include { BLAST_DATABASE } from './modules/blast_database.nf'
include { GET_OUTGROUP } from './modules/get_outgroup.nf'
include { OUTGROUP_READS } from './modules/outgroup_reads.nf'
include { OUTGROUP_ALIGNMENT } from './modules/outgroup_alignment.nf'
include { OUTGROUP_CONSENSUS } from './modules/outgroup_consensus.nf'
include { ALIGNMENT } from './modules/alignment.nf'
include { ALIGNMENT_SUMMARY } from './modules/alignment_summary.nf'
include { NORMALIZATION } from './modules/normalization.nf'
include { UPDATE_NORMALIZATION } from './modules/update_normalization.nf'
include { BCFTOOLS_CONSENSUS } from './modules/bcftools_consensus.nf'
include { GATK_CONSENSUS } from './modules/gatk_consensus.nf'
include { PLOT_COVERAGE_COMPLETENESS } from './modules/plot_coverage_completeness.nf.nf'
include { COVERAGE_BOUNDS } from './modules/coverage_bounds.nf'
include { UPDATE_MATRIX } from './modules/update_matrix.nf'
include { HEATMAP } from './modules/heatmap.nf'
include { UPDATE_PLOT_COVERAGE_COMPLETENESS } from './modules/update_plot_coverage_completeness.nf.nf'
include {  } from './modules/.nf'
include {  } from './modules/.nf'
include {  } from './modules/.nf'
include {  } from './modules/.nf'
include {  } from './modules/.nf'
include {  } from './modules/.nf'
include {  } from './modules/.nf'
include {  } from './modules/.nf'
include {  } from './modules/.nf'
include {  } from './modules/.nf'



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

	MAKE_PANGENOME(ANNOTATE.out.prokka_gff, params.panaroo_pangenome_mode, params.pangenome_identity_threshold, params.panaroo_pangenome_threads)

	FORMATTING_PANGENOME(MAKE_PANGENOME.out.panSequence)

	BLAST_DATABASE(FORMATTING_PANGENOME.out.map { pangenome_reference, pangenome_dict, pangenome_index -> pangenome_reference})

        GET_OUTGROUP(params.outgroup_tax_id)

        OUTGROUP_READS(GET_OUTGROUP.out.outgroupFasta)

        OUTGROUP_ALIGNMENT(OUTGROUP_READS.out.outgroupReads, FORMATTING_PANGENOME.out.map { pangenome_reference, pangenome_dict, pangenome_index -> pangenome_reference},
				params.outgroup_alignment_threads)

        OUTGROUP_CONSENSUS(OUTGROUP_ALIGNMENT.out.outgroupFastaPostAlignment, 
				FORMATTING_PANGENOME.out.map { pangenome_reference, pangenome_dict, pangenome_index -> pangenome_reference})

	ALIGNMENT(params.data, FORMATTING_PANGENOME.out.map { pangenome_reference, pangenome_dict, pangenome_index -> pangenome_reference}, params.config, 
		tuple(params.alignment_threads, params.missing_prob, params.seed, params.gap_fraction, params.min_read_length, params.max_read_length, params.alignment_parallel))

	ALIGNMENT_SUMMARY(configFile, alignment.out.postAlignedBams, params.alignment_parallel)

	NORMALIZATION(ALIGNMENT_SUMMARY.out.refLenght, ALIGNMENT_SUMMARY.out.rawCoverage, params.alignment_parallel)

	UPDATE_NORMALIZATION(NORMALIZATION.out.geneNormalizedSummary, ALIGNMENT_SUMMARY.out.completenessSummary)


	if (params.genotyper == "gatk") {
		GATK_CONSENSUS(FORMATTING_PANGENOME.out.map { pangenome_reference, pangenome_dict, pangenome_index ->  pangenome_reference, pangenome_dict, pangenome_index}, 
				ALIGNMENT_SUMMARY.out.postAlignmentFiles, params.alignment_parallel)
		extractedSequencesFasta = GATK_CONSENSUS.out.gatkConsensusSequences
		vcfFile = GATK_CONSENSUS.out.gatkGenotypes

	} else if (params.genotyper == "bcftools") {
		BCFTOOLS_CONSENSUS(FORMATTING_PANGENOME.out.map { pangenome_reference, pangenome_dict, pangenome_index -> pangenome_reference}, 
					ALIGNMENT_SUMMARY.out.postAlignmentFiles, params.alignment_parallel)
		extractedSequencesFasta = BCFTOOLS_CONSENSUS.out.consensusSequences

	} else {
		error "Invalid option for --genotyper. Please choose 'gatk' or 'bcftools'."
	}

	PLOT_COVERAGE_COMPLETENESS(UPDATE_NORMALIZATION.out.geneNormalizedUpdated, params.gene_completeness, params.lower_coverage_bound)

        COVERAGE_BOUNDS(UPDATE_NORMALIZATION.out.geneNormalizedUpdated,  params.lower_coverage_bound, params.upper_coverage_bound, params.gene_completeness)

	UPDATE_MATRIX(MAKE_PANGENOME.out.initialMatrix , NORMALIZATION.out.globalMeanCoverage, COVERAGE_BOUNDS.out.geneNormalizedUpdatedFiltered, 
		params.gene_completeness, params.lower_coverage_bound, params.upper_coverage_bound)

	HEATMAP(UPDATE_MATRIX.out.finalCsv, UPDATE_MATRIX.out.index ,UPDATE_MATRIX.out.matrix, UPDATE_MATRIX.out.sampleNames)

	UPDATE_PLOT_COVERAGE_COMPLETENESS(COVERAGE_BOUNDS.out.geneNormalizedUpdatedFiltered, params.gene_completeness, params.lower_coverage_bound)

	FILTER_GENE_ALIGNMENTS(MAKE_PANGENOME.out.alignedGenesSeqs, extractedSequencesFasta, REMOVE_REDUNDANCY.out.nonRedundant_files.map { fasta, gb -> fasta }, 
			params.genomes, OUTGROUP_CONSENSUS.out.extractedSequencesOutgroupFasta, HEATMAP.out.blackListed)

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


