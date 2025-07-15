#!/usr/bin/env nextflow

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
This is the main document. To modify CPU usage and parameters fine tuning please go to nextflow.config
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
include { FILTER_GENE_ALIGNMENTS } from './modules/filter_gene_alignments.nf'
include { BUILD_MSA } from './modules/build_msa.nf'
include { TREE_THRESHOLD } from './modules/tree_threshold.nf'
include { TREE_CORE } from './modules/tree_core.nf'
include { TREE_ACCESSORY } from './modules/tree_accessory.nf'
include { TREE_ANCIENT } from './modules/tree_ancient.nf'
include { GET_RESULTS } from './modules/get_results.nf'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Print pipeline metadata: Version and Help
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
                GET_DATA(params.genomes, params.tax_id, params.get_data_parallel)
                fastaFiles = GET_DATA.out.fasta_files
                gffFiles = GET_DATA.out.gbk_files

        } else {
                fastaFiles = Channel.of(files("${params.trusted_data}/*fasta"))
                gffFiles = Channel.of(files("${params.trusted_data}/*gb"))

        }

	PARSE_GENBANK(gffFiles, fastaFiles)

	REMOVE_REDUNDANCY(PARSE_GENBANK.out.validFiles.map { fasta, gb -> fasta , gb}, params.remove_redundancy_parallel)

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

	ALIGNMENT_SUMMARY(params.config, ALIGNMENT.out.postAlignedBams, params.alignment_parallel)

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
			params.genomes, OUTGROUP_CONSENSUS.out.extractedSequencesOutgroupFasta, HEATMAP.out.blackListed, params.filter_gene_alignments_parallel)

	BUILD_MSA(FILTER_GENE_ALIGNMENTS.out.genesAlnSeq, HEATMAP.out.maskedMatrixGenesNoUbiquitous, HEATMAP.out.maskedMatrixGenesOnlyAncient, 
		HEATMAP.out.maskedMatrixGenesUbiquitous, HEATMAP.out.genesAbovePercentSeries, FILTER_GENE_ALIGNMENTS.out.sampleNames)

	TREE_THRESHOLD(BUILD_MSA.out.genesAbovePercentMSA, params.tree_threads)

	TREE_CORE(BUILD_MSA.out.maskedMatrixGenesUbiquitousMSA, params.tree_threads)

	TREE_ACCESSORY(BUILD_MSA.out.maskedMatrixGenesNoUbiquitousMSA, params.tree_threads)

	TREE_ANCIENT(BUILD_MSA.out.maskedMatrixGenesOnlyAncientMSA, params.tree_threads)

	// GET_RESULTS(
	// resultsDir, fastaDatabase.out.validFasta , fastaDatabase.out.validGff , fastaDatabase.out.fastaDatabaseLogFile , fastaDatabase.out.theFastaDatabase, 
	// clustering.out.clusteredDatabase, clustering.out.clusteringLog, prokkaMakeAnnotations.out.prokkaGFF, prokkaMakeAnnotations.out.prokkaLogfile,  makePangenome.out.panarooLog,
	// filterGeneAlignments.out.genesAlnSeq, formattingPangenome.out.panGenomeReference, updateNormalization.out.geneNormalizedUpdated, normalizationFunction.out.globalMeanCoverage,
	// alignmentSummary.out.postAlignmentFiles, alignmentSummary.out.refLenght, alignmentSummary.out.rawCoverage, alignmentSummary.out.completenessSummary, buildHeatmap.out.finalMatrix,
	// buildHeatmap.out.presenceAbsence, buildHeatmap.out.maskedMatrixGenesOnlyAncient, buildHeatmap.out.maskedMatrixGenesUbiquitous, buildHeatmap.out.maskedMatrixGenesNoUbiquitous,
	// buildHeatmap.out.genesAbovePercentSeries, treeThreshold.out.genesAbovePercentMSAIqtree ,treeThreshold.out.genesAbovePercentMSALog , treeThreshold.out.genesAbovePercentMSATreefile,
	// treeUbiquitous.out.maskedMatrixGenesUbiquitousMSAIqtree, treeUbiquitous.out.maskedMatrixGenesUbiquitousMSALog, treeUbiquitous.out.maskedMatrixGenesUbiquitousMSATreefile, 
	// treeNoUbiquitous.out.maskedMatrixGenesNoUbiquitousMSAIqtree , treeNoUbiquitous.out.maskedMatrixGenesNoUbiquitousMSALog , treeNoUbiquitous.out.maskedMatrixGenesNoUbiquitousMSATreefile,
	// treeAncient.out.maskedMatrixGenesOnlyAncientMSAIqtree , treeAncient.out.maskedMatrixGenesOnlyAncientMSALog , treeAncient.out.maskedMatrixGenesOnlyAncientMSATreefile
	// )

	// pMauve(fastaDatabase.out.validFasta)
	// xmfaToFasta(pMauve.out.pMauveCoreGenome)
	// filterMauveFasta(xmfaToFasta.out.pMauveFastaMSA)
	// startingTree(filterMauveFasta.out.concatenatedSeqtkMauveFastaMSA)
	// findRecombinationSpots(filterMauveFasta.out.concatenatedSeqtkMauveFastaMSA, startingTree.out.startingTreeMauveFasta, startingTree.out.kappa)
	// mapRecombinantsToGenes(findRecombinationSpots.out.recombinationMap, filterMauveFasta.out.concatenatedSeqtkMauveFastaMSA, blastMe.out.panGenomeReferenceDB, prokkaMakeAnnotations.out.prokkaGFF)
}
