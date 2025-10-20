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
include { EXTEND_SEQUENCES } from './modules/extend_sequences.nf'
include { FORMATTING_PANGENOME } from './modules/formatting_pangenome.nf'
include { BLAST_DATABASE } from './modules/blast_database.nf'
include { GET_OUTGROUP } from './modules/get_outgroup.nf'
include { OUTGROUP_READS } from './modules/outgroup_reads.nf'
include { OUTGROUP_ALIGNMENT } from './modules/outgroup_alignment.nf'
include { OUTGROUP_CONSENSUS } from './modules/outgroup_consensus.nf'
include { ALIGNMENT } from './modules/alignment.nf'
include { ALIGNMENT_SUMMARY } from './modules/alignment_summary.nf'
include { NORMALIZE } from './modules/normalize.nf'
include { UPDATE_NORMALIZATION } from './modules/update_normalization.nf'
include { BCFTOOLS_CONSENSUS } from './modules/bcftools_consensus.nf'
include { GATK_CONSENSUS } from './modules/gatk_consensus.nf'
include { PLOT_COVERAGE_COMPLETENESS } from './modules/plot_coverage_completeness.nf'
include { COVERAGE_BOUNDS } from './modules/coverage_bounds.nf'
include { UPDATE_MATRIX } from './modules/update_matrix.nf'
include { HEATMAP } from './modules/heatmap.nf'
include { UPDATE_PLOT_COVERAGE_COMPLETENESS } from './modules/update_plot_coverage_completeness.nf'
include { FILTER_GENE_ALIGNMENTS } from './modules/filter_gene_alignments.nf'
include { REALIGN_GENE_ALIGNMENTS } from './modules/realign.nf'
include { BUILD_MSA } from './modules/build_msa.nf'
include { TREE_THRESHOLD } from './modules/tree_threshold.nf'
include { TREE_CORE } from './modules/tree_core.nf'
include { TREE_ACCESSORY } from './modules/tree_accessory.nf'
include { TREE_ANCIENT } from './modules/tree_ancient.nf'
include { MAPDAMAGE } from './modules/mapdamage.nf'

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

// Running check in case we messed up something.

println "\n\033[1;33mCHECKING PARAMETERS\033[0m"
println "======================="
println "\033[1;31mBASIC OPTIONS\033[0m"
println "\033[1;37mData\033[0m: ${params.data}"
println "\033[1;37mOutput\033[0m: ${params.output}"
println "\033[1;37mGene completeness\033[0m: ${params.gene_completeness}"
println "\033[1;37mGenomes to download\033[0m: ${params.genomes}"
println "\033[1;37mConfig\033[0m: ${params.config}"
println "\033[1;37mParallel\033[0m: ${params.parallel}"
println "\033[1;37mTrusted data\033[0m: ${params.trusted_data}"
println "======================="
println "\033[1;31mProcess: get_data.nf\033[0m"
println "\033[1;37mParallel\033[0m: ${params.get_data_parallel}"
println "======================="
println "\033[1;31mProcess: remove_redundancy.nf\033[0m"
println "\033[1;37mParallel\033[0m: ${params.remove_redundancy_parallel}"
println "\033[1;37m[fastANI] threads\033[0m: ${params.fastANI_threads}"
println "======================="
println "\033[1;31mProcess: gene_clustering.nf\033[0m"
println "\033[1;37mGene identity threshold\033[0m: ${params.gene_identity_clustering}"
println "\033[1;37m[cd-hit] threads\033[0m: ${params.cd_hit_threads}"
println "======================="
println "\033[1;31mProcess: annotate.nf\033[0m"
println "\033[1;37m[prokka] Parallel\033[0m: ${params.prokka_annotate_parallel}"
println "\033[1;37m[prokka] Threads\033[0m: ${params.prokka_annotate_threads}"
println "======================="
println "\033[1;31mProcess: make_pangenome.nf\033[0m"
println "\033[1;37m[panaroo] Identity threshold\033[0m: ${params.pangenome_identity_threshold}"
println "\033[1;37m[panaroo] Pangeome mode\033[0m: ${params.panaroo_pangenome_mode}"
println "\033[1;37m[panaroo] Threads\033[0m: ${params.panaroo_pangenome_threads}"
println "\033[1;37m[panaroo] Alignment type (core,pan)\033[0m: ${params.panaroo_alignment_type}"
println "======================="
println "\033[1;31mProcess: outgroup_alignment.nf\033[0m"
println "\033[1;37m[bwa/samtools] Threads\033[0m: ${params.outgroup_alignment_threads}"
println "======================="
println "\033[1;31mProcess: extend_sequences.nf\033[0m"
println "\033[1;37m[bedtools] Sequence extension\033[0m: ${params.bedtools_slop}"
println "\033[1;37mParallel\033[0m: ${params.extend_sequences_parallel}"
println "======================="
println "\033[1;31mProcess: alignment.nf\033[0m"
println "\033[1;37m[bwa/samtools] Threads\033[0m: ${params.alignment_threads}"
println "\033[1;37m[bwa/samtools] Missing prob.\033[0m: ${params.missing_prob}"
println "\033[1;37m[bwa/samtools] Gap fraction\033[0m: ${params.gap_fraction}"
println "\033[1;37m[bwa/samtools] Seed\033[0m: ${params.seed}"
println "\033[1;37m[bwa/samtools] Min. read length\033[0m: ${params.min_read_length}"
println "\033[1;37m[bwa/samtools] Max. read length\033[0m: ${params.max_read_length}"
println "\033[1;37m[bwa/samtools] Mapping quality\033[0m: ${params.mapping_quality}"
println "\033[1;37m[bwa/samtools] Parallel\033[0m: ${params.alignment_parallel}"
println "======================="
println "\033[1;31mGenotyping\033[0m"
println "\033[1;37mGenotyper: \033[0m: ${params.genotyper}"
println "======================="
println "\033[1;31mProcess: filter_gene_alignments.nf\033[0m"
println "\033[1;37mParallel \033[0m: ${params.filter_gene_alignments_parallel}"
println "======================="
println "\033[1;31mProcess: mapdamage.nf\033[0m"
println "\033[1;37mParallel \033[0m: ${params.mapdamage_parallel}"
println "======================="
println "\033[1;31mProcess: realign.nf\033[0m"
println "\033[1;37mParallel \033[0m: ${params.realign_parallel}"
println "\033[1;37m[mafft] Threads\033[0m: ${params.mafft_threads}"
println "======================="
println "\033[1;31mPhylogeny\033[0m"
println "\033[1;37m[iq-tree] Threads\033[0m: ${params.tree_threads}"
println "=======================\n"


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

	REMOVE_REDUNDANCY(PARSE_GENBANK.out.validFiles.map { fasta, gb -> tuple(fasta , gb)}, params.remove_redundancy_parallel, params.fastANI_threads)

	GENE_FASTA_DATABASE(REMOVE_REDUNDANCY.out.nonRedundant_files.map { fasta, gb -> gb})

	GENE_CLUSTERING(GENE_FASTA_DATABASE.out.fastaDatabase, params.gene_identity_clustering, params.cd_hit_threads)

	ANNOTATE(GENE_CLUSTERING.out.clusteredDatabase, params.prokka_annotate_threads, REMOVE_REDUNDANCY.out.nonRedundant_files.map { fasta, gb -> tuple(fasta, gb)},
		params.prokka_annotate_parallel)

	MAKE_PANGENOME(ANNOTATE.out.prokka_gff, params.panaroo_pangenome_mode, params.pangenome_identity_threshold, params.panaroo_pangenome_threads, params.panaroo_alignment_type)

	EXTEND_SEQUENCES(ANNOTATE.out.prokka_gff, REMOVE_REDUNDANCY.out.nonRedundant_files.map { fasta, gb -> fasta }, MAKE_PANGENOME.out.pangenome_metadata.map { graph, gene_data -> tuple(graph, gene_data)},
			params.bedtools_slop, params.extend_sequences_parallel)

	FORMATTING_PANGENOME(EXTEND_SEQUENCES.out.extended_reference, MAKE_PANGENOME.out.panSequence)

	BLAST_DATABASE(FORMATTING_PANGENOME.out.indexed_pangenome.map { pangenome_reference, pangenome_dict, pangenome_index -> pangenome_reference})

	GET_OUTGROUP(params.outgroup_tax_id)

	OUTGROUP_READS(GET_OUTGROUP.out.outgroupFasta)

	OUTGROUP_ALIGNMENT(OUTGROUP_READS.out.outgroupReads, FORMATTING_PANGENOME.out.originalPangenomeReference,
				params.outgroup_alignment_threads)

	OUTGROUP_CONSENSUS(OUTGROUP_ALIGNMENT.out.outgroupFastaPostAlignment, 
				FORMATTING_PANGENOME.out.originalPangenomeReference)

	ALIGNMENT(params.data, FORMATTING_PANGENOME.out.indexed_pangenome.map { pangenome_reference, pangenome_dict, pangenome_index -> pangenome_reference}, params.config, 
		tuple(params.alignment_threads, params.missing_prob, params.seed, params.gap_fraction, params.min_read_length, params.max_read_length, 
		params.alignment_parallel, params.mapping_quality))

	MAPDAMAGE(FORMATTING_PANGENOME.out.indexed_pangenome.map { pangenome_reference, pangenome_dict, pangenome_index -> pangenome_reference}, ALIGNMENT.out.pan_index, ALIGNMENT.out.bam_mapdamage,
		params.mapdamage_parallel)

	ALIGNMENT_SUMMARY(params.config, ALIGNMENT.out.postAlignedBams, params.alignment_parallel, params.bedtools_slop)

	NORMALIZE(ALIGNMENT_SUMMARY.out.refLength, ALIGNMENT_SUMMARY.out.rawCoverage, params.alignment_parallel)

	UPDATE_NORMALIZATION(NORMALIZE.out.geneNormalizedSummary, ALIGNMENT_SUMMARY.out.completenessSummary, params.update_normalization_parallel)

	COVERAGE_BOUNDS(UPDATE_NORMALIZATION.out.geneNormalizedUpdated,  params.lower_coverage_bound, params.upper_coverage_bound, params.gene_completeness)


	if (params.genotyper == "gatk") {
		GATK_CONSENSUS(FORMATTING_PANGENOME.out.indexed_pangenome.map { pangenome_reference, pangenome_dict, pangenome_index ->  tuple(pangenome_reference, pangenome_dict, pangenome_index)}, 
				ALIGNMENT_SUMMARY.out.postAlignmentFiles, params.alignment_parallel)
		extractedSequencesFasta = GATK_CONSENSUS.out.gatkConsensusSequences
		vcfFile = GATK_CONSENSUS.out.gatkGenotypes

	} else if (params.genotyper == "bcftools") {
		BCFTOOLS_CONSENSUS(FORMATTING_PANGENOME.out.indexed_pangenome.map { pangenome_reference, pangenome_dict, pangenome_index -> pangenome_reference}, 
					ALIGNMENT_SUMMARY.out.postAlignmentFiles, params.alignment_parallel, tuple(params.bcftools_map_quality , params.bcftools_base_quality), 
					params.bedtools_slop)
		extractedSequencesFasta = BCFTOOLS_CONSENSUS.out.consensusSequences

	} else {
		error "Invalid option for --genotyper. Please choose 'gatk' or 'bcftools'."
	}

	PLOT_COVERAGE_COMPLETENESS(UPDATE_NORMALIZATION.out.geneNormalizedUpdated, params.gene_completeness, params.lower_coverage_bound)

	UPDATE_MATRIX(MAKE_PANGENOME.out.initialMatrix , NORMALIZE.out.globalMeanCoverage, COVERAGE_BOUNDS.out.geneNormalizedUpdatedFiltered, 
		params.gene_completeness, params.lower_coverage_bound, params.upper_coverage_bound)

	HEATMAP(UPDATE_MATRIX.out.finalCsv, UPDATE_MATRIX.out.index ,UPDATE_MATRIX.out.matrix, UPDATE_MATRIX.out.sampleNames)

	UPDATE_PLOT_COVERAGE_COMPLETENESS(COVERAGE_BOUNDS.out.geneNormalizedUpdatedFiltered, params.gene_completeness, params.lower_coverage_bound)

	FILTER_GENE_ALIGNMENTS(MAKE_PANGENOME.out.alignedGenesSeqs, extractedSequencesFasta, REMOVE_REDUNDANCY.out.nonRedundant_files.map { fasta, gb -> fasta }, 
			params.genomes, OUTGROUP_CONSENSUS.out.extractedSequencesOutgroupFasta, HEATMAP.out.blackListed, params.filter_gene_alignments_parallel)

	REALIGN_GENE_ALIGNMENTS(FILTER_GENE_ALIGNMENTS.out.genesAlnSeq, params.realign_parallel, params.mafft_threads)

	BUILD_MSA(REALIGN_GENE_ALIGNMENTS.out.re_aligned, HEATMAP.out.maskedMatrixGenesNoUbiquitous, HEATMAP.out.maskedMatrixGenesOnlyAncient, 
		HEATMAP.out.maskedMatrixGenesUbiquitous, HEATMAP.out.genesAbovePercentSeries, FILTER_GENE_ALIGNMENTS.out.sampleNames)

	TREE_THRESHOLD(BUILD_MSA.out.genesAbovePercentMSA, params.tree_threads)

	TREE_CORE(BUILD_MSA.out.maskedMatrixGenesUbiquitousMSA, params.tree_threads)

	TREE_ACCESSORY(BUILD_MSA.out.maskedMatrixGenesNoUbiquitousMSA, params.tree_threads)

	TREE_ANCIENT(BUILD_MSA.out.maskedMatrixGenesOnlyAncientMSA, params.tree_threads)

	// pMauve(fastaDatabase.out.validFasta)
	// xmfaToFasta(pMauve.out.pMauveCoreGenome)
	// filterMauveFasta(xmfaToFasta.out.pMauveFastaMSA)
	// startingTree(filterMauveFasta.out.concatenatedSeqtkMauveFastaMSA)
	// findRecombinationSpots(filterMauveFasta.out.concatenatedSeqtkMauveFastaMSA, startingTree.out.startingTreeMauveFasta, startingTree.out.kappa)
	// mapRecombinantsToGenes(findRecombinationSpots.out.recombinationMap, filterMauveFasta.out.concatenatedSeqtkMauveFastaMSA, blastMe.out.panGenomeReferenceDB, prokkaMakeAnnotations.out.prokkaGFF)
}
