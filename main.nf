#!/usr/bin/env nextflow

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
This is the main document. To modify CPU usage and parameters fine tuning please go to nextflow.config. 
Do not modify anything here unless you know what you are doing.
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
include { GENOTYPING } from './modules/genotyping.nf'
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
include { TEST } from './modules/test.nf'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Print pipeline metadata: Version and Help
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

	
def print_help() {
	 println "\n\033[1;31mSYNOPSIS\033[0m"	

	 println "\n\033[1;33mUSAGE\033[0m"
	 println "\nnextflow run main.nf --data <PATH> --output <PATH> --tax_id <INT> --config <PATH/FILE> [..OPTIONS..]"
	 
	 println "\n\033[1;33mMANDATORY\033[0m"
	 println "  --data                               <PATH>  data directory (fastq, fastq.gz)"
	 println "  --output                             <PATH>  output directory"
	 println "  --tax_id                              <INT>  taxonomical ID value"
	 println "  --outgroup_tax_id                     <INT>  outgroup taxonomical ID value"
	 println "  --config                        <PATH/FILE>  config file PATH"
	 
	 println "\n\033[1;33mOPTIONS\033[0m"
	 println "  --gene_completeness             <INT/FLOAT>  gene completeness/breadth of coverage threshold (Current value: ${params.gene_completeness})"
	 println "  --upper_coverage_bound                <INT>  maximum normalised coverage threshold (Current value: ${params.upper_coverage_bound})"
	 println "  --lower_coverage_bound                <INT>  mininum normalised coverage threshold (Current value: ${params.lower_coverage_bound})"
	 println "  --parallel                            <INT>  global parallel computing value (Current value: ${params.parallel})"
	 println "  --trusted_data                        <INT>  user curated input data PATH (*FASTA,*gb)"
	 println "  --genomes                             <INT>  number of genomes to download (Current value: ${params.genomes})"
	 println "  --get_data_parallel                   <INT>  number of samples to be downloaded in parallel. Do not use a value bigger than 3. (Current value: ${params.get_data_parallel})"
	 println "  --remove_redundancy_parallel          <INT>  parallel computing value for remove_redundancy.nf module [single core] (Current value: ${params.remove_redundancy_parallel})"
	 println "  --fastANI_threads                     <INT>  thread usage for fastANI (Current value: ${params.fastANI_threads})"
	 println "  --gene_identity_clustering            <INT>  gene clustering threshold for CD-HIT (Current value: ${params.gene_identity_clustering})"
	 println "  --cd_hit_threads                      <INT>  thread usage for CD-HIT (Current value: ${params.cd_hit_threads})"
	 println "  --prokka_annotate_parallel            <INT>  parallel computing value for gene annotation [multi core] (Current value: ${params.prokka_annotate_parallel})"
	 println "  --prokka_annotate_threads             <INT>  thread usage for prokka (Current value: ${params.prokka_annotate_threads})"
	 println "  --pangenome_identity_threshold        <INT>  gene clustering threshold for Panaroo (Current value: ${params.pangenome_identity_threshold})"
	 println "  --panaroo_pangenome_mode           <STRING>  Panaroo pangenone mode [strict,moderate,sensitive] (Current value: ${params.panaroo_pangenome_mode})"
	 println "  --panaroo_pangenome_threads           <INT>  thread usage for Panaroo pangenome building (Current value: ${params.panaroo_pangenome_threads})"
	 println "  --panaroo_alignment_type           <STRING>  Panaroo gene alignment mode [core,pan] (Current value: ${params.panaroo_alignment_type})"
	 println "  --outgroup_alignment_threads          <INT>  thread usage for outgroup alignment (Current value: ${params.outgroup_alignment_threads})"
	 println "  --alignment_threads                   <INT>  thread usage for alignment (Current value: ${params.alignment_threads})"
	 println "  --alignment_parallel                  <INT>  parallel computing value for alignment [multi core] (Current value: ${params.alignment_parallel})"
	 println "  --missing_prob                  <INT/FLOAT>  missing prob value for [bwa aln] (Current value: ${params.missing_prob})"
	 println "  --gap_fraction                  <INT/FLOAT>  maximum number or fraction of gap opens [bwa aln] (Current value: ${params.gap_fraction})"
	 println "  --seed                          <INT/FLOAT>  seed length [bwa aln] (Current value: ${params.seed})"
	 println "  --min_read_length                     <INT>  minimum read lenght cutoff after alignment (Current value: ${params.min_read_length})"
	 println "  --max_read_length                     <INT>  maximum read length cutoff after alignment (Current value: ${params.max_read_length})"
	 println "  --mapping_quality                     <INT>  minimum mapping quality cutoff after alignment (Current value: ${params.mapping_quality})"
	 println "  --update_normalization_parallel       <INT>  parallel computing value for update_normalization.nf module [single core] (Current value: ${params.update_normalization_parallel})"
	 println "  --bedtools_slop                       <INT>  extending sequences value (bp) (Current value: ${params.bedtools_slop})"    
	 println "  --extend_sequences_parallel           <INT>  parallel computing value for extend_sequences.nf module [single core] (Current value: ${params.extend_sequences_parallel})"
	 println "  --mapdamage_parallel                  <INT>  parallel computing value for mapDamage (Current value: ${params.mapdamage_parallel})"
	 println "  --realign_parallel                    <INT>  parallel computing value for realign.nf module [multi core] (Current value: ${params.realign_parallel})"
	 println "  --mafft_threads                       <INT>  thread usage for mafft (Current value: ${params.mafft_threads})"
	 println "  --bcftools_map_quality                <INT>  minimum mapping quality after genotyping (Current value: ${params.bcftools_map_quality})"
	 println "  --bcftools_base_quality               <INT>  minimum base quality after genotyping (Current value: ${params.bcftools_base_quality})"
	 println "  --variant_call_quality                <INT>  minimum variant call quality for genotyping (Current value: ${params.variant_call_quality})"
	 println "  --force_homozigosity                  <1/0>  exclude any site having both ALT and REF alleles. Set 1 to activate it (Current value: ${params.force_homozygosity})"
	 println "  --filter_gene_alignments_parallel     <INT>  parallel computing thread usage for filter_gene_alignment.nf module [single core] (Current value: ${params.filter_gene_alignments_parallel})"
	 println "  --tree_threads                        <INT>  thread usage for IQ-TREE. WARNING: This value will be x4 as there will be 4 phylogenetic runs in parallel. (Current value: ${params.tree_threads})"
	 println "  --n_samples_heatmap                   <INT>  maximum number of samples to include per heatmap (Current value: ${params.n_samples_heatmap})"
	 println "  --threshold_value_heatmap           <FLOAT>  custom gene set cutoff value. Include a gene if it is present in N percent of samples (Current value: ${params.threshold_value_heatmap})"
	 println "  --version                                    print version and exit"
	 println "  --help                                       print this page and exit"

    exit 0
}

def version() {
	println "You are using panchronos version 1.0"
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
println "\033[1;37mForce homozygosity: \033[0m: ${params.force_homozygosity}"
println "\033[1;37mVariant call quality: \033[0m: ${params.variant_call_quality}"
println "\033[1;37mBase quality: \033[0m: ${params.bcftools_base_quality}"
println "\033[1;37mMapping quality: \033[0m: ${params.bcftools_map_quality}"
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
                        params.bedtools_slop, params.extend_sequences_parallel, MAKE_PANGENOME.out.gene_list)

        FORMATTING_PANGENOME(EXTEND_SEQUENCES.out.extended_reference, MAKE_PANGENOME.out.panSequence)

        BLAST_DATABASE(EXTEND_SEQUENCES.out.unextended_reference)

        GET_OUTGROUP(params.outgroup_tax_id)

        OUTGROUP_READS(GET_OUTGROUP.out.outgroupFasta)

        OUTGROUP_ALIGNMENT(OUTGROUP_READS.out.outgroupReads, FORMATTING_PANGENOME.out.originalPangenomeReference,
                                params.outgroup_alignment_threads)

        OUTGROUP_CONSENSUS(OUTGROUP_ALIGNMENT.out.outgroupFastaPostAlignment,
                                FORMATTING_PANGENOME.out.originalPangenomeReference)

        if (params.test) {

            TEST(REMOVE_REDUNDANCY.out.nonRedundant_files.map { fasta, gb -> fasta})

                ALIGNMENT(TEST.out.test_data, FORMATTING_PANGENOME.out.indexed_pangenome.map { pangenome_reference, pangenome_dict, pangenome_index -> pangenome_reference}, params.config,
                        tuple(params.alignment_threads, params.missing_prob, params.seed, params.gap_fraction, params.min_read_length, params.max_read_length,
                        params.alignment_parallel, params.mapping_quality))
        } else {

                ALIGNMENT(params.data, FORMATTING_PANGENOME.out.indexed_pangenome.map { pangenome_reference, pangenome_dict, pangenome_index -> pangenome_reference}, params.config,
                        tuple(params.alignment_threads, params.missing_prob, params.seed, params.gap_fraction, params.min_read_length, params.max_read_length,
                        params.alignment_parallel, params.mapping_quality))
        }

        MAPDAMAGE(FORMATTING_PANGENOME.out.indexed_pangenome.map { pangenome_reference, pangenome_dict, pangenome_index -> pangenome_reference}, ALIGNMENT.out.pan_index, ALIGNMENT.out.bam_mapdamage,
                params.mapdamage_parallel)

        ALIGNMENT_SUMMARY(params.config, ALIGNMENT.out.postAlignedBams, params.alignment_parallel, params.bedtools_slop)

        NORMALIZE(EXTEND_SEQUENCES.out.pangenome_length, ALIGNMENT_SUMMARY.out.rawCoverage, params.alignment_parallel)

        UPDATE_NORMALIZATION(NORMALIZE.out.geneNormalizedSummary, ALIGNMENT_SUMMARY.out.completenessSummary, params.update_normalization_parallel)


        GENOTYPING(FORMATTING_PANGENOME.out.indexed_pangenome.map { pangenome_reference, pangenome_dict, pangenome_index -> pangenome_reference},
                           ALIGNMENT_SUMMARY.out.postAlignmentFiles, params.alignment_parallel, tuple(params.bcftools_map_quality , params.bcftools_base_quality, params.variant_call_quality),
                           params.bedtools_slop, params.force_homozygosity)
        extractedSequencesFasta = GENOTYPING.out.consensusSequences


        COVERAGE_BOUNDS(UPDATE_NORMALIZATION.out.geneNormalizedUpdated,  params.lower_coverage_bound, params.upper_coverage_bound, params.gene_completeness)

        PLOT_COVERAGE_COMPLETENESS(UPDATE_NORMALIZATION.out.geneNormalizedUpdated, params.gene_completeness, params.lower_coverage_bound, params.upper_coverage_bound,
                                                                params.normalised_coverage_boundary_plot)

        UPDATE_MATRIX(MAKE_PANGENOME.out.initialMatrix , NORMALIZE.out.globalMeanCoverage, COVERAGE_BOUNDS.out.geneNormalizedUpdatedFiltered,
                params.gene_completeness, params.lower_coverage_bound, params.upper_coverage_bound, EXTEND_SEQUENCES.out.final_list_genes)

        HEATMAP(UPDATE_MATRIX.out.finalCsv, UPDATE_MATRIX.out.index ,UPDATE_MATRIX.out.matrix, UPDATE_MATRIX.out.sampleNames, params.threshold_value_heatmap, params.n_samples_heatmap)

        UPDATE_PLOT_COVERAGE_COMPLETENESS(COVERAGE_BOUNDS.out.geneNormalizedUpdatedFiltered, params.gene_completeness, params.lower_coverage_bound, params.upper_coverage_bound,
                                                                                params.normalised_coverage_boundary_plot)


        FILTER_GENE_ALIGNMENTS(MAKE_PANGENOME.out.alignedGenesSeqs, extractedSequencesFasta, REMOVE_REDUNDANCY.out.nonRedundant_files.map { fasta, gb -> fasta },
                        params.genomes, OUTGROUP_CONSENSUS.out.extractedSequencesOutgroupFasta, HEATMAP.out.blackListed, params.filter_gene_alignments_parallel, EXTEND_SEQUENCES.out.final_list_genes, HEATMAP.out.genesIndex)

        REALIGN_GENE_ALIGNMENTS(FILTER_GENE_ALIGNMENTS.out.genesAlnSeq, params.realign_parallel, params.mafft_threads)

        BUILD_MSA(REALIGN_GENE_ALIGNMENTS.out.re_aligned, HEATMAP.out.maskedMatrixGenesNoUbiquitous, HEATMAP.out.maskedMatrixGenesOnlyAncient,
                HEATMAP.out.maskedMatrixGenesUbiquitous, HEATMAP.out.genesAbovePercentSeries, FILTER_GENE_ALIGNMENTS.out.sampleNames, params.parallel)

        TREE_THRESHOLD(BUILD_MSA.out.genesAbovePercentMSA, params.tree_threads)

        TREE_CORE(BUILD_MSA.out.maskedMatrixGenesUbiquitousMSA, params.tree_threads)

        TREE_ACCESSORY(BUILD_MSA.out.maskedMatrixGenesNoUbiquitousMSA, params.tree_threads)

        TREE_ANCIENT(BUILD_MSA.out.maskedMatrixGenesOnlyAncientMSA, params.tree_threads)
}
