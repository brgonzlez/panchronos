/*
 * MAKE_PANGENOME{} will build a gene pangenome based on prokka gff output.
 */
 
process MAKE_PANGENOME {
	conda "${projectDir}/envs/panaroo.yaml"
	
	input:
	path prokka_gff
	val pangenomeMode
	val pangenomeThreshold
	val threads

	output:
	path 'pan_genome_reference.fa' , emit: panSequence
	path 'gene_presence_absence.Rtab' , emit: initialMatrix
	path 'aligned_gene_sequences/*' , emit: alignedGenesSeqs

	script:
	"""
	panaroo -i *.gff -o ./ --clean-mode $pangenomeMode -a core --core_threshold $pangenomeThreshold -t $threads

	cat .command.out >> makePangenome.log
	"""
}
