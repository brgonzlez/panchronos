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
	val alignment

	output:
	path 'pan_genome_reference.fa' , emit: panSequence
	path 'gene_presence_absence.Rtab' , emit: initialMatrix
	path 'aligned_gene_sequences/*' , emit: alignedGenesSeqs
	tuple path('final_graph.gml'), path('gene_data.csv'), emit: pangenome_metadata
	path 'pangenome_length.txt', emit: pangenomeLength

	script:
	"""
	#!/bin/bash

	#make pangenome
	panaroo -i *.gff -o ./ --clean-mode $pangenomeMode -a $alignment --core_threshold $pangenomeThreshold -t $threads

	#get pangenome length
	seqtk seq pan_genome_reference.fa | awk '!/^>/ {line_length += length(\$0)} END {print line_length}' > pangenome_length.txt
	
	cat .command.out >> makePangenome.log
	"""
}
