/*
 * REALIGN_GENE_ALIGNMENTS{} will re-align each gene MSA.
 */


process REALIGN_GENE_ALIGNMENTS {
	conda "${projectDir}/envs/realign.yaml"

	input:
	path gene_msa	
	val parallel
	
	output:
	path 'realigned/*fasta', emit: re_aligned

	script:
	"""
	#!/bin/bash

	mkdir -p realigned

	realign() {
	file=\$1
	name=\$(basename "\${file%.fasta}")
		mafft --auto 
	}
	export -f realign
	find ./ -name "*fasta" | parallel -j $parallel realign
	"""
}
