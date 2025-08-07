/*
 * REALIGN_GENE_ALIGNMENTS{} will re-align each gene MSA.
 */


process REALIGN_GENE_ALIGNMENTS {
	conda "${projectDir}/envs/realign.yaml"

	input:
	path gene_msa	
	val parallel
	val threads

	output:
	path 'realigned/*fasta', emit: re_aligned

	script:
	"""
	#!/bin/bash

	mkdir -p ${params.output}/GENE_MSA
	mkdir -p realigned

	realign() {
	file=\$1
	name=\$(basename "\${file%.fasta}")
		mafft --auto --thread $threads "\$file" | seqtk seq > realigned/"\${name}".fasta
	}
	export -f realign
	find ./ -name "*fasta" | parallel -j $parallel realign

	cp -r realigned/* ${params.output}/GENE_MSA
	"""
}
