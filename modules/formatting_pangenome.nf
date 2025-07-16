/*
 * FORMATTING_PANGENOME{} will change format and index pangenome.
 */

process  FORMATTING_PANGENOME {
	conda "${projectDir}/envs/formatting_pangenome.yaml"

	input:
	path panGenomeReference

	output:
	tuple path('panGenomeReference.fasta'), path('panGenomeReference.dict'), path('panGenomeReference.fasta.fai'), emit: indexed_pangenome

	script:
	"""
	#!/bin/bash

	mkdir -p ${params.output}/PANGENOME

	seqtk seq $panGenomeReference > panGenomeReference.fasta
	picard CreateSequenceDictionary -R panGenomeReference.fasta
	samtools faidx panGenomeReference.fasta

	cp panGenomeReference.fasta ${params.output}/PANGENOME
	"""
}
