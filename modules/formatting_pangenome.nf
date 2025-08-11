/*
 * FORMATTING_PANGENOME{} will change format and index pangenome.
 */

process  FORMATTING_PANGENOME {
	conda "${projectDir}/envs/formatting_pangenome.yaml"

	input:
	path panGenomeReference_extended
	path panRef

	output:
	tuple path('panGenomeReference.fasta'), path('panGenomeReference.dict'), path('panGenomeReference.fasta.fai'), emit: indexed_pangenome

	script:
	"""
	#!/bin/bash

	mkdir -p ${params.output}/PANGENOME

	seqtk seq $panGenomeReference_extended > panGenomeReference_extended.fasta
	picard CreateSequenceDictionary -R panGenomeReference_extended.fasta
	samtools faidx panGenomeReference_extended.fasta

	seqtk seq $panRef >  panGenomeReference.fasta 

	cp panGenomeReference_extended.fasta ${params.output}/PANGENOME
	cp panGenomeReference.fasta ${params.output}/PANGENOME
	"""
}
