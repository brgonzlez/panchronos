/*
 * FORMATTING_PANGENOME{} will change format and index pangenome.
 */

process  FORMATTING_PANGENOME {
	conda "${projectDir}/envs/formatting_pangenome.yaml"

	input:
	path panGenomeReference_extended
	path panRef

	output:
	tuple path('panGenomeReference_extended.fasta'), path('panGenomeReference_extended.dict'), path('panGenomeReference_extended.fasta.fai'), emit: indexed_pangenome
	path 'panGenomeReference.fasta', emit: originalPangenomeReference

	script:
	"""
	#!/bin/bash

	mkdir -p ${params.output}/PANGENOME

	seqtk seq $panGenomeReference_extended > extended_pangenome_reference.fasta
	picard CreateSequenceDictionary -R extended_pangenome_reference.fasta
	samtools faidx extended_pangenome_reference.fasta

	seqtk seq $panRef >  panGenomeReference_panaroo.fasta 

	cp extended_pangenome_reference.fasta ${params.output}/PANGENOME
	cp panGenomeReference_panaroo.fasta ${params.output}/PANGENOME
	"""
}
