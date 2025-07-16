/*
 * BLAST_DATABASE{} will build a blast database for pangenome reference sequence.
 */

process BLAST_DATABASE {
	conda "${projectDir}/envs/blast.yaml"

	input:
	path panSeq

	output:
	path 'pangenome_blast_database*', emit: blast_database

	script:
	"""
	#!/bin/bash

	mkdir -p ${params.output}/PANGENOME

	makeblastdb -in $panSeq -dbtype nucl -out pangenome_blast_database

	cp pangenome_blast_database* ${params.output}/PANGENOME
	"""
}
