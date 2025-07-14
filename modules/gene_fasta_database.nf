/*
 * GENE_FASTA_DATABASE{} process will build a gene fasta database from non-redundant samples.
 */

process GENE_FASTA_DATABASE {

	conda "${projectDir}/envs/biopython.yaml"

	input:
	path gb_files

	output:
	path 'clustered_sequences.fasta' , emit: fastaDatabase

	script:
	"""
	#!/bin/bash

	parsing_and_contatenating.py
	"""
}
