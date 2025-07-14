/*
 * GENE_CLUSTERING{} process will cluster gene sequences based on identity threshold.
 */

process GENE_CLUSTERING {

	conda "${projectDir}/envs/cdhit.yaml"

	input:
	path fastaDB
	val clustering
	val threads

	output:
	path "clustered_non_redundant_genes.fasta", emit: clusteredDatabase

	script:
	"""
	#!/bin/bash

	cd-hit-est -i $fastaDB -o clustered_non_redundant_genes.fasta -c $clustering -T $threads -d 0 -g 1 -M 0

	cat .command.out >> clustering.log
	"""
}
