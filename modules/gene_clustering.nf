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

	mkdir -p ${params.output}/GENE_DATABASE/

	awk '/^>/ {keep = (\$0 !~ /hypothetical/)} keep' $fastaDB > clustered_sequences_no_hypotheticals.fasta

	cd-hit-est -i clustered_sequences_no_hypotheticals.fasta -o clustered_non_redundant_genes.fasta -c $clustering -T $threads -d 0 -g 1 -M 0
	
	cp clustered_non_redundant_genes.fasta ${params.output}/GENE_DATABASE/

	cat .command.out >> clustering.log
	"""
}
