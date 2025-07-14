/*
 * OUTGROUP_READS{} will generate artificial reads for outgroup.
 */

process OUTGROUP_READS {
	conda "${projectDir}/envs/art.yaml"

	input:
	path outgroupFasta
	
	output:
	path 'outgroupReads.fq', emit: outgroupReads


	script:
	"""
	#!/bin/bash

	art_illumina -i $outgroupFasta -l 150 -f 100 -o outgroupReads
	"""
}
