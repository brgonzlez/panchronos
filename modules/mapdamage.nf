/*
 * MAPDAMAGE{} will run mapdamage on aligned bam.
 */

process MAPDAMAGE {
	conda "${projectDir}/envs/mapdamage.yaml"

  	input:
	path panRef 
	path index
	path bams
	val parallel

	output:
	stdout


	script:
	"""
	#!/bin/bash

	mkdir -p ${params.output}/MAPDAMAGE

	run_mapdamage() {
	file=\$1
	name=\$(basename "\${file%_sorted_mappedreads.bam}")
	
		mapDamage -i "\${file}" -r $panRef -d mapdamage_"\${name}"
	}
	export -f run_mapdamage
	find ./* -name "*.bam" | parallel -j $parallel run_mapdamage

	cp -r mapdamage_* ${params.output}/MAPDAMAGE
	"""
}
	
