/*
 * PLOT_COVERAGE_COMPLETENESS{} will plot normalized gene coverage vs gene completeness (breadth of coverage).
 */

process PLOT_COVERAGE_COMPLETENESS {
	conda "${projectDir}/envs/plot.yaml"
	
	input:
	path geneNormalizedUpdated
	val completeness
	val coverage_lower
	val coverage_upper

	output:
	path 'plotCoverage_vs_Completeness.png', emit: plotCoverage_vs_Completeness
	
	script:
	"""
 	#!/bin/bash

	mkdir -p ${params.output}/PLOTS

	plot_cvg_vs_completeness.py $geneNormalizedUpdated $completeness $coverage_lower $coverage_upper

	cp *png ${params.output}/PLOTS
	"""
}
