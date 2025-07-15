/*
 * MAPPING{} will index a graph and map reads against it.
 */

process UPDATE_PLOT_COVERAGE_COMPLETENESS {
	conda "${projectDir}/envs/plot.yaml"
	
	input:
	path geneNormalizedUpdatedFiltered
	val completeness
	val coverage

	output:
	path 'plotCoverageVsCompletenessFiltered.png', emit: plotCoverageVsCompletenessFiltered
	
	script:
	"""
	plot_cvg_vs_completeness.py $geneNormalizedUpdatedFiltered $completeness $coverage
	mv plotCoverage_vs_Completeness.png ./plotCoverageVsCompletenessFiltered.png
	"""
}
