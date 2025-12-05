/*
 * MAPPING{} will index a graph and map reads against it.
 */

process UPDATE_PLOT_COVERAGE_COMPLETENESS {
	conda "${projectDir}/envs/plot.yaml"
	
	input:
	path geneNormalizedUpdatedFiltered
	val completeness
	val coverage_lower
	val coverage_upper
	val normalised_boundary_plot

	output:
	path 'plotCoverage_vs_Completeness_filtered*', emit: plotCoverageVsCompletenessFiltered
	
	script:
	"""
	#!/bin/bash

	mkdir -p ${params.output}/PLOTS

	#clean-up the names
	sed -i -e 's/postPangenomeAlignment_//g' $geneNormalizedUpdatedFiltered

	#we make one file per sample
	awk 'NR>1 {print \$1}' $geneNormalizedUpdatedFiltered | sort | uniq > samples.txt

	while read -r sample;do
		name=\$(echo "\${sample#postPangenomeAlignment_}")
		awk 'NR==1{print \$0}' $geneNormalizedUpdatedFiltered > "\$name"_individual_normalised.tab
		grep -w "\$sample" $geneNormalizedUpdatedFiltered >> "\$name"_individual_normalised.tab
	done < samples.txt

	plot_cov() {
	tab_file=\$1

	name=\$(basename "\${tab_file%_individual_normalised.tab}")
	plot_cvg_vs_completeness.py "\$tab_file" $completeness $coverage_lower $coverage_upper $normalised_boundary_plot
	mv plotCoverage_vs_Completeness_"\${name}".png plotCoverage_vs_Completeness_filtered_"\${name}".png
	}
	export -f plot_cov
	find ./ -name "*_individual_normalised.tab" | parallel -j $task.cpus plot_cov

	

	cp *png ${params.output}/PLOTS
	"""
}
