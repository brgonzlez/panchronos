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

	#we make one file per sample
	awk 'NR>1 {print \$1}' $geneNormalizedUpdated | sort | uniq > samples.txt


	while read -r sample;do
		name=\$(echo "\${sample#postPangenomeAlignment_}")
		awk 'NR==1{print \$0}' $geneNormalizedUpdated > "\$name"_individual_normalised.tab
		grep -w "\$sample" $geneNormalizedUpdated >> "\$name"_individual_normalised.tab
	done < samples.txt

	plot_cov() {
	tab_file=\$1

	plot_cvg_vs_completeness.py "\$tab_file" $completeness $coverage_lower $coverage_upper
	}
	export -f plot_cov
	find ./ -name "*_individual_normalised.tab" | parallel -j $task.cpus plot_cov

	plot_cvg_vs_completeness.py $geneNormalizedUpdated $completeness $coverage_lower $coverage_upper

	cp *png ${params.output}/PLOTS
	"""
}
