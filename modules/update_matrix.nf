/*
 * UPDATE_MATRIX{} will .
 */


process UPDATE_MATRIX {
	conda "${projectDir}/envs/plot.yaml"

	input:
	path pangenomeRtab
	path globalMeanCoverage
	path normalized
	val gCompleteness
	val lowerBound
	val upperBound


	output:
	path 'matrix.tab', emit: matrix
	path '*_final.csv', emit: finalCsv
	path 'sample_names', emit: sampleNames
	path 'INDEX', emit: INDEX

	script:
	"""
	awk 'NR==1{print \$0}' $pangenomeRtab > matrix.tab
	awk 'NR>1 {print \$0}' $pangenomeRtab | sort -k 1 -t \$'\t' >> matrix.tab
	awk 'NR>1 {print \$1}' $pangenomeRtab | sort -k 1 -t \$'\t' > INDEX

	awk 'NR>1 {print \$1}' $globalMeanCoverage > sample_names

	while read -r name; do

		echo -e "Gene\tnormalizedCoverage\tcompleteness" > "\${name}"_index.tmp
		grep -w "\$name" $normalized | awk '{print \$2, \$3, \$NF}' >> "\${name}"_index.tmp

	done < sample_names


	for i in *_index.tmp; do

		sed -i -e 's/ /\t/g' "\$i"
		lambda.py "\$i" $gCompleteness $lowerBound $upperBound

	done
	"""
}
