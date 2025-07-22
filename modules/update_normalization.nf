/*
 * UPDATE_NORMALIZATION{} will .
 */

process UPDATE_NORMALIZATION {

	input:
	path normalized
	path completeness
	val parallel

	output:
	path 'geneNormalizedUpdated.tab', emit: geneNormalizedUpdated

	script:
	"""
	#!/bin/bash

	mkdir -p ${params.output}/STATS

	echo -e "sampleID\tgene\tnormalizedGeneSimple\tnormalizedGeneScaled\tnormalizedGenomeSimple\tnormalizedGenomeScaled\tgeneCompleteness" > geneNormalizedUpdated.tab
	sed -i -e 's/~/_/g' geneNormalizedSummary.txt

	awk '{print \$1}' $completeness | uniq > samples.txt

	while read -r sample; do
		echo "\$sample" > "\$sample".map
	done < samples.txt
	
	normalize_updating() {
	file=\$1
	name=\$(basename "\${file%.map}")

		grep -w "\$name" geneNormalizedSummary.txt > "\$name"_geneNormalizedSummary
		grep -w "\$name" completenessSummary.tab > "\$name"_completenessSummary

		awk 'NR>1{print \$1"XYZ"\$2, \$3, \$4, \$5, \$6}' "\$name"_geneNormalizedSummary > "\$name"_TMP1

		awk '{print \$1"XYZ"\$2, \$3}' "\$name"_completenessSummary > "\$name"_TMP2

		while read -r ID completeness;do

			if grep -wq "\${ID}" "\$name"_TMP1; then
				oldLine=\$(grep -w "\${ID}" "\$name"_TMP1)
				specificCompleteness=\$(grep -w "\${ID}" "\$name"_TMP2 | awk '{print \$NF}')
				echo -e "\${oldLine}\t\${specificCompleteness}" >> "\$name"_geneNormalizedUpdated.tab
			fi
	
		done < "\$name"_TMP2
	
		sed -i -e 's/XYZ/\t/g' "\$name"_geneNormalizedUpdated.tab
		sed -i -e 's/ /\t/g' "\$name"_geneNormalizedUpdated.tab
	}
	export -f normalize_updating
	find ./ -name "*.map" | parallel -j $parallel normalize_updating

	cat *_geneNormalizedUpdated.tab >> geneNormalizedUpdated.tab

	cp geneNormalizedUpdated.tab ${params.output}/STATS
	"""
}
