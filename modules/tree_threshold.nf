

process TREE_THRESHOLD {
	conda "${projectDir}/envs/iqtree.yaml"
	
	input:
	path genesMSA, stageAs: 'genesAbovePercentMSA.fasta'
	val threads

	output:
	path 'genesAbovePercentMSA.contree', emit: genesAbovePercentMSAContree , optional: true
	path 'genesAbovePercentMSA.iqtree', emit: genesAbovePercentMSAIqtree , optional: true
	path 'genesAbovePercentMSA.log', emit: genesAbovePercentMSALog , optional: true
	path 'genesAbovePercentMSA.treefile', emit: genesAbovePercentMSATreefile , optional: true

	script:
	"""
	#!/bin/bash

	mkdir -p ${params.output}/TREE/


	awk '/^>/ {
    		header = \$0
    		getline seq
    		if (seq ~ /^[nN-]+\$/) {
			print header
    		}
	}' genesAbovePercentMSA.fasta  > to_remove

	awk 'NR==FNR {remove[\$0]; next}
     		/^>/ {keep = !(\$0 in remove)}
	keep' to_remove genesAbovePercentMSA.fasta  > threshold.fasta

	iqtree -s threshold.fasta --prefix threshold -T $threads -B 1000 -m MFP 

	cp threshold* ${params.output}/TREE/
	"""
}
