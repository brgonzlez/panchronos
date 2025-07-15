

process TREE_THRESHOLD {
	conda "${projectDir}/envs/iqtree.yaml"
	
	input:
	path genesMSA, stageAs: 'genesAbovePercentMSA.fasta'
	val threads

	output:
	path 'genesAbovePercentMSA.contree', emit: genesAbovePercentMSAContree
	path 'genesAbovePercentMSA.iqtree', emit: genesAbovePercentMSAIqtree
	path 'genesAbovePercentMSA.log', emit: genesAbovePercentMSALog
	path 'genesAbovePercentMSA.treefile', emit: genesAbovePercentMSATreefile

	script:
	"""
	iqtree -s genesAbovePercentMSA.fasta --prefix genesAbovePercentMSA -T $threads -B 1000 -m MFP 
	"""
}
