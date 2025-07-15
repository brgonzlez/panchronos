
process TREE_ANCIENT {
        conda "${projectDir}/envs/iqtree.yaml"
        
        input:
        path ancientMSA, stageAs: 'maskedMatrixGenesOnlyAncientMSA.fasta'
	val threads

        output:
        path 'maskedMatrixGenesOnlyAncientMSA.contree', emit: maskedMatrixGenesOnlyAncientMSAConqtree
	path 'maskedMatrixGenesOnlyAncientMSA.iqtree', emit: maskedMatrixGenesOnlyAncientMSAIqtree
	path 'maskedMatrixGenesOnlyAncientMSA.log', emit: maskedMatrixGenesOnlyAncientMSALog
	path 'maskedMatrixGenesOnlyAncientMSA.treefile', emit: maskedMatrixGenesOnlyAncientMSATreefile

        script:
        """
	#!/bin/bash

        iqtree -s maskedMatrixGenesOnlyAncientMSA.fasta --prefix maskedMatrixGenesOnlyAncientMSA -T $threads -B 1000 -m MFP
        """
}
