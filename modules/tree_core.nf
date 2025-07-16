

process TREE_CORE {
        conda "${projectDir}/envs/iqtree.yaml"
        
        input:
        path ubiquitousMSA, stageAs: 'maskedMatrixGenesUbiquitousMSA.fasta'
	val threads

        output:
        path 'maskedMatrixGenesUbiquitousMSA.contree', optional: true
	path 'maskedMatrixGenesUbiquitousMSA.iqtree', optional: true
	path 'maskedMatrixGenesUbiquitousMSA.log', optional: true
	path 'maskedMatrixGenesUbiquitousMSA.treefile', optional: true

        script:
        """
	#!/bin/bash

	mkdir -p ${params.output}/TREE/

        iqtree -s maskedMatrixGenesUbiquitousMSA.fasta --prefix maskedMatrixGenesUbiquitousMSA -T $threads -B 1000 -m MFP

	cp maskedMatrixGenesUbiquitousMSA* ${params.output}/TREE/
        """
}
