

process TREE_CORE {
        conda "${projectDir}/envs/iqtree.yaml"
        
        input:
        path ubiquitousMSA, stageAs: 'maskedMatrixGenesUbiquitousMSA.fasta'
	val threads

        output:
        path 'maskedMatrixGenesUbiquitousMSA.contree', emit: maskedMatrixGenesUbiquitousMSAContree
	path 'maskedMatrixGenesUbiquitousMSA.iqtree', emit: maskedMatrixGenesUbiquitousMSAIqtree
	path 'maskedMatrixGenesUbiquitousMSA.log', emit: maskedMatrixGenesUbiquitousMSALog
	path 'maskedMatrixGenesUbiquitousMSA.treefile', emit: maskedMatrixGenesUbiquitousMSATreefile

        script:
        """
        iqtree -s maskedMatrixGenesUbiquitousMSA.fasta --prefix maskedMatrixGenesUbiquitousMSA -T $threads -B 1000 -m MFP
        """
}
