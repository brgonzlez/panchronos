
process TREE_ACCESSORY {
        conda "${projectDir}/envs/iqtree.yaml"
        
        input:
        path noUbiquitousMSA, stageAs: 'maskedMatrixGenesNoUbiquitousMSA.fasta'
	val threads

        output:
        path 'maskedMatrixGenesNoUbiquitousMSA.contree', emit: maskedMatrixGenesNoUbiquitousMSAContree
	path 'maskedMatrixGenesNoUbiquitousMSA.iqtree', emit: maskedMatrixGenesNoUbiquitousMSAIqtree
	path 'maskedMatrixGenesNoUbiquitousMSA.log', emit: maskedMatrixGenesNoUbiquitousMSALog
	path 'maskedMatrixGenesNoUbiquitousMSA.treefile', emit: maskedMatrixGenesNoUbiquitousMSATreefile

        script:
        """
        iqtree -s maskedMatrixGenesNoUbiquitousMSA.fasta --prefix maskedMatrixGenesNoUbiquitousMSA -T $threads -B 1000 -m MFP
        """
}
