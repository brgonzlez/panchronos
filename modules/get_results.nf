
process GET_RESULTS {

	input:
	path makeDir
	path checkedFastas, stageAs: 'checkedFasta/*'
	path checkedGffs, stageAs: 'checkedGff/*'
	path fastaDatabaseLog, stageAs: 'fastaDatabase.log'
	path fastaDatabaseSeqs, stageAs: 'clusteredSequences.fasta'
	path clusteredSequences, stageAs: 'clusteredNonRedundantGenes.fasta'
	path clusteringLogFile, stageAs: 'clustering.log'
	path prokkaGff, stageAs: 'prokkaGff/*'
	path prokkaLog, stageAs: 'prokka.log'
	path panarooLog, stageAs: 'makePangenome.log'
	path genesMSA, stageAs: 'geneMSA/*'
	path panrefG, stageAs: 'pangenomeReferenceGenome.fasta'
	path geneNormalizedUpdated, stageAs: 'geneNormalizedUpdated.tab'
	path globalMeanCoverage, stageAs: 'globalMeanCoverage.tab'
	path PostAlignmentFiles, stageAs: '*'
	path RefLenghts, stageAs: '*'
	path RawCoverage, stageAs: '*'
	path CompletenessSummary, stageAs: 'completenessSummary.tab'
	path FinalMatrix, stageAs: 'finalMatrix.tab'
	path presenceAbsenceplots, stageAs: '*'
	path MaskedMatrixGenesOnlyAncient, stageAs: 'maskedMatrixGenesOnlyAncient.txt'
	path MaskedMatrixGenesUbiquitous, stageAs: 'maskedMatrixGenesUbiquitous.txt'
	path MaskedMatrixGenesNoUbiquitous, stageAs: 'maskedMatrixGenesNoUbiquitous.txt'
	path GenesAbovePercentSeries, stageAs: 'genesAbovePercentSeries.txt'
	path GenesAbovePercentMSAIqtree, stageAs: 'genesAbovePercentMSA.iqtree'
	path GenesAbovePercentMSALog, stageAs: 'genesAbovePercentMSA.log'
	path GenesAbovePercentMSATreefile, stageAs: 'genesAbovePercentMSA.treefile'
	path MaskedMatrixGenesUbiquitousMSAIqtree, stageAs: 'maskedMatrixGenesUbiquitousMSA.iqtree'
	path MaskedMatrixGenesUbiquitousMSALog, stageAs: 'maskedMatrixGenesUbiquitousMSA.log'
	path MaskedMatrixGenesUbiquitousMSATreefile, stageAs: 'maskedMatrixGenesUbiquitousMSA.treefile'
	path MaskedMatrixGenesNoUbiquitousMSAIqtree, stageAs: 'maskedMatrixGenesNoUbiquitousMSA.iqtree'
	path MaskedMatrixGenesNoUbiquitousMSALog, stageAs: 'maskedMatrixGenesNoUbiquitousMSA.log'
	path MaskedMatrixGenesNoUbiquitousMSATreefile, stageAs: 'maskedMatrixGenesNoUbiquitousMSA.treefile'
	path MaskedMatrixGenesOnlyAncientMSAIqtree, stageAs: 'maskedMatrixGenesOnlyAncientMSA.iqtree'
	path MaskedMatrixGenesOnlyAncientMSALog, stageAs: 'maskedMatrixGenesOnlyAncientMSA.log'
	path MaskedMatrixGenesOnlyAncientMSATreefile, stageAs: 'maskedMatrixGenesOnlyAncientMSA.treefile'

	output:
	stdout

	script:
	"""
	#!/bin/bash

	mkdir -p "${makeDir}/modernData/"

	mv checkedFasta/* "${makeDir}/modernData/"
	mv checkedGff/* "${makeDir}/modernData/"
	mv fastaDatabase.log "${makeDir}/modernData/"

	mkdir -p "${makeDir}/clusteredSequences"
	mv clusteredSequences.fasta "${makeDir}/clusteredSequences/"
	mv clusteredNonRedundantGenes.fasta "${makeDir}/clusteredSequences/"
	mv clustering.log "${makeDir}/clusteredSequences/"

	mkdir -p "${makeDir}/prokkaResults/"

	mv prokkaGff/* "${makeDir}/prokkaResults/"
	mv prokka.log "${makeDir}/prokkaResults/"

	mkdir -p "${makeDir}/pangenomeFiles"

	mv makePangenome.log "${makeDir}/pangenomeFiles/"
	mv geneMSA/* "${makeDir}/pangenomeFiles/"
	mv pangenomeReferenceGenome.fasta "${makeDir}/pangenomeFiles/"

	mkdir -p "${makeDir}/alignmentResults"

	mv geneNormalizedUpdated.tab "${makeDir}/alignmentResults/"
	mv globalMeanCoverage.tab "${makeDir}/alignmentResults/"
	mv *_refLength.txt "${makeDir}/alignmentResults/"
	mv *Coverage.txt "${makeDir}/alignmentResults/"
	mv *bam "${makeDir}/alignmentResults/"
	mv completenessSummary.tab "${makeDir}/alignmentResults/"

	mkdir -p "${makeDir}/matrixResults"
	mv finalMatrix.tab "${makeDir}/matrixResults"

	mkdir -p "${makeDir}/plotsResults"
	mv *png "${makeDir}/plotsResults"

	mkdir -p "${makeDir}/MSAs"
	mv maskedMatrixGenesUbiquitous.txt "${makeDir}/MSAs/"
	mv maskedMatrixGenesOnlyAncient.txt "${makeDir}/MSAs/"
	mv maskedMatrixGenesNoUbiquitous.txt "${makeDir}/MSAs/"
	mv genesAbovePercentSeries.txt "${makeDir}/MSAs/"
	mv *MSA* "${makeDir}/MSAs/"

	"""
}
