/*
 * BUILD_MSA{} will index a graph and map reads against it.
 */

process BUILD_MSA {

	input:
	path genesAlnSeq, stageAs: 'genes/*'
	path maskedMatrixGenesNoUbiquitous, stageAs: 'maskedMatrixGenesNoUbiquitous.txt'
	path maskedMatrixGenesOnlyAncient, stageAs: 'maskedMatrixGenesOnlyAncient.txt'
	path maskedMatrixGenesUbiquitous, stageAs: 'maskedMatrixGenesUbiquitous.txt'
	path genesAbovePercentSeries, stageAs: 'genesAbovePercentSeries.txt'
	path sampleNames, stageAs: 'sampleNames.txt'

	output:
	path 'genesAbovePercentMSA.fasta', emit: genesAbovePercentMSA
	path 'maskedMatrixGenesNoUbiquitousMSA.fasta', emit: maskedMatrixGenesNoUbiquitousMSA
	path 'maskedMatrixGenesOnlyAncientMSA.fasta', emit: maskedMatrixGenesOnlyAncientMSA
	path 'maskedMatrixGenesUbiquitousMSA.fasta', emit: maskedMatrixGenesUbiquitousMSA
	

	script:
	"""
	#!/bin/bash

	touch genesAbovePercentMSA.fasta
	touch maskedMatrixGenesUbiquitousMSA.fasta
	touch maskedMatrixGenesOnlyAncientMSA.fasta
	touch maskedMatrixGenesNoUbiquitousMSA.fasta

	sed -i -e 's/~/_/g' genesAbovePercentSeries.txt
	sed -i -e 's/~/_/g' maskedMatrixGenesNoUbiquitous.txt
	sed -i -e 's/~/_/g' maskedMatrixGenesOnlyAncient.txt
	sed -i -e 's/~/_/g' maskedMatrixGenesUbiquitous.txt



	while read -r gene; do
		file="genes/\${gene}.fasta"
		if [[ -f "\$file" ]] ; then
			paste genesAbovePercentMSA.fasta "\$file" > TMP; mv TMP genesAbovePercentMSA.fasta
		else
			echo "Warning: File \$file not found, skipping..."
		fi
	done < genesAbovePercentSeries.txt
	sed -i -e 's/\t//g' genesAbovePercentMSA.fasta
	awk -F'>' '/^>/ {print ">" \$2} !/^>/' genesAbovePercentMSA.fasta > TMPg; mv TMPg genesAbovePercentMSA.fasta



	while read -r gene; do
		file="genes/\${gene}.fasta"
		if [[ -f "\$file" ]] ; then
			paste maskedMatrixGenesNoUbiquitousMSA.fasta "\$file" > TMP; mv TMP maskedMatrixGenesNoUbiquitousMSA.fasta
		else
			echo "Warning: File \$file not found, skipping..."
		fi
	done < maskedMatrixGenesNoUbiquitous.txt
        sed -i -e 's/\t//g' maskedMatrixGenesNoUbiquitousMSA.fasta
        awk -F'>' '/^>/ {print ">" \$2} !/^>/' maskedMatrixGenesNoUbiquitousMSA.fasta > TMPnU; mv TMPnU maskedMatrixGenesNoUbiquitousMSA.fasta



	while read -r gene; do
		file="genes/\${gene}.fasta"
		if [[ -f "\$file" ]] ; then
			paste maskedMatrixGenesOnlyAncientMSA.fasta "\$file" > TMP; mv TMP maskedMatrixGenesOnlyAncientMSA.fasta
		else
			echo "Warning: File \$file not found, skipping..."
		fi
	done < maskedMatrixGenesOnlyAncient.txt
        sed -i -e 's/\t//g' maskedMatrixGenesOnlyAncientMSA.fasta
        awk -F'>' '/^>/ {print ">" \$2} !/^>/' maskedMatrixGenesOnlyAncientMSA.fasta > TMPo2; mv TMPo2 maskedMatrixGenesOnlyAncientMSA.fasta

	while read -r gene; do
		file="genes/\${gene}.fasta"
		if [[ -f "\$file" ]] ; then
			paste maskedMatrixGenesUbiquitousMSA.fasta "\$file" > TMP; mv TMP maskedMatrixGenesUbiquitousMSA.fasta
		else
			echo "Warning: File \$file not found, skipping..."
		fi
	done < maskedMatrixGenesUbiquitous.txt

        sed -i -e 's/\t//g' maskedMatrixGenesUbiquitousMSA.fasta
        awk -F'>' '/^>/ {print ">" \$2} !/^>/' maskedMatrixGenesUbiquitousMSA.fasta > TMPU2; mv TMPU2 maskedMatrixGenesUbiquitousMSA.fasta

	"""
}
