/*
 * BUILD_MSA{} will concatenate individual gene MSA.
 */

process BUILD_MSA {

	input:
	path genesAlnSeq, stageAs: 'genes/*'
	path maskedMatrixGenesNoUbiquitous, stageAs: 'filesList_maskedMatrixGenesNoUbiquitous.txt'
	path maskedMatrixGenesOnlyAncient, stageAs: 'filesList_maskedMatrixGenesOnlyAncient.txt'
	path maskedMatrixGenesUbiquitous, stageAs: 'filesList_maskedMatrixGenesUbiquitous.txt'
	path genesAbovePercentSeries, stageAs: 'filesList_genesAbovePercentSeries.txt'
	path sampleNames, stageAs: 'sampleNames.txt'
	val parallel

	output:
	path 'genesAbovePercentMSA.fasta', emit: genesAbovePercentMSA
	path 'maskedMatrixGenesNoUbiquitousMSA.fasta', emit: maskedMatrixGenesNoUbiquitousMSA
	path 'maskedMatrixGenesOnlyAncientMSA.fasta', emit: maskedMatrixGenesOnlyAncientMSA
	path 'maskedMatrixGenesUbiquitousMSA.fasta', emit: maskedMatrixGenesUbiquitousMSA
	
	script:
	"""
	#!/bin/bash

	touch maskedMatrixGenesNoUbiquitous.fasta
	touch maskedMatrixGenesOnlyAncient.fasta
	touch maskedMatrixGenesUbiquitous.fasta
	touch genesAbovePercentSeries.fasta

	sed -i -e 's/~/_/g' $genesAbovePercentSeries
	sed -i -e 's/~/_/g' $maskedMatrixGenesNoUbiquitous
	sed -i -e 's/~/_/g' $maskedMatrixGenesOnlyAncient
	sed -i -e 's/~/_/g' $maskedMatrixGenesUbiquitous

	rename() {
	file=\$1

		mv "\${file}" "\${file%_trimmed.fasta}.fasta"
	}
	export -f rename
	find ./genes/ -name "*.fasta" | parallel -j $parallel rename

	build_msa() {
	inputfile=\$1
	filename=\$(basename "\${inputfile}")
	filename="\${filename#filesList_}"
	filename="\${filename%.txt}"

		# Build a list of existing files
		files=()
		while read -r gene; do
    		file="genes/\${gene}.fasta"
    		if [[ -f "\$file" ]]; then
        		files+=("\$file")
    		else
        		echo "Warning: File \$file not found, skipping..."
    		fi
		done < "\$inputfile"

		#append genes array to MSA into array
		files=("\$filename".fasta "\${files[@]}")

		#paste everything inside array
		paste "\${files[@]}" | tr -d '\t' > TMP_"\${filename}"
		mv TMP_"\${filename}" "\${filename}".fasta

		#Clean headers (if needed)
		awk -F'>' '/^>/ {print ">" \$2} !/^>/' "\${filename}".fasta > TMPg_"\${filename}"
		mv TMPg_"\${filename}" "\${filename}".fasta

	}
	export -f build_msa
	find ./ -name "filesList_*" | parallel -j $parallel build_msa

	#rename for now
	mv maskedMatrixGenesNoUbiquitous.fasta maskedMatrixGenesNoUbiquitousMSA.fasta
	mv maskedMatrixGenesOnlyAncient.fasta maskedMatrixGenesOnlyAncientMSA.fasta
	mv maskedMatrixGenesUbiquitous.fasta maskedMatrixGenesUbiquitousMSA.fasta
	mv genesAbovePercentSeries.fasta genesAbovePercentMSA.fasta
	"""
}
