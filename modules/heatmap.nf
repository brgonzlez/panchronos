/*
 * MAPPING{} will index a graph and map reads against it.
 */

process HEATMAP {
	conda "${projectDir}/envs/heatmap.yaml"

	input:
	path fCSV, stageAs: 'fCSV/*'
	path INDEX, stageAs: 'INDEX/*'
	path matrix, stageAs: 'matrix/*'
	path names, stageAs: 'names/*'

	output:
	path 'final_matrix.tab', emit: finalMatrix
	path 'presenceAbsence*.png', emit: presenceAbsence
	path 'maskedMatrixGenesOnlyAncient.txt', emit: maskedMatrixGenesOnlyAncient
	path 'maskedMatrixGenesUbiquitous.txt', emit: maskedMatrixGenesUbiquitous
	path 'maskedMatrixGenesNoUbiquitous.txt', emit: maskedMatrixGenesNoUbiquitous
	path 'sampleOrdernoUbiquitous.txt', emit: sampleOrdernoUbiquitous 
	path 'sampleOrderonlyAncient.txt', emit: sampleOrderonlyAncient
	path 'genesAbovePercentSeries.txt', emit: genesAbovePercentSeries
	path 'blackListedQualityChecked.txt', emit: blackListed

	script:
	"""
  	#!/bin/bash

	mkdir ${params.output}/MATRIX
	mkdir ${params.output}/PLOTS

	for i in fCSV/*_final.csv; do

		name=\$(basename "\$i")
		sed  -e 's/,/\t/g' "\$i" | awk 'NR>1{print \$0}' > "\${name%_index.tmp_final.csv}"_INDEX.Z

	done


	for i in *_INDEX.Z; do
		name=\$(basename "\$i")
		#Create the FINAL_INDEX file
		echo "\${name%_INDEX.Z}" > "\${name}"_FINAL_INDEX
    
		#Process the INDEX file
		while read -r gene; do
			toprint=\$(echo "\$gene 0")
			if grep -wq "\$gene" "\$i"; then
				grep -w "\$gene" "\$i" >> "\${name}"_FINAL_INDEX
			else
				echo "\$toprint" >> "\${name}"_FINAL_INDEX
			fi
		done < INDEX/INDEX
    
		#Extract the last column
		awk '{print \$NF}' "\${name}"_FINAL_INDEX > "\${name}"_FINALCOLUMN

	done
        # I need to add a quality control step right here. User samples can be false positives sometimes or just super low quality and have 0 genes after filtering
        # Then, black list unwanted samples and exclude them from the final_matrix.tab document.
        for checkSample in *_FINALCOLUMN; do
                sampleName=\$(basename "\${checkSample%_INDEX.Z_FINALCOLUMN}")
                counts=\$(grep -c "0" "\$checkSample")
                totalLines=\$(wc -l "\$checkSample" | awk '{print \$1 - 1}')
                proportion=\$(awk -v absence="\$counts" -v record="\$totalLines" 'BEGIN { print (absence / record ) }')         

                if  (( \$(awk -v p="\$proportion" 'BEGIN { print (p > 0.95) }' ) )); then
                        echo "\$sampleName" >> blackListedQualityChecked.txt
                fi
        done
	
	# Do this only if blacklisted
	if [[ -s blackListedQualityChecked.txt ]]; then
		while read -r removeMe; do
			mv "\${removeMe}_INDEX.Z_FINALCOLUMN" "\${removeMe}LowQualitySample"    
			grep -v "\${removeMe}" names/sample_names > names/sample_names.tmp
			mv names/sample_names.tmp names/sample_names
		done < blackListedQualityChecked.txt    
	else
		touch blackListedQualityChecked.txt
	
	fi

	paste matrix/matrix.tab *_FINALCOLUMN > final_matrix.tab
	tr '\n' ' ' < names/sample_names > names_heatmap

	heatmap.py final_matrix.tab names_heatmap

	cp final_matrix.tab ${params.output}/MATRIX
	cp *png ${params.output}/PLOTS
	"""
}

