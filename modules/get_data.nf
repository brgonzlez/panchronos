
/*
 * GET_DATA() will download FASTA and GenBank files from NCBI as long as a taxonomical ID is provided.
 */

process GET_DATA {
	conda "${projectDir}/envs/get_data.yaml"

	input:
	val genomes
	val tax_id
	val parallel

	output:
	path '*fna', emit: fasta_files
	path '*gb', emit: gbk_files

	script:
	"""
	#!/bin/bash

	#get accessions list file
    	query="txid\${tax_id}[Organism] AND complete genome[Title]" 
        esearch -db nuccore -query "\$query" | efetch -format uid | head -n $genomes | efetch -db nuccore -format acc >> accessions.txt
    		
	mkdir -p accessions

	while read -r accession; do

		touch accessions/"\$accession".map

	done < accessions.txt


	get_data() {
	sample=\$1
	ACCESSION=\$(echo "\${sample%.map}")

		efetch -db nuccore -id "\$ACCESSION" -format gbwithparts -mode text > "\${ACCESSION}.gb"

		efetch -db nucleotide -id "\$ACCESSION" -format fasta > "\${ACCESSION}.fasta"

		# TEst integrity without actually compressing them
		if gzip -c "\${ACCESSION}.gb" > /dev/null && gzip -c "\${ACCESSION}.fasta" > /dev/null; then
			echo "\${ACCESSION} files passed integrity check"
		else
			echo "Corrupted files detected for \${ACCESSION}. Retrying..."
			rm -f "\${ACCESSION}.gb" "\${ACCESSION}.fasta"
			sleep 3
		fi
	}
	export -f getData
	find ./accessions/ -name "*.map" | parallel -j $parallel get_data

	mv ./accessions/*fasta ./
	mv ./accessions/*gb ./

	cat .command.out >> get_data.log
	"""
}
