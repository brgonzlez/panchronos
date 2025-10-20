/*
 * REALIGN_GENE_ALIGNMENTS{} will re-align each gene MSA and perform mandatory and optional QC.
 */


process REALIGN_GENE_ALIGNMENTS {
	conda "${projectDir}/envs/realign.yaml"

	input:
	path gene_msa	
	val parallel
	val threads
	val trim
	path indexes
	val mask_seqs

	output:
	path 'realigned_qc/*fasta', emit: re_aligned

	script:
	"""
	#!/bin/bash

	mkdir -p ${params.output}/GENE_MSA
	mkdir -p realigned
	mkdir -p realigned_qc

	realign() {
	file=\$1
	name=\$(basename "\${file%.fasta}")
		mafft --auto --thread $threads "\$file" | seqtk seq > realigned/"\${name}".fasta
	}
	export -f realign
	find ./ -name "*fasta" | parallel -j $parallel realign


	trimming_artifacts() {
	fasta_file=\$1

		samples_n=\$(grep -c "^>" "\$fasta_file")
		to_trim=\$(awk '!/^>/ { 
								last_character = substr(\$0, length($0), 1)
								if (last_character == "n" || last_character == "-")
									non_nucleotides++
								}

								END {
									print (non_nucleotides ? non_nucleotides : 0)}' "\$fasta_file")
		
	    if [[ "\$to_trim" -eq "\$samples_n" ]]; then
	        awk '{

	        	if (\$0 ~ /^>/) {
	        		print \$0
	        	}

	        	else {
	        		print substr(\$0, 1, length(\$0) - 1)
	        	}
	        }' "\$fasta_file" > trimmed_"\$fasta_file"

	        mv trimmed_"\$fasta_file" ./"\$fasta_file"

	    else
	    	echo "No trimming needed for \$fasta_file"
	    	mv "\$fasta_file" ./realigned_qc/
	    fi

		}
	export -f trimming_artifacts
	find ./realigned/ -name "*fasta" | parallel -j $parallel trimming_artifacts 


	mask_sequences() {
	

	}
	export -f mask_sequences

	if [[ $mask_seqs -eq 1 ]]; then
	    find ./realigned/ -name "*fasta" | parallel -j $parallel mask_sequences
	else
	    echo "--mask_sequences set to 0. No masking will take place."
	fi

	cp -r realigned_qc/* ${params.output}/GENE_MSA
	"""
}
