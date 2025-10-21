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

	#trim columns if they only have either n or -.
	trimming_artifacts() {
    fasta_file="\$1"

    	awk '
    	/^>/ {
	        headers[++count] = \$0
        	next
    	}
	
	    {
        	seqs[count] = \$0
        	len = length(\$0)
        	ftrim = 0
        	btrim = 0
	
        	#scan backwards from end
        	for (bpos = len; bpos >= 1; bpos--) {
	            c = substr(\$0, bpos, 1)
            	if (c == "n" || c == "N" || c == "-")
	                btrim++
            	else
	                break
        	}
	
        	#scan forward from start
        	for (fpos = 1; fpos <= len; fpos++) {
	            c = substr(\$0, fpos, 1)
            	if (c == "n" || c == "N" || c == "-")
	                ftrim++
            	else
	                break
        	}
	
        	forward_trim[count] = ftrim
        	backward_trim[count] = btrim
    	}
	
	    END {
        	#find smallest backward trim
        	back_min = backward_trim[1]
        	for (i = 2; i <= count; i++)
	            if (backward_trim[i] < back_min)
                	back_min = backward_trim[i]
	
        	#find smallest forward trim
        	fwd_min = forward_trim[1]
        	for (i = 2; i <= count; i++)
	            if (forward_trim[i] < fwd_min)
                	fwd_min = forward_trim[i]
	
        	if (back_min > 0)
	            print "Backwards trimming", back_min, "characters from", FILENAME > "/dev/stderr"
        	else
	            print "No backwards trimming needed for", FILENAME > "/dev/stderr"
	
        	if (fwd_min > 0)
	            print "Forward trimming", fwd_min, "characters from", FILENAME > "/dev/stderr"
        	else
	            print "No forward trimming needed for", FILENAME > "/dev/stderr"
	
        	#trim the sequences
        	for (i = 1; i <= count; i++) {
	            print headers[i]
            	print substr(seqs[i], 1 + fwd_min, length(seqs[i]) - fwd_min - back_min)
        	}
    	}' "\$fasta_file" > "\${fasta_file%.fasta}_trimmed.fasta"
	}
	export -f trimming_artifacts
	find ./ -name "*.fasta" | parallel -j $parallel trimming_artifacts


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
