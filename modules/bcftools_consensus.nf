/*
 * BCFTOOLS_CONSENSUS{} will generate consensus sequences for each gene from aligned data using bcftools.
 */
 
process BCFTOOLS_CONSENSUS {
	conda "${projectDir}/envs/consensus.yaml"

	input:
	path panGenomeRef
	path bamFiles
	val parallel

	output:
	path 'extractedSequences*.fasta', emit: consensusSequences

	script:
	"""
	#!/bin/bash

  	bcfconsensus() {
   	bam_file=\$1
    
		basename=\$(basename "\${bam_file%.bam}")
		bcftools mpileup -f $panGenomeRef "\$bam_file" | bcftools call -c | vcfutils.pl vcf2fq > extractedSequences"\${basename}".fq
		seqtk seq -a extractedSequences"\${basename}".fq > extractedSequences"\${basename}".fasta
	}
 	export -f bcfconsensus
  	find ./ -name "*.bam" | parallel -j $parallel bcfconsensus
	"""
}
