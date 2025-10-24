/*
 * OUTGROUP_CONSENSUS{} will generate outgroup consensus sequences for MSA concatenating.
 */


process OUTGROUP_CONSENSUS {
	conda "${projectDir}/envs/consensus.yaml"

	input:
	path outgroupFastaPostAlignment
	path panGenomeRef

	output:
	path 'outgroup_genotyped.fasta', emit: extractedSequencesOutgroupFasta

	script:
	"""
	#!/bin/bash

	mkdir -p ${params.output}/GENOTYPING

	outgroupConsensus() {
        bam_file=\$1
        basename=\$(basename "\${bam_file%.bam}")

                bcftools mpileup -f $panGenomeRef -q 30 -Q 20 "\$bam_file" > "\${basename}"_mpileup_file
                bcftools call -c "\${basename}"_mpileup_file > "\${basename}".vcf
                bgzip -i -c "\${basename}".vcf > "\${basename}".vcf.gz
                bcftools index "\${basename}".vcf.gz
                bcftools consensus -a - -f $panGenomeRef "\${basename}".vcf.gz > extractedSequences"\${basename}".fq
                seqtk seq -a extractedSequences"\${basename}".fq > outgroup_genotyped.fasta
                rm -f extractedSequences"\${basename}".fq

		mv "\${basename}".vcf.gz ./outgroup.vcf.gz
	}
 	export -f outgroupConsensus
  	find ./ -name "*.bam" | parallel -j 1 outgroupConsensus

	cp outgroup.vcf.gz ${params.output}/GENOTYPING
	cp outgroup_genotyped.fasta ${params.output}/GENOTYPING
	"""
}

