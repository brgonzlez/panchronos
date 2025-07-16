/*
 * GATK_CONSENSUS{} will generate consensus sequences for aligned data using gatk.
 */

process GATK_CONSENSUS {
	conda "${projectDir}/envs/gatk.yaml"

	input:
	tuple path(panGenomeRef), path(panGenomeRefDictionary),	path(panGenomeReferenceIndex)
	path bamFiles
	val parallel

	output:
	path 'extractedSequences*.fasta', emit: gatkConsensusSequences
	path '*GenotypedNormalized.vcf.gz', emit: gatkGenotypes

	script:
	"""
	#!/bin/bash

	mkdir -p ${params.output}/GENOTYPING

	gatkconsensus() {
	bam_file=\$1

		name=\$(basename "\${bam_file%.bam}")
		gatk3 -T UnifiedGenotyper --min_base_quality_score 30 \
			--genotype_likelihoods_model BOTH --annotateNDA \
			--genotyping_mode DISCOVERY --output_mode EMIT_ALL_SITES \
			-I "\$b" -R $panGenomeRef \
			-o "\${name}"Genotyped.vcf 

		bcftools norm -f $panGenomeRef "\${name}"Genotyped.vcf > "\${name}"GenotypedNormalized.vcf
		bgzip -i "\${name}"GenotypedNormalized.vcf
		bcftools index "\${name}"GenotypedNormalized.vcf.gz
		bcftools consensus -a N -M N -f $panGenomeRef "\${name}"GenotypedNormalized.vcf.gz -o "\${name}"GenotypedNormalizedConsensus.fasta
		seqtk seq "\${name}"GenotypedNormalizedConsensus.fasta > extractedSequences"\${name}".fasta
	}
	export -f gatkconsensus
	find ./ -name "*.bam" | parallel -j $parallel gatkconsensus

	cp *GenotypedNormalized.vcf.gz ${params.output}/GENOTYPING
	cp extractedSequences* ${params.output}/GENOTYPING
	"""
}
