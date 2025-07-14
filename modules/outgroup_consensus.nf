/*
 * OUTGROUP_CONSENSUS{} will generate outgroup consensus sequences for MSA concatenating.
 */


process OUTGROUP_CONSENSUS {
	conda "${projectDir}/envs/consensusOutgroup.yaml"

	input:
	path outgroupFastaPostAlignment
	path panGenomeRef

	output:
	path 'extractedSequencesOutgroup.fasta', emit: extractedSequencesOutgroupFasta
	path 'extractedSequencesOutgroup.fq', emit: extractedSequencesOutgroupFastq

	script:
	"""
	bcftools mpileup -f $panGenomeRef $outgroupFastaPostAlignment | bcftools call -c | vcfutils.pl vcf2fq > extractedSequencesOutgroup.fq
	seqtk seq -a extractedSequencesOutgroup.fq > extractedSequencesOutgroup.fasta

	"""
}

