/*
 * OUTGROUP_ALIGNMENT{} will align artificial outgroup reads against pangenome reference sequences.
 */

process OUTGROUP_ALIGNMENT {
	conda "${projectDir}/envs/alignment.yaml"
	
	input:
	path outgroupReads
	path panGenomeRef
	val threads

	output:
	path 'outgroup_aligned.bam', emit: outgroupFastaPostAlignment

	script:
	"""
	#!/bin/bash

	bwa index $panGenomeRef
	bwa mem -B 1 -E 1 $panGenomeRef $outgroupReads -t $threads > outgroupFasta.sam
	samtools view -bS outgroupFasta.sam > outgroupFasta.bam
	samtools quickcheck outgroupFasta.bam
	samtools sort -o outgroupFastaSorted.bam -O bam -@ $threads outgroupFasta.bam
	samtools index outgroupFastaSorted.bam
	samtools view -b -@ 10 -F 4 outgroupFastaSorted.bam > outgroupFastaSortedMappedreads.bam
	samtools index outgroupFastaSortedMappedreads.bam
	samtools sort -o outgroup_aligned.bam -O bam -@ $threads outgroupFastaSortedMappedreads.bam
	"""
}
