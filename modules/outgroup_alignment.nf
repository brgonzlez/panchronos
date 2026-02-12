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

	mkdir -p ${params.output}/ALIGNMENT

	bwa index $panGenomeRef
	bwa mem -B 1 -E 1 $panGenomeRef $outgroupReads -t $threads > outgroupFasta.sam
	samtools view -@ $threads -bS outgroupFasta.sam > outgroupFasta.bam
	samtools quickcheck outgroupFasta.bam
	samtools sort -o outgroupFastaSorted.bam -O bam -@ $threads outgroupFasta.bam
	samtools index outgroupFastaSorted.bam
	samtools view -b -@ $threads -F 4 outgroupFastaSorted.bam > outgroupFastaSortedMappedreads.bam
	samtools index -@ $threads outgroupFastaSortedMappedreads.bam
	samtools sort -o outgroup_aligned.bam -O bam -@ $threads outgroupFastaSortedMappedreads.bam

	cp outgroup_aligned.bam ${params.output}/ALIGNMENT
	"""
}
