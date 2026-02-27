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
        samtools view -bS outgroupFasta.sam > outgroupFasta.bam
        samtools quickcheck outgroupFasta.bam
        samtools sort -o outgroupFastaSorted.bam -O bam -@ $threads outgroupFasta.bam
        picard MarkDuplicates OPTICAL_DUPLICATE_PIXEL_DISTANCE=100 REMOVE_DUPLICATES=TRUE I=outgroupFastaSorted.bam  O=outgroup_deduped.bam M=outgroup_deduped.stats
        samtools index outgroup_deduped.bam
        samtools view -b -@ 10 -F 4 outgroup_deduped.bam > outgroupFastaSortedMappedreads.bam
        samtools index outgroupFastaSortedMappedreads.bam
        samtools sort -o outgroup_aligned.bam -O bam -@ $threads outgroupFastaSortedMappedreads.bam

        cp outgroup_aligned.bam ${params.output}/ALIGNMENT
        """
}
