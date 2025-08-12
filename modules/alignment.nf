/*
 * ALIGNMENT{} will align user data agianst pangenome reference sequences.
 */

process ALIGNMENT {
	conda "${projectDir}/envs/alignment.yaml"

	input:
	path data
	path panRef
	path configFile
	tuple val(threadsGlobal), val(missingProb), val(seedAlignment), val(gapFraction),val(minReadLength),val(maxReadLength),val(parallel)

	output:
	path '*_aligned.bam', emit: postAlignedBams
	path '*_final.fastq', emit: postAlignedReads
	path '*_sorted_mappedreads.bam', emit: bam_mapdamage
	path 'panGenomeReference_extended.fasta.*' , emit: pan_index

	script:
	"""
	#!/bin/bash

	mkdir -p ${params.output}/ALIGNMENT

	bwa index $panRef

	#align() will align the data and perform several post-alignment computations
	align() {
	sample=\$1
    	
		name=\$(basename "\${sample%.fastq*}")
    		
    		softClip=\$(grep "\$name" $configFile | awk '{print \$2}')
       		# Making read groups
    		rg_id="\${name}"  # sample name as id
    		rg_sm="\${name}" # sample name again
    		rg_pl="illumina"        # I dont think this is very important for this pipeline so its going to be just illumina because why not
    		rg_lb="lib1"            # group id
    		rg_pu="unit1"           # not sure what Ill put here

    		bwa aln -l $seedAlignment -n $missingProb -o $gapFraction -t $threadsGlobal $panRef "\$sample" > "\${name}.sai"
    		bwa samse -r "@RG\\tID:\$rg_id\\tSM:\$rg_sm\\tPL:\$rg_pl\\tLB:\$rg_lb\\tPU:\$rg_pu" \
		$panRef "\${name}.sai" "\$sample" > "\${name}.sam"
    
    		samtools view -bS "\${name}.sam" > "\${name}.bam"
    		samtools quickcheck "\${name}.bam"
    		samtools sort -o "\${name}_sorted.bam" -O bam -@ $threadsGlobal "\${name}.bam"
    		samtools index "\${name}_sorted.bam"
    		rm "\${name}.bam"
    		samtools view -b -@ $threadsGlobal -F 4 "\${name}_sorted.bam" > "\${name}_sorted_mappedreads.bam"
    		samtools index "\${name}_sorted_mappedreads.bam"
    		bam trimBam "\${name}_sorted_mappedreads.bam" "\${name}_softclipped.bam" -L "\$softClip" -R "\$softClip" --clip
    		samtools view -q 30 -o "\${name}_qc.bam" "\${name}_softclipped.bam"
    		samtools view -e 'length(seq)>$minReadLength && length(seq)<$maxReadLength' -O BAM -o "\${name}_lg.bam" "\${name}_qc.bam"
    		samtools sort -o "\${name}_aligned.bam" -O bam -@ $threadsGlobal "\${name}_lg.bam"
    		samtools coverage "\${name}_aligned.bam" > "\${name}_genomicsMetrics.txt"
    		samtools fastq -@ $threadsGlobal "\${name}_aligned.bam" > "\${name}_final.fastq"
	}
	export -f align
	find $data/* -name "*.fastq*" | parallel -j $parallel align


	rm *sam *sai *_lg.bam *_qc.bam 

	cp *_aligned.bam ${params.output}/ALIGNMENT
	cp *_final.fastq ${params.output}/ALIGNMENT

	cat .command.out >> alignment.log
	"""
}
