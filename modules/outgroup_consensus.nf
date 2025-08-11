/*
 * OUTGROUP_CONSENSUS{} will generate outgroup consensus sequences for MSA concatenating.
 */


process OUTGROUP_CONSENSUS {
	conda "${projectDir}/envs/consensus.yaml"

	input:
	path outgroupFastaPostAlignment
	path panGenomeRef
	val extension

	output:
	path 'extractedSequencesOutgroup.fasta', emit: extractedSequencesOutgroupFasta
	path 'extractedSequencesOutgroup.fq', emit: extractedSequencesOutgroupFastq

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
                bcftools consensus -a N -f $panGenomeRef "\${basename}".vcf.gz > extractedSequences"\${basename}".fq
                seqtk seq -a extractedSequences"\${basename}".fq > extractedSequences"\${basename}".fasta
                rm -f extractedSequences"\${basename}".fq

		#now we need to do padding to make sure consensus sequences lengths are the same as in the reference genome
		awk '
			FNR==NR {
				if(/^>/) {
					header = \$0
					getline seq
					gene[header] = length(seq)
				}
				next
			}
			{
				if (/^>/) {
					print
        				current_header = \$0
        				next
    				} else {
        				ref_len = gene[current_header]
        				while (length(\$0) < ref_len) {
            					\$0 = \$0 "n"
        				}
        				print
    				}
			}
		' $panGenomeRef extractedSequences"\${basename}".fasta > "\${basename}"_padded && mv "\${basename}"_padded ./extractedSequences"\${basename}".fasta

		#now remove the extended sequences
		awk -v trim=$extension '
			/^>/ { 
				print 
				next	
			}

			!/^>/ { 	
				len = length(\$0)
    				print substr(\$0, trim + 1, len - 2 * trim)
			}' extractedSequences"\${basename}".fasta > trimmed_"\${basename}" && mv trimmed_"\${basename}" ./outgroup_genotyped.fasta

		mv "\${basename}".vcf.gz ./outgroup.vcf.gz
	}
 	export -f outgroupConsensus
  	find ./ -name "*.bam" | parallel -j 1 outgroupConsensus

	cp outgroup.vcf.gz ${params.output}/GENOTYPING
	cp outgroup_genotyped.fasta ${params.output}/GENOTYPING
	"""
}

