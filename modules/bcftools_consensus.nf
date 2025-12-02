/*
 * BCFTOOLS_CONSENSUS{} will generate consensus sequences for each gene from aligned data using bcftools.
 */
 
process BCFTOOLS_CONSENSUS {
	conda "${projectDir}/envs/consensus.yaml"

	input:
	path panGenomeRef
	path bamFiles
	val parallel
	tuple val(mapq) , val(baseq) , val(call_qual)
	val extension
	path genes_2_mask
	
	output:
	path 'extractedSequences*.fasta', emit: consensusSequences

	script:
	"""
	#!/bin/bash

	mkdir -p ${params.output}/GENOTYPING
    bcfconsensus() {
    bam_file=\$1
    basename=\$(basename "\${bam_file%.bam}")

                bcftools mpileup -f $panGenomeRef -q $mapq -Q $baseq "\$bam_file" > "\${basename}"_mpileup_file
                bcftools call --ploidy 1 -m "\${basename}"_mpileup_file > "\${basename}"_raw.vcf
				bcftools filter -i 'QUAL>$call_qual' "\${basename}"_raw.vcf > "\${basename}".vcf
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
			}' extractedSequences"\${basename}".fasta > trimmed_"\${basename}" && mv trimmed_"\${basename}" ./extractedSequences"\${basename}".fasta

	}
    export -f bcfconsensus
    find ./ -name "*.bam" | parallel -j $parallel bcfconsensus


	#Now we mask those sus genes
	

	cp extractedSequences* ${params.output}/GENOTYPING
	cp *.vcf.gz* ${params.output}/GENOTYPING
	"""
}
