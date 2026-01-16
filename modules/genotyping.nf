/*
 * GENOTYPING{} will generate consensus sequences for each gene from aligned data using bcftools.
 */

process GENOTYPING {
        conda "${projectDir}/envs/consensus.yaml"

        input:
        path panGenomeRef
        path bamFiles
        val parallel
        tuple val(mapq) , val(baseq) , val(call_qual)
        val extension
		val force_homozygot

        output:
        path 'extractedSequences*.fasta', emit: consensusSequences
		path '*_per_site.txt', emit: per_site
		path '*_per_gene_and_global.txt', emit: per_gene_and_global

        script:
        """
        #!/bin/bash

        mkdir -p ${params.output}/GENOTYPING
        bcfconsensus() {
        bam_file=\$1
        basename=\$(basename "\${bam_file%.bam}")

                bcftools mpileup -f $panGenomeRef -q $mapq -Q $baseq "\$bam_file" > "\${basename}"_mpileup_file
                bcftools call --ploidy 1 -m "\${basename}"_mpileup_file > "\${basename}"_raw.vcf

                if [[ $force_homozygot -eq 1 ]]; then
                        bcftools filter -i 'QUAL>$call_qual && ((DP4[2]+DP4[3]==0) || (DP4[0]+DP4[1]==0))' "\${basename}"_raw.vcf > "\${basename}".vcf
                else
                        bcftools filter -i 'QUAL>$call_qual' "\${basename}"_raw.vcf > "\${basename}".vcf
                fi

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

        #Calculate genome-wide heterozygosity per sample as a measure of contamination, called heteroplasmy analysis

        #remove intermediate files
        rm *_raw.vcf

        heteroplasmy() {
        file=\$1
        name=\$(basename "\${file%.vcf}")

        # first I need to process the VCF file and only get the fields of interest. I need to include the pos field to remove extended sequences.

         awk -F";" 'BEGIN {OFS="\t"} !/^#/ {print \$0}' "\$file" | sed -e 's/;/\t/g' -e 's/DP4=//g' | awk '{print \$1, \$2, \$(NF-3)}' | sed -e 's/,/\t/g' -e 's/ /\\t/g' > "\$name"_vcf_first_filter

        #note to myself. condition ? value_if_true : value_if_false

        awk -v to_trim=$extension '
        {
            gene = \$1
            pos  = \$2
            maxpos[gene] = (pos > maxpos[gene] ? pos : maxpos[gene])
            lines[NR] = \$0
            genes[NR] = gene
            poss[NR]  = pos
        }
        END {
            for (i = 1; i <= NR; i++) {
                g = genes[i]
                p = poss[i]
                if (p > to_trim && p <= maxpos[g] - to_trim) {
                    print lines[i]
                }
            }
        }
        ' "\$name"_vcf_first_filter > "\$name"_vcf_first_filter_no_extd


        #Now we can run the AWK script to calculate genome-wide heteroplasmy.


        awk 'BEGIN{OFS="\t"}
        {
            #per site
            gene = \$1
            pos = \$2
            geneList[gene] = 1
            geneLength[gene] += 1   #simply sum up +1 for every line for that gene
            RefSite = \$3 + \$4
            Altsite = \$5 + \$6
            heteroplasmySite = RefSite - Altsite
            if (heteroplasmySite < 0) heteroplasmySite = -heteroplasmySite
            allreadsSite = RefSite + Altsite

            hetSite = (heteroplasmySite / allreadsSite) * 100
            print gene, pos, hetSite

            #per gene
            geneCumulativePersite[gene] = geneCumulativePersite[gene] + hetSite

            #global, sum up everything
            globalHetsites +=  hetSite
        }

        END {

            globalHet = (globalHetsites/NR)
            print "Global heteroplasmy:", globalHet

            #per gene
            for(gene in geneList) {

                geneHet[gene] = geneCumulativePersite[gene]/geneLength[gene]

                print gene, geneHet[gene]

            }
        }' "\$name"_vcf_first_filter_no_extd > "\${name}"_heteroplasmy_per_sample

        awk -v sampleName="\$name" '
        /Global heteroplasmy:/ {
            out = sampleName "_per_gene_and_global.txt"
        }
        {
            print > out
        }
        BEGIN {
            out = sampleName "_per_site.txt"
        }
        ' "\${name}"_heteroplasmy_per_sample

        }
        export -f heteroplasmy
        find ./ -name "*vcf" | parallel -j $parallel heteroplasmy



        cp extractedSequences* ${params.output}/GENOTYPING
        cp *.vcf.gz* ${params.output}/GENOTYPING
        """
}
