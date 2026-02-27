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
        val allelic_site
        path refLength
        path rawCoverage
        val dp_threshold
        path config

        output:
        path 'extractedSequences*.fasta', emit: consensusSequences
        path '*_per_site.txt', emit: per_site
        path '*_per_gene_and_global.txt', emit: per_gene_and_global

        script:
        """
        #!/bin/bash

        mkdir -p ${params.output}/GENOTYPING

        awk '\$4=="M"{ print > "modern_config.tab"} \$4=="A" { print > "ancient_config.tab"}' $config

        #lets get global mean first

        global_mean_depth() {
        file=\$1
        name=\$(basename "\${file%_rawCoverage.txt}")
        name="\${name#postPangenomeAlignment_}"

                refCount=\$(cat $refLength)
                globalMean=\$(awk -v count="\$refCount" '{sum += \$3} END {if (count > 0) print sum / count; else print "Something went wrong, check log file"}' "\$file")
                echo -e "\$name\t\$globalMean" > "\${name}"_globalMeanCoverage.txt
        }
        export -f global_mean_depth
        find ./ -name "*_rawCoverage.txt" | parallel -j $parallel global_mean_depth



        bcfconsensus() {
        bam_file=\$1
        file_name=\$(basename "\${bam_file%.bam}")
        file_name="\${file_name#postPangenomeAlignment_}"

                # flow contrl because mpileup seems to regularly produce corrupted files if server is too busy
                max_attempts=5
                attempt=1

                while [ "\$attempt" -le "\$max_attempts" ]; do
                        bcftools mpileup -f $panGenomeRef -q $mapq -Q $baseq "\$bam_file" > "\${file_name}"_mpileup_file

                        if awk '!/^#/ && \$4 != "N" {found=1; exit} END {exit !found}' "\${file_name}"_mpileup_file; then
                                echo -e "mpileup file for \$file_name looks fine. Moving on"
                                break
                        else
                                echo -e "mpileup has 1 or more N for \$file_name. Looks corrupted. Retrying . . ."
                                rm -rf "\${file_name}"_mpileup_file
                        fi

                                ((attempt++))
                done


                                if (( attempt > max_attempts )); then
                                echo "ERROR: mpileup failed after \${max_attempts} attempts" >&2
                                exit 1
                                fi

                bcftools call --ploidy 1 -m "\${file_name}"_mpileup_file > "\${file_name}"_raw.vcf

                if [[ $force_homozygot -eq 1 ]]; then
                        bcftools filter -i 'QUAL>$call_qual && ((DP4[2]+DP4[3]==0) || (DP4[0]+DP4[1]==0))' "\${file_name}"_raw.vcf > "\${file_name}".vcf
                else
                        bcftools filter -i 'QUAL>$call_qual' "\${file_name}"_raw.vcf > "\${file_name}".vcf
                fi


                #Calculate genome-wide heterozygosity per sample as a measure of contamination, called heteroplasmy analysis
                #remove intermediate files
                rm "\${file_name}"_raw.vcf

                # first I need to process the VCF file and only get the fields of interest. I need to include the pos field to remove extended sequences.
                awk -F";" 'BEGIN {OFS="\t"} !/^#/ {print \$0}' "\${file_name}".vcf | sed -e 's/;/\t/g' -e 's/DP4=//g' -e 's/DP=//g' |  awk '{print \$1,\$2,\$8,\$(NF-3)}' | sed -e 's/,/\t/g' -e 's/ /\t/g' > "\$file_name"_vcf_first_filter

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
                ' "\$file_name"_vcf_first_filter > "\$file_name"_vcf_first_filter_no_extd

                #split file into two: one for heteroplasmy and one for DP analysis
                awk 'BEGIN {OFS="\t"} {print \$1, \$2, \$3}' "\$file_name"_vcf_first_filter_no_extd > "\$file_name"_vcf_first_filter_by_DP_no_extd
                awk 'BEGIN {OFS="\t"} {print \$1, \$2, \$4, \$5, \$6, \$7}' "\$file_name"_vcf_first_filter_no_extd > "\$file_name"_vcf_first_filter_no_extd_TMP && mv "\$file_name"_vcf_first_filter_no_extd_TMP  ./"\$file_name"_vcf_first_filter_no_extd


                #DP analysis. DP value cannot be bigger than N times global mean, save output to lines to exclude

                gM=\$(awk '{print \$NF}' "\${file_name}"_globalMeanCoverage.txt)
                gfactor=\$(awk -v cutoff=$dp_threshold -v globalMean="\$gM" 'BEGIN {print (cutoff * globalMean)}')
                echo -e "gM value: \$gM, gfactor value \$gfactor"

                awk -v dp_threshold="\$gfactor" 'BEGIN {OFS="\t"}
                        \$3 > dp_threshold {
                        print \$1 , \$2
                }' "\$file_name"_vcf_first_filter_by_DP_no_extd >> "\${file_name}"_sites_to_exclude


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
                }' "\$file_name"_vcf_first_filter_no_extd > "\${file_name}"_heteroplasmy_per_sample

                awk -v sampleName="\$file_name" '
                /Global heteroplasmy:/ {
                out = sampleName "_per_gene_and_global.txt"
                }
                {
                print > out
                }
                BEGIN {
                out = sampleName "_per_site.txt"
                }
                ' "\${file_name}"_heteroplasmy_per_sample

                #now apply per site allelic balance threshold on vcf file.
                awk -v allelic_cutoff=$allelic_site 'BEGIN {OFS="\t"} \$3 < allelic_cutoff { print \$1, \$2 }' "\${file_name}"_per_site.txt >> "\${file_name}"_sites_to_exclude

                awk 'BEGIN{
                    OFS="\t"
                    #exclude into an array
                    while ( (getline < "'\${file_name}_sites_to_exclude'") > 0 ) {
                        sites[\$1","\$2] = 1
                    }
                }

                /^#/ { print; next }

                #clean columns
                {
                    g = \$1
                    p = \$2
                    gsub(/^[ \\t]+|[ \\t]+\$/, "", g)
                    gsub(/^[ \\t]+|[ \\t]+\$/, "", p)
                }

                #akip if site is in exclusion list
                !(g","p in sites) { print }
                ' "\${file_name}.vcf" > "\${file_name}.filtered.vcf" && mv "\${file_name}.filtered.vcf" ./"\${file_name}.vcf"


                bgzip -i -c "\${file_name}".vcf > "\${file_name}".vcf.gz
                bcftools index "\${file_name}".vcf.gz

                file_name_group="\${file_name#mergedGroup}"
                if grep -q "\${file_name_group}" ancient_config.tab; then
                        bcftools consensus -a N -f $panGenomeRef "\${file_name}".vcf.gz > extractedSequences"\${file_name}".fq
                else
                        bcftools consensus -a - -f $panGenomeRef "\${file_name}".vcf.gz > extractedSequences"\${file_name}".fq
                fi

                seqtk seq -a extractedSequences"\${file_name}".fq > extractedSequences"\${file_name}".fasta
                rm -f extractedSequences"\${file_name}".fq

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
                ' $panGenomeRef extractedSequences"\${file_name}".fasta > "\${file_name}"_padded && mv "\${file_name}"_padded ./extractedSequences"\${file_name}".fasta



                #now remove the extended sequences
                awk -v trim=$extension '
                        /^>/ {
                                print
                                next
                        }

                        !/^>/ {
                                len = length(\$0)
                                print substr(\$0, trim + 1, len - 2 * trim)
                        }' extractedSequences"\${file_name}".fasta > trimmed_"\${file_name}" && mv trimmed_"\${file_name}" ./extractedSequences"\${file_name}".fasta

        }
        export -f bcfconsensus
        find ./ -name "*.bam" | parallel -j $parallel bcfconsensus

        cp extractedSequences* ${params.output}/GENOTYPING
        cp *.vcf.gz* ${params.output}/GENOTYPING
        """
}
