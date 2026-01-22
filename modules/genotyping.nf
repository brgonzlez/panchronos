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

        output:
        path 'extractedSequences*.fasta', emit: consensusSequences
        path '*_per_site.txt', emit: per_site
        path '*_per_gene_and_global.txt', emit: per_gene_and_global

        script:
        """
        #!/bin/bash

        mkdir -p ${params.output}/GENOTYPING

        #lets get global mean first


        global_mean_depth() {
        file=\$1
        name=\$(basename "\${file%_rawCoverage.txt}")

                refCount=\$(cat $refLength)
                globalMean=\$(awk -v count="\$refCount" '{sum += \$3} END {if (count > 0) print sum / count; else print "Something went wrong, check log file"}' "\$file")
                echo -e "\$name\t\$globalMean" > "\${name}"_globalMeanCoverage.txt
        }
        export -f global_mean_depth
        find ./ -name "*_rawCoverage.txt" | parallel -j $parallel global_mean_depth



        bcfconsensus() {
        bam_file=\$1
        basename=\$(basename "\${bam_file%.bam}")

                # flow contrl because mpileup seems to regularly produce corrupted files if server is too busy
                max_attempts=5
                attempt=1

                while [ "\$attempt" -le "\$max_attempts" ]; do
                        bcftools mpileup -f $panGenomeRef -q $mapq -Q $baseq "\$bam_file" > "\${basename}"_mpileup_file

                        if awk '!/^#/ && \$4 != "N" {found=1; exit} END {exit !found}' "\${basename}"_mpileup_file; then
                                echo -e "mpileup file for \$basename looks fine. Moving on"
                                break
                        else
                                echo -e "mpileup has 1 or more N for \$basename. Looks corrupted. Retrying . . ."
                                rm -rf "\${basename}"_mpileup_file
                        fi

                                ((attempt++))
                done


                                if (( attempt > max_attempts )); then
                                echo "ERROR: mpileup failed after \${max_attempts} attempts" >&2
                                exit 1
                                fi

                bcftools call --ploidy 1 -m "\${basename}"_mpileup_file > "\${basename}"_raw.vcf

                if [[ $force_homozygot -eq 1 ]]; then
                        bcftools filter -i 'QUAL>$call_qual && ((DP4[2]+DP4[3]==0) || (DP4[0]+DP4[1]==0))' "\${basename}"_raw.vcf > "\${basename}".vcf
                else
                        bcftools filter -i 'QUAL>$call_qual' "\${basename}"_raw.vcf > "\${basename}".vcf
                fi


                #Calculate genome-wide heterozygosity per sample as a measure of contamination, called heteroplasmy analysis
                #remove intermediate files
                rm *_raw.vcf

                cp "\${basename}".vcf ./"\${basename}"COPY.vcf
                # first I need to process the VCF file and only get the fields of interest. I need to include the pos field to remove extended sequences.

                awk -F";" 'BEGIN {OFS="\t"} !/^#/ {print \$0}' "\${basename}".vcf | sed -e 's/;/\t/g' -e 's/DP4=//g' -e 's/DP=//g' |  awk '{print \$1,\$2,\$8,\$(NF-3)}' | sed -e 's/,/\t/g' -e 's/ /\t/g' > "\$basename"_vcf_first_filter

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
                ' "\$basename"_vcf_first_filter > "\$basename"_vcf_first_filter_no_extd

                #split file into two: one for heteroplasmy and one for DP analysis
                awk 'BEGIN {OFS="\t"} {print \$1, \$2, \$3}' "\$basename"_vcf_first_filter_no_extd > "\$basename"_vcf_first_filter_by_DP_no_extd
                awk 'BEGIN {OFS="\t"} {print \$1, \$2, \$4, \$5, \$6, \$7}' "\$basename"_vcf_first_filter_no_extd > "\$basename"_vcf_first_filter_no_extd_TMP && mv "\$basename"_vcf_first_filter_no_extd_TMP  ./"\$basename"_vcf_first_filter_no_extd


                #DP analysis. DP value cannot be bigger than N times global mean, save output to lines to exclude

                gM=\$(awk '{print \$NF}' "\${basename}"_globalMeanCoverage.txt)
                gfactor=\$(awk -v cutoff=$dp_threshold -v globalMean="\$gM" 'BEGIN {print (cutoff * globalMean)}')
                echo -e "gM value: \$gM, gfactor value \$gfactor"

                awk -v dp_threshold="\$gfactor" 'BEGIN {OFS="\t"}
                        \$3 > dp_threshold {
                        print \$1 , \$2
                }' "\$basename"_vcf_first_filter_by_DP_no_extd >> "\${basename}"_sites_to_exclude

                #debug
                cp "\${basename}"_sites_to_exclude ./"\${basename}"_sites_to_exclude_TMP


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
                }' "\$basename"_vcf_first_filter_no_extd > "\${basename}"_heteroplasmy_per_sample

                awk -v sampleName="\$basename" '
                /Global heteroplasmy:/ {
                out = sampleName "_per_gene_and_global.txt"
                }
                {
                print > out
                }
                BEGIN {
                out = sampleName "_per_site.txt"
                }
                ' "\${basename}"_heteroplasmy_per_sample


                #now apply per site allelic balance threshold on vcf file.
                awk -v allelic_cutoff=$allelic_site 'BEGIN {OFS="\t"} \$3 < allelic_cutoff { print \$1, \$2 }' "\${basename}"_per_site.txt >> "\${basename}"_sites_to_exclude

                #Now we need to apply a threshdold for DP values based on global mean. DP shouldnt be higher than N times global mean.
                globalMean=\$(cat "\${basename}"_globalMeanCoverage.txt  | awk '{print \$2}')
                DP_threshold=\$(awk -v gM="\$globalMean" -v dp=$dp_threshold 'BEGIN{ print gM * dp}')


                awk 'BEGIN{
                    OFS="\t"
                    # Read sites to exclude into an array
                    while ( (getline < "'\${basename}_sites_to_exclude'") > 0 ) {
                        sites[\$1","\$2] = 1
                    }
                }
                # Print headers as is
                /^#/ { print; next }

                # Trim whitespace from first two columns
                {
                    g = \$1
                    p = \$2
                    gsub(/^[ \\t]+|[ \\t]+\$/, "", g)
                    gsub(/^[ \\t]+|[ \\t]+\$/, "", p)
                }

                # Skip if site is in exclusion list
                !(g","p in sites) { print }
                ' "\${basename}.vcf" > "\${basename}.filtered.vcf" && mv "\${basename}.filtered.vcf" ./"\${basename}.vcf"




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

        cp extractedSequences* ${params.output}/GENOTYPING
        cp *.vcf.gz* ${params.output}/GENOTYPING
        """
}
