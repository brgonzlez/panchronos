/*
 * Synthetic reads sub-workflow will generate artificial FASTQ files for input FASTA used in panaroo, align them against extended pangenome, genotype and extract gene sequences to replace them in MSA.
 */


process SYNTHETIC_READS {
        conda "${projectDir}/envs/art.yaml"

        input:
        path fasta_files
        val depth
        val length
        val parallel
        path key

        output:
        path '*.fq' , emit: synthetic_reads_files

        script:
        """
        #!/bin/bash


        make_reads() {
        fasta_file=\$1

                name=\$(basename "\${fasta_file%.fasta}")
                art_illumina -i "\${fasta_file}" -l $length -f $depth -o "\${name}"

        }
        export -f make_reads
        find ./ -name "*fasta" | parallel -j $parallel make_reads
        """
}


process SYNTHETIC_READS_ALIGNMENT {
        conda "${projectDir}/envs/alignment.yaml"

        input:
        path synthetic_reads
        val parallel
        val threads
        path pangenome

        output:
        path '*_aligned.bam', emit: synthetic_bam

        script:
        """
        #!/bin/bash


        bwa index $pangenome

        align_synthetic() {
        synthetic_reads=\$1

                name=\$(basename "\${synthetic_reads%.fq}")

                bwa mem -B 1 -E 1 $pangenome "\${synthetic_reads}" -t $threads > "\${name}".sam
                samtools view -bS "\${name}".sam > "\${name}".bam
                samtools quickcheck "\${name}".bam
                samtools sort -o "\${name}"_sorted.bam -O bam -@ $threads "\${name}".bam
                picard MarkDuplicates OPTICAL_DUPLICATE_PIXEL_DISTANCE=100 REMOVE_DUPLICATES=TRUE I="\${name}"_sorted.bam  O="\${name}_deduped.bam" M="\${name}_deduped.stats"
                samtools index "\${name}_deduped.bam"
                samtools view -b -@ $threads -F 4 "\${name}_deduped.bam" > "\${name}"_sortedMappedreads.bam
                samtools index "\${name}"_sortedMappedreads.bam
                samtools sort -o "\${name}"_aligned.bam -O bam -@ $threads "\${name}"_sortedMappedreads.bam
        }
        export -f align_synthetic
        find ./ -name "*fq" | parallel -j $parallel align_synthetic
        """
}

process SYNTHETIC_READS_GENOTYPING {
        conda "${projectDir}/envs/consensus.yaml"

        input:
        path synthetic_bam
        val parallel
        val extended
        path panGenomeRef
        tuple val(mapq) , val(baseq) , val(call_qual)
        val force_homozygot
        val dp_threshold
        val allelic_site
        path refLength
        path raw_coverage

        output:
        path 'extractedSequences*', emit: synthetic_reads_sequences
        path '*_per_gene_and_global.txt', emit: synthetic_per_gene_and_global

        script:
        """
        #!/bin/bash


        global_mean_depth() {
        file=\$1
        name=\$(basename "\${file%_rawCoverage.txt}")

                refCount=\$(cat $refLength)
                globalMean=\$(awk -v count="\$refCount" '{sum += \$3} END {if (count > 0) print sum / count; else print "Something went wrong, check log file"}' "\$file")
                echo -e "\$name\t\$globalMean" > "\${name}"_globalMeanCoverage.txt
        }
        export -f global_mean_depth
        find ./ -name "*_rawCoverage.txt" | parallel -j $parallel global_mean_depth


        synthetic_consensus() {
        bam_file=\$1

                name=\$(basename "\${bam_file%.bam}")


                # flow contrl because mpileup seems to regularly produce corrupted files if server is too busy
                max_attempts=5
                attempt=1

                while [ "\$attempt" -le "\$max_attempts" ]; do
                        bcftools mpileup -f $panGenomeRef -q $mapq -Q $baseq "\$bam_file" > "\${name}"_mpileup_file

                        if awk '!/^#/ && \$4 != "N" {found=1; exit} END {exit !found}' "\${name}"_mpileup_file; then
                                echo -e "mpileup file for \$name looks fine. Moving on"
                                break
                        else
                                echo -e "mpileup has 1 or more N for \$name. Looks corrupted. Retrying . . ."
                                rm -rf "\${name}"_mpileup_file
                        fi

                                ((attempt++))
                done


                                if (( attempt > max_attempts )); then
                                        echo "ERROR: mpileup failed after \${max_attempts} attempts" >&2
                                        exit 1
                                fi


                bcftools call --ploidy 1 -m "\${name}"_mpileup_file > "\${name}"_raw.vcf


                if [[ $force_homozygot -eq 1 ]]; then
                        bcftools filter -i 'QUAL>$call_qual && ((DP4[2]+DP4[3]==0) || (DP4[0]+DP4[1]==0))' "\${name}"_raw.vcf > "\${name}".vcf
                else
                        bcftools filter -i 'QUAL>$call_qual' "\${name}"_raw.vcf > "\${name}".vcf
                fi


                # first I need to process the VCF file and only get the fields of interest. I need to include the pos field to remove extended sequences.
                awk -F";" 'BEGIN {OFS="\t"} !/^#/ {print \$0}' "\${name}".vcf | sed -e 's/;/\t/g' -e 's/DP4=//g' -e 's/DP=//g' |  awk '{print \$1,\$2,\$8,\$(NF-3)}' | sed -e 's/,/\t/g' -e 's/ /\t/g' > "\$name"_vcf_first_filter

                #note to myself. condition ? value_if_true : value_if_false

                awk -v to_trim=$extended '
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

                #split file into two: one for heteroplasmy and one for DP analysis
                awk 'BEGIN {OFS="\t"} {print \$1, \$2, \$3}' "\$name"_vcf_first_filter_no_extd > "\$name"_vcf_first_filter_by_DP_no_extd
                awk 'BEGIN {OFS="\t"} {print \$1, \$2, \$4, \$5, \$6, \$7}' "\$name"_vcf_first_filter_no_extd > "\$name"_vcf_first_filter_no_extd_TMP && mv "\$name"_vcf_first_filter_no_extd_TMP  ./"\$name"_vcf_first_filter_no_extd


                #DP analysis. DP value cannot be bigger than N times global mean, save output to lines to exclude

                gM=\$(awk '{print \$NF}' "\${name}"_globalMeanCoverage.txt)
                gfactor=\$(awk -v cutoff=$dp_threshold -v globalMean="\$gM" 'BEGIN {print (cutoff * globalMean)}')
                echo -e "gM value: \$gM, gfactor value \$gfactor"

                awk -v dp_threshold="\$gfactor" 'BEGIN {OFS="\t"}
                        \$3 > dp_threshold {
                        print \$1 , \$2
                }' "\$name"_vcf_first_filter_by_DP_no_extd >> "\${name}"_sites_to_exclude


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

                #now apply per site allelic balance threshold on vcf file.
                awk -v allelic_cutoff=$allelic_site 'BEGIN {OFS="\t"} \$3 < allelic_cutoff { print \$1, \$2 }' "\${name}"_per_site.txt >> "\${name}"_sites_to_exclude

                awk 'BEGIN{
                    OFS="\t"
                    #exclude into an array
                    while ( (getline < "'\${name}_sites_to_exclude'") > 0 ) {
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
                ' "\${name}.vcf" > "\${name}.filtered.vcf" && mv "\${name}.filtered.vcf" ./"\${name}.vcf"

                bgzip -i -c "\${name}".vcf > "\${name}".vcf.gz
                bcftools index "\${name}".vcf.gz
                bcftools consensus -a - -f $panGenomeRef "\${name}".vcf.gz > extractedSequences"\${name}".fq


                seqtk seq -a extractedSequences"\${name}".fq > extractedSequences"\${name}".fasta
                rm -f extractedSequences"\${name}".fq

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
                ' $panGenomeRef extractedSequences"\${name}".fasta > "\${name}"_padded && mv "\${name}"_padded ./extractedSequences"\${name}".fasta



                #now remove the extended sequences
                awk -v trim=$extended '
                        /^>/ {
                                print
                                next
                        }

                        !/^>/ {
                                len = length(\$0)
                                print substr(\$0, trim + 1, len - 2 * trim)
                        }' extractedSequences"\${name}".fasta > trimmed_"\${name}" && mv trimmed_"\${name}" ./extractedSequences"\${name}".fasta


        }
        export -f synthetic_consensus
        find ./ -name "*.bam" | parallel -j $parallel synthetic_consensus
        """
}

/*
 * SYNTHETIC_READS_ALIGNMENT_SUMMARY{} will generate alignment summary files for synthetic reads.
 */


process SYNTHETIC_READS_ALIGNMENT_SUMMARY {
        conda "${projectDir}/envs/alignment.yaml"
        label 'demand_3'

        input:
        path bamfiles
        val parallel
        val extend

        output:
        path 'synthetic_completeness_summary.tab', emit: synthetic_completeness_summary
        path '*_rawCoverage.txt' , emit: synthetic_raw_coverage


        script:
        """
        #!/bin/bash

        mkdir -p ${params.output}/STATS

        #Getting some stats AND Re group the data so genotyping is simpler
        rawStats() {
        bam_file=\$1

                samplename=\$(basename "\${bam_file%.bam}")
                samtools index "\$bam_file"
                samtools depth -a "\$bam_file" > "\${samplename}_PreRawCoverage.txt"

                                #Remove the extended positions from raw coverage
                                awk -v ext=$extend '
                                {
                                        last[\$1] = \$2   # remember last position seen for each gene
                                        data[NR] = \$0   # store all lines
                                        id[NR] = \$1
                                        pos[NR] = \$2
                                }
                                END {
                                        for (g in last) {
                                        start[g] = ext + 1
                                        stop[g]  = last[g] - ext
                                        }

                                        for (i = 1; i <= NR; i++) {
                                        g = id[i]
                                        p = pos[i]
                                        if (p >= start[g] && p <= stop[g])
                                                print data[i]
                                        }
                                }'  "\${samplename}_PreRawCoverage.txt" > "\${samplename}_rawCoverage.txt"

                                samtools coverage "\$bam_file" | awk -v samplename="\$samplename" 'NR>1 {print samplename, \$1, \$6}' | sed -e 's/~/_/g' | sed -e 's/ /\t/g' | sort -k 1 -t \$'\t' >> "\${samplename}"_completenessSummary.tab
                                mv "\$bam_file" ./"\${samplename}_TMP.bam"
                                mv "\$bam_file".bai ./"\${samplename}_TMP.bam.bai"
                                picard AddOrReplaceReadGroups I="\${samplename}_TMP.bam" O="\${samplename}.bam" RGLB="\${samplename}" RGSM="\${samplename}" RGPU=Illumina RGPL=ILLUMINA RGID="\${samplename}" RGDS="\${samplename}"
                                samtools index "\${samplename}.bam"
        }
        export -f rawStats
        find ./ -name "*bam" | parallel -j $parallel rawStats

        rm *TMP.bam*
        cat *_completenessSummary.tab > synthetic_completeness_summary.tab
        cat .command.out >> synthetic_alignment_summary.log
        """
}

/*
 * SYNTHETIC_READS_NORMALIZATION{} will apply coverage normalization on aligned data.
 */


process SYNTHETIC_READS_NORMALIZATION {

        input:
        path refLength
        path rawCoverage
        val parallel
        path heteroplasmy

        output:
        path 'panchronos_synthetic_reads_per_gene_statistics.tab', emit: synthetic_reads_gene_normalized_summary
        path 'panchronos_synthetic_reads_global_statistics.tab' , emit: synthetic_reads_global_mean_coverage

        script:
        """
        #!/bin/bash

        mkdir -p ${params.output}/STATS

        echo -e "sampleID\tgene\tnormalizedGeneSimple\tnormalizedGeneScaled\tnormalizedGenomeSimple\tnormalizedGenomeScaled\tallelicBalance" > panchronos_synthetic_reads_per_gene_statistics.tab
        echo -e "sampleID \t sampleCoverage \t refCount \t globalMean \t allelicBalance"  > panchronos_synthetic_reads_global_statistics.tab

        normalization() {
        file=\$1

                name=\$(basename "\${file%_rawCoverage.txt}")
                name="\${name#postPangenomeAlignment_}"

                #Get allelic balance results per sample
                grep "Global heteroplasmy" "\$name"_per_gene_and_global.txt | awk -v sample="\$name" 'BEGIN {OFS="\t"} {print sample, \$NF}' > "\$name"_global_allelic_balance.txt
                grep -v "Global heteroplasmy" "\$name"_per_gene_and_global.txt | awk -v sample="\$name" 'BEGIN {OFS="\t"} {print sample, \$0}' > "\$name"_per_gene_allelic_balance.txt


                #Compute global mean coverage
                refCount=\$(cat $refLength)
                globalMean=\$(awk -v count="\$refCount" '{sum += \$3} END {if (count > 0) print sum / count; else print "Something went wrong, check log file"}' "\$file")
                finalCount=\$(awk '\$3 > 0{fcount++} END {print fcount}' "\$file")
                echo -e "\$name\t\$finalCount\t\$refCount\t\$globalMean" > "\${name}"_globalMeanCoverage.txt

                #Normalize coverage per gene
                awk -v globalMean="\$globalMean" -v name="\$name" -v sampleCoverage="\$finalCount" -v refCount="\$refCount" '
                        {
                        if (\$2 > geneLength[\$1]) {
                                geneLength[\$1] = \$2
                        }
                        sumgene[\$1] += \$3
                        countgene[\$1]++
                        }
                END {
                        for (gene in geneLength) {
                                geneMean = sumgene[gene] / countgene[gene]
                                normalizedGeneScaled = (geneMean / globalMean) * geneLength[gene]
                                normalizedGeneSimple = (geneMean / globalMean)
                                normalizedGenomeSimple = (geneMean / globalMean) * (geneLength[gene] / sampleCoverage)
                                normalizedGenomeScaled = (geneMean / globalMean) * (geneLength[gene] / refCount)
                                print name"\t"gene"\t"normalizedGeneSimple"\t"normalizedGeneScaled"\t"normalizedGenomeSimple"\t"normalizedGenomeScaled
                        }
                }
                ' "\$file" > "\${name}"_geneNormalizedSummary.txt
        }
        export -f normalization
        find ./ -name "*_rawCoverage.txt" | parallel -j $parallel normalization

        cat *_geneNormalizedSummary.txt > geneNormalizedSummary.txt
        cat *_globalMeanCoverage.txt > globalMeanCoverage.txt

        #update the files with allelic balance results

        cat *_global_allelic_balance.txt > global_allelic_balance_summary.txt
        cat *_per_gene_allelic_balance.txt > per_gene_allelic_balance_summary.txt

        #global
        while read -r sample value;do
                 awk -v s="\$sample" -v v="\$value" '\$1 == s { print \$0, v }' globalMeanCoverage.txt >> panchronos_synthetic_reads_global_statistics.tab
        done < global_allelic_balance_summary.txt


        #per gene
        awk '

        function trim_whitespaces(w) {
                 gsub(/^[[:space:]]+|[[:space:]]+\$/, "", w)
                 return w
        }

        BEGIN {OFS="\t"}

        FNR==NR {

                first_pair = trim_whitespaces(\$1) OFS trim_whitespaces(\$2)
                value_first_pair[first_pair] = trim_whitespaces(\$3)
                next

         }

         {

                second_pair = trim_whitespaces(\$1) OFS trim_whitespaces(\$2)
                whole_line = \$0
                if(second_pair in value_first_pair) {

                        print whole_line, value_first_pair[second_pair]

                }

         }' per_gene_allelic_balance_summary.txt geneNormalizedSummary.txt >> panchronos_synthetic_reads_per_gene_statistics.tab

        cp panchronos_synthetic_reads_global_statistics.tab ${params.output}/STATS/
        """
}

/*
 * SYNTHETIC_READS_UPDATE_NORMALIZATION{} will append completeness/breadth of coverage per gene into the summary tab file.
 */

process SYNTHETIC_READS_UPDATE_NORMALIZATION {

        input:
        path normalized
        path completeness
        val parallel

        output:
        path 'panchronos_synthetic_reads_per_gene_statistics_threshold.tab', emit: synthetic_gene_normalized_updated

        script:
        """
        #!/bin/bash

        mkdir -p ${params.output}/STATS

        echo -e "sampleID\tgene\tnormalizedGeneSimple\tnormalizedGeneScaled\tnormalizedGenomeSimple\tnormalizedGenomeScaled\tallelicBalance\tgeneCompleteness" > panchronos_synthetic_reads_per_gene_statistics_threshold.tab

        sed -i -e 's/~/_/g' $normalized

        sed -i -e 's/postPangenomeAlignment_//g' $completeness

        awk '{print \$1}' $completeness | uniq > samples.txt

        while read -r sample; do
                echo "\$sample" > "\$sample".map
        done < samples.txt

        normalize_updating() {
        file=\$1
        name=\$(basename "\${file%.map}")

                grep -w "\$name" $normalized > "\$name"_geneNormalizedSummary
                grep -w "\$name" $completeness > "\$name"_completenessSummary

                awk '{print \$1"XYZ"\$2, \$3, \$4, \$5, \$6, \$7}' "\$name"_geneNormalizedSummary > "\$name"_TMP1

                awk '{print \$1"XYZ"\$2, \$3}' "\$name"_completenessSummary > "\$name"_TMP2

                while read -r ID completeness;do

                        if grep -wq "\${ID}" "\$name"_TMP1; then
                                oldLine=\$(grep -w "\${ID}" "\$name"_TMP1)
                                specificCompleteness=\$(grep -w "\${ID}" "\$name"_TMP2 | awk '{print \$NF}')
                                echo -e "\${oldLine}\t\${specificCompleteness}" >> "\$name"_geneNormalizedUpdated.tab
                        fi

                done < "\$name"_TMP2

                sed -i -e 's/XYZ/\t/g' "\$name"_geneNormalizedUpdated.tab
                sed -i -e 's/ /\t/g' "\$name"_geneNormalizedUpdated.tab
        }
        export -f normalize_updating
        find ./ -name "*.map" | parallel -j $parallel normalize_updating

        cat *_geneNormalizedUpdated.tab >> panchronos_synthetic_reads_per_gene_statistics_threshold.tab

        rm -f *TMP1 *TMP2

        cp panchronos_synthetic_reads_per_gene_statistics_threshold.tab ${params.output}/STATS/panchronos_synthetic_reads_per_gene_statistics_threshold.tab
        """
}
