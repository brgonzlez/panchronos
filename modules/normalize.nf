/*
 * NORMALIZE{} will apply coverage normalization on aligned data.
 */


process NORMALIZE {

        input:
        path refLength
        path rawCoverage
        val parallel
        path heteroplasmy

        output:
        path 'panchronos_per_gene_statistics.txt', emit: geneNormalizedSummary
        path 'panchronos_global_statistics.txt' , emit: globalMeanCoverage

        script:
        """
        #!/bin/bash

        mkdir -p ${params.output}/STATS

        echo -e "sampleID\tgene\tnormalizedGeneSimple\tnormalizedGeneScaled\tnormalizedGenomeSimple\tnormalizedGenomeScaled\tallelicBalance" > panchronos_per_gene_statistics.txt
        echo -e "sampleID \t sampleCoverage \t refCount \t globalMean \t allelicBalance"  > panchronos_global_statistics.txt

        normalization() {
        file=\$1

                name=\$(basename "\${file%_rawCoverage.txt}")

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
                 awk -v s="\$sample" -v v="\$value" '\$1 == s { print \$0, v }' globalMeanCoverage.txt >> panchronos_global_statistics.txt
        done < global_allelic_balance_summary.txt

        #per gene
        while read -r sample gene value; do
                awk -v s="\$sample" -v g="\$gene" -v v="\$value" '\$1 == s && \$2 == g { print \$0, v }' geneNormalizedSummary.txt >> panchronos_per_gene_statistics.txt
        done < per_gene_allelic_balance_summary.txt

        cp panchronos_global_statistics.txt ${params.output}/STATS/
        """
}
