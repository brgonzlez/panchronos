/*
 * COVERAGE_BOUNDS{} will apply coverage bounds thresholds.
 */


process COVERAGE_BOUNDS {

        input:
        path geneNormalizedUpdated
        val normalizedCoverageDown
        val normalizedCoverageUp
        val completenessBound
        path final_list_genes

        output:
        path 'panchronos_normalisation_summary_filtered.tab', emit: geneNormalizedUpdatedFiltered

        script:
        """
        #!/bin/bash

        # if synthetic reads file (panchronos_synthetic_reads_per_gene_statistics_threshold.tab), then

        if [[ -f panchronos_synthetic_reads_per_gene_statistics_threshold.tab ]]; then
                awk 'NR>1{print \$0}' panchronos_synthetic_reads_per_gene_statistics_threshold.tab > panchronos_synthetic_reads_per_gene_statistics_threshold_no_header.tab
                cat panchronos_per_gene_statistics.tab panchronos_synthetic_reads_per_gene_statistics_threshold_no_header.tab > panchronos_mixed_per_gene_statistics_threshold.tab
                awk 'NR==1{print \$0}' panchronos_mixed_per_gene_statistics_threshold.tab > header
                awk -v UpBound=$normalizedCoverageUp '\$3 < UpBound {print \$0}' panchronos_mixed_per_gene_statistics_threshold.tab > TMP1
                awk -v DownBound=$normalizedCoverageDown '\$3 > DownBound {print \$0}' TMP1 > TMP2
                awk -v completenessBound=$completenessBound '\$NF > completenessBound {print \$0}' TMP2 > TMP3
                cat header TMP3 > panchronos_normalisation_summary_filtered.tab

                #make raw completeness matrix
                awk '
                BEGIN {FS=OFS="\t"}
                        NR>1 {
                                if (!(\$1 in sample)) {
                                        sample[\$1] = 1
                                        sample_order[++item_sample]=\$1   #items order index for printing
                                }
                                gene[\$2] = 1
                                samples_gene_completeness[\$1,\$2]=\$NF
                        }
                END {
                        header="Gene"
                        for(i=1;i<=item_sample;i++) {
                                header=header OFS sample_order[i]
                        }
                        print header
                        for(g in gene) {
                                gene_row=""
                                for(i=1;i<=item_sample;i++) {
                                        s=sample_order[i]
                                        values_per_record=((s,g) in samples_gene_completeness ? samples_gene_completeness[s,g] : 0)
                                        gene_row=gene_row OFS values_per_record
                                }
                                print g gene_row
                        }
                }' panchronos_mixed_per_gene_statistics_threshold.tab > raw_gene_completeness_matrix.tab

        else
                awk 'NR==1{print \$0}' panchronos_per_gene_statistics.tab > header
                awk -v UpBound=$normalizedCoverageUp '\$3 < UpBound {print \$0}' panchronos_per_gene_statistics.tab > TMP1
                awk -v DownBound=$normalizedCoverageDown '\$3 > DownBound {print \$0}' TMP1 > TMP2
                awk -v completenessBound=$completenessBound '\$NF > completenessBound {print \$0}' TMP2 > TMP3
                cat header TMP3 > panchronos_normalisation_summary_filtered.tab

                #make raw completeness matrix
                awk '
                BEGIN {FS=OFS="\t"}
                        NR>1 {
                                if (!(\$1 in sample)) {
                                        sample[\$1] = 1
                                        sample_order[++item_sample]=\$1   #items order index for printing
                                }
                                gene[\$2] = 1
                                samples_gene_completeness[\$1,\$2]=\$NF
                        }
                END {
                        header="Gene"
                        for(i=1;i<=item_sample;i++) {
                                header=header OFS sample_order[i]
                        }
                        print header
                        for(g in gene) {
                                gene_row=""
                                for(i=1;i<=item_sample;i++) {
                                        s=sample_order[i]
                                        values_per_record=((s,g) in samples_gene_completeness ? samples_gene_completeness[s,g] : 0)
                                        gene_row=gene_row OFS values_per_record
                                }
                                print g gene_row
                        }
                }' panchronos_per_gene_statistics.tab > raw_gene_completeness_matrix.tab
        fi

        rm TMP1 TMP2 TMP3 header

        #Check if raw_gene_completeness_matrix.tab contains every gene listed in final_list_genes.txt
        awk 'BEGIN{OFS="\t"}
        FNR==NR && NR == 1{
                print \$0
                number_of_fields=(NF-1)
                next
        }
        FNR==NR && NR > 1 {
                matrix_rows[\$1] = \$0
                next
        }
        {
                list_of_genes[\$1]=1
                next
        }
        END {
                zeroes_list=0
                
                while(z_count <= (number_of_fields-2)) {
                        zeroes_list = zeroes_list OFS "0"
                        z_count++
                        }

                for(gene in list_of_genes) {
                        if(gene in matrix_rows) {
                                print matrix_rows[gene]
                        } else {
                                print gene, zeroes_list
                        }
                }
        }' raw_gene_completeness_matrix.tab final_list_genes.txt  > tmp_raw_gene_completeness_matrix.tab && mv tmp_raw_gene_completeness_matrix.tab raw_gene_completeness_matrix.tab

        cp raw_gene_completeness_matrix.tab ${params.output}/STATS/panchronos_raw_gene_completeness_matrix.tab
        cp panchronos_normalisation_summary_filtered.tab ${params.output}/STATS/panchronos_per_gene_statistics_after_thresholds.tab
        """
}
