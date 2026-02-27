/*
 * COVERAGE_BOUNDS{} will apply coverage bounds thresholds.
 */


process COVERAGE_BOUNDS {

        input:
        path geneNormalizedUpdated
        val normalizedCoverageDown
        val normalizedCoverageUp
        val completenessBound

        output:
        path 'panchronos_normalisation_summary_filtered.tab', emit: geneNormalizedUpdatedFiltered

        script:
        """
        #!/bin/bash

        # if synthetic reads file (panchronos_synthetic_reads_per_gene_statistics_threshold.tab), then

        if [[ -f panchronos_synthetic_reads_per_gene_statistics_threshold.tab ]]; then
                awk 'NR>1{print \$0}' panchronos_synthetic_reads_per_gene_statistics_threshold.tab > panchronos_synthetic_reads_per_gene_statistics_threshold_no_header.tab
                cat panchronos_per_gene_statistics_threshold.tab panchronos_synthetic_reads_per_gene_statistics_threshold_no_header.tab > panchronos_mixed_per_gene_statistics_threshold.tab
                awk 'NR==1{print \$0}' panchronos_mixed_per_gene_statistics_threshold.tab > header
                awk -v UpBound=$normalizedCoverageUp '\$3 < UpBound {print \$0}' panchronos_mixed_per_gene_statistics_threshold.tab > TMP1
                awk -v DownBound=$normalizedCoverageDown '\$3 > DownBound {print \$0}' TMP1 > TMP2
                awk -v completenessBound=$completenessBound '\$NF > completenessBound {print \$0}' TMP2 > TMP3
                cat header TMP3 > panchronos_normalisation_summary_filtered.tab
        else
                awk 'NR==1{print \$0}' $geneNormalizedUpdated > header
                awk -v UpBound=$normalizedCoverageUp '\$3 < UpBound {print \$0}' $geneNormalizedUpdated > TMP1
                awk -v DownBound=$normalizedCoverageDown '\$3 > DownBound {print \$0}' TMP1 > TMP2
                awk -v completenessBound=$completenessBound '\$NF > completenessBound {print \$0}' TMP2 > TMP3
                cat header TMP3 > panchronos_normalisation_summary_filtered.tab
        fi

        cp panchronos_normalisation_summary_filtered.tab ${params.output}/STATS/panchronos_normalisation_summary_filtered.tab
        """
}
