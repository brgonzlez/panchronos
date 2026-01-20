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
        awk 'NR==1{print \$0}' $geneNormalizedUpdated > header
        awk -v UpBound=$normalizedCoverageUp '\$3 < UpBound {print \$0}' $geneNormalizedUpdated > TMP1
        awk -v DownBound=$normalizedCoverageDown '\$3 > DownBound {print \$0}' TMP1 > TMP2
        awk -v completenessBound=$completenessBound '\$NF> completenessBound {print \$0}' TMP2 > TMP3
        cat header TMP3 > panchronos_normalisation_summary_filtered.tab

        cp panchronos_normalisation_summary_filtered.tab ${params.output}/STATS/panchronos_normalisation_summary_filtered.tab
        """
}
