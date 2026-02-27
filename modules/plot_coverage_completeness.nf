/*
 * PLOT_COVERAGE_COMPLETENESS{} will plot normalized gene coverage vs gene completeness (breadth of coverage).
 */

process PLOT_COVERAGE_COMPLETENESS {
        conda "${projectDir}/envs/plot.yaml"

        input:
        path geneNormalizedUpdated
        val completeness
        val coverage_lower
        val coverage_upper
        val normalised_boundary_plot

        output:
        path 'plotCoverage_vs_Completeness*', emit: plotCoverage_vs_Completeness

        script:
        """
        #!/bin/bash

        mkdir -p ${params.output}/PLOTS


        if [[ -f panchronos_synthetic_reads_per_gene_statistics_threshold.tab ]]; then
                awk 'NR>1{print \$0}' panchronos_synthetic_reads_per_gene_statistics_threshold.tab > panchronos_synthetic_reads_per_gene_statistics_threshold_no_header.tab
                cat panchronos_per_gene_statistics_threshold.tab panchronos_synthetic_reads_per_gene_statistics_threshold_no_header.tab > panchronos_mixed_per_gene_statistics_threshold.tab

                #clean-up the names
                sed -i -e 's/postPangenomeAlignment_//g' panchronos_mixed_per_gene_statistics_threshold.tab
                #we make one file per sample
                awk 'NR>1 {print \$1}' panchronos_mixed_per_gene_statistics_threshold.tab | sort | uniq > samples.txt

                while read -r sample;do
                        name=\$(echo "\${sample}")
                        awk 'NR==1{print \$0}' panchronos_mixed_per_gene_statistics_threshold.tab > "\$name"_individual_normalised.tab
                        grep -w "\$sample" panchronos_mixed_per_gene_statistics_threshold.tab >> "\$name"_individual_normalised.tab
                done < samples.txt
        else
                #clean-up the names
                sed -i -e 's/postPangenomeAlignment_//g' $geneNormalizedUpdated

                #we make one file per sample
                awk 'NR>1 {print \$1}' $geneNormalizedUpdated | sort | uniq > samples.txt


                while read -r sample;do
                        name=\$(echo "\${sample#postPangenomeAlignment_}")
                        awk 'NR==1{print \$0}' $geneNormalizedUpdated > "\$name"_individual_normalised.tab
                        grep -w "\$sample" $geneNormalizedUpdated >> "\$name"_individual_normalised.tab
                done < samples.txt
        fi

        plot_cov() {
        tab_file=\$1

        plot_cvg_vs_completeness.py "\$tab_file" $completeness $coverage_lower $coverage_upper $normalised_boundary_plot
        }
        export -f plot_cov
        find ./ -name "*_individual_normalised.tab" | parallel -j $task.cpus plot_cov


        cp *png ${params.output}/PLOTS
        """
}
