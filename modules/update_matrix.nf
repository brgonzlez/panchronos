/*
 * UPDATE_MATRIX{} will update presence/absence matrix with user data.
 */


process UPDATE_MATRIX {
        conda "${projectDir}/envs/plot.yaml"

        input:
        path pangenomeRtab
        path globalMeanCoverage
        path normalized
        val gCompleteness
        val lowerBound
        val upperBound
        path gene_list

        output:
        path 'matrix.tab', emit: matrix
        path '*_final.csv', emit: finalCsv
        path 'sample_names', emit: sampleNames
        path 'INDEX', emit: index

        script:
        """
        #!/bin/bash

        awk 'NR==1{print \$0}' $pangenomeRtab > matrix.tab

        #now get only the genes from gene_list
        while read -r gene;do
                awk -v gene="\$gene" 'BEGIN {FS=OFS="\t"} NR > 1 && \$1 == gene' $pangenomeRtab >> filtered.Rtab
        done < $gene_list

        awk 'NR>1 {print \$0}' filtered.Rtab | sort -k 1 -t \$'\t' >> matrix.tab
        awk 'NR>1 {print \$1}' filtered.Rtab | sort -k 1 -t \$'\t' > INDEX


        if [[ -f panchronos_synthetic_reads_global_statistics.tab ]]; then
                awk 'NR>1{print \$0}' panchronos_synthetic_reads_global_statistics.tab > panchronos_synthetic_reads_global_statistics_no_header.tab
                cat panchronos_global_statistics.txt panchronos_synthetic_reads_global_statistics_no_header.tab > panchronos_mixed_global_statistics.tab
                awk 'NR>1 {print \$1}' panchronos_mixed_global_statistics.tab > sample_names

                while read -r name; do
                        echo -e "Gene\tnormalizedCoverage\tcompleteness" > "\${name}"_index.tmp
                        awk -v name="\$name" 'BEGIN {FS=OFS="\t"} \$1 == name {print \$2, \$3, \$NF}' $normalized >> "\${name}_index.tmp"
                done < sample_names

        else
                awk 'NR>1 {print \$1}' $globalMeanCoverage > sample_names

                while read -r name; do

                        echo -e "Gene\tnormalizedCoverage\tcompleteness" > "\${name}"_index.tmp
                        awk -v name="\$name" 'BEGIN {FS=OFS="\t"} \$1 == name {print \$2, \$3, \$NF}' $normalized >> "\${name}_index.tmp"

                done < sample_names
        fi

        for i in *_index.tmp; do
                sed -i -e 's/ /\t/g' "\$i"
                lambda.py "\$i" $gCompleteness $lowerBound $upperBound
        done
        """
}
