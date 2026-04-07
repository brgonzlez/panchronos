/*
 * UPDATE_NORMALIZATION{} will append completeness/breadth of coverage per gene into the summary tab file.
 */

process UPDATE_NORMALIZATION {

        input:
        path normalized
        path completeness
        val parallel

        output:
        path 'panchronos_per_gene_statistics_threshold.tab', emit: geneNormalizedUpdated
        path 'pass_the_key', emit: key_to_synth

        script:
        """
        #!/bin/bash

        mkdir -p ${params.output}/STATS

        echo -e "sampleID\tgene\tnormalizedGeneSimple\tnormalizedGeneScaled\tnormalizedGenomeSimple\tnormalizedGenomeScaled\tallelicDominance\tgeneCompleteness" > panchronos_per_gene_statistics_threshold.tab
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

        cat *_geneNormalizedUpdated.tab >> panchronos_per_gene_statistics_threshold.tab

        rm -f *TMP1 *TMP2

        cp panchronos_per_gene_statistics_threshold.tab ${params.output}/STATS/panchronos_per_gene_statistics_threshold.tab

        touch pass_the_key
        """
}
