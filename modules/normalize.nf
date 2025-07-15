/*
 * NORMALIZE{} will apply coverage normalization on aligned data.
 */


process NORMALIZE {

	input:
	path refLength
	path rawCoverage
	val parallel

	output:
	path 'geneNormalizedSummary.txt', emit: geneNormalizedSummary
	path 'globalMeanCoverage.txt' , emit: globalMeanCoverage

	script:
	"""
	#!/bin/bash

	echo -e "sampleID\tgene\tnormalizedGeneSimple\tnormalizedGeneScaled\tnormalizedGenomeSimple\tnormalizedGenomeScaled" > geneNormalizedSummary.txt
	echo -e "sampleID \t sampleCoverage \t refCount \t globalMean"  > globalMeanCoverage.txt

	normalization() {
	file=\$1

		name=\$(basename "\${file%_rawCoverage.txt}")

		#Compute global mean coverage
		globalMean=\$(awk -v name="\$name" '{sum += \$3; count++} END {if (count > 0) print sum / count; else print "Something went wrong, check log file"}' "\$file")
		finalCount=\$(awk -v name="\$name" '{fcount++} END {print fcount}' "\$file")
		refCount=\$(cat "\${name}_refLength.txt")
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

	cat *_geneNormalizedSummary.txt >> geneNormalizedSummary.txt
	cat *_globalMeanCoverage.txt >> globalMeanCoverage.txt
	"""
}
