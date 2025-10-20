/*
 * ALIGNMENT_SUMMARY{} will .
 */


process ALIGNMENT_SUMMARY {
	conda "${projectDir}/envs/alignment_summary.yaml"
	
	input:
	path configFile
	path bamfiles
	val parallel
	val trim

	output:
    path 'postPangenomeAlignment*bam' , emit: postAlignmentFiles
	path 'completenessSummary.tab', emit: completenessSummary
	path '*_refLength.txt', emit: refLength
	path '*_rawCoverage.txt' , emit: rawCoverage


	script:
	"""
	#!/bin/bash

	mkdir -p ${params.output}/STATS

	awk '{print \$NF}' $configFile | uniq > groups.txt

	while read -r groupID; do
		groupName=\$(echo "\$groupID")
		grep -w "\$groupID" $configFile | awk '{print \$1}' > "\$groupName"ID
	done < groups.txt

	merging() {
	file=\$1
                IDs=\$(basename "\${file%ID}")

                bamFiles=()

                while read -r sampleName; do
                        bamFiles+=("\${sampleName%.fastq*}_aligned.bam")
                        echo "Adding BAM file: \${sampleName%.fastq*}_aligned.bam"
		done < "\$file"

                if [ \${#bamFiles[@]} -eq 1 ]; then
                        cp "\${bamFiles[0]}" postPangenomeAlignment_"\${IDs}".bam
                elif [ \${#bamFiles[@]} -gt 1 ]; then
                        samtools merge postPangenomeAlignment_mergedGroup"\${IDs}".bam "\${bamFiles[@]}"
                fi
	}
	export -f merging
	find ./ -name "*ID" | parallel -j $parallel merging



	#Getting some stats AND Re group the data so genotyping is simpler
	rawStats() {
	bam_file=\$1

                samplename=\$(basename "\${bam_file%.bam}")
                samtools index "\$bam_file"
				bam trimBam "\$bam_file" trimmed_"\$bam_file" -L $trim -R $trim --clip
				samtools index trimmed_"\$bam_file"
                samtools depth -a trimmed_"\$bam_file" > "\${samplename}_rawCoverage.txt"
                samtools idxstats trimmed_"\$bam_file" | awk '{sum += \$2} END {print sum}' > "\${samplename}_refLength.txt"
                samtools coverage trimmed_"\$bam_file" | awk -v samplename="\$samplename" 'NR>1 {print samplename, \$1, \$6}' | sed -e 's/~/_/g' | sed -e 's/ /\t/g' | sort -k 1 -t \$'\t' >> "\${samplename}"_completenessSummary.tab
				mv trimmed_"\$bam_file" ./"\${samplename}_TMP.bam"
				mv trimmed_"\$bam_file".bai ./"\${samplename}_TMP.bam.bai"
				picard AddOrReplaceReadGroups I="\${samplename}_TMP.bam" O="\${samplename}.bam" RGLB="\${samplename}" RGSM="\${samplename}" RGPU=Illumina RGPL=ILLUMINA RGID="\${samplename}" RGDS="\${samplename}"
				samtools index "\${samplename}.bam"
	}
	export -f rawStats
	find ./ -name "postPangenomeAlignment*bam" | parallel -j $parallel rawStats

	rm *TMP.bam*

	cat *_completenessSummary.tab > completenessSummary.tab

	cp completenessSummary.tab ${params.output}/STATS

	cat .command.out >> alignmentSummary.log
	"""
}
