/*
 * GET_OUTGROUP{} will download 1 genome for specified tax id
 */

process GET_OUTGROUP {
	conda "${projectDir}/envs/get_data.yaml"

	input:
	val outgroup_ID
	
	output:
	path '*fna', emit: outgroupFasta

	script:
	"""
	#!/bin/bash
	counter=0
	esearch -db assembly -query "txid${outgroup_ID}[Organism] AND (latest[filter] AND (complete genome[filter] OR chromosome level[filter]))" | esummary | xtract -pattern DocumentSummary -element FtpPath_RefSeq | while read url; do
        
		if [ "\$counter" -ge 1 ]; then
			break
		fi
	
		if [ -z "\$url" ]; then
			continue
		fi
        
		fname="\$(basename "\$url")"
		downloadMeAndCheck() {
			wget "\$url/\${fname}_genomic.gbff.gz"
	        	wget "\$url/\${fname}_genomic.fna.gz"

			gunzip -f "\${fname}_genomic.gbff.gz"
			gunzip -f "\${fname}_genomic.fna.gz"


			if [ -f "\${fname}_genomic.gbff.gz" ] || [ -f "\${fname}_genomic.fna.gz" ]; then	
				rm -f "\${fname}_genomic.gbff.gz" "\${fname}_genomic.fna.gz"
				return 1
			fi

			return 0
		}

		while ! downloadMeAndCheck; do
			echo "Files were corrupted. Retrying"
			sleep 3
		done
		# Adding to the counter
	        counter="\$((counter + 1))"     
                fna_count=\$(ls -1 *.fna 2>/dev/null | wc -l)
                gbff_count=\$(ls -1 *.gbff 2>/dev/null | wc -l)

                # check that there is already one file
                if [ "\$fna_count" -lt 1 ] || [ "\$gbff_count" -lt 1 ]; then
                        echo "It seems a file was already downloaded. Moving on"
                        break
                else
                        continue
                fi

	done
	"""
}
