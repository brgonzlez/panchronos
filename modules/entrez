
/*
 * entrez() will download FASTA and GenBank files from NCBI as long as a taxonomical ID was provided (mandatory value)
 */

process entrez {
	conda "${projectDir}/envs/entrez.yaml"
        publishDir "${makeDir}/results/NCBI", mode: 'copy', overwrite: true

	input:
	val gE
	val txid
	path makeDir

	output:

	path '*fna', emit: fastaFiles
	path '*gbff', emit: gffFiles

	script:
	"""
	#!/bin/bash
	counter=0
	esearch -db assembly -query "txid${txid}[Organism] AND (latest[filter] AND (complete genome[filter] OR chromosome level[filter]))" | esummary | xtract -pattern DocumentSummary -element FtpPath_RefSeq | while read url; do
        
	if [ "\$counter" -ge "${gE}" ]; then
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
	
	# IF AFTER DOING gunzip -f WE STILL FIND THE EXTENSION *gz, then remove both files with the same {fname} (even if is just "\${fname}_genomic.gbff.gz" or "\${fname}_genomic.fna.gz" or both)
	# If we remove any file, then:
	#	1. We don't add +1 to the counter (it doesn't make sense to add +1 if the files were corrupted)
	#	2. We download both *gbff.gz and *fna.gz again and we proceed to test the integrity again by checking if there is any *gz extension after gunzip.

	counter="\$((counter + 1))"	
	fna_count=\$(ls -1 *.fna 2>/dev/null | wc -l)
	gbff_count=\$(ls -1 *.gbff 2>/dev/null | wc -l)

	# Continue downloading until the count matches gE
	if [ "\$fna_count" -lt "\$gE" ] || [ "\$gbff_count" -lt "\$gE" ]; then
		echo "Still need more files. Current count: \$fna_count .fna files, \$gbff_count .gbff files."
		continue
	fi


	done
	cat .command.out >> entrez.log
	"""
}
