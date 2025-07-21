
/*
 * GET_DATA() will download FASTA and GenBank files from NCBI as long as a taxonomical ID is provided.
 */

process GET_DATA {
	conda "${projectDir}/envs/get_data.yaml"

	input:
	val genomes
	val tax_id
	val parallel

	output:
	path '*fasta', emit: fasta_files
	path '*gb', emit: gbk_files

	script:
	"""
	#!/bin/bash
        #get downloading links
        esearch -db assembly -query "txid${tax_id}[Organism] AND (latest[filter] AND (complete genome[filter] OR chromosome level[filter]))" | esummary \ 
	| xtract -pattern DocumentSummary -element FtpPath_RefSeq FtpPath_GenBank | head -n $genomes >> links.txt

        #split those links to use parallel

        mkdir -p links

        while read -r link; do

                touch links/"\\$link".map

        done < links.txt


        download_and_check() {
        file=\$1
        link=\$(cat "\$file")
        filename=\$(basename "\$link")

        max_attempts=5
        attempt=1

                while [ "\$attempt" -le "\$max_attempts" ]; do

                        wget -q "\$link/\${filename}_genomic.gbff.gz"
                        wget -q "\$link/\${filename}_genomic.fna.gz"

                        gunzip -f "\${filename}_genomic.gbff.gz" 2>/dev/null
                        gunzip -f "\${filename}_genomic.fna.gz" 2>/dev/null

                        #check if .gz files are still present after decompressing them
                        if [ -f "\${filename}_genomic.gbff.gz" ] || [ -f "\${filename}_genomic.fna.gz" ]; then
                                echo "Corrupted files detected. Cleaning up and retrying..."
                                rm -f "\${filename}_genomic.gbff.gz" "\${filename}_genomic.fna.gz"
                                sleep 3
                                attempt=\$((attempt + 1))
                        else
                                #change file extension and name
                                name=\$(head -n 1 "\${filename}_genomic.fna" | awk -F',' '{print \$1}' | sed -e 's/^>//' -e 's/[ -/().+]/_/g')

                                mv "\${filename}_genomic.gbff" ./"\${name}.gb"
                                mv "\${filename}_genomic.fna" ./"\${name}.fasta"

                                return 0
                        fi
                done

                echo "Failed to download \$filename after \$max_attempts attempts. Skipping."
                return 1


        }
        export -f download_and_check
        find ./links/ -name "*.map" | parallel -j $parallel download_and_check


	cat .command.out >> get_data.log
	"""
}
