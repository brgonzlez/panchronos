/*
 * FILTER_GENE_ALIGNMENTS{} will complete each gene MSA and parse them to find artifacts.
 */


process FILTER_GENE_ALIGNMENTS {
	conda "${projectDir}/envs/filterGeneAlignments.yaml"

	input:
	path genesAln, stageAs: 'panaroo_genes/*'
 	path extractedSequencesFasta, stageAs: 'user_genes/*'
	path fFiles, stageAs: 'modern_data/*'
	val genomes
	path outgroupSeq, stageAs: 'outgroup'
	path blackListed, stageAs: 'blackListed.txt'
	val parallel
	path final_list_genes

	output:
	path 'sorted/*.fasta', emit: genesAlnSeq
	path 'sampleNames.txt', emit: sampleNames

	script:
	"""
	#!/bin/bash

	shopt -s nullglob  #Prevent literal interpretation of wildcards if there are no matchings
	#sanity check before starting the process

	if ! find panaroo_genes/ -name "*.aln.fas" | grep -q .; then
    		echo -e "No files with .aln.fas extension were found in panaroo_genes/. Check previous process. Stopping the pipeline."
    		exit 1
	else
		echo -e "Fasta files with expected extension .aln.fas were found. Proceeding with the process."
	fi

	if [[ -e outgroup ]]; then
	        sed -i -e 's/~/_/g' outgroup
	else
        	echo -e "outgroup file was not found. Stopping the pipeline."
        	exit 1
	fi

	if ! find modern_data/ -name "*fasta" | grep -q .; then
    		echo -e "No modern sequences with .fasta extension were found in modern_data/. Check previous process. Stopping the pipeline."
    		exit 1
	else
        	echo -e "Fasta files with expected extension .fasta were found. Proceeding with the process."
	fi

	#Remove duplicated genes from the dataset
	mkdir ./redundant_genes

	#valid genes for matching
	grep -F -f $final_list_genes /dev/null >/dev/null 2>&1

	#move all files NOT in the list
	for file in panaroo_genes/*; do
	    base=\$(basename "\$file")
	    gene=\${base%%.*}
	
	    # If gene not in list, move it
    	if ! grep -Fxq "\$gene" $final_list_genes; then
	    	   	mv "\$file" ./redundant_genes/
    	fi
	done

	# modern_samples_list() will create a text file with modern genomes names
	modern_samples_list() {
	fasta=\$1
        	name=\$(basename "\${fasta%.fasta}")
        	echo "\${name}" > "\${name}"_modern_TMP
	}
	export -f modern_samples_list
	find modern_data/ -name "*.fasta" | parallel -j $parallel modern_samples_list

	cat *_modern_TMP >> modernSampleNames.txt
	rm *_modern_TMP


	# Adding the outgroup to this as it is modern too
	echo outgroup >> modernSampleNames.txt
	echo -e "Standardise fasta suffixes"

	panaroo_fasta_suffix() {
	fasta_file=\$1
	filename=\$(basename "\${fasta_file%.aln.fas}")

		mv "\${fasta_file}" panaroo_genes/"\${filename}.fasta"
	}
	export -f panaroo_fasta_suffix
	find panaroo_genes/ -name "*.aln.fas" | parallel -j 10 panaroo_fasta_suffix

	echo -e "Done"

	echo -e "Fixing FASTA headers and sequences formatting with seqtk in existing gene alignments"

	mkdir -p panaroo_parsed
	parsing_panaroo_msa() {
	fasta_file=\$1
                name=\$(basename "\${fasta_file%.fasta}" | sed -e 's/~/_/g') # Also replace  ~ characters with _, with double % to remove two dots
                seqtk seq "\${fasta_file}" | awk '/^>/ {sub(/;.*/, "", \$0)} {print}' > panaroo_parsed/"\${name}_parsing_panaroo.fasta"
	}
	export -f parsing_panaroo_msa
	find panaroo_genes/ -name "*.fasta" | parallel -j $parallel parsing_panaroo_msa

	echo -e "Done"



	mkdir -p blacklisted
	echo -e "Checking if there are low quality samples"
	if [[ -s blackListed.txt ]]; then
	        while read -r removeMe; do
        	        mv user_genes/*"\${removeMe}.fasta" blacklisted/"\${removeMe}"
                	echo "\${removeMe} has been removed from analysis due to low quality."
        	done < blackListed.txt
	else
	        echo -e "Every sample passed quality checks."
	fi

	#index_and_formatting() send user samplenames to userSampleNames.txt and replaces ~ from gene names with _
	index_and_formatting() {
	sample=\$1
        	sampleName=\$(basename "\${sample%.fasta}")
        	echo "\${sampleName}" >> userSampleNames.txt
		sed -i -e 's/~/_/g' "\${sample}"
	}
	export -f index_and_formatting
	find user_genes/ -name "*.fasta" | parallel -j $parallel index_and_formatting

	echo -e "Done"


	# Make a file with every sample combined
	cat modernSampleNames.txt userSampleNames.txt > sampleNames.txt


	echo -e "Adding user sample genes sequences to each particular gene MSA and replace gene name with sample name"

	#add_user_sample_sequences() adds user sample gene sequences into each panaroo msa and replaces gene name with user sample name
	add_user_sample_sequences() {
	fasta_file=\$1
	        if [[ -e "\$fasta_file" ]]; then
                	name=\$(basename "\${fasta_file%_parsing_panaroo.fasta}")

                	while read -r sampleName; do
                                grep -w -A 1 "\$name" "user_genes/\${sampleName}.fasta" | awk -v newHeader="\$sampleName" '/^>/ {sub(/^>.*/, ">" newHeader, \$0)} {print}' >> "\${fasta_file}"
                	done < userSampleNames.txt
        	else
                	echo -e "There are no files with .fasta extension in panaroo_parsed/ folder. Stopping the process.\n"
			exit 1
        	fi
	}

	export -f add_user_sample_sequences
	find panaroo_parsed/ -name "*.fasta" | parallel -j $parallel add_user_sample_sequences
	echo -e "Done"


	echo -e "Adding outgroup gene sequences into each gene MSA file"
	#similarly as add_user_sample_sequences(), add_outgroup_sequences() will add outgroup aligned gene sequences into each gene panaroo msa file
	add_outgroup_sequences() {
	fasta_file=\$1

                name=\$(basename "\${fasta_file%_parsing_panaroo.fasta}")
                grep -w -A 1 "\$name" outgroup | awk -v outgroup="outgroup" '/^>/ {sub(/^>.*/, ">" outgroup, \$0)} {print}' >> "\${fasta_file}"
	}
	export -f add_outgroup_sequences
	find panaroo_parsed/ -name "*_parsing_panaroo.fasta" | parallel -j $parallel add_outgroup_sequences

	echo -e "Done"


	echo -e "Padding incomplete sequences to the right"

	panaroo_parsed_padding() {
	MSA=\$1
        	name=\$(basename "\${MSA}")
        	seqLength=\$(awk '!/^>/ { if (length > max) max = length } END { print max }' "\$MSA")

        	awk -v seq_length="\$seqLength" '{
	                if (\$0 ~ /^>/) {
                        	print
                	} else {
                	while ( length(\$0) < seq_length) {
	                        \$0 = \$0 "n"
                	}
                	print
                	}
        	}' "\${MSA}" > tmp_"\${name}" && mv tmp_"\${name}" "\${MSA}"
	}
	export -f panaroo_parsed_padding
	find panaroo_parsed/ -name "*parsing_panaroo.fasta" | parallel -j $parallel panaroo_parsed_padding

	echo -e "Done"



	echo -e "Checking if there are Panaroo headers artifacts"

	check_panaroo_headers_artifacts() {
	file_msa=\$1

	        while read -r sampleName; do
        	        header_name=">\$sampleName"
			matched_header=\$(grep "\$sampleName" "\$file_msa")

                	if [[ -n "\$matched_header" ]]; then

                        	while IFS= read -r matched; do

                                	if [[ "\$matched" != "\$header_name" ]]; then
                                        	sed -i -e "s/\${matched}/\${header_name}/g" "\$file_msa"
                                        	echo -e "\$matched name was found in \$file_msa but it should have been \$header_name instead. Fixed"
                                	fi
                        	done <<< "\$matched_header"
                	fi

        	done < modernSampleNames.txt
	}
	export -f check_panaroo_headers_artifacts
	find panaroo_parsed/ -name "*_parsing_panaroo.fasta" | parallel -j $parallel check_panaroo_headers_artifacts
	echo -e "Done"



	#modern_append_gaps() will check if there are modern strains missing in panaroo_parsed alignments. If yes, append the sample and padding with -
	echo -e "Padding missing samples. Gaps for modern strains and n for ancient samples"

	modern_append_gaps() {
	fasta_file=\$1

		seq_length=\$(awk '!/^>/ { if (length > max) max = length } END { print max }' "\$fasta_file")

                while read -r strain; do
                        if ! grep -wq "\$strain" "\$fasta_file"; then
                                echo ">\$strain" >> "\$fasta_file"
                                gaps_length=\$(printf '%*s' "\$((seq_length))" | tr ' ' '-')
                                echo "\$gaps_length" >> "\$fasta_file"
                        fi
                done < modernSampleNames.txt
	}
	export -f modern_append_gaps
	find panaroo_parsed/ -name "*_parsing_panaroo.fasta" | parallel -j $parallel modern_append_gaps



	#ancient_append_missingness() will check if there are ancient strains missing in panaroo_parsed alignments. If yes, append the sample and fill it with n as we don't know the ancient status for that gene.
	ancient_append_missingness() {
	fasta_file=\$1

	        seq_length=\$(awk '!/^>/ { if (length > max) max = length } END { print max }' "\$fasta_file")

		while read -r strain; do
                	if ! grep -wq "\$strain" "\$fasta_file"; then
                        	echo ">\$strain" >> "\$fasta_file"
                                fakeSeq=\$(printf '%*s' "\$((seq_length))" | tr ' ' 'n')
                                echo "\$fakeSeq" >> "\$fasta_file"
                        fi
                done < userSampleNames.txt
	}

	export -f ancient_append_missingness
	find panaroo_parsed/ -name "*_parsing_panaroo.fasta" | parallel -j $parallel ancient_append_missingness

	echo -e "Done"


	# At this point there will be many msa that are finished: number of headers == number of samples AND headers names == all samples names. I need to find them based on this logic and send them
	# to the final folder. There will be some msa with header duplication problems and they need to be send to special_cases.

	mkdir -p sanitised_msa

	sanitised_msa() {
	msa_file=\$1

		name=\$(basename "\${msa_file%_parsing_panaroo.fasta}.fasta")
	        sample_count=\$(grep -c '^>' "\$msa_file")
	        total_sample_count=\$(wc -l < sampleNames.txt)

		if [[ "\$sample_count" -ne "\$total_sample_count" ]]; then  #if they are equal then keep going. if not skip current file.
			return
		fi

		missing_samples=0
		while read -r strain; do
			if ! grep -wq "\$strain" "\$msa_file"; then
				missing_samples=1
				break
			fi
		done < sampleNames.txt

		if [[ "\$missing_samples" -eq 0 ]]; then
			mv "\$msa_file" "sanitised_msa/\${name}"
		fi

	}
	export -f sanitised_msa
	find panaroo_parsed/ -name "*_parsing_panaroo.fasta" | parallel -j $parallel sanitised_msa

	# Dealing with fragmented Panaroo gene alignments multi-entries
	mkdir -p special_cases
	
	echo -e "Checking if there are leftovers after sanitising"

	unsanitised=\$(ls panaroo_parsed/*fasta | wc -l )
	if (( unsanitised != 0 )); then
		echo -e "There are \$unsanitised samples that have problems in panaroo_parsed/ directory. Inspecting them. ."
	else
		echo -e "It seems every gene msa was succesfully cleaned. Moving on."
	fi

	fix_duplicated_entries() {
	fasta_file=\$1

               	if [[ -e "\$fasta_file" ]] ; then
                       	name=\$(basename "\${fasta_file%_parsing_panaroo.fasta}")
               		# Identifying repeated headers and save them to .dup
                       	awk '/^>/ {count[\$0]++} END {for (header in count) if (count[header] > 1) print substr(header, 2)}' "\$fasta_file" > special_cases/"\${name}.dup"
	
               		# Extracting sequences for each repeated entry into temporary files
                	while read -r entry; do
                        	awk -v sampleName="\$entry" '
                        	\$0 ~ "^>" sampleName {print_header=1; next}
                        	/^>/ {print_header=0}
                        	print_header {print}' "\$fasta_file" > special_cases/"\${name}_duplicated_\${entry}"
                       done < special_cases/"\${name}.dup"
               # Remove repeated entries and their sequences from the original FASTA file
                       awk -v dpdFile=special_cases/"\${name}.dup" '
                       BEGIN {
                       while (getline < dpdFile) {
                               repeated[\$0] = 1
                               }
                       }
                       /^>/ { header = substr(\$0, 2)
                       if (repeated[header]) {
                               skip = 1
                       } else {
                               skip = 0}
                       }
                       !skip' "\$fasta_file" > special_cases/"\${name}.fasta"

                       for indexSeqs in special_cases/"\${name}_duplicated_"*; do
                               geneName=\$(basename "\${indexSeqs%_duplicated_*}")
                               sampleName=\$(basename "\${indexSeqs##*_duplicated_}")

                              sequence=\$(awk '
                                {
                                    seq_length = length(\$0)
                                    line[NR] = \$0
                                    num_lines = NR
                                }

                                END {
                                    final = ""
                                    for (pos = 1; pos <= seq_length; pos++) {
                                        merged_char = "-"
                                        for (l = 1; l <= num_lines; l++) {
                                            c = substr(line[l], pos, 1)
                                            if (c ~ /[aAtTgGcC]/) {
                                                merged_char = c
                                                break
                                            }
                                        }
                                        final = final merged_char
                                    }

                                    # Remove dashes
                                    gsub(/-/, "", final)
                                    print final
                                }' "\$indexSeqs")



                       # Finally add the selected sequence back to the cleaned original FASTA file with the header as well
                               echo ">\${sampleName}" >> special_cases/"\${name}.fasta"
                               echo "\$sequence" >> special_cases/"\${name}.fasta"
                       done

               else
                       echo -e "No files found with fasta extension in special_cases folder"    
       		# Cleaning temporary files
               fi
       }

	export -f fix_duplicated_entries

	if (( unsanitised != 0 )); then
		find panaroo_parsed/ -name "*_parsing_panaroo.fasta" | parallel -j $parallel fix_duplicated_entries
	fi

	#finally checking that every seq has the same length

	checking_alignment_lengths() {
	msa_file=\$1
	name=\$(basename "\${msa_file}")

               seq_length=\$(awk 'NR%2 == 0 && length > max { max = length } END { print max }' "\$msa_file")

               awk -v numCols="\$seq_length" '{
                       if (\$0 ~ /^>/) {
                               print
                       } else {
                       while ( length(\$0) < numCols) {
                               \$0 = \$0 "-"
                       }
                       print
                       }
               }' "\$msa_file" > special_cases/tmp_"\${name}" && mv special_cases/tmp_"\${name}" "\${msa_file}"
       }

	#export -f checking_alignment_lengths
	#find special_cases/ -name "*.fasta" | parallel -j $parallel checking_alignment_lengths
	#echo -e "Done"
	#NOTE: We may not need this function anymore as we are re-aligning everything with mafft.

	#now check if these cleaned fasta passed the sanitised_msa() test.
	sanitised_msa_special_cases() {
	msa_file=\$1

        	name=\$(basename "\${msa_file}")
        	sample_count=\$(grep -c '^>' "\$msa_file")
        	total_sample_count=\$(wc -l < sampleNames.txt)

        	if [[ "\$sample_count" -ne "\$total_sample_count" ]]; then  #if they are equal then keep going. if not skip current file.
	                return
	        fi

	        missing_samples=0
	        while read -r strain; do
                	if ! grep -wq "\$strain" "\$msa_file"; then
                        	missing_samples=1
                        	break
                	fi
        	done < sampleNames.txt

        	if [[ "\$missing_samples" -eq 0 ]]; then
	                mv "\$msa_file" "sanitised_msa/\${name}"
	        fi

	}
	export -f sanitised_msa_special_cases
	find special_cases/ -name "*.fasta" | parallel -j $parallel sanitised_msa_special_cases

	# final steps: sort everything.

	# Make INDEX file so we can sort every file in the same order before building MSA based on first file.
	firstFile=\$(ls -1 sanitised_msa/ | awk 'NR==1 {print \$0}') && awk '/^>/ {print \$0}' sanitised_msa/"\$firstFile" > sorting_index

	mkdir -p sorted

	echo -e "Sorting MSAs"

	sortingMSA() {
	fasta=\$1

	        name=\$(basename "\${fasta%.fasta}")

	        while read -r header; do
                	awk -v headerName="\$header" '
	                        \$0 == headerName {
                        	print \$0     # Print the matched header
                        	getline      # Get the sequence line after the header and store it in line
                        	print \$0     # Print the sequene
                        	}
                	' "\$fasta" >> sorted/"\${name}.fasta"
        	done < sorting_index
	}

	export -f sortingMSA
	find sanitised_msa/ -name "*.fasta" | parallel -j $parallel sortingMSA

	echo -e "Done"

	cat .command.out >> filterGeneAlignments.log
	"""
}
