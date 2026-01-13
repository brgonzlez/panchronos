/*
 * REMOVE_REDUNDANCY{} process will find and remove redundant plasmid or chromosome entries and output a nonRedundant.fasta file.
 */

process REMOVE_REDUNDANCY {

	conda "${projectDir}/envs/remove_redundancy.yaml"

	input:
	tuple path(fasta_files), path(gb_files)
	val threads
	val parallel

	output:
	tuple path("cleaned/*fasta"), path("cleaned/*gb"), emit: nonRedundant_files

	script:
	"""
	#!/bin/bash


	echo -e "Concatenating, formatting and replacing unusual characters in FASTA sequences.\n"
	# First I need to clean up sequences by removing unusual characters
	sanitising_nucleotides() {
	file=\$1

		seqtk seq "\${file}" | awk '/^>/ {
	   	print \$0
	   	next
		} 
		{ 
	    	# Replace unusual characters with N
	    	gsub(/[^AaCcGgTtNn]/, "N")
	    	print
		}' > "\${file%.fasta}_cleaned.fasta" && mv "\${file%.fasta}_cleaned.fasta" ./"\${file}"
	}
	export -f sanitising_nucleotides
	find ./ -name '*fasta' | parallel -j $parallel sanitising_nucleotides


	echo -e  "Clustering sequences\n"

	find ./ -name "*.fasta" > genomesList.txt
	fastANI --ql genomesList.txt --rl genomesList.txt --threads $threads -o fastaniOutput.txt


	# first step: Ignore identical pairs, find samples with 100% identity and only pick one representative.

	awk '\$3 == 100 && \$1 != \$2 {

        if(!(\$1 in groups)) {
                groups[\$1] = \$1
        }

        if(index(groups[\$1], \$2) == 0) {
                groups[\$1] = groups[\$1]","\$2
        }

        if(!(\$2 in groups)) {
                groups[\$2] = \$2
        }

        if(index(groups[\$2], \$1) == 0) {
                groups[\$2] = groups[\$2]","\$1
        }
	}

	END {

        if (length(groups) == 0) {
        	print "There were no cluster of samples with 100% identity"
        } else {
        	for(id in groups) {
                	print "Initial cluster:" groups[id] "\n"
                	n = split(groups[id], arr, ",")
                        	asort(arr)  # sort the array
                        	normalized = ""
                        	for(i=1; i<=n; i++) normalized = normalized (i==1?"":",") arr[i]
                        	unique_groups[normalized] = 1  # store only unique sets
                	}
	
        	for(cluster in unique_groups) {
                	print "Non-redundant cluster:" cluster
	
	
                	#find if there are any samples that starts with NZ_ or NC_ and print it, if not, print the first one from the cluster.
                	n = split(cluster, nonRedArr, ",")
                	rep = nonRedArr[1]
	
                	for(i=1;i<=n;i++) {
                        	if(nonRedArr[i] ~ /^NZ_/ || nonRedArr[i] ~ /^NC_/) {
                                	rep = nonRedArr[i]
                                	break
                        	}
                	}
        	print "Representative sample for cluster: " rep "\n"
        	}
    	}
	}' fastaniOutput.txt > remove_redundancy_summary.txt


	if grep -q "There were no cluster of samples with 100% identity" remove_redundancy_summary.txt; then
        	echo "There were no clusters of samples with 100% identity"
	else
        	grep "Non-redundant cluster" remove_redundancy_summary.txt | awk -F"," '{
        	# Remove prefix from first field
        	gsub(/^Non-redundant cluster:/, "", $1)
        	# Print all fields
        	for(i=1;i<=NF;i++) print $i
        	}' | sed 's|^\./||' > samples_from_cluster.txt
	
        	grep "Representative sample for cluster" remove_redundancy_summary.txt | awk '{print $NF}' | sed 's|^\./||' > to_keep.txt
        	grep -Fv -f to_keep.txt samples_from_cluster.txt > redundants.txt    #remove representatives from redundants

			echo -e "Removing redundant sequences\n"

			while read -r ID; do
				rm "\${ID}" "\${ID%.fasta}.gb" 
			done < redundants.txt 
	fi
	
	mkdir -p cleaned

	mv *fasta cleaned/
	mv *gb cleaned/

	mkdir -p ${params.output}/DOWNLOADED

	cp -r cleaned/* ${params.output}/DOWNLOADED

	cat .command.out >> CLUSTERING.log
	"""
}
