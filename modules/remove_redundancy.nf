/*
 * REMOVE_REDUNDANCY{} process will find and remove redundant plasmid or chromosome entries and output a nonRedundant.fasta file.
 */

process REMOVE_REDUNDANCY {

	conda "${projectDir}/envs/remove_redundancy.yaml"

	input:
	tuple path(fasta_files), path(gb_files)
	val parallel

	output:
	tuple path("*fasta"), path("*gb"), emit: nonRedundant_files

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
	fastANI --ql genomesList.txt --rl genomesList.txt --threads $task.cpus -o fastaniOutput.txt

	# awk '\$3 == 100 && \$1 != \$2 fastaniOutput.txt | is a two field input with multiple lines, where lines will have repeated strings on same of different field. Example:
	# ./NZ_CP009712.1.fasta ./CP009712.1.fasta # line 1
	# ./CP009712.1.fasta ./NZ_CP009712.1.fasta # line 2
	# ./NZ_CP009999.1.fasta ./CP009999.1.fasta # line 3
	# ./CP009999.1.fasta ./NZ_CP009999.1.fasta # line 4
	# I need to make an array with KEY only names. KEY will be \$1 field. Give priority to RefSeq entries NZ and NC. First line:
	# nonRedundantSample[\$1]
	# But then if second line has \$1 OR \$2 that is already in the array, then move to the next line and DON'T append anything to the array.
	# The final array in this example should only contain two names: NZ_CP009712.1.fasta AND NZ_CP009999.1.fasta.
	# Output only redundant entries so we can exclude them.



	awk '\$3 == 100 && \$1 != \$2 {
		if (!(seenField[\$1] || seenField[\$2])) {
			if (\$1 ~ /\\/NZ_/) preferredName = \$1
			else if (\$2 ~ /\\/NZ_/) preferredName = \$2
			else if (\$1 ~ /\\/NC_/) preferredName = \$1
			else if (\$2 ~ /\\/NC_/) preferredName = \$2
			else preferredName = \$1

				seenField[\$1] = 1
				seenField[\$2] = 1
        			nonRedundant[preferredName] = 1
		}
	}	
	END {
		for (sampleId in seenField) {
			if (!(sampleId in nonRedundant))
				redundant[sampleId] = 1
		}
		for (id in redundant) print id
	}' fastaniOutput.txt > redundants.txt 

	echo -e "Removing redundant sequences\n"
	sed -i 's/\\.\\///g' redundants.txt

	while read -r ID; do
		rm "\${ID}" "\${ID%.fasta}.gb" 
	done < redundants.txt 


	cat .command.out >> CLUSTERING.log
	"""
}
