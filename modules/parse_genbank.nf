/*
 * PARSE_GENBANK() will perform genbank format quality check.
 */

 
process PARSE_GENBANK {
	conda "${projectDir}/envs/biopython.yaml"
	
	input:
	path gffFiles
	path fastaFiles
		
	output:
	tuple path('*fasta'), path('*gb'), emit: validFiles

	script:
	"""
	parseTest.py > parseTest.txt
	grep "is not a valid GenBank file" parseTest.txt | awk '{print \$1}' > blackListed.txt

	if [[ -s blackListed.txt ]] ;then

		while read -r removeMe; do
			rm "\${removeMe%gb}fasta" *"\$removeMe"	
  			echo -e "\${removeMe%.gb} sample has been removed due to format problems."
		done < blackListed.txt
	else
		echo -e "Every file passed the test. Moving on."
	fi

	cat .command.out .command.err >> fastaDatabase.log
	"""
}
