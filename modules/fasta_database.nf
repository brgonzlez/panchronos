/*
 * FASTA_DATABASE() will perform genbank format quality check and build a gene FASTA database for prokka proteins guide.
 */

 
process FASTA_DATABASE {
	conda "${projectDir}/envs/biopython.yaml"
	
	input:
	path gffFiles
	path fastaFiles
 	val parallel
		
	output:
	path 'clustered_sequences.fasta' , emit: fastaDatabase
	tuple path('*fasta'), path('*gb'), emit: validFiles

	script:
	"""
	parseTest.py > parseTest.txt
	grep "is not a valid GenBank file" parseTest.txt | awk '{print \$1}' | sed -e 's/gff\\///g' > blackListed.txt
	
	while read -r removeMe; do
		rm "\${removeMe%gb}fasta" *"\$removeMe"	
  		echo -e "\${removeMe%.gb} sample has been removed due to format problems."
	done < blackListed.txt

	parsing_and_contatenating.py

	cat .command.out .command.err >> fastaDatabase.log
	"""
}
