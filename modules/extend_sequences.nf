/*
 * EXTEND_sEQUENCES{} process will extend gene sequences from pangenome reference.
 */


process EXTEND_SEQUENCES {

	conda "${projectDir}/envs/extend_sequences.yaml"

	input:
	path gff
	path fasta
	tuple path(final_graph), path(gene_data)
	val extend
	val parallel
	

	output:
	path 'extended_pangenome_reference_sequence.fasta' , emit: extended_reference

	script:
	"""
	#!/bin/bash

	#fetch data from graph
	awk '
		BEGIN {OFS= "\t"}

		/^[[:space:]]*node/ {      #if node, activate the script and re-set values
	
			activator = 1
			seq_ID = ""
			gene_name = ""
			next
	
		}
	
		activator {
	
			if (\$1 == "seqIDs" && seq_ID == "" && \$2 ~ /^"[0-9]/ && \$2 !~ /refound/) {  #id has to start with number
	
				seq_ID = \$2
	
			}
	
				else if (\$1 == "name") {
	
				gene_name = \$2
				print seq_ID, gene_name
				activator = 0
	
			}
		}' $final_graph > seq_ids_per_node
	
	sed -i -e 's/"//g' -e 's/ /\t/g' seq_ids_per_node
	
	
	while read -r seqID gene; do
	
		grep -w "\$seqID" $gene_data | awk -v gene_name=\$gene -F',' '{print \$1, \$3, \$4, gene_name}' >> updated_seq_ids_per_node
		
	done < seq_ids_per_node
	
	#now we read updated_seq_ids_per_node to generate new gff files

	while read -r sample seqID gene_tag gene_name; do

		grep -w "\$gene_tag" "\${sample}".gff >> "\${sample}"_selected_lines.gff

	done < updated_seq_ids_per_node


	#now get the bed files from _selected_lines.gff

	make_bed() {
	file=\$1
	name=\$(basename "\${file%_selected_lines.gff}")

		awk '\$3 == "CDS" || \$3 == "gene" {
    		sample = \$1
    		start = \$4 - 1  #GFF is 1-based, BED is 0-based
    		end = \$5
    		strand = \$7
    		attr = \$9

    		name = "."
    		if (attr ~ /ID=/) {
        		match(attr, /ID=[^;]+/)
        		name = substr(attr, RSTART+3, RLENGTH-3)
    		}

    		print sample, start, end, name, ".", strand
		}' OFS="\t" "\$file" > "$name".bed

	}
	export -f make_bed
	find ./ -name "*_selected_lines.gff" | parallel -j $parallel make_bed


	#finally faidx fasta files and then do bedtools slop

	index_fasta() {
	file=\$1

		samtools faidx "\$file"
	}
	export -f index_fasta
	find ./ -name "*.fasta" | parallel -j $parallel index_fasta

	#now extend the sequences
	extend_sequences() {
	file=\$1
	name=\$(basename "\${file%.bed}")

		bedtools slop -i "\$file" -g "\$name".fasta.fai -b $extend > "\$name"_extended_sequences.bed
		bedtools getfasta -bed "\$name"_extended_sequences.bed -fi "\$name".fasta -name > "\$name"_extended_sequences.fasta
	}
	export -f extend_sequences
	find ./ -name "*.bed" | parallel -j $parallel extend_sequences

	#now we need the gene names back!


	rename_extended_sequences() {
	file=\$1
	name=\$(basename "\${file%_extended_sequences.fasta}")

		awk '

			FNR==NR {

				gene_tag = \$3
				gene_name = \$4
				pair[gene_tag] = gene_name
				next
			}
			
			/^>/ { 
		
				split(substr(\$0,2), subStrings, "::")
				gene_tag_header = subStrings[1]
	
				if (gene_tag_header in pair) {
	
					print ">" pair[gene_tag_header]
	
				} 
	
				next
	
			}
	
			{
	
				print
	
			}
			'	updated_seq_ids_per_node "\$file" > "\$name"_extended_reference.fasta
	
	}
	export -f rename_extended_sequences
	find ./ -name "*_extended_sequences.fasta" | parallel -j $parallel rename_extended_sequences
	

	#now we just concatenate them

	cat *_extended_reference.fasta > extended_pangenome_reference_sequence.fasta
	"""
}
