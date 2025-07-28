# extend_sequences

'''
What do we need?
1. GFF files from prokka
2. FASTA files (with faidx)
3. final_graph.gml
4. gene_data.csv



3/ final_graph.gml

This file contains 
node [
seqIDs
]
Which we need to fetch: one seqID for each node. Then we go to gene_data.csv and get the gene TAG ID paired with seqID and we fetch the GFF line.
For each sample ID make a subset of GFF lines and output to gff. Convert that to BED. Then extend sequences. So we also need the gene name from final_graph.gml

Fields we need from final_graph:

seqIDs (hopefully only the first one. I say hopefully because what if genomic architecture is different when comparing a specific modern strain against ancient? That can be a bit cumbersome to solve)

'''



File final_graph has this structure:

  node [
    id 0
    label "251"
    size 4
    centroid "0_0_0"
    maxLenId 0
    members 0
    members 1
    members 2
    members 3
    seqIDs "3_0_0"
    seqIDs "2_0_0"
    seqIDs "0_0_0"
    seqIDs "1_0_0"
    hasEnd 1
    protein "MFVPHAKKPEIYENQRDTSLADDLSLGFTTVWNAVVSELNGESNTDDEATNDSTLVTPLTPQQRAWLNLVQPLTIIEGFALLSVPSSFVQNEIERHLRTPITDALSRRLGQQIQLGVRIAPPSTDHIDDNSSSADVLLTDDCGTDTDENYGEPLTGEYQGLPTYFTERPHHTESTVTGGTSLNRRYTFETFVIGASNRFAHAAALAIAEAPARAYNPLFIWGESGLGKTHLLHAAGNYAQRLFPGMRVKYVSTEEFTNDFINSLRDDRKVAFKRSYRDVDVLLVDDIQFIEGKEGIQEEFFHTFNTLHNANKQIVISSDRPPKQLATLEDRLRTRFEWGLITDVQPPELETRIAILRKKAQMERLAVPGDVLELIASSIERNIRELEGALIRVTAFASLNKTAIDKALAEIVLRDLIADASTMQISAATIMTATAEYFDTTIEELRGPGKTRALAQSRQIAMYLCRELTDLSLPKIGQAFGRDHTTVMYAQRKILSEMAERREVFDHVKELTTRIRQRSKR"
    dna "ATGTTTGTACCGCACGCCAAAAAGCCCGAAATTTACGAGAACCAGAGAGATACGTCGTTGGCCGATGACCTTAGTCTAGGTTTCACCACGGTTTGGAACGCAGTCGTCTCCGAACTCAACGGCGAATCCAACACAGACGACGAAGCCACCAACGACAGCACCCTAGTCACTCCGCTAACTCCTCAGCAAAGAGCATGGCTAAATCTGGTTCAACCACTCACCATCATCGAGGGATTTGCTCTTTTATCGGTGCCCAGCAGCTTTGTCCAAAATGAAATTGAACGTCATCTACGAACGCCAATCACCGATGCACTCAGCCGTCGACTCGGACAACAGATACAGCTCGGAGTCCGTATCGCACCGCCCTCTACCGACCATATTGACGACAATTCCTCGTCAGCCGACGTCCTTCTAACCGACGATTGCGGCACAGATACAGACGAAAATTACGGGGAGCCTCTTACAGGCGAGTACCAGGGTTTGCCAACCTACTTCACCGAACGTCCGCACCATACCGAATCAACCGTCACGGGAGGTACCAGCCTTAATCGCCGTTACACCTTCGAAACGTTCGTTATTGGCGCGTCGAATCGGTTCGCGCATGCTGCCGCGCTAGCGATAGCCGAAGCACCGGCCCGAGCCTACAACCCCCTTTTCATTTGGGGCGAGTCAGGTCTTGGCAAAACCCACCTATTGCACGCCGCCGGGAACTACGCACAACGACTGTTTCCCGGCATGCGGGTCAAATATGTCTCCACAGAAGAATTCACCAACGACTTCATCAACTCGCTCCGTGACGACCGCAAAGTAGCGTTCAAACGCAGCTACCGCGACGTCGATGTGCTACTGGTCGATGACATCCAATTCATCGAAGGAAAAGAAGGTATACAAGAAGAGTTCTTCCATACCTTTAATACCTTACATAACGCCAACAAGCAAATCGTCATCTCTTCTGACCGCCCACCGAAACAACTCGCCACCCTCGAAGACCGACTAAGGACCCGGTTCGAGTGGGGGCTTATTACCGACGTACAACCCCCTGAACTAGAAACCCGCATCGCCATCTTGCGTAAGAAAGCCCAGATGGAACGCCTAGCGGTGCCTGGCGATGTCCTCGAACTCATCGCCAGCAGTATCGAACGTAACATCCGTGAACTCGAGGGAGCTCTCATCAGAGTCACCGCGTTTGCTTCGCTCAACAAGACTGCAATCGACAAAGCATTAGCGGAAATCGTACTGCGTGACCTGATCGCAGACGCCAGCACGATGCAAATCAGTGCGGCAACCATAATGACAGCCACCGCCGAATACTTCGATACCACCATCGAAGAACTCCGTGGGCCAGGCAAAACCCGAGCACTGGCCCAGTCACGCCAGATCGCGATGTATTTGTGTCGTGAACTCACCGACCTCTCGCTACCCAAGATCGGCCAAGCATTCGGCCGTGACCACACCACGGTTATGTACGCACAACGAAAAATCTTGTCCGAGATGGCTGAACGTCGCGAAGTGTTCGACCACGTCAAGGAACTCACCACTCGCATCCGGCAACGGTCTAAGCGCTGA"
    annotation "dnaA"
    description "Chromosomal replication initiator protein DnaA"
    lengths 1566
    lengths 1566
    lengths 1566
    lengths 1566
    longCentroidID 1566
    longCentroidID "0_0_0"
    paralog 0
    mergedDNA 0
    genomeIDs "0;1;2;3"
    geneIDs "3_0_0;2_0_0;0_0_0;1_0_0"
    degrees 1
    name "dnaA"
  ]
  node [
    id 1
    label "483"
    size 4
    centroid "0_0_1"
    maxLenId 0
    members 0
    members 1
    members 2
    members 3
    seqIDs "3_0_1"
    seqIDs "2_0_1"
    seqIDs "1_0_1"
    seqIDs "0_0_1"
    hasEnd 0
    protein "MDLAKTNVGCSDLKFCLARESFASAVSWVAKYLPTRPTVPVLSGVLLTGSDSGLTISGFDYEVSAEVQVAAEIASSGSVLVSGRLLSDITRALPNKPVHFYVDGNRVALTCGSARFSLPTMAVEDYPTLPTLPDETGTLPSDVFAEAIGQVAIAAGRDYTLPMLTGIRIEISGDTVVLAATDRFRLAVRELKWSVLSSDFEASVLVPAKTLVEVAKAGTDGSGVCLSLGAGVGVGKDGLFGISGGGKRSTTRLLDAEFPKFRQLLPAEHTAVATIDVAELTEAIKLVALVADRGAQVRMEFGDGILRLSAGADDVGRAEEDLAVAFTGEPLTIAFNPNYLTDGLASVHSERVSFGFTTPSKPALLRPTSNDDVHPTHDGPFPALPTDYVYLLMPVRLPG"
    dna "ATGGACCTGGCCAAAACCAATGTTGGTTGCAGCGATTTAAAATTTTGTTTAGCACGTGAGTCGTTCGCCAGCGCGGTGTCATGGGTAGCGAAGTATCTCCCCACTAGACCAACGGTACCGGTGCTATCCGGCGTGCTGCTGACCGGTTCAGATAGCGGCCTAACGATCTCCGGATTCGATTACGAAGTTTCCGCAGAGGTGCAGGTTGCTGCCGAGATAGCTTCTTCTGGAAGCGTTTTGGTATCTGGGAGGTTGTTGTCAGATATTACTCGGGCGCTTCCGAACAAACCTGTTCACTTTTATGTGGACGGTAATCGGGTTGCATTGACTTGCGGAAGTGCGAGGTTCTCGTTGCCGACGATGGCGGTTGAAGACTACCCTACACTGCCTACTTTGCCGGATGAAACTGGCACGTTGCCATCGGATGTGTTCGCTGAGGCGATAGGTCAGGTTGCTATTGCGGCCGGTCGTGACTATACCTTGCCCATGCTTACCGGGATCCGGATCGAGATCTCGGGTGACACGGTAGTCTTGGCCGCCACCGATAGGTTCCGTTTGGCGGTCCGCGAGTTAAAATGGTCGGTGTTGTCATCGGATTTCGAGGCATCGGTGCTGGTGCCGGCCAAAACTTTGGTGGAAGTTGCTAAAGCCGGTACCGACGGCTCTGGTGTTTGTCTGTCATTGGGCGCTGGAGTCGGCGTTGGAAAAGATGGGCTTTTCGGCATTAGCGGTGGCGGGAAGCGCAGTACTACTCGACTTCTTGACGCTGAGTTTCCGAAATTCAGGCAGTTATTGCCGGCTGAGCACACCGCAGTGGCCACCATCGACGTGGCTGAGTTGACCGAGGCGATCAAACTGGTGGCGTTGGTAGCTGACCGTGGCGCGCAAGTGCGCATGGAATTCGGTGACGGTATATTACGGCTTTCCGCCGGCGCCGATGATGTGGGCCGAGCCGAAGAAGATCTTGCTGTTGCTTTTACTGGTGAACCGTTGACTATAGCGTTTAACCCGAATTATCTGACTGATGGACTTGCATCGGTGCATTCAGAGCGGGTGTCATTTGGTTTCACGACACCGAGTAAGCCAGCATTGCTGCGTCCAACGTCTAATGATGATGTCCATCCGACGCATGATGGCCCCTTTCCCGCACTGCCAACTGACTATGTGTATTTACTGATGCCAGTTCGGTTGCCTGGATAA"
    annotation "dnaN"
    description "Beta sliding clamp"
    lengths 1200
    lengths 1200
    lengths 1200
    lengths 1200
    longCentroidID 1200
    longCentroidID "0_0_1"
    paralog 0
    mergedDNA 0
    genomeIDs "0;1;2;3"
    geneIDs "3_0_1;2_0_1;1_0_1;0_0_1"
    degrees 2
    name "dnaN"
  ]


Given the structure of final_graph.gml, we need to read by node blocks and store relevant data. We could print each node in one line:

3_0_0 dnaA
3_0_1 dnaN


awk '
	BEGIN {OFS= "\t"}

	/^[[:space:]]*node/ {      #if node, activate the script and re-set values

		activator = 1
		seq_ID = ""
		gene_name = ""
		next

	}

	activator {

		if ($1 == "seqIDs" && seq_ID == "" && $2 ~ /^"[0-9]/) {  #id has to start with number

			seq_ID = $2

		}

		else if ($1 == "name") {

			gene_name = $2
			print seq_ID, gene_name
			activator = 0

		}
	}' final_graph.gml > seq_ids_per_node

sed -i -e 's/"//g' -e 's/ /\t/g' seq_ids_per_node


while read -r seqID gene; do

	grep -w "$seqID" gene_data.csv | awk -v gene_name=$gene -F',' '{print $1, $3, $4, gene_name}' >> updated_seq_ids_per_node

done < seq_ids_per_node

#now we read updated_seq_ids_per_node to generate new gff files

while read -r sample seqID gene_tag gene_name; do

	grep -w "$gene_tag" "${sample}".gff >> "${sample}"_selected_lines.gff

done < updated_seq_ids_per_node


#now get the bed files from _selected_lines.gff

make_bed() {
file=$1
name=$(basename "${file%_selected_lines.gff}")

awk '$3 == "CDS" || $3 == "gene" {
    sample = $1
    start = $4 - 1  #GFF is 1-based, BED is 0-based
    end = $5
    strand = $7
    attr = $9

    name = "."
    if (attr ~ /ID=/) {
        match(attr, /ID=[^;]+/)
        name = substr(attr, RSTART+3, RLENGTH-3)
    }

    print sample, start, end, name, ".", strand
}' OFS="\t" "$file" > "$name".bed

}
export -f make_bed
find ./ -name "*_selected_lines.gff" | parallel -j 10 make_bed


#finally faidx fasta files and then do bedtools slop

index_fasta() {
file=$1

	samtools faidx "$file"
}
export -f index_fasta
find ./ -name "*.fasta" | parallel -j 10 index_fasta

#now extend the sequences

extend_sequences() {
file=$1
name=$(basename "${file%.bed}")

	bedtools slop -i "$file" -g "$name".fasta.fai -b 150 > "$name"_extended_sequences.bed
	bedtools getfasta -bed "$name"_extended_sequences.bed -fi "$name".fasta -name > "$name"_extended_sequences.fasta
}
export -f extend_sequences
find ./ -name "*.bed" | parallel -j 10 extend_sequences

#now we need the gene names back!


rename_extended_sequences() {
file=$1
name=$(basename "${file%_extended_sequences.fasta}")

	awk '

		FNR==NR {

			gene_tag = $3
			gene_name = $4
			pair[gene_tag] = gene_name
			next
		}
		
		/^>/ { 
	
			split(substr($0,2), subStrings, "::")
			gene_tag_header = subStrings[1]

			if (gene_tag_header in pair) {

				print ">" pair[gene_tag_header]

			} 

			next

		}

		{

			print

		}
		'	updated_seq_ids_per_node "$file" > "$name"_extended_reference.fasta

}
export -f rename_extended_sequences
find ./ -name "*_extended_sequences.fasta" | parallel -j 10 rename_extended_sequences


#now we just concatenate them

cat *_extended_reference.fasta > extended_pangenome_reference_sequence.fasta
