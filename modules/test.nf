/*
 * TEST{} will generate ancient DNA simulations.
 */

process TEST {
	conda "${projectDir}/envs/test.yaml"
	
	input:
	path fasta_genomes

	output:
	path 'test.fastq.gz', emit: test_data, optional: true

	script:
	"""
	#!/bin/bash

	export LD_LIBRARY_PATH=\$CONDA_PREFIX/lib:\$LD_LIBRARY_PATH
	export PATH=\$CONDA_PREFIX/bin:\$PATH

	mkdir -p bact
	mkdir -p endo
	mkdir -p cont

	to_move=\$(ls *fasta | awk 'NR == 1{print \$0}')
	mv "\$to_move" ./endo
	rm *fasta

	samtools faidx ./endo/*
	
	cp ${projectDir}/test/misincorporation.txt ./
  cp ${projectDir}/test/frag_freq.txt ./

	gargammel -mapdamage misincorporation.txt double -f frag_freq.txt -se -c 120 -rl 151  -ss MSv1 -o ./ ./

	mv _s.fq.gz ./test.fastq.gz
	"""
}
