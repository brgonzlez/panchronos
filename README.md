# panchronos workflow

`panchronos` is a `nextflow` pipeline. This pipeline was developed to perform the main computation processes for the paper: "" . `panchronos` aims to perform microbial phylogenetic reconstruction based on pangenome building and it's designed to handle low quality data (typical of ancient DNA datasets).

# 1/ Installation

First you need to have `nextflow` and `conda/mamba` installed. Please visit: https://www.nextflow.io/docs/latest/install.html and follow the instructions. In the case of conda, if you don't have any conda version installed I would recomend miniforge (https://github.com/conda-forge/miniforge).

### Then git clone this repository:

>`$ git clone https://github.com/brgonzlez/panchronos/`

After downloading this repository you should see the folders `bin/` `config/` `envs/` and the files `main.nf` and `nextflow.config`. To test `nextflow` availability, you can try with `nextflow run main.nf --help`. 

# 2/ First steps

`panchronos` workflow assumes you already know what bacteria/virus is in your metagenomic data. If you don't know what's inside your data, you need to do metagenomic screening first.

Once you have an organism in mind, you need to search for its taxonomic ID. Additionally, do the same for another specie that you want to use as outgroup for phylogenetic reconstruction. 

# 3/ Input/Ouput

Esentially, the workflow needs these things to run: {1} *.fastq files, {2} a config.tab file, {3} taxonomic IDs and {4} PATHS for your data and output folders. 

{1} The fastq files can be either compressed (.gz) or not, and they must be located in the same folder (if you have multiple datasets). The PATH for this folder has to be specified in the `--data` parameter (See `configuration` section). `panchronos` does not support paired-end data, but you can include unmerged libraries for the same individual (Read about `config.tab` below). **Do not** include files that you don't want to include in the analysis in the same folder.

{2} config.tab file is a tab-separated text file with 3 fields: $1 Sample name, $2 Soft-clipping value and $3 Group ID. The pipeline will read this file and will perform sample-specific soft-clipping and merge aligned data by group ID (if you have different libraries that belongs to the same individual, you can specify same group ID for each). You can see an example of this file in the config/ folder.

{3} There are two taxonomic IDs that you need to collect: One for the species you want to use for pangenome construction and one for the outgroup (this is used for the sole purpose of rooting the phylogenetic tree). The taxonomic ID for your specie of interest will be handled by the workflow to download the files that will be needed for pangenome building. You can control the number of samples you want to download with the `--genomes` parameter. If you already have both `FASTA` and `GenBank` files that you want to use for pangenome building, you can specify the PATH in `--trusted_data` to let the workflow use your curated dataset instead of downloading samples for you. Beware that the `FASTA` and `GenBank` filenames must coincide and have the appropiate extensions.










# 4/ Parameters

The workflow's behavior can be controlled by modifying the `nextflow.config` file.









### Making the `config.tab` file

The pipeline reads your fastq files through the config.tab file. An example of the file's format can be found in the `config/` folder. The first column correspond to the sample name, second column to the amount of softclipping you want to apply on each end for every read and the third column is a group identification label. 

globalMeanCoverage.txt contains:

sampleID

sampleCoverage = Total number of bases of the reference genome genes where there is at least 1 read covering.

refCount = Total number of bases of the reference genome, including genes that are not being covered by reads.

globalMean = Mean depth of coverage of every gene where there is at least 1 read covering.

geneNormalizedSummary.txt contains:

sampleID
gene = Gene name

normalizedGeneSimple = (geneMeanDepth / globalMean)
-	geneMeanDepth -> Mean Depth of coverage of a particular gene.

normalizedGeneScaled = (geneMeanDepth / globalMean) * geneLength[gene]
-	geneLength[gene] -> Length of a particular gene.

normalizedGenomeSimple = (geneMeanDepth / globalMean) * (geneLength[gene] / sampleCoverage)

normalizedGenomeScaled = (geneMeanDepth / globalMean) * (geneLength[gene] / refCount)




# 4/ TODO


Add paired end reads options.


Add heteroplasmy process before getting genotypes and update heterozygosis values into the updatedNormalized table. 


Make documentation and improve --help message.


Add a small test run with tiny dataset.


Fix python heatmap labels related to number of samples. Make labels a bit bigger.


Improve the diagram by making letters bigger.


IMPORTANT: In COVERAGE_BOUNDS we are indeed removing the genes that do not satisfy the completeness threshold but we are currently NOT masking those sequences after genotyping. This has to be fixed.
