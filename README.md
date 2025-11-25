# panchronos workflow

`panchronos` is a `Nextflow` pipeline designed to perform the core computational analyses for the paper: _"The genomic identity of early smallpox in South America."_  It provides an end-to-end workflow for microbial phylogenetic reconstruction based on pangenome building, with a particular focus on handling low quality and fragmented data (typical of ancient DNA datasets).

![Workflow overview](https://github.com/brgonzlez/panchronos/blob/main/.paper-diagram-v2.drawio.png)
_Pipeline's overview_


# 1/ Installation

First you need to have `Nextflow` and `Conda/Mamba` installed. To install `Nextflow`, please follow the instructions: https://www.nextflow.io/docs/latest/install.html. For `Conda`, if you don't already have any version installed, we recomend using Miniforge: https://github.com/conda-forge/miniforge

### Then git clone this repository:

>`$ git clone https://github.com/brgonzlez/panchronos/`

After downloading this repository you should see the folders `bin/` `config/` `envs/` and the files `main.nf` and `nextflow.config`. To test `nextflow` availability, you can try with `nextflow info` on the terminal. 

# 2/ First steps


`panchronos` workflow assumes you already know what bacteria/virus is in your metagenomic data. If you don't know what's inside your data, you need to do metagenomic screening first.

Once you have an organism in mind, you need to search for its taxonomic ID. Additionally, do the same for another specie that you want to use as outgroup for phylogenetic reconstruction. 

Before you run the workflow with your data, you can try and do a test run (this is optional) to check that the pipeline intended behaviors are working well. To do so, you can try:

>`nextflow run main.nf --test`

This will install the neccessary conda environments that will be use whenever you run the workflow. This is not mandatory as they will be created on any first run you do.

# 3/ Input

Esentially, the workflow needs these things to run: {1} *.fastq files, {2} a config.tab file, {3} taxonomic IDs and {4} PATHS for your data and output folders. 

{1} The fastq files can be either compressed (.gz) or not, and they must be located in the same folder (if you have multiple datasets). `panchronos` does not support paired-end data, which means you need to collapse them first. You can include multiple unmerged single-end/collapsed libraries for the same individual (Read about `config.tab` below).

{2} config.tab file is a tab-separated text file with 3 fields: $1 Sample name, $2 Soft-clipping value and $3 Group ID. The pipeline will read this file and will perform sample-specific soft-clipping and merge aligned data by group ID (if you have different libraries that belongs to the same individual, you can specify same group ID for each). You can see an example of this file in the config/ folder.

{3} There are two taxonomic IDs that you need to collect: One for the species you want to use for pangenome construction and one for the outgroup (this is used for the sole purpose of rooting the phylogenetic tree). The taxonomic ID for your specie of interest will be handled by the workflow to download the files that will be needed for pangenome building. You can control the number of samples you want to download with the `--genomes` parameter. If you already have both `FASTA` and `GenBank` files that you want to use for pangenome building, you can specify the PATH in `--trusted_data` to let the workflow use your curated dataset instead of downloading samples for you. Beware that the `FASTA` and `GenBank` filenames must coincide and have the appropiate extensions (*.fasta and *.gb).

{4} You need to have two folders ready before running the workflow: One for your fastq/.gz data and one for storing the outputs. **Only store the data that will be included in the analysis in the data folder and nothing else.**
You can tell the pipeline where to locate these folder with `--data` and `--output` (See Configuration section below).

# 4/ Configuration

The workflow's behavior can be controlled by modifying the `nextflow.config` file. In this file, you will find every parameter that is available for fine tuning. You can directly specify your analysis settings by replacing the default values. 

Parameters are often linked to a specific module and you should pay special attention to the CPU usage. Besides each software thread usage the pipeline may also have threads allocated to parallel computing. 
For example, if you set up `--alignment_threads 10`, `--alignment_parallel 5` and you have 5 samples in the config.tab file, the pipeline will take 50 threads to execute that particular process.


**NOTE: Don't assign more than 3 threads to --get_data_parallel as NCBI will deny a request if it is >3.**

# 5/ Running the pipeline


If you used `nextflow.config` file to specify your settings (including PATHs), then you can simply exectue the workflow with: 


>`nextflow run main.nf -resume` 

I would recommend including the option `-resume` on every run as it will make good use of nextflow's cache system. 
You can also directly set up a parameter value through the terminal by adding double dashes before the parameter name. Note that these values will overwrite whatever is in `nextflow.config` file:

>`nextflow run main.nf --data /new/path/ --panaroo_alignment_type core -resume`

In this case, the values "/new/path/" and "core" will be used by the workflow instead of the ones that are specified in nextflow.config file for those particular parameters. 
Aditionally, when you execute the workflow it will print out every parameter that will be use in the current run, so you can double check your settings.


# 6/ Output

Once the workflow is running, important outputs will start populating the `--output` directory as soon as those files are generated and they will be organised by folders:

| Folder | Description | Format |
| :--- | :--- | :---: | 
|ALIGNMENT| Aligned data against pangenome reference sequences | .bam , .fastq |
|DOWNLOADED| Input data for pangenome building. If `--trusted_data` was used it will store your own curated datasets | .fasta , .gb |
|GENE_DATABASE| Guide file for prokka. Includes clustered genes from `GenBank` files prior to annotation | .fasta |
|GENE_MSA| Each individual gene multiple sequence alignment | .fasta |
|GENOTYPING| Results from variant calling, including final consensus sequences per sample per gene plus vcf files. | .fasta, .vcf* |
|MAPDAMAGE| Results from damage pattern assessment using `mapDamage` | Multiple files |
|MATRIX| Final matrix of presence/absence of genes | .tab |
|PANGENOME| Outputs from pangenome building and extending. Included are the original Panaroo reference genome plus the extended/unextended versions. Additionally, a BLAST database from the extended version of the pangenome is included | .fasta, BLAST database |
|PLOTS| Presence/absence heatmap and coverage vs completeness plots | .png |
|STATS| Several files with basic statistics, particularly showcasing normalisation coverage and completeness per gene per sample | .tab , .txt |
|TREE| Results from IQTREE. The concatenated MSA files that were used to reconstruct the phylogenies are also included | .fasta , .treefile |



# 7/ Documentation


If unaware of the workflow's settings, try out:

>`nextflow run main.nf --help`

Or you can follow this link (Shigeki's website) for a more comprehensive documentation.
