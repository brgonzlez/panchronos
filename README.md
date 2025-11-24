# panchronos workflow

`panchronos` is a `nextflow` pipeline. This pipeline was developed to perform the main computation processes for the paper: "" . `panchronos` aims to perform microbial phylogenetic reconstruction based on pangenome building and it's designed to handle low quality data (typical of ancient DNA datasets).

![Workflow overview](https://github.com/brgonzlez/panchronos/blob/main/.paper-diagram-v2.drawio.png)
_Pipeline's overview_

![Workflow overview](https://github.com/brgonzlez/panchronos/blob/main/.paper-diagram-v2.drawio.svg)

# 1/ Installation

First you need to have `nextflow` and `conda/mamba` installed. Please visit: https://www.nextflow.io/docs/latest/install.html and follow the instructions. In the case of conda, if you don't have any conda version installed I would recomend miniforge (https://github.com/conda-forge/miniforge).

### Then git clone this repository:

>`$ git clone https://github.com/brgonzlez/panchronos/`

After downloading this repository you should see the folders `bin/` `config/` `envs/` and the files `main.nf` and `nextflow.config`. To test `nextflow` availability, you can try with `nextflow run main.nf --help`. 

# 2/ First steps


`panchronos` workflow assumes you already know what bacteria/virus is in your metagenomic data. If you don't know what's inside your data, you need to do metagenomic screening first.

Once you have an organism in mind, you need to search for its taxonomic ID. Additionally, do the same for another specie that you want to use as outgroup for phylogenetic reconstruction. 

# 3/ Input

Esentially, the workflow needs these things to run: {1} *.fastq files, {2} a config.tab file, {3} taxonomic IDs and {4} PATHS for your data and output folders. 

{1} The fastq files can be either compressed (.gz) or not, and they must be located in the same folder (if you have multiple datasets). `panchronos` does not support paired-end data, which means you need to collapse them first. You can include multiple unmerged single-end/collapsed libraries for the same individual (Read about `config.tab` below).

{2} config.tab file is a tab-separated text file with 3 fields: $1 Sample name, $2 Soft-clipping value and $3 Group ID. The pipeline will read this file and will perform sample-specific soft-clipping and merge aligned data by group ID (if you have different libraries that belongs to the same individual, you can specify same group ID for each). You can see an example of this file in the config/ folder.

{3} There are two taxonomic IDs that you need to collect: One for the species you want to use for pangenome construction and one for the outgroup (this is used for the sole purpose of rooting the phylogenetic tree). The taxonomic ID for your specie of interest will be handled by the workflow to download the files that will be needed for pangenome building. You can control the number of samples you want to download with the `--genomes` parameter. If you already have both `FASTA` and `GenBank` files that you want to use for pangenome building, you can specify the PATH in `--trusted_data` to let the workflow use your curated dataset instead of downloading samples for you. Beware that the `FASTA` and `GenBank` filenames must coincide and have the appropiate extensions (*.fasta and *.gb).

{4} You need to have two folders ready before running the workflow: One for your fastq/.gz data and one for storing the outputs. **Only store the data that will be included in the analysis in the data folder and nothing else.**
You can tell the pipeline where to locate these folder with `--data` and `--output` (See `configuration` section).

# 4/ Configuration

The workflow's behavior can be controlled by modifying the `nextflow.config` file. In this file, you will find every parameter that is available for fine tuning and you can directly specify your analysis settings by replacing the default values. Parameters are often linked to a specific module and you should pay special attention to the CPU usage. Besides each software thread usage the pipeline may also have threads allocated to parallel computing. 
For example, if you set up --alignment threads 10, --alignment parallel 5 and you have 5 samples in the config.tab file, the pipeline will take 50 threads to execute that particular process.


# 5/ Running the pipeline


These values can be overwritten when you call a parameter directly from the terminal:
>'nextflow run main.nf --data /new/path' 
In this case /new/path will be used by the workflow instead of the path that is specified in nextflow.config file.


# 6/ Output


# 7/ Documentation

To read a more comprehensive documentation you can follow this link (Shigeki's website).

