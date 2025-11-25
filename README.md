# panchronos workflow

`panchronos` is a `Nextflow` pipeline designed to perform the core computational analyses for the paper: _"The genomic identity of early smallpox in South America."_  It provides an end-to-end workflow for microbial phylogenetic reconstruction based on pangenome building, with a particular focus on handling low quality and fragmented data (typical of ancient DNA datasets).

![Workflow overview](https://github.com/brgonzlez/panchronos/blob/main/.paper-diagram-v2.drawio.png)
_Pipeline's overview_

# 1/ Installation

First you need to have `Nextflow` and `Conda/Mamba` installed. To install `Nextflow`, please follow the instructions: https://www.nextflow.io/docs/latest/install.html. For `Conda`, if you don't already have any version installed, we recomend using Miniforge: https://github.com/conda-forge/miniforge

### Then git clone this repository:

>`$ git clone https://github.com/brgonzlez/panchronos/`

After cloning the repositor, you should see the directories `bin/` `config/` `envs/`, along with the files `main.nf` and `nextflow.config`. To veryfy that `Nextflow` is properly installed and availabile in your environment, please run: `nextflow info` 

# 2/ First steps

The `panchronos` workflow assumes that you already know which bacteria or virus is present in your metagenomic dataset. If you are unsure, you should perform metagenomic screening beforehand to identify candidate organisms.

Once you have selected the organism of interest, look up its NCBI Taxonomic ID. You will also need the Taxonomic ID of a closely related species to use as an outgroup for phylogenetic reconstruction. 

Before you run the workflow with your own data, you can perform an optional test run to ensure that the pipeline behaves as expected:

>`nextflow run main.nf --test`

This test run will install the required `Conda` environments, which will be reused for future executions. Running the test is optional—these environments will also be created automatically the first time you run the pipeline with real data.

# 3/ Input

To run `panchronos`, the workflow requires the following four components:
(1) `.fastq` files
(2) `config.tab` file
(3) Taxonomic IDs
(4) Paths to your data and output directories

# 1. FASTQ files
- The workflow accepts compressed (`.gz`) or uncompressed `.fastq` files.
- All input FASTQ files should be placed in the same directory if you are analysing multiple datasets.
- `panchronos` does not support paired-end data. If you have paired-end reads, you must collapse/merge them beforehand.
- Multiple single-end/collapsed libraries from the same individual can be included—see `config.tab` below for grouping instructions.

# 2. `config.tab` file
`config.tab` is a tab-separated text file with three fields:
| Sample name | Soft-clipping value | Group ID |
|----------|----------|----------|
| Sample_A.fastq | 2 | Neolithic_pestis  |

The workflow uses this file to:
- Apply sample-specific soft clipping, and
- Merge aligned data belonging to the same individual (i.e., those sharing the same Group ID).
You can see an example file in the `config/` directory of the repository.

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
