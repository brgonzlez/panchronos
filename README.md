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

# (1). FASTQ files
- The workflow accepts compressed (`.gz`) or uncompressed `.fastq` files.
- All input FASTQ files should be placed in the same directory if you are analysing multiple datasets.
- `panchronos` does not support paired-end data. If you have paired-end reads, you must collapse/merge them beforehand.
- Multiple single-end/collapsed libraries from the same individual can be included—see `config.tab` below for grouping instructions.

# (2). `config.tab` file
`config.tab` is a tab-separated text file with three fields:
| Sample name | Soft-clipping value | Group ID |
|----------|----------|----------|
| Sample_A.fastq | 2 | Neolithic_pestis  |
| Sample_B.fastq | 2 | Neolithic_pestis  |
| Sample_C.fastq | 5 | Neolithic_pestis  |

The workflow uses this file to:
- Apply sample-specific soft clipping, and
- Merge aligned data belonging to the same individual (i.e., those sharing the same Group ID).
You can see an example file in the `config/` directory of the repository.

# (3). Required taxonomic IDs
You need to provide two NCBI taxonomic IDs:
1. Target species taxonomic ID — used for pangenome construction. The workflow will automatically download the necessary genomic data based on this ID.
2. Outgroup taxonomic ID — used only for rooting the phylogenetic tree.
You can control how many genomes to download for pangenome building using the `--genomes` parameter.

If you already have your own curated dataset (FASTA + GenBank files), you can provide its path using `--trusted_data`.
Important requirements:
- FASTA and GenBank filenames must match exactly (aside from extensions).
- Required extensions: `*.fasta` and `*.gb`

# (4). Data and output directories
Before running the workflow, prepare two directories:
- A folder containing only your FASTQ/FASTQ.gz (i.e. input) files. **Only store the data that will be included in the analysis.**
- A separate folder where pipeline outputs will be written.

Provide their paths with the parameters:
- `--data` (input FASTQ directory)
- `--output` (output directory)

See the Configuration section below.

# 4/ Configuration
The workflow's configuration is controlled through the `nextflow.config` file. This file contains all available parameters, allowing you to fine-tune the pipeline according to your analysis needs. You can adjust the default values directly to match your setup and computational resources.

Many parameters are associated with specific workflow modules, and CPU/thread usage is particularly important to configure carefully. In addition to threads used by individual tools, the pipeline may spawn parallel jobs, multiplying the total number of threads in use.
For example:
- `--alignment_threads 10`
- `--alignment_parallel 5`
- 5 samples in the `config.tab`

This combination means the alignment process will use 10 × 5 = 50 threads simultaneously.

**Important NOTE:**
**Don't assign more than 3 threads to `--get_data_parallel`.**

NCBI will reject download requests if the number of parallel queries exceeds 3.

# 5/ Running the pipeline
Once you have configured your settings in the `nextflow.config` file (including all required paths), you can exectue the workflow with: 

>`nextflow run main.nf -resume` 

It is recommended to include the `-resume` flag for every run, as this enables Nextflow’s caching system and prevents unnecessary recomputation.

**Overriding parameters from the command line**

You can override any parameter directly from the terminal by prefixing it with `--`.
Values provided this way take precedence over those defined in the `nextflow.config` file:

>`nextflow run main.nf --data /new/path/ --panaroo_alignment_type core -resume`

In this example, the values `/new/path/` and `core` will be used instead of the corresponding entries in `nextflow.config`.

Whenever the workflow starts, Nextflow prints a summary of all parameters used for the current run.
This allows you to quickly verify that your settings are being applied correctly.

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
