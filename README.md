### Panchronos workflow
`panchronos` is a `Nextflow` pipeline designed to perform the core computational analyses for the paper: "_The genomic identity of early smallpox in South America_."  It provides an end-to-end workflow for microbial phylogenetic reconstruction based on pangenome building, with a particular focus on handling low-quality and fragmented data (typical of ancient DNA data).

![Workflow overview](https://github.com/brgonzlez/panchronos/blob/main/.paper-diagram-v2.drawio.png)
_Workflow's overview_

# 1/ Installation
First you need to have `Nextflow` and `Conda/Mamba` installed. To install `Nextflow`, please follow the instructions: https://www.nextflow.io/docs/latest/install.html. For `Conda`, if you don't already have any version installed, we recomend using Miniforge: https://github.com/conda-forge/miniforge

To verify that Nextflow is properly installed and available in your environment, please run: 

>`nextflow info`


### Then git clone this repository:

```
git clone https://github.com/brgonzlez/panchronos/
```

After cloning the repository, you should see the directories bin/ config/ envs/, along with the files `main.nf` and `nextflow.config`.

### Perform test run

To perform test run, please go to the main directory where you cloned `panchronos` (where the `main.nf` file is located) and execute:

**IMPORTANT: When `nextflow` is installing an environment, do not cancel it as it will produce a broken installation.**

```
dir=$(pwd)
nextflow run main.nf --test --config "$dir"/test/config_test.tab
```

This test run will install the required `Conda` environments, which will be reused for future executions. Running the test is optional—these environments will also be created automatically the first time you run the pipeline with real data.


# 3/ Input

The `panchronos` workflow assumes that you already know which bacteria or virus is present in your metagenomic dataset. If you are unsure, you should perform metagenomic screening beforehand to identify candidate organisms.

To run `panchronos`, the workflow requires the following four components:
(1) `.fastq` files
(2) `config.tab` file
(3) Taxonomic IDs
(4) Paths to your data and output directories

# (1) FASTQ files
- The workflow accepts compressed (`.gz`) or uncompressed `.fastq` or `.fq` files.
- All input FASTQ files should be placed in the same directory if you are analysing multiple datasets.
- `panchronos` does not support paired-end data. If you have paired-end reads, you must collapse/merge them beforehand.
- Multiple single-end/collapsed libraries from the same individual can be included. Please see `config.tab` below for grouping instructions.

# (2) `config.tab` file
`config.tab` is a tab-separated text file with three fields:
| Sample name | Soft-clipping value | Group ID |
|----------|----------|----------|
| Sample_A.fastq | 2 | Individual_1  |
| Sample_B.fastq | 2 | Individual_1  |
| Sample_C.fastq | 5 | Individual_1  |

The workflow uses this file to:
- Apply sample-specific soft clipping, and
- Merge aligned data belonging to the same individual (i.e., those sharing the same Group ID) after alignment.
You can see an example file in the `config/` directory of the repository.

# (3) Required taxonomic IDs
You need to provide two NCBI taxonomic IDs:
- Target species taxonomic ID — used for pangenome construction. The workflow will automatically download the necessary genomic data based on this ID.
- Outgroup taxonomic ID — used only for rooting the phylogenetic tree.
You can control how many genomes to download for pangenome building using the `--genomes` parameter.

If you already have your own curated dataset (FASTA + GenBank files), you can provide its path using `--trusted_data`.

Important requirements:
- FASTA and GenBank filenames must match exactly (aside from extensions).
- Required extensions: `.fasta` and `.gb`

# (4) Data and output directories
Before running the workflow, prepare two directories:
- A folder containing only your `FASTQ/FASTQ.gz` files. **Only store the data that will be included in the analysis.**
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

Whenever the workflow starts, `panchronos` prints a summary of all parameters used for the current run.
This allows you to quickly verify that your settings are being applied correctly.

# 6/ Output
As the workflow runs, the `--output` directory will begin to populate with results. Files are written as soon as they are generated and are automatically organized into structured subdirectories, making it easy to track each stage of the analysis.

| Folder | Description | Format |
| :--- | :--- | :---: | 
|**ALIGNMENT**| Aligned reads against the pangenome reference sequences | `.bam`, `.fastq` |
|**DOWNLOADED**| Input genomes used for pangenome construction. If `--trusted_data` is provided, this folder stores the user-supplied curated datasets | `.fasta`, `.gb` |
|**GENE_DATABASE**| Guide file for `Prokka`, containing clustered genes extracted from `GenBank` files prior to annotation | `.fasta` |
|**GENE_MSA**| Multiple sequence alignments for each individual gene | `.fasta` |
|**GENOTYPING**| Variant-calling results, including consensus sequences per sample per gene, and associated VCF files | `.fasta`, `.vcf*` |
|**MAPDAMAGE**| Output from DNA damage pattern assessment using `mapDamage` | Multiple files |
|**MATRIX**| Final gene presence/absence matrix | `.tab` |
|**PANGENOME**| Outputs from pangenome construction and extension. Includes the original Panaroo reference genome, extended/unextended versions, and a BLAST database | `.fasta`, BLAST database |
|**PLOTS**| Visualizations including presence/absence heatmaps and coverage vs. completeness plots | `.png` |
|**STATS**| Summary statistics, including gene-level coverage and completeness normalization per sample | `.tab`, `.txt` |
|**TREE**| Results from IQ-TREE, including the concatenated MSAs used for reconstruction | `.fasta`, `.treefile`, other outputs from IQ-TREE |

# 7/ Tools installed by `panchronos`
The following tools will be automatically installed:
- Entrez Direct
- Biopython v1.85
- fastANI v1.34
- CD-HIT v4.8.1
- Prokka v1.14.6
- Panaroo v1.5.2 (`panchronos` has been only tested with this version of Panaroo, and compatibility with newer releases is not guaranteed.)
- seqtk v1.5
- BWA v0.7.18
- bamUtil v1.0.15
- SAMtools v1.19
- mapDamage v2.2.3
- art_illumina v2.5.8
- bcftools v1.21
- MAFFT v7.525
- IQ-TREE v3.0.1

# 8/ Documentation
If you are unsure about the available workflow settings or parameters, you can display all options directly from the terminal using:

>`nextflow run main.nf --help`

For more extensive and detailed documentation, you can also refer to the materials available on: https://shigekinakagomelab.com/Software/.
