#!/usr/bin/env nextflow

params.data = ""
reads = Channel.of(params.data)

params.output = ""
resultsDir = Channel.of(params.output)

params.completeness = 50
geneCompleteness = Channel.of(params.completeness)

params.coverageDown = 0.5
normalizedCoverageDown = Channel.of(params.coverageDown)

params.coverageUp = 3
normalizedCoverageUp = Channel.of(params.coverageUp)

params.threads = 10
threadsGlobal = Channel.of(params.threads)

params.tax_id = 0
taxID = Channel.of(params.tax_id)

params.genomes = 100
downloadGenomes = Channel.of(params.genomes)

params.clustering = 0.95
cdHitCluster = Channel.of(params.clustering)

params.core = 0.95
pangenomeThreshold = Channel.of(params.core)

params.clean = "strict"
pangenomeMode = Channel.of(params.clean)

params.config = ""
configFile = Channel.of(params.config)

params.outgroup = ""
outTax = Channel.of(params.outgroup)

params.trustedGenomes = false
def trustedDataChannel = params.trustedGenomes ? Channel.fromPath(params.trustedGenomes) : Channel.empty()

params.missing = 0.01
missingProb = Channel.of(params.missing)

params.gaps = 2
gapFraction = Channel.of(params.gaps)

params.seed = 16500
seedAlignment = Channel.of(params.seed)

params.mapq = 30
mappingQuality = Channel.of(params.mapq)

params.minlength = 34
minReadLength = Channel.of(params.minlength)

params.maxlength = 300
maxReadLength = Channel.of(params.maxlength)

params.genotyper = "gatk"
genotypeMethod = Channel.of(params.genotyper)


params.help = false


// Enable DSL2
nextflow.enable.dsl=2

def printHelp {
	println "\n\033[1;31mSYNOPSIS\033[0m"

	println "\n\033[1;33mUSAGE\033[0m"
	println "\nnextflow run aBPA.nf --data <PATH> --output <PATH> --tax_id <INT> --config <PATH/FILE> [..OPTIONS..]"
	
	println "\n\033[1;33mMANDATORY\033[0m"
	println "  --data <PATH>		Set data file PATH"
	println "  --output <PATH>		Set output directory PATH"
	println "  --tax_id <INT>		Set taxonomical ID value <INT>"
	println "  --config <PATH>		Set config file PATH"
	
	println "\n\033[1;33mOPTIONS\033[0m"
	println "  --threads <INT>		Set number of threads (default: 10)"
	println "  --completeness <INT/FLOAT>	Set gene completeness/breadth of coverage threshold (default: 50)"
	println "  --coverage <INT/FLOAT>	Set mean depth of coverage threshold (default: 0.5)"
	println "  --genomes <INT>		Set number of genomes to download (default: 100)"
	println "  --clustering <INT/FLOAT>	Set clustering threshold (default 0.95)"
	println "  --core-threshold <FLOAT>	Set core genome threshold (default: 0.01)"
	println "  --clean-mode <STRING>	Set pangenome mode (default: strict)"
	println "  --help			Print help page and exit"
	
	println "\n\033[1;31mDESCRIPTION\033[0m"
	println "\n\033[1;33m--data <PATH>\033[0m"
	println "Please specify the full PATH of your data. Example: /home/user/mydata/data"
	
	println "\n\033[1;33m--output <PATH>\033[0m"
	println "Please specify the full PATH of your output folder. You need to make the folder first before running the program."

	println "\n\033[1;33m--tax_id <INT>\033[0m"
	println "Please specify the taxonomical ID for your bacteria. It should be a discrete and unique number."
	
	println "\n\033[1;33m--config <PATH>\033[0m"
	println "\nPlease set file PATH of your config.tab file. Example: /home/user/me/aBPA/config/config.tab"
	println "config.tab file should contain 3 fields separated by tab. First field should have the sample name, second field softclipping value <INT> and third field group ID."
	println "\nExample: \n	SAMPLE1	5	NONUDG\n	SAMPLE2	2	UDG"

	println "\n\033[1;33m--threads <INT>\033[0m"
	println "Set amount of threads to be used globally."

	println "\n\033[1;33m--completeness <INT>\033[0m"
	println "Set gene breadth of coverage threshold as percentage. Genes that have a value less than <INT> will be considered absent."

	println "\n\033[1;33m--coverage <INT/FLOAT>\033[0m"
	println "Set gene normalized coverage threshold. Currently aBPA is using the simplest statistic for normalization: (Gene mean depth/Global mean depth)."

        println "\n\033[1;33m--genomes <INT/FLOAT>\033[0m"
	println "Set amount of FASTA/GENBANK files to be downloaded. Bear in mind disk space."

        println "\n\033[1;33m--clustering <INT/FLOAT>\033[0m"
	println "Set clustering threshold <INT/FLOAT> for FASTA database.\nA value of 0.9 means any group of sequences with identity values equal or bigger than 0.9 will be clustered together and a consensus representative sequence will be produced."

        println "\n\033[1;33m--core-threshold <INT/FLOAT>\033[0m"
	println "Set threshold for core genome building. Similarly as clustering flag but during pangenome step."

        println "\n\033[1;33m--clean-mode <INT/FLOAT>\033[0m"
	println "Set behaviour of pangenome building. Visit Panaroo documentation to know more about this.\n\n"

    exit 0
}


if (params.help) {
    printHelp
}


workflow {


}
