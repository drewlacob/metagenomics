# Metagenomics Project

### This repository contains the code for a metagenomics project. The goal is to take a reads file and map it to the genomes provided. The difficulty lies in that the reads could come from any of the 1000 genomes given. Thus one must first solve the read mapping problem for a single genome and then apply that solution across all of the genomes efficiently.

## Inputs

- reads.fasta - reads file with 20,000 reads, each numbered
- genome\_###.fasta - 1,000 reference genomes, where ### is the number of the genome

## Outputs

- ans.txt - answer file that maps each read to the genome that it belongs to in the format '>read*## Genome_Number*###'. Note that the desired output is to have a genome for every read so all reads must be mapped to a genome (answer file should be 20,000 reads).

## Usage

Can be run using `python3 main.py` after changing desired constants for use case and data. Sample reads file and genomes files are provided in the data subfolder.

## Results

- Successfully identifies present genomes based on custom threshold.
- Performs metagenomics analysis to map each read to a genome, account for SNPs in the form of substitutions.
- Successfully mapped 76% of reads to the correct genome in case of 20,000 reads with 1,000 possible genomes each of length 10,000 in under a minute.

## Modules

### `main.py`

The main module is the entry point for the script and solves the metagenomics problem using the meta module. The basic steps are to:

- Load the reads and genomes
- Determine which genomes are actually present in the reads since only a small subset of the 1000 total will be present
- Align the reads to their best match based on the present genomes
- Map the read numbers to their genome numbers, randomly choosing a present genome for any unmapped reads

### `meta.py`

The metagenomics module performs the metagenomics analysis for the reads by calling the aligner. The following functions are defined:

- load_reads - loads the reads
- determine_genomes_present - uses the reads and the aligner submodule to determine which genomes are present in the reads
- find_read_to_genome_matches - match each read to each present genome
- find_read_num_to_genome_greedy - map the read number back to the genome, guessing randomly from the present genomes when a read was not aligned to be in any genome
- make_result_file - creates the output ans.txt file in the desired format

### `aligner.py`

The aligner module performs the actual aligning of a the reads to a single genome. It uses a hashing approach where reads are partitioned and mapped to a 'phone book' in order to significantly reduce the time taken to align the reads. The basic steps are to:

- Load the reads and the genome
- Create hash table of kmers of thirds of the length of the reads of the genome known as the genome phone book
- For each read:
  - split the read into thirds
  - map each third to the phone book
  - one of the the three thirds will match
  - for the match, use the full read and full reference genome slice to check the hamming distance and determine if the read actually maps to that location
- Check for insertions or deletions for each read

### Todo / Future Improvements

- Account for indels during alignment by finding a more efficient way of solving for them. Use this to further increase accuracy.
- Apply further optimizations such as bloom filters, minimizers, or minimum perfect hashing to enable scaling for even larger datasets.
- Implement a consensus algorithm in cases where there are duplicated regions between genomes so a read may
  map to more than one genome.
