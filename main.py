import meta

GENOME_FILE_BASE = './data/genome_'
READS_FILE = './data/reads.fasta'
NUM_GENOMES = 1000 

#load the reads
reads = meta.load_reads(READS_FILE)
ALIGNED_READ_THRESHOLD = len(reads) / 25

#check which genomes are present in the reads
genomes_present = meta.determine_genomes_present(NUM_GENOMES, GENOME_FILE_BASE, reads, ALIGNED_READ_THRESHOLD)
print(genomes_present)
print(len(genomes_present), 'genomes present with threshold of', ALIGNED_READ_THRESHOLD)

#align the reads to their best match based on the present genomes
reads_to_genome = meta.find_read_to_genome_matches(reads, genomes_present, GENOME_FILE_BASE)
print('Reads aligned:', len(reads_to_genome))

#map the read numbers to their genome numbers, randomly choosing a present genome for any unmapped reads
read_nums_to_genome_num = meta.find_read_num_to_genome_greedy(reads, reads_to_genome, genomes_present)
print('Reads aligned after greedy random choices:', len(read_nums_to_genome_num))

#create the result file, prints to ans.txt
meta.make_result_file(read_nums_to_genome_num)