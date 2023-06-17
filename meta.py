import aligner
import random
import time

#load the reads and return them as a list
def load_reads(reads_file):
    reads_f = open(reads_file, 'r')
    reads = reads_f.read().splitlines()
    reads = [line for i, line in enumerate(reads) if i%2 == 1]
    return reads

#determine which of the possible genomes are actually present in the sample reads and return them as a list
def determine_genomes_present(num_genomes, genome_file_base, reads, aligned_read_threshold = 0):
    genomes_present = []
    print('Determining which genomes are present...')
    start_time = time.time()

    #for each genome
    for i in range(num_genomes):
        #call the aligner on the reads for the current genome
        genomeFile = genome_file_base + str(i) + '.fasta'
        alignedReads = aligner.call_aligner(genomeFile, reads)

        #determine if the current genome was present in the reads based on aligned read length
        if len(alignedReads) > aligned_read_threshold:
            genomes_present.append(i)

        #progess report every 100 genomes searched
        if (i + 1) % 100 == 0:
            elapsed_time = time.time() - start_time
            iterations_remaining = num_genomes - (i + 1)
            estimated_time_remaining = (elapsed_time / (i + 1)) * iterations_remaining
            print(f"Iteration: {i + 1}, Elapsed Time: {elapsed_time:.2f}s, Estimated Time Remaining: {estimated_time_remaining:.2f}s")

    elapsed_time = time.time() - start_time
    print(f"Finished determing present genomes in {elapsed_time:.2f}s")
    return genomes_present

#match each read to each present genome
def find_read_to_genome_matches(reads, genomes_present, genome_file_base):
    read_to_genome = {}
    print('Finding read to present genome matches...')
    start_time = time.time()

    #for each present genome
    for i in genomes_present:
        #align the reads to the genome
        genomeFile = genome_file_base + str(i) + '.fasta'
        alignedReads = aligner.call_aligner(genomeFile, reads)
        
        #if the read has not already been mapped to a genome, map it
        for read in alignedReads:
            if read not in read_to_genome:
                read_to_genome[read] = i

    elapsed_time = time.time() - start_time
    print(f"Finished finding read to present genome matches in {elapsed_time:.2f}s")
    return read_to_genome

#map the read number back to the genome, guessing randomly from the present genomes when a read was not aligned to be in any genome
def find_read_num_to_genome_greedy(reads, reads_to_genome, genomes_present):
    read_num_to_genome_number = {}
    #for each read
    for i, read in enumerate(reads):
        #if we have mapped the read to a genome, then take that genome number
        if read in reads_to_genome: 
            read_num_to_genome_number[i] = reads_to_genome[read]
        else: #else choose a random genome to map it to from the genome list
            random_present_genome = random.choice(genomes_present)
            read_num_to_genome_number[i] = random_present_genome

    return read_num_to_genome_number

def make_result_file(read_number_to_genome_number):
    #loop through all the reads and print them to an output file with their matched genome
    ansF = open('./ans.txt', 'w')
    for i in range(len(read_number_to_genome_number)):
        ansF.write('>read_' + str(i) + ' Genome_Number_' + str(read_number_to_genome_number[i]) + '\n')
    ansF.close()
    