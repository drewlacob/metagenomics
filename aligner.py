DIFFERENCE_THRESHOLD = 2
READ_LENGTH = 50

#return the hamming distance between two strings
#defined as the number of mismatches between two strings of equal length
def find_hamming_distance(s1, s2):
    if len(s1) != len(s2): return READ_LENGTH
    hd = 0
    for i, c in enumerate(s1):
        if c != s2[i]:
            hd += 1
    return hd

#check if there is an insertion or deletion between a string and a reference string
def check_indel(string, reference):
    #if strings of unequal lengths return
    if len(string) != len(reference):
        return -3, -3

    #if the strings are the same return
    if string == reference:
        return -2, -2 

    #find the first difference between the two strings
    first_difference = -1
    for i, c in enumerate(string):
        if c != reference[i]:
            first_difference = i
            break

    #check if an insertion was found by first assuming there was and then disproving it if not
    insertion_found = True
    for i in range(first_difference, len(string)):
        if i + 1 < len(string):
            if string[i+1] != reference[i]:
                insertion_found = False
                break

    #if no insertion was found, now check for a deletion
    deletion_found = True
    if not insertion_found:
        c = 0
        for i in range(first_difference+1, len(reference)-1):
            if c > 4: break
            if i + 1 < len(reference):
                if reference[i] != string[i-1]:
                    deletion_found = False
                    break
            c += 1

    #return codes based on insertion, deletion, or subsitution found
    if insertion_found:
        return first_difference, 1
    elif deletion_found:
        return first_difference, 2
    else:
        return -1, -1

#align reads to the genome using the phone book and split hash method
#by splitting our reads of length 50 into 3 parts, we can guarrantee
#based on coverage in the reads that one section of the 50
#will not have errors, thus can be mapped using a hash map
def align_reads(reads, genome, GPB):
    read_to_position_map = dict()

    #for each read
    for read in reads:
        #find the length of each third of the read based on the read length
        partition_length = int(READ_LENGTH / 3)

        #partition the read into thirds
        first_third_read = read[0:partition_length]
        middle_third_read = read[partition_length : 2 * partition_length]
        last_third_read = read[2 * partition_length + 2 :]

        #create lists for the potential indices where each third may map based on the GPB
        potential_indexes_first_third = []
        potential_indexes_second_third = []
        potential_indexes_last_third = []

        #if first third mapped in the hash map
        if first_third_read in GPB:
            potential_indexes_first_third = GPB[first_third_read]
            #for each potential index
            for index in potential_indexes_first_third:
                #find the reference string for the genome
                genome_string = genome[index:index + READ_LENGTH]
                #check the hamming distance between the full read and the reference
                #if it is below the threshold, we have found a correct mapping of the full read
                if find_hamming_distance(read, genome_string) <= DIFFERENCE_THRESHOLD:
                    read_to_position_map[read] = index
        #else check the middle third of the read
        elif middle_third_read in GPB:
            potential_indexes_second_third = GPB[middle_third_read]
            #for each potential index
            for index in potential_indexes_second_third:
                #find the reference string for the genome
                genome_string = genome[index-partition_length:index + (2*partition_length)+2]
                #check the hamming distance between the full read and the reference
                #if it is below the threshold, we have found a correct mapping of the full read
                if find_hamming_distance(read, genome_string) <= DIFFERENCE_THRESHOLD:
                    read_to_position_map[read] = index - partition_length 
        #else check the last third of the read
        elif last_third_read in GPB:
            potential_indexes_last_third = GPB[last_third_read]
            #for each potential index
            for index in potential_indexes_last_third:
                #find the reference string for the genome
                genome_string = genome[index-(2*partition_length)-2:index+partition_length]
                #check the hamming distance between the full read and the reference
                #if it is below the threshold, we have found a correct mapping of the full read
                if find_hamming_distance(read, genome_string) <= DIFFERENCE_THRESHOLD: 
                    read_to_position_map[read] = index - (2*partition_length) - 2
        
    return read_to_position_map

#TODO: improve this function to identify indels, need to also edit the main alignment step to determine the indexes of potential indels
#indentify any indels in the list of potential indels given the current read map and the genome
def indentify_indels(read_to_position_map, genome, potential_indels):
    insertions_seen_before = set()
    deletions_seen_before = set()
    
    for read in potential_indels:
        genome_position = potential_indels[read]
        genome_string = genome[genome_position:genome_position+len(read)]
        indel_position, indel_type = check_indel(read, genome_string)
        if indel_position or indel_position == 0:
            if indel_type == 1:
                mutation = '>I' + str(genome_position+indel_position-1) + ' ' + read[indel_position]
                if mutation in insertions_seen_before:
                    read_to_position_map[read] = genome_position+indel_position
                else:
                    insertions_seen_before.add(mutation)
            elif indel_type == 2:
                mutation = '>D' + str(genome_position+indel_position-2) + ' ' +  genome_string[indel_position]
                if mutation in deletions_seen_before:
                    read_to_position_map[read] = genome_position+indel_position
                else:
                    deletions_seen_before.add(mutation)

    return read_to_position_map

def create_genome_phone_book(genome):
    GPB = dict()
    key_length = int(READ_LENGTH / 3)
    for i in range(0, len(genome)-key_length+1):
        genome_string = genome[i:i+key_length]
        if genome_string in GPB:
            GPB[genome_string].append(i)
        else:
            GPB[genome_string] = [i]
    return GPB

def call_aligner(genome_file, reads):
    #load the reference genome
    ref_file = open(genome_file, 'r')
    genome = ref_file.read().splitlines()
    genome = ''.join(genome[1:])

    #create a phone book of kmers for the genome
    GPB = create_genome_phone_book(genome)

    # align the reads and obtain a mapping from each read to the position it starts in the genome
    read_to_position = align_reads(reads, genome, GPB)

    #TODO: identify indels during mapping and account for that when aligning using indentify_indels function

    return read_to_position
