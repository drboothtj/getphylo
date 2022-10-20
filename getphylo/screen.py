'''
Screen fasta files and blast databases for singletons and extract sequences

Functions:
'''
import os
import glob
from collections import Counter
from random import shuffle
from getphylo import console, diamond, io, parser

def get_locus(file, locus):
    '''returns a sequence from a fasta file with the provided locus name'''
    fasta = io.read_file(file)
    line_number = 0
    for line in fasta:
        line_number += 1
        if locus in line:
            sequence = fasta[line_number]
    return sequence

def get_unique_hits_from_tsv(file):
    '''finds unique hits from the provided dmnd result file'''
    hits = []
    unique_hits = []
    lines = io.read_tsv(file)
    for line in lines:
        hits.append(line[0])
    counter = Counter(hits)
    for hit in counter:
        if counter[hit] == 1:
            unique_hits.append(hit)
    return unique_hits

def get_seed_paths(seed):
    '''
    ensures a seed is defined and provides file names for with
    .fasta, .dmnd and .tsv extensions
    '''
    seed_fasta = io.change_extension(seed, "fasta")
    seed_fasta = "fasta/" + seed_fasta
    seed_dmnd = io.change_extension(seed, "dmnd")
    seed_dmnd = "dmnd/" + seed_dmnd
    seed_tsv = io.change_extension(seed, "tsv")
    seed_tsv = "tsv/" + seed_tsv
    return seed_fasta, seed_dmnd, seed_tsv

def get_singletons_from_seed(seed):
    '''uses diamond to identify singletons in the seed genome'''
    try:
        os.mkdir('tsv')
    except OSError as error:
        if error.errno == 17:
            console.print_to_system("ALERT: The directory './tsv/' already exists. Exiting!")
            exit()
        else:
            raise
    else:
        console.print_to_system("Identifying singletons in seed genome...")
        seed_fasta, seed_dmnd, seed_tsv = get_seed_paths(seed)
        diamond.run_diamond_search(seed_fasta, seed_dmnd, seed_tsv)
        unique_loci = get_unique_hits_from_tsv(seed_tsv)
        console.print_to_system("Found " + str(len(unique_loci)) + (" unique loci."))
        shuffle(unique_loci)
        loci = 0
        loci_fasta = []
        candidate_loci = []
        loci_to_find, loci_min_length, loci_max_length = parser.get_seed_thresholds()
        for locus in unique_loci:
            if loci_to_find < 0 or loci < loci_to_find:
                sequence = get_locus(seed_fasta, locus)
                if loci_max_length > len(sequence) > loci_min_length:
                    candidate_loci.append(locus)
                    loci_fasta.append(">" + locus)
                    loci_fasta.append(sequence)
                    loci += 1
            else:
                break
        console.print_to_system(str(loci) + " loci found. Writing to tsv/loci.fasta")
        io.write_to_file("tsv/seed_loci.fasta", loci_fasta)
        io.write_to_file("tsv/seed_loci.txt", candidate_loci)
        return candidate_loci

def get_loci_from_file(file):
    '''gets a list of loci from the provided text file'''
    loci = io.read_file(file)
    loci = [locus[:len(locus) -1] for locus in loci] #remove newline (fix, unsafe!)
    return loci

def search_candidates():
    '''searches for candidates in othe genomes'''
    try:
        os.mkdir('tsvs')
    except OSError as error:
        if error.errno == 17:
            console.print_to_system("ALERT: The directory './tsvs/' already exists. Exiting!")
            exit()
        else:
            raise
    else:
        console.print_to_system("Screening candidate loci against other genomes...")
        diamond_databases = glob.glob('dmnd/*.dmnd')
        for database in diamond_databases:
            output = database.split("/")[1]
            output = "tsvs/" + io.change_extension(output, "tsv")
            diamond.run_diamond_search("tsv/seed_loci.fasta", database, output)
        #allow fiddling with dmnd options

def score_locus(locus, files):
    '''scores the locus'''
    presence_counter = 0
    unique_flag = True
    for file in files:
        counter = 0
        lines = io.read_file(file)
        for line in lines:
            if locus in line:
                counter += 1
        if counter > 0:
            presence_counter += 1
        if counter > 1:
            unique_flag = False
    return presence_counter, unique_flag

def threshold_loci(target_loci):
    '''Score loci for singleton status and presence in dataset'''
    presence_threshold, minimum_loci = parser.get_loci_tresholds()
    if os.path.exists('final_loci.txt'):
        console.print_to_system('ALERT: final_loci.txt already exists. Exiting!')
        exit() 
        #either all return statements in a function should return an expression, or none of them should.
    else:
        console.print_to_system("Thresholding targets...")
        final_loci = []
        thresholding_data = ["locus;" + "presence;" + "unique"]
        for locus in target_loci:
            presence, unique = score_locus(locus, glob.glob("tsvs/*.tsv"))
            presence_percent = (presence / len(glob.glob("tsvs/*.tsv"))) * 100
            thresholding_string = str(locus) + ";" + str(presence_percent) + ";" + str(unique)
            thresholding_data.append(thresholding_string)
            if presence_percent >= presence_threshold and unique is True:
                final_loci.append(locus)
        io.write_to_file("thresholding_data.txt", thresholding_data)
        number_of_loci = len(final_loci)
        console.print_to_system(str(number_of_loci) + " loci selected for MLST...")
        if number_of_loci < minimum_loci:
            console.print_to_system("Number of loci below defined threshold. Exiting...")
            exit()
            #either all return statements in a function should return an expression, or none of them should.
        io.write_to_file("final_loci.txt", final_loci)
        return final_loci

def get_target_proteins(checkpoint, seed):
    '''main routine for screen.py'''
    if checkpoint < 3:
        console.print_to_system("CHECKPOINT 2: Identifying singletons...")
        candidate_loci = get_singletons_from_seed(seed)
    if checkpoint > 2:
        #since candidate loci is created before checkpoint 3
        #all later checkpoints require this list to be read in
        candidate_loci = get_loci_from_file('tsv/seed_loci.txt')
    if checkpoint < 4:
        console.print_to_system("CHECKPOINT 3: Searching for singletons in other genomes...")
        search_candidates()
    if checkpoint < 5:
        console.print_to_system("CHECKPOINT 4: Applying loci thresholds...")
        final_loci = threshold_loci(candidate_loci)
        return final_loci
