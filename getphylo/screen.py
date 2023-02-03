'''
Screen fasta files and blast databases for singletons and extract sequences

Functions:
'''
import os
import glob
from collections import Counter
from random import shuffle
from getphylo import console, diamond, io

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

def get_seed_paths(seed, output):
    '''
    ensures a seed is defined and provides file names for with
    .fasta, .dmnd and .tsv extensions
    '''
    seed_fasta = io.change_extension(seed, "fasta")
    seed_fasta = f'{output}/fasta/{seed_fasta}'
    seed_dmnd = f'{output}/dmnd/{io.change_extension(seed, "dmnd")}'
    seed_tsv = f'{output}/tsv/{io.change_extension(seed, "tsv")}'
    return seed_fasta, seed_dmnd, seed_tsv

def get_singletons_from_seed(seed, output, thresholds):
    #too many local variables
    '''uses diamond to identify singletons in the seed genome'''
    io.make_folder(f'{output}/tsv')
    console.print_to_system("Identifying singletons in seed genome...")
    seed_fasta, seed_dmnd, seed_tsv = get_seed_paths(seed, output)
    diamond.run_diamond_search(seed_fasta, seed_dmnd, seed_tsv)
    unique_loci = get_unique_hits_from_tsv(seed_tsv)
    console.print_to_system("Found " + str(len(unique_loci)) + (" unique loci."))
    shuffle(unique_loci)
    loci = 0
    loci_fasta = []
    candidate_loci = []
    loci_to_find, loci_min_length, loci_max_length, _, _ = thresholds
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
    console.print_to_system(f'{loci} loci found. Writing to {output}/tsv/loci.fasta')
    io.write_to_file(f'{output}/tsv/seed_loci.fasta', loci_fasta)
    io.write_to_file(f'{output}/tsv/seed_loci.txt', candidate_loci)
    return candidate_loci

def get_loci_from_file(file):
    '''gets a list of loci from the provided text file'''
    loci = io.read_file(file)
    loci = [locus[:len(locus) -1] for locus in loci] #remove newline (fix, unsafe!)
    return loci

def search_candidates(output):
    '''searches for candidates in othe genomes'''
    io.make_folder(f'{output}/tsvs')
    console.print_to_system("Screening candidate loci against other genomes...")
    diamond_databases = glob.glob(f'{output}/dmnd/*.dmnd')
    for database in diamond_databases:
        tsv_name = io.change_extension(database.split("/")[2], 'tsv')
        tsv_name = f'{output}/tsvs/{tsv_name}'
        diamond.run_diamond_search(f'{output}/tsv/seed_loci.fasta', database, tsv_name)
    #allow fiddling with dmnd options

def score_locus(locus, files):
    '''scores the locus'''
    pa_data = []
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
        pa_data.append(counter)
    return presence_counter, unique_flag, pa_data

def threshold_loci(target_loci, thresholds, output):
    #Too many local variables
    '''Score loci for singleton status and presence in dataset'''
    _, _, _, presence_threshold, minimum_loci = thresholds
    if os.path.exists(f'{output}/final_loci.txt'):
        console.print_to_system('ALERT: final_loci.txt already exists. Exiting!')
        exit()
        #either all return statements in a function should return an expression
        #or none of them should.
    else:
        console.print_to_system("Thresholding targets...")
        final_loci = []
        pa_table = []
        files = glob.glob(f'{output}/tsvs/*.tsv')
        pa_table.append([file.split('/')[2].split('.')[0] for file in files])
        thresholding_data = ["locus;" + "presence;" + "unique"]
        for locus in target_loci:
            presence, unique, pa_data = score_locus(locus, files)
            pa_table.append(pa_data)
            presence_percent = (presence / len(glob.glob(f'{output}/tsvs/*.tsv'))) * 100
            thresholding_string = str(locus) + ";" + str(presence_percent) + ";" + str(unique)
            thresholding_data.append(thresholding_string)
            if presence_percent >= presence_threshold and unique:
                final_loci.append(locus)
        io.write_to_file(f'{output}/thresholding_data.txt', thresholding_data)
        write_pa_table(pa_table, target_loci, output)
        number_of_loci = len(final_loci)
        console.print_to_system(str(number_of_loci) + " loci selected for MLST...")
        if number_of_loci < minimum_loci:
            console.print_to_system("Number of loci below defined threshold. Exiting...")
            exit()
            #either all return statements in a function should return an expression
            #or none of them should.
        io.write_to_file(f'{output}/final_loci.txt', final_loci)
        return final_loci

def write_pa_table(pa_table, loci, output):
    '''transposes and writes the presence absence table to the specified folder'''
    new_table = []
    header = ['strain']
    header.extend(loci)
    new_table.append(";".join(header))
    transposed_data = list(zip(*pa_table))
    for data in transposed_data:
        data_string = ';'.join([str(datum) for datum in data])
        new_table.append(data_string)
    io.write_to_file(f'{output}/presence_absence_table.csv', new_table)
    #check_table

def get_target_proteins(checkpoint, output, seed, thresholds):
    '''main routine for screen.py'''
    if checkpoint < 3:
        console.print_to_system("CHECKPOINT 2: Identifying singletons...")
        candidate_loci = get_singletons_from_seed(seed, output, thresholds)
    if checkpoint > 2:
        #since candidate loci is created before checkpoint 3
        #all later checkpoints require this list to be read in
        candidate_loci = get_loci_from_file(f'{output}/tsv/seed_loci.txt')
    if checkpoint < 4:
        console.print_to_system("CHECKPOINT 3: Searching for singletons in other genomes...")
        search_candidates(output)
    if checkpoint < 5:
        console.print_to_system("CHECKPOINT 4: Applying loci thresholds...")
        final_loci = threshold_loci(candidate_loci, thresholds, output)
        return final_loci
