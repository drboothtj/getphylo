'''
Screen fasta files and blast databases for singletons and extract sequences

Functions:
'''
import os
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
        target_loci = []
        loci_to_find, loci_min_length, loci_max_length = parser.get_seed_thresholds()
        for locus in unique_loci:
            if loci_to_find < 0 or loci < loci_to_find:
                sequence = get_locus(seed_fasta, locus)
                if loci_max_length > len(sequence) > loci_min_length:
                    target_loci.append(locus)
                    loci_fasta.append(">" + locus)
                    loci_fasta.append(sequence)
                    loci += 1
            else:
                break
        console.print_to_system(str(loci) + " loci found. Writing to tsv/loci.fasta")
        io.write_to_file("tsv/seed_loci.fasta", loci_fasta)
        io.write_to_file("tsv/seed_loci.txt", target_loci)
        return target_loci


def get_target_proteins(seed):
    '''main routine for screen.py'''
    checkpoint = parser.get_checkpoint()
    if checkpoint < 3:
        target_loci = get_singletons_from_seed(seed)
