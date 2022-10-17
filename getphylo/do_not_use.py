'''
Main routine for fastMLST

Functions:

'''
import os
import glob
from collections import Counter
from random import shuffle
from fastMLST import console, diamond, extract, fasttree, io, muscle, parser


#clean into more modules!

def concatenate_alignments():
    concatenated_fasta = []
    for file in glob.glob("./aligned_fasta/*.fasta"):
        lines = io.read_file(file)
        concatenated_fasta.append(lines)
    io.write_to_file("./aligned_fasta/concatenated_alignment.fasta", lines)

def do_alignments():
    try:
        os.mkdir('aligned_fasta')
    except OSError as e:
        if e.errno == 17:
            console.print_to_system("The directory './aligned_fasta/' already exists. Skipping alignments...")
        else:
            raise
    else:    
        console.print_to_system("Building alignments...")
        for file in glob.glob("./unaligned_fasta/*.fasta"):
            outfile = "aligned_fasta/" + file.split("/")[2]
            muscle.run_muscle(file, outfile)

def get_loci_from_file(file):
    loci = io.read_file(file)
    loci = [locus[:len(locus) -2] for locus in loci] #remove newline
    return loci   

def get_fasta_for_alignment(loci_list):
    try:
        os.mkdir('unaligned_fasta')
    except OSError as e:
        if e.errno == 17:
            console.print_to_system("The directory './unaligned_fasta/' already exists. Skipping fasta extraction...")
        else:
            raise
    else:
        console.print_to_system("Extracting targets for alignment...")
        if loci_list is None:
            loci_list = get_loci_from_file('final_loci.txt') #catch error if missing

        files = glob.glob("fasta/*.fasta")

        for locus in loci_list:
            write_lines = []
            for file in files:
                tsv_name = "tsvs/" + io.change_extension(file.split("/")[1], "tsv")
                tsv = io.read_tsv(tsv_name)
                for line in tsv:
                    if locus in line[0]:
                        sequence = get_locus(file, line[1])
                        write_lines.append(">" + file + line[1])
                        write_lines.append(sequence)
                        break
            io.write_to_file("unaligned_fasta/" + locus + ".fasta", write_lines)

def get_locus(file, locus):
    fasta = io.read_file(file)
    line_number = 0
    for line in fasta:
        line_number += 1
        if locus in line:
            sequence = fasta[line_number]
    return sequence

def get_unique_hits_from_tsv(file):
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

def identify_target_genes(seed, loci_to_find, loci_min_length):
    try: 
        os.mkdir('tsv')
    except OSError as e:
        if e.errno == 17:
            console.print_to_system("The directory './tsv/' already exists. Skipping singleton identification...")
        else:
            raise
    else:
        console.print_to_system("Identifying target genes from seed genome...")
        seed_fasta = io.change_extension(seed, "fasta")
        seed_fasta = "fasta/" + seed_fasta
        seed_dmnd = io.change_extension(seed, "dmnd")
        seed_dmnd = "dmnd/" + seed_dmnd
        seed_tsv = io.change_extension(seed, "tsv")
        seed_tsv = "tsv/" + seed_tsv
        diamond.run_diamond_search(seed_fasta, seed_dmnd, seed_tsv)
        unique_loci = get_unique_hits_from_tsv(seed_tsv)
        console.print_to_system("Found " + str(len(unique_loci)) + (" unique loci."))
        shuffle(unique_loci)

        loci = 0
        loci_fasta = []
        target_loci = []
        for locus in unique_loci:
            if loci < loci_to_find:
                sequence = get_locus(seed_fasta, locus)
                if len(sequence) > loci_min_length:
                    target_loci.append(locus)
                    loci_fasta.append(">" + locus)
                    loci_fasta.append(sequence)
                    loci += 1
            else:
                break
        console.print_to_system(str(loci) + " loci found. Writing to loci.fasta")
        io.write_to_file("seed_loci.fasta", loci_fasta) ###WIPE BEOFRE WRITING
        io.write_to_file("seed_loci.txt", target_loci) ###WIPE BEFORE WRITING
        
        return target_loci

def make_trees():
    try:
        os.mkdir('trees')
    except OSError as e:
        if e.errno == 17:
            console.print_to_system("The directory './trees/' already exists. Skipping tree building...")
        else:
            raise
    else:
        console.print_to_system("Building trees...")
        for file in glob.glob("aligned_fasta/*.fasta"):
            outfile = "trees/" + io.change_extension(file, "tree").split("/")[1]
            fasttree.run_fasttree(file, outfile)
    


def screen_target_loci(target_loci):
    try: os.mkdir('tsvs')
    except OSError as e:
        if e.errno == 17:
            console.print_to_system("The directory './tsvs/' already exists. Skipping locus screening...")
        else:
            raise
    else:
        console.print_to_system("Screening target loci against other genomes...")
        if target_loci is None:
            target_loci = get_loci_from_file('seed_loci.txt') #catch error if missing
        
        diamond_databases = glob.glob('dmnd/*.dmnd')

        for database in diamond_databases:
            output = database.split("/")[1]
            output = "tsvs/" + io.change_extension(output, "tsv")
            diamond.run_diamond_search("seed_loci.fasta", database, output)
        #unique_hits.append(get_unique_hits_from_tsv(io.change_extension(database, "tsv")))
        #allow fiddling with dmnd options

def score_locus(locus, files):
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



def threshold_loci(target_loci, presence_threshold, minimum_loci):
    if os.path.exists('final_loci.txt'):
        console.print_to_system('final_loci.txt already exists. Skipping thresholding.')
    else:
        if target_loci is None:
            target_loci = get_loci_from_file('seed_loci.txt') #catch error if missing
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
        
        #add a flag to write thesholding results...
        io.write_to_file("thresholding_data.txt", thresholding_data)
        
        ###new function!###
        number_of_loci = len(final_loci)
        console.print_to_system(str(number_of_loci) + " loci selected for MLST...")
        if number_of_loci < minimum_loci:
            console.print_to_system("Number of loci below defined threshold. Exiting...")
            exit()

        io.write_to_file("final_loci.txt", final_loci)
        return final_loci

def main():
 
    ###WORKFLOW###
    seed = set_seed(seed)
    extract.extract_data(list_of_genbanks)
    target_loci = identify_target_genes(seed, loci_to_find, loci_min_length) #get singletons and write to files
    screen_target_loci(target_loci)
    final_loci = threshold_loci(target_loci, presence_threshold, minimum_loci)
    get_fasta_for_alignment(final_loci)
    do_alignments()
    #trim_alignments()
    concatenate_alignments()
    make_trees()

    ###END PRINT###
    console.print_to_system("Done!")

#concatenate alignments
#check for tree congruence?

#clean into modules
#run pylint
#add doc strings
#check comments
#clean output -> use folders

if __name__ == "__main__":
    main()
