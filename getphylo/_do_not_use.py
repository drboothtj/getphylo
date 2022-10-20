'''
Main routine for fastMLST

Functions:

'''
import os
import glob
from collections import Counter
from random import shuffle
from fastMLST import console, diamond, extract, fasttree, io, muscle, parser

def concatenate_alignments():
    concatenated_fasta = []
    for file in glob.glob("./aligned_fasta/*.fasta"):
        lines = io.read_file(file)
        concatenated_fasta.append(lines)
    io.write_to_file("./aligned_fasta/concatenated_alignment.fasta", lines)





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
    
def parse_args():
    '''get the arguments from the console via the parser'''
    arg_parser = parser.get_parser()
    args = arg_parser.parse_args()
    return args




def set_seed(seed):
    if seed is None:
        gbks = glob.glob("*.gbk")
        shuffle(gbks)
        seed = io.change_extension(gbks[0], "fasta") 
    else:
        seed = io.change_extension(seed, "fasta")
    
    console.print_to_system("The seed genome is " + seed + "...")
    return seed



def main():
    ###PRINT START###
    console.print_to_system("Running fastMLST version 0.0.1")
    ###PARSE ARGS### FIX ALL THIS!!!
    args = parse_args()
    ##make arguments get_args()
    seed = args.seed
    loci_to_find = 20
    loci_min_length = 500
    #add option to use all avalibale loci 
    #add extra dmnd options
    presence_threshold = 100
    minimum_loci = 3
    genbank_path = args.gbk
    #add max threshold

    ###WORKFLOW###
    seed = set_seed(seed)
    extract.extract_data(list_of_genbanks)
    target_loci = identify_target_genes(seed, loci_to_find, loci_min_length) #get singletons and write to files
    screen_target_loci(target_loci)
    
    
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
