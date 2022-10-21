'''
Build trees from a directory containing .fasta alignments

Functions:


'''
import os
import glob
from getphylo import console, fasttree, io

def make_trees():
    '''main routine for trees'''
    try:
        os.mkdir('trees')
    except OSError as e:
        if e.errno == 17:
            console.print_to_system("The directory './trees/' already exists. Exiting")
            exit()
        else:
            raise
    else:
        console.print_to_system("Building trees...")
        for filename in glob.glob("aligned_fasta/*.fasta"):
            outfile = "trees/" + io.change_extension(filename, "tree").split("/")[1]
            fasttree.run_fasttree(filename, outfile)


