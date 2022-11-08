'''
Build trees from a directory containing .fasta alignments

Functions:


'''
import os
import glob
from getphylo import console, fasttree, io

def make_trees(output):
    '''main routine for trees'''
    io.make_folder(f'{output}/trees')
    console.print_to_system("Building trees...")
    for filename in glob.glob(f'{output}/aligned_fasta/*.fasta'):
        outfile = f'{output}/trees/{io.change_extension(filename, "tree").split("/")[2]}'
        fasttree.run_fasttree(filename, outfile)
