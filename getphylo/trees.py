'''
Build trees from a directory containing .fasta alignments

Functions:
    def make_trees(output: str) -> None
'''
import os
import glob
import logging
from getphylo.utils import io
from getphylo.ext import fasttree

def make_trees(output: str) -> None:
    '''Main routine for trees.
        Arguments:
            output: path to the output directory
        Returns:
            None
    '''
    tree_directory = os.path.join(output, 'trees')
    io.make_folder(tree_directory)
    logging.info("Building trees...")
    for filename in glob.glob(os.path.join(output, 'aligned_fasta/*.fasta')):
        outfile = os.path.join(
            tree_directory, os.path.basename(io.change_extension(filename, "tree"))
            )
        fasttree.run_fasttree(filename, outfile)
