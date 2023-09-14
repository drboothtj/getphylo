'''
Build trees from a directory containing .fasta alignments

Functions:
    def make_trees(output: str, build_all: bool) -> None
'''
import os
import glob
import logging
from getphylo.utils import io
from getphylo.ext import fasttree

def make_trees(output: str, build_all: bool, cpus: int) -> None:
    '''Main routine for trees.
        Arguments:
            output: path to the output directory
        Returns:
            None
    '''
    tree_directory = os.path.join(output, 'trees')
    io.make_folder(tree_directory)
    logging.info("Building trees...")
    if build_all is True:
        args_list = []
        for filename in glob.glob(os.path.join(output, 'aligned_fasta/*.fasta')):
            outfile = os.path.join(
                tree_directory, os.path.basename(io.change_extension(filename, "tree"))
                )
            args_list.append([filename, outfile])
            io.run_in_parallel(fasttree.run_fasttree, args_list, cpus)
    else:
        filename = os.path.join(output, 'aligned_fasta/combined_alignment.fasta')
        outfile = os.path.join(tree_directory, 'combined_alignment.tree')
        fasttree.run_fasttree(filename, outfile)
    logging.info("CHECKPOINT: TREES_BUILT")
