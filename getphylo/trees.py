'''
Build trees from a directory containing .fasta alignments

Functions:
    build_all_trees(files: List, cpus: int, method: str, output: str) -> None
    make_trees(output: str, build_all: bool, method: str, cpus: int) -> None
'''
import os
import glob
import logging
from getphylo.utils import io
from getphylo.ext import fasttree, iqtree
from getphylo.utils.errors import GetphyloError
from typing import List


def build_all_trees(files: List, cpus: int, method: str, tree_directory: str, output: str) -> None:
    '''
    builds all trees in from a list of files
        Arguments:
            files: list of alignment files to be processed
            cpus: number of cpus for parallelisation
            method: pyhlogenetic method (e.g. fasttree)
        Returns:
            None 
    '''
    args_list = []
    if method == 'fasttree':
        for filename in files:
            outfile = os.path.join(
                tree_directory, os.path.basename(io.change_extension(filename, "tree"))
                )
            args_list.append([filename, outfile])
        io.run_in_parallel(fasttree.run_fasttree, args_list, cpus)
    elif method == 'iqtree':
        partition = os.path.join(output, 'partition.txt')
        for filename in files:
            outfile = os.path.join(
                tree_directory, os.path.basename(os.path.splitext(filename)[0])
                )
            args_list.append([filename, outfile])
        io.run_in_parallel(iqtree.run_iqtree, args_list, cpus)
    else:
        raise GetphyloError(method + ' is not a phylogenetic tool.')

def make_trees(output: str, build_all: bool, method: str, cpus: int) -> None:
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
        files = glob.glob(os.path.join(output, 'aligned_fasta/*.fasta'))
        build_all_trees(files, cpus, method, tree_directory, output)
    else:
        filename = os.path.join(output, 'aligned_fasta/combined_alignment.fasta')
        if method == 'fasttree':
            output = os.path.join(tree_directory, 'combined_alignment.tree')
            fasttree.run_fasttree(filename, output)
        elif method == 'iqtree':
            partition = os.path.join(output, 'partition.txt')
            output = os.path.join(tree_directory, 'combined_alignment')
            iqtree.run_iqtree(filename, output, partition)
        else:
            raise GetphyloError(method + ' is not a phylogenetic tool.')
    logging.info("CHECKPOINT: TREES_BUILT")
