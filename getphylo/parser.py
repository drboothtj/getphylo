'''
Create an argument parser using argparse

Functions:
    get_parser() -> Parser
    parse_args() -> List
'''

import argparse
from argparse import RawTextHelpFormatter

import logging
from getphylo.utils.checkpoint import Checkpoint

def get_parser():
    ''''Create a parser object specific to getphylo'''
    parser = argparse.ArgumentParser(
        "getphylo",
        description="""
        getphylo: a python package to produce heuristic phylogentic trees from genomic data\n
        example usage: getphylo -g 'input/*.gb' -c 4
        """,
        epilog="Written by Dr. Thom Booth, 2022.",
        formatter_class=RawTextHelpFormatter
        )
    parser.add_argument(
        '-b',
        '--build-all',
        action='store_true',
        help=(
            'build phylogenetic trees for all loci, not just concatenated alignment '
            '(default: %(default)s)'
        )
        )
    parser.add_argument(
        '-cp',
        '--checkpoint',
        default='START',
        type=str,
        choices=[cp.name for cp in Checkpoint],
        help=(
            'string indicating the checkpoint to start from '
            'START = default\n '
            'FASTA_EXTRACTED = Skip extracting fasta sequences from genbank files\n '
            'DIAMOND_BUILT = Skip building diamond databases\n '
            'SINGLETONS_IDENTIFIED = Skip identifying singletons from the seed genome\n '
            'SINGLETONS_SEARCHED = Skip searching singletons against other genomes\n '
            'SINGLETONS_THRESHOLDED = Skip thresholding of singletons\n '
            'SINGLETONS_EXTRACTED = Skip extract fasta sequences for alignments\n '
            'SINGLETONS_ALIGNED = Skip individual protein alignments\n '
            'ALIGNMENTS_COMBINED = Skip combining alignments\n '
            'TREES_BUILT = Skip building trees\n '
            'DONE = Done\n '
            '(default: %(default)s)'
        )
        )
    parser.add_argument(
        '-c',
        '--cpus',
        default=1,
        type=int,
        help=(
            'The number of cpus to use for paralleslisation '
            '(default: %(default)s)'
        )
        )
    parser.add_argument(
        '-f',
        '--find',
        default=-1,
        type=int,
        help=(
            'integer indicating the number of loci to find in the seed genome '
            '(default: %(default)s)'
        )
        )
    parser.add_argument(
        '-g',
        '--gbks',
        default="*.gbk",
        type=str,
        help='string indicating the genbank files to use in the phylogeny (default: %(default)s)'
        )
    parser.add_argument(
        '-ia',
        '--ignore-bad-annotations',
        action='store_true',
        help='ignore missing annotations - NOT RECCOMMENDED (default: %(default)s)'
        )
    parser.add_argument(
        '-ir',
        '--ignore-bad-records',
        action='store_true',
        help='ignore poorly formated records - NOT RECCOMMENDED (default: %(default)s)'
        )
    parser.add_argument(
        '-l',
        '--logging',
        default='ERROR',
        choices=list(logging._nameToLevel.keys()),
        help='set the logging level (default: %(default)s)'
    )
    parser.add_argument(
        '-m',
        '--method',
        default='fasttree',
        choices=['fasttree', 'iqtree'],
        help=(
            'choose the phylogenetic method '
            '- NOTE: Using iqtree will test individual gene modeles but will exponentially increase the run time '
            '(default: %(default)s)'
        )
    )
    parser.add_argument(
        '-max',
        '--maxlength',
        default=2000,
        type=int,
        help=(
            'interger indicating the minimum length of loci to be included in the analysis '
            '(default: %(default)s)'
        )
        )
    parser.add_argument(
        '-min',
        '--minlength',
        default=200,
        type=int,
        help=(
            'interger indicating the minimum length of loci to be included in the analysis '
            '(default: %(default)s)'
        )
        )
    parser.add_argument(
        '-minl',
        '--minloci',
        default=1,
        type=int,
        help=(
            'minimum number of loci required to continue to alignment and tree building steps '
            '(default: %(default)s)'
        )
        )
    parser.add_argument(
        '-maxl',
        '--maxloci',
        default=1000,
        type=int,
        help=(
            'maximum number of loci required to continue to alignment and tree building steps '
            '(default: %(default)s)'
        )
        )
    parser.add_argument(
        '-o',
        '--output',
        default='output',
        type=str,
        help=(
            'a string designating the name of the folder to output the results'
            '(default: %(default)s)'
        )
        )
    parser.add_argument(
        '-p',
        '--presence',
        default=100,
        type=float,
        help=(
            'interger indicating the percentage of genomes each loci must be present in '
            '(default: %(default)s)'
        )
        )
    parser.add_argument(
        '-r',
        '--random-seed-number',
        default=None,
        type=int,
        help=(
            'interger to be used as a seed for randomising loci selection, random if left as None'
            '(default: None)'
        )
        )
    parser.add_argument(
        '-s',
        '--seed',
        default=None,
        type=str,
        help='path to a genbankfile with for the target organism (default: %(default)s)'
        )
    parser.add_argument(
        '-t',
        '--tag',
        default="locus_tag",
        type=str,
        help='string indicating the GenBank annotations to extract (default: %(default)s)'
        )
    return parser

def parse_args():
    '''get the arguments from the console via the parser'''
    arg_parser = get_parser()
    args = arg_parser.parse_args()
    return args
