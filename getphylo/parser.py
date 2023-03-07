'''
Create an argument parser using argparse

Functions:
    get_parser() -> Parser
    parse_args() -> List
'''

import argparse
import logging
from getphylo.utils.checkpoint import Checkpoint

def get_parser():
    ''''Create a parser object specific to getphylo'''
    parser = argparse.ArgumentParser(
        "getphylo",
        description=
        "getphylo: a python package to produce heuristic phylogentic trees from genomic data",
        epilog="Written by Dr. Thom Booth, 2022."
        )
    parser.add_argument(
        '-c',
        '--checkpoint',
        default='START',
        type=str,
        choices=[cp.name for cp in Checkpoint],
        help=(
            'string indicating the checkpoint to start from '
            'START = default '
            'FASTA_EXTRACTED = Skip extracting fasta sequences from genbank files '
            'DIAMOND_BUILT = Skip building diamond databases '
            'SINGLETONS_IDENTIFIED = Skip identifying singletons from the seed genome '
            'SINGLETONS_SEARCHED = Skip searching singletons against other genomes '
            'SINGLETONS_THRESHOLDED = Skip thresholding of singletons '
            'SINGLETONS_EXTRACTED = Skip extract fasta sequences for alignments '
            'SINGLETONS_ALIGNED = Skip individual protein alignments '
            'ALIGNMENTS_COMBINED = Skip combining alignments '
            'TREES_BUILT = Skip building trees '
            'DONE = Done '
            'default: %(default)s'
        )
        )
    parser.add_argument(
        '-f',
        '--find',
        default=-1,
        type=int,
        help=(
            'integer indicating the number of loci to find in the seed genome '
            'default: %(default)s)'
        )
        )
    parser.add_argument(
        '-g',
        '--gbks',
        default="*.gbk",
        type=str,
        help='string indicating the genbank files to use in the phylogeny. default: %(default)s'
        )
    parser.add_argument(
        '-ia',
        '--ignore-bad-annotations',
        default=False,
        type=bool,
        help='ignore missing annotations - NOT RECCOMMENDED (default: %(default)s)'
        )
    parser.add_argument(
        '-ir',
        '--ignore-bad-records',
        default=False,
        type=bool,
        help='ignore poorly formated records - NOT RECCOMMENDED (default: %(default)s)'
        )
    parser.add_argument(
        '-l',
        '--logging',
        default='WARNING',
        choices=list(logging._nameToLevel.keys()),
        help='set the logging level (default: %(default)s)'
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
        default=1,
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
    
    #add extra dmnd options
    #add upper threshold
    