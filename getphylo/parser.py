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
    return parser

def get_config_parser(arg_parser):
    '''
    Create an argument group for basic config details.
        Arguments:
            arg_parser: the basic argument parser
        Returns:
            arg_parser: the argument parser with arguments added
    '''
    config_parser = arg_parser.add_argument_group(
            'basic configuration', 'basic configuration of getphylo'
            )
    config_parser.add_argument(
        '-c',
        '--cpus',
        default=1,
        type=int,
        help=(
            'The number of cpus to use for paralleslisation\n'
            '(default: %(default)s)'
        )
        )
    config_parser.add_argument(
        '-cp',
        '--checkpoint',
        default='START',
        type=str,
        choices=[cp.name for cp in Checkpoint],
        help=(
            'string indicating the checkpoint to start from:\n'
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
    config_parser.add_argument(
        '-l',
        '--logging',
        default='ERROR',
        choices=[
            logging.getLevelName(level) for level in [logging.DEBUG, logging.INFO, logging.WARNING]
            ],
        help='set the logging level\n'
        '(default: %(default)s)'
    )
    return arg_parser

def get_phylo_parser(arg_parser):
    '''
    Create an argument group for phylogenetic analysis
        Arguments:
            arg_parser: the basic argument parser
        Returns:
            arg_parser: the argument parser with arguments added
    '''
    phylo_parser = arg_parser.add_argument_group(
        'phylogenetics', 'parameters for building phylogentic trees'
        )
    phylo_parser.add_argument(
    '-b',
    '--build-all',
    action='store_true',
    help=(
        'build phylogenetic trees for all loci, not just concatenated alignment\n'
        'NOTE: not recommended to use in combination with iqtree\n'
        '(default: %(default)s)'
    )
    )
    phylo_parser.add_argument(
        '-m',
        '--method',
        default='fasttree',
        choices=['fasttree', 'iqtree'],
        help=(
            'choose the phylogenetic method\n'
            'NOTE: Using iqtree will test individual gene models\n'
            'but will exponentially increase the run time\n'
            'NOTE: Not recommended to use in combination with --build-all\n'
            '(default: %(default)s)'
        )
    )
    return arg_parser

def get_records_parser(arg_parser):
    '''
    Create an argument group for ignoring malformed record and annotations.
        Arguments:
            arg_parser: the basic argument parser
        Returns:
            arg_parser: the argument parser with arguments added
    '''
    record_parser = arg_parser.add_argument_group(
            'records and annotations - NOT RECOMMENDED', 
            'deal with malformed records and annotations \n' +
            'its is NOT RECOMMENDED to use these options unless you are confident with their usage'
            )
    record_parser.add_argument(
        '-ia',
        '--ignore-bad-annotations',
        action='store_true',
        help='ignore missing annotations - NOT RECOMMENDED\n (default: %(default)s)'
        )
    record_parser.add_argument(
        '-ir',
        '--ignore-bad-records',
        action='store_true',
        help='ignore poorly formated records - NOT RECOMMENDED\n (default: %(default)s)'
        )
    return arg_parser

def get_search_parser(arg_parser):
    '''
    Create an argument group for changing search parameters
        Arguments:
            arg_parser: the basic argument parser
        Returns:
            arg_parser: the argument parser with arguments added
    '''
    search_parser = arg_parser.add_argument_group(
            'search filtering', 
            'alter the parameters of the orthologue search'
            )
    search_parser.add_argument(
        '-f',
        '--find',
        default=-1,
        type=int,
        help=(
            'integer indicating the number of loci to find in the seed genome\n'
            '(default: %(default)s)'
        )
        )
    search_parser.add_argument(
        '-max',
        '--maxlength',
        default=2000,
        type=int,
        help=(
            'interger indicating the minimum length of loci to be included in the analysis\n'
            '(default: %(default)s)'
        )
        )
    search_parser.add_argument(
        '-min',
        '--minlength',
        default=200,
        type=int,
        help=(
            'interger indicating the minimum length of loci to be included in the analysis\n'
            '(default: %(default)s)'
        )
        )
    search_parser.add_argument(
        '-minl',
        '--minloci',
        default=1,
        type=int,
        help=(
            'minimum number of loci required to continue to alignment and tree building steps\n'
            '(default: %(default)s)'
        )
        )
    search_parser.add_argument(
        '-maxl',
        '--maxloci',
        default=1000,
        type=int,
        help=(
            'maximum number of loci required to continue to alignment and tree building steps\n'
            '(default: %(default)s)'
        )
        )
    search_parser.add_argument(
        '-p',
        '--presence',
        default=100,
        type=float,
        help=(
            'interger indicating the percentage of genomes each loci must be present in\n'
            '(default: %(default)s)'
        )
        )
    return arg_parser

def get_seed_parser(arg_parser):
    '''
    Create an argument group for changing search parameters
        Arguments:
            arg_parser: the basic argument parser
        Returns:
            arg_parser: the argument parser with arguments added
    '''
    seed_parser = arg_parser.add_argument_group(
            'seed parameters', 
            'seed options for the genome or loci seed'
            )
    seed_parser.add_argument(
        '-r',
        '--random-seed-number',
        default=None,
        type=int,
        help=(
            'interger to be used as a seed for randomising loci selection\n'
            'random if left as None\n'
            '(default: None)'
        )
        )
    seed_parser.add_argument(
        '-s',
        '--seed',
        default=None,
        type=str,
        help='path to a genbankfile with for the target organism\n'
        'first in glob if left as None\n'
        'NOTE: using the smallest genome will generally result in lower runtimes\n'
        'NOTE: this will only effect the results if -p is used\n'
        '(default: %(default)s)'
        )
    return arg_parser

def get_io_parser(arg_parser):
    '''
    Create an argument group for changing input and output parameters
        Arguments:
            arg_parser: the basic argument parser
        Returns:
            arg_parser: the argument parser with arguments added
    '''
    io_parser = arg_parser.add_argument_group(
            'basic input and output', 
            'set options for input and output'
            )
    io_parser.add_argument(
        '-g',
        '--gbks',
        default="*.gbk",
        type=str,
        help='string indicating the genbank files to use in the phylogeny\n'
        '(default: %(default)s)'
        )
    io_parser.add_argument(
        '-o',
        '--output',
        default='output',
        type=str,
        help=(
            'a string designating the name of the folder to output the results\n'
            '(default: %(default)s)'
        )
    )
    io_parser.add_argument(
        '-t',
        '--tag',
        default="locus_tag",
        type=str,
        help='string indicating the GenBank annotations to extract\n'
        '(default: %(default)s)'
        )
    return arg_parser

def get_exe_parser(arg_parser):
    '''
    Create an argument group for providing custom executable paths
        Arguments:
            arg_parser: the basic argument parser
        Returns:
            arg_parser: the argument parser with arguments added
    '''
    exe_parser = arg_parser.add_argument_group(
            'executables'
            'define custom executable paths, give the full path or just the command if set in PATH', 
            )
    exe_parser.add_argument(
        '-d',
        '--diamond',
        default="diamond",
        type=str,
        help=(
        'path to DIAMOND executable \n'
        '(default: %(default)s)'
        )
    )
    exe_parser.add_argument(
        '-mu',
        '--muscle',
        default="muscle",
        type=str,
        help=(
        'path to MUSCLE executable \n'
        '(default: %(default)s)'
        )
    )
    exe_parser.add_argument(
        '-ft',
        '--fasttree',
        default="fasttree",
        type=str,
        help=(
        'path to fasttree2 executable \n'
        '(default: %(default)s)'
        )
    )
    exe_parser.add_argument(
        '-iq',
        '--iqtree',
        default="iqtree",
        type=str,
        help=(
        'path to IQ-TREE executable \n'
        '(default: %(default)s)'
        )
    )
    return arg_parser

def get_arguments(arg_parser):
    '''
    Add arguments and argument groups to the parser
        Arguments:
            arg_parser: the basic argument parser
        Returns:
            arg_parser: the argument parser with arguments added
    '''
    arg_parser = get_io_parser(arg_parser)
    arg_parser = get_config_parser(arg_parser)
    arg_parser = get_search_parser(arg_parser)
    arg_parser = get_phylo_parser(arg_parser)
    arg_parser = get_seed_parser(arg_parser)
    arg_parser = get_records_parser(arg_parser)
    arg_parser = get_exe_parser(arg_parser)
    return arg_parser

def parse_args():
    '''
    get the arguments from the console via the parser
    '''
    arg_parser = get_parser()
    arg_parser = get_arguments(arg_parser)
    args = arg_parser.parse_args()
    return args
