'''
Create an argument parser using argparse

Functions:
    get_seed_thresholds -> loci_to_find, loci_min_length, loci_max_length
    get_checkpoint() -> checkpoint
    get_gbks() -> gbks
    get_parser() -> parser
    get_seed() -> seed
    parse_args() -> args
'''

import argparse

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
        default=0,
        type=int,
        help=(
            'interger indicating the checkpoint from which to continue the analysis (default:0)'
            '(checkpoint info: '
            '1: skip extraction of CDSs from genbank; '
            '2: skip creation of diamond databases; '
            '3: skip identifying singletons from seed genome'
            '4: skip confirming singletons in other genomes'
            '5: skip thresholding and identifying final loci for alignments'
        )
        )
    parser.add_argument(
        '-f',
        '--find',
        default=-1,
        type=int,
        help=(
            'integer indicating the number of loci to find in the seed genome '
            '(default: all suitable loci (-1))'
        )
        )
    parser.add_argument(
        '-g',
        '--gbks',
        default="*.gbk",
        type=str,
        help='string indicating the genbank files to use in the phylogeny (default: *.gbk)'
        )
    parser.add_argument(
        '-max',
        '--maxlength',
        default=2000,
        type=int,
        help=(
            'interger indicating the minimum length of loci to be included in the analysis '
            '(default: 2000)'
        )
        )
    parser.add_argument(
        '-min',
        '--minlength',
        default=200,
        type=int,
        help=(
            'interger indicating the minimum length of loci to be included in the analysis '
            '(default: 200)'
        )
        )
    parser.add_argument(
        '-ml',
        '--minloci',
        default=1, #set types for args!!
        type=int,
        help=(
            'minimum number of loci required to continue to alignment and tree building steps '
            '(default: 1)'
        )
        )
    parser.add_argument(
        '-o',
        '--output',
        default='output',
        type=str,
        help=(
            'a string designating the name of the folder to output the results'
            '(default: output)'
        )
        )
    parser.add_argument(
        '-p',
        '--presence',
        default=100,
        type=float,
        help=(
            'interger indicating the percentage of genomes each loci must be present in '
            '(default: 100)'
        )
        )
    parser.add_argument(
        '-s',
        '--seed',
        default=None,
        type=str,
        help='path to a genbankfile with for the target organism (default: random)'
        )
    return parser

def parse_args():
    '''get the arguments from the console via the parser'''
    arg_parser = get_parser()
    args = arg_parser.parse_args()
    return args

    #add extra dmnd options
    #add max threshold?
    #set types for args!!
    