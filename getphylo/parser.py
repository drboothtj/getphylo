'''
Create an argument parser using argparse

Functions:
    get_checkpoint() -> checkpoint
    get_gbks() -> gbks
    get_parser() -> parser
    get_seed() -> seed
    parse_args() -> args
'''

import argparse

def get_checkpoint():
    '''get checkpoint from the parser'''
    args = parse_args()
    checkpoint = int(args.checkpoint)
    return checkpoint

def get_gbks():
    '''get genbank path from the parser'''
    args = parse_args()
    gbks = args.gbks
    return gbks

def get_parser():
    ''''Create a parser object specific to getphylo'''
    parser = argparse.ArgumentParser(
        "getphylo",
        description=
        "getphylo: a python package to produce heuristic phylogentic trees from genomic data",
        epilog="Written by Dr. Thom Booth, 2022."
        )
    parser.add_argument(
        '-s',
        '--seed',
        default=None,
        help='path to a genbankfile with for the target organism (default: random)'
        )
    parser.add_argument(
        '-g',
        '--gbks',
        default="*.gbk",
        help='list of genbank files to use in the phylogeny (default: *.gbk)'
        )
    parser.add_argument(
        '-f',
        '--find',
        default=None,
        help=(
            'integer indicating the number of loci to find in the seed genome '
            '(default: all suitable loci)'
        )
        )
    parser.add_argument(
        '-l',
        '--length',
        default=200,
        help=(
            'interger indicating the minimum length of loci to be included in the analysis '
            '(default: 200)'
        )  
        )
    parser.add_argument(
        '-p',
        '--presence',
        default=100,
        help=(
            'interger indicating the percentage of genomes each loci must be present in '
            '(default: 100)'
        )
        )
    parser.add_argument(
        '-c',
        '--checkpoint',
        default=0,
        help=(
            'interger indicating the checkpoint from which to continue the analysis (default:0)'
            '(checkpoint info: '
            '1: skip extraction of CDSs from genbank; '
            '2: skip creation of diamond databases; '
        )
        )
    return parser

def get_seed():
    '''get seed information from the parser'''
    args = parse_args()
    seed = args.seed
    return seed

def parse_args():
    '''get the arguments from the console via the parser'''
    arg_parser = get_parser()
    args = arg_parser.parse_args()
    return args

    #add extra dmnd options
    #add max threshold?
