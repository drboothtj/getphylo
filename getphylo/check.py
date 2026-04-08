'''
check inputs for getphylo
    functions:
        check_executables(args) -> None
        check_seed(checkpoint: Checkpoint, gbk_search_string: str) -> str
        check_gbks(gbks: str) -> None
'''
import logging
import glob
import os

from shutil import copy as cp
from typing import Tuple

from getphylo.utils import io
from getphylo.utils.checkpoint import Checkpoint
from getphylo.utils.errors import (
    BadInputError, #too generic
    BadMethodError,
    BadSeedError, #too generic
    FolderExistsError,
    )

def check_executables(args) -> None:
    '''
    check excutables are defined and break early if not
        arguments:
            args: the args from the args parsers
        returns:
            None
    '''
    logging.debug("Checking diamond...")
    io.run_in_command_line([args.diamond, 'help'])
    logging.debug("Checking muscle...")
    io.run_in_command_line([args.muscle])
    if args.method =="fasttree":
        logging.debug("Checking fasttree...")
        io.run_in_command_line([args.fasttree])
    elif args.method =="iqtree":
        logging.debug("Checking iqtree...")
        io.run_in_command_line([args.iqtree])
    else:
        raise BadMethodError(args.method)
    logging.debug("Executables checked successfully.")

def check_input(gbks: str) -> None:
    '''
    check at least three files are found by the provided search string
    otherwise, raise BadInputError
        arguments:
            gbks: search string from the parser
        returns:
            None
    '''
    gbk_count = len(glob.glob(gbks))
    if os.path.isdir(gbks):
        raise BadInputError(
            gbks + ' is a directory. Please provide a search string (e.g. \'my_dir/*.gbk\').'
            )
    if gbk_count < 3:
        raise BadInputError(
            'getphylo requires at least 3 input sequences. '
            f'{gbk_count} provided. '
            'Please, check input search sting parameter (-g/-f) and try again.'
            )

def check_seed(checkpoint: Checkpoint, gbk_search_string: str) -> str:
    '''
    Set a seed for a new analysis and raise an error if continuing an old analysis.
        Arguments:
            checkpoint: the checkpoint supplied by the user
            gbk_search_string: the string used to filter the glob (e.g. '*.gbk')
        Returns:
            seed: the filename of the selected seed genome
    '''
    if checkpoint > 0:
        raise BadSeedError('A checkpoint has been set! Please ensure the seed is defined.')
    gbks = glob.glob(gbk_search_string)
    if not gbks:
        raise BadSeedError(f'No files found in {gbk_search_string}.')
    seed = gbks[0]
    logging.warning(
        'No seed defined. Using first file in glob (%s) as seed.', seed
        )
    return seed

def check_fastas(files: str, seed: str, output: str, checkpoint: Checkpoint) -> str:
    '''
    handles fasta file inputs by creating copies within the correct file structure
    also sets the checkpoint if appropriate
        arguments:
            files: the glob string for the fasta files (e.g. *.fasta)
            seed: the name of the seed file, if provided, else None
            output: path to output dir
        returns:
            checkpoint: checkpoint adjusted to FASTA_EXTRACTED if necissary
            seed: selected seed from fasta input
    '''
    check_input(files)
    logging.warning(
        'fasta has been set, ignoring genbank files and using %s as input', files
    )
    fastas = glob.glob(files)
    if seed is None:
        seed = fastas[0]
    fasta_path = os.path.join(output, 'fasta')
    io.make_folder(fasta_path)
    for file in fastas:
        cp(file, os.path.join(fasta_path, os.path.splitext(os.path.basename(file))[0] + '.fasta'))
    checkpoint = max(checkpoint, Checkpoint.FASTA_EXTRACTED)
    return checkpoint, files, seed

def initialise_analysis(args) -> Tuple[Checkpoint, str, str]:
    '''
    perform initialisation checks
        arguments:
            args: arguments object from argparse
        returns:
            checkpoint: getphylo checkpoint
            files: glob string for the input files
            seed: seed genome for the analysis
    '''
    logging.debug('Performing initialisation checks...')
    output = os.path.abspath(args.output)
    checkpoint = Checkpoint[args.checkpoint]
    try:
        io.make_folder(output)
    except FolderExistsError:
        logging.warning(
            '%s already exists. Continuing analysis in that directory.', output
            )
    seed = args.seed
    if args.fasta:
        checkpoint, files, seed = check_fastas(
            args.fasta, seed, output, Checkpoint[args.checkpoint.upper()]
            )
    else:
        files = args.gbks
        check_input(files)
    if seed is None:
        seed = check_seed(checkpoint, files)
    logging.info('The seed genome is %s!', seed)
    check_executables(args)
    logging.debug('Initialisation checks passed!')
    return checkpoint, files, seed
