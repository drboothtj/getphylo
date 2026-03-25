'''
check inputs for getphylo
    functions:
        !!!
'''
import logging
import glob
import os

from getphylo.utils import io
from getphylo.utils.checkpoint import Checkpoint
from getphylo.utils.errors import (
    BadInputError, #too generic
    BadMethodError,
    BadSeedError, #too generic
    FolderExistsError,
    NoFinalLociError
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

def check_gbks(gbks: str) -> None:
    '''
    check at least three files are  found by the provided search string
    otherwise, raise BadInputError
        arguments:
            gbks: search string from the parser
        returns:
            None
    '''
    gbk_count = len(glob.glob(gbks))
    if gbk_count < 3:
        raise BadInputError(
            'getphylo requires at least 3 input sequences. '
            f'{gbk_count} provided. '
            'Please, check input search sting parameter (-g) and try again.'
            )
    if os.path.isdir(gbks):
        raise BadInputError(
            gbks + ' is a directory. Please provide a search string (e.g. \'my_dir/*.gbk\').'
            )
