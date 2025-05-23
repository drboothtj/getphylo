'''
Main routine for getphylo

Functions:
    initialize_logging() -> None
    check_seed(checkpoint: Checkpoint, gbk_search_string: str) -> str
    check_gbks(gbks: str) -> None
    main()
'''
import glob
import logging
import os
from getphylo import align, extract, parser, screen, trees
from getphylo.utils import io
from getphylo.utils.errors import (
    BadInputError,
    BadMethodError,
    BadSeedError,
    FolderExistsError,
    NoFinalLociError
    )
from getphylo.utils.checkpoint import Checkpoint

def initialize_logging() -> None:
    '''Set up and configure logging.
        Arguments: None
        Returns: None'''
    logging_level = logging.DEBUG
    logging.basicConfig(
        #filename='getphylo.log',
        level=logging_level,
        format='[%(asctime)s] %(levelname)-10s: %(message)s',
        datefmt='%H:%M:%S')
    logging.info("Running getphylo version 1.0.0.")

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

def main():
    '''main routine for getphylo
        Arguments: None
        Returns: None'''
    args = parser.parse_args()
    logging.getLogger().setLevel(args.logging)

    gbks = args.gbks
    check_gbks(gbks)
    checkpoint = Checkpoint[args.checkpoint.upper()]
    seed = args.seed
    output = os.path.abspath(args.output)
    diamond_args = (args.diamond, args.identity, args.query_coverage, args.subject_coverage)

    if os.path.isdir(gbks):
        raise BadInputError(
            gbks + ' is a directory. Please provide a search string (e.g. \'my_dir/*.gbk\').'
            )
    if seed is None:
        seed = check_seed(checkpoint, gbks)
    logging.info('The seed genome is %s!', seed)

    ### Begin main workflow
    ### extract.py
    if checkpoint < Checkpoint.DIAMOND_BUILT:
        try:
            io.make_folder(output)
        except FolderExistsError:
            logging.warning(
                'ALERT: %s already exists. Continuing analysis in that directory.', output
                )
        extract.extract_data(
            checkpoint, output, gbks, args.tag, args.ignore_bad_annotations,
            args.ignore_bad_records, args.cpus, args.diamond
            )
    ### screen.py
    final_loci = None
    if checkpoint < Checkpoint.SINGLETONS_THRESHOLDED:
        thresholds = [
            args.find, args.minlength, args.maxlength, args.presence, args.minloci, args.maxloci,
            ]
        final_loci = screen.get_target_proteins(
            checkpoint, output, seed, thresholds, args.cpus, args.random_seed_number, diamond_args
            )
    ### before continuing check final loci is defined, otherwise read from file
    try:
        assert final_loci
    except AssertionError:
        logging.info(
            'Final loci not detected. This is normal if restarting from a later checkpoint.'
            )
        try:
            final_loci_path = os.path.join(output, 'final_loci.txt')
            logging.info('Attempting to read %s', final_loci_path)
            final_loci = screen.get_loci_from_file(final_loci_path)
        except:
            raise NoFinalLociError(
                'Final loci could not be read from final_loci.txt.'
                'If restarting from a checkpoint ensure there is a final_loci.txt file'
                'in the specified output folder.'
                )

    ### align.py
    if checkpoint < Checkpoint.ALIGNMENTS_COMBINED:
        align.make_alignments(checkpoint, output, final_loci, gbks, args.cpus, args.muscle)

    ### trees.py
    if checkpoint < Checkpoint.TREES_BUILT:
        build_all = args.build_all
        if args.method == 'fasttree':
            tree_builder = args.fasttree
        elif args.method == 'iqtree':
            tree_builder = args.iqtree
        else:
            raise BadMethodError(
                'Neither fasttree or iqtree was selected.'
                'It should not be possible for you to generate this error - please report!')
        trees.make_trees(output, build_all, args.method, args.cpus, tree_builder)
    logging.info("CHECKPOINT: DONE")
    logging.info("Analysis complete. Thank you for using getphylo!")
