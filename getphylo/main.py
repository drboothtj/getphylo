'''
Main routine for getphylo

Functions:
    initialize_logging() -> None
    main()
'''
import logging
import os
from getphylo import align, check, extract, parser, screen, trees
from getphylo.utils.errors import (
    BadMethodError,
    NoFinalLociError
    ) #move all to check in some way
from getphylo.utils.checkpoint import Checkpoint

def initialize_logging() -> None: #move call to main.py so we can add other arguments
    '''Set up and configure logging.
        Arguments: None
        Returns: None
        '''
    logging_level = logging.DEBUG
    logging.basicConfig(
        level=logging_level,
        format='[%(asctime)s] %(levelname)-10s: %(message)s',
        datefmt='%H:%M:%S',
        handlers=[
            logging.StreamHandler(),
            logging.FileHandler('getphylo.log') #make customisable
        ])
    logging.info("Running getphylo version 1.1.0.")

def main():
    '''
    main routine for getphylo
        Arguments: None
        Returns: None
    '''
    args = parser.parse_args()

    logging.getLogger().setLevel(args.logging) #default to upper (in parser)!
    ###ALWAYS SET LOGGING LEVEL FIRST!

    output = os.path.abspath(args.output)
    diamond_args = (args.diamond, args.identity, args.query_coverage, args.subject_coverage)

    checkpoint, files, seed = check.initialise_analysis(args)

    ### Begin main workflow
    ### extract.py
    if checkpoint < Checkpoint.DIAMOND_BUILT:
        extract.extract_data(
            checkpoint, output, files, args.tag, args.ignore_bad_annotations,
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
        align.make_alignments(checkpoint, output, final_loci, files, args.cpus, args.muscle)

    ### trees.py
    if checkpoint < Checkpoint.TREES_BUILT:
        build_all = args.build_all
        if args.method == 'fasttree':
            tree_builder = args.fasttree
        elif args.method == 'iqtree':
            tree_builder = args.iqtree
        else:
            raise BadMethodError(args.method)
        trees.make_trees(output, build_all, args.method, args.cpus, tree_builder)
    logging.info("CHECKPOINT: DONE")
    logging.info("Analysis complete. Thank you for using getphylo!")
