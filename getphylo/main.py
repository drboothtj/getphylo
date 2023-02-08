'''
Main routine for getphylo

Functions:
main()
set_seed()

'''
import glob
import logging
import os
from typing import List
from getphylo import align, extract, parser, screen, trees
from getphylo.utils.errors import BadSeedError, NoFinalLociError
from getphylo.utils.checkpoint import Checkpoint

def initialize_logging():
    '''Set up and configure logging.
        Arguments: None
        Returns: None'''
    logging_level = logging.DEBUG
    logging.basicConfig(
        #filename='getphylo.log',
        level=logging_level,
        format='[%(asctime)s] %(levelname)-10s: %(message)s',
        datefmt='%H:%M:%S') #maybe change to run time?
    logging.info("Running getphylo version 0.0.1")

def check_seed(checkpoint: Checkpoint, gbk_search_string: str) -> str:
    '''Set a seed for a new analysis and raise an error if continuing an old analysis.
        Arguments: 
            checkpoint: the checkpoint supplied by the user
            gbk_search_string: the string used to filter the glob (e.g. *.gbk)
        Returns:
            seed: the filename of the selceted seed genome'''
    if checkpoint > 0:
        raise BadSeedError('A checkpoint has been set! Please ensure the seed is defined.')
    gbks = glob.glob(gbk_search_string)
    if not gbks:
        raise BadSeedError(f'No files found in {gbk_search_string}.')
    seed = gbks[0]
    logging.info(
        f'No seed defined. Using first file in glob ({seed}) as seed.'
        )
    return seed

def main():
    '''main routine for getphylo
        Arguments: None
        Returns: None'''
    args = parser.parse_args()
    logging.getLogger().setLevel(args.logging) 
    
    gbks = args.gbks
    checkpoint = Checkpoint[args.checkpoint.upper()]
    seed = args.seed
    output = os.path.abspath(args.output)
    
    if seed is None: check_seed(checkpoint, gbks)
    logging.info(f'The seed genome is {seed}.')
    
    ### Begin main workflow
    ### extract.py
    if checkpoint < DIAMOND_BUILT:
        try:
            io.make_folder(output)
        except:
            logging.warning(
                f'ALERT: {output} already exists. Continuing analysis in that directory.'
                )
        tag_args = [args.tag, args.ignore]
        extract.extract_data(checkpoint, output, gbks, tag_args)
    ### screen.pys
    if checkpoint < SINGLETONS_THRESHOLDED:
        thresholds = [args.find, args.minlength, args.maxlength, args.presence, args.minloci]
        final_loci = screen.get_target_proteins(checkpoint, output, seed, thresholds)

    ### before continuing check final loci is defined, otherwise read from file
    try:
        assert final_loci
    except AssertionError:
        logging.info('Final loci not detected. This is normal if restarting from a later checkpoint.')
        try:
            final_loci_path = os.path.join(output, 'final_loci.txt')
            logging.info(f'Attempting to read {final_loci_path}')
            final_loci = screen.get_loci_from_file(final_loci_path)
        except:
            raise NoFinalLociError(
                'Final loci could not be read from final_loci.txt.'
                'If restarting from a checkpoint ensure there is a final_loci.txt file in the specified output folder.'
                )

    ### align.py
    if checkpoint < Checkpoint.ALIGNMENTS_COMBINED:
        align.make_alignments(checkpoint, output, final_loci, gbks)
        
    ### trees.py
    if checkpoint < Checkpoint.TREES_BUILT:
        trees.make_trees(output)
    logging.info("Analysis complete. Thank you for using getphylo!")
    