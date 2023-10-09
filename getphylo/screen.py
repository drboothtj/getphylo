'''
Screen fasta files and blast databases for singletons and extract sequences

Functions:
    get_unique_hits_from_tsv(file: str) -> List
    get_seed_paths(seed: str, output: str) -> Tuple[str, str, str]
    get_singletons_from_seed(seed, output, thresholds, random_seed_number)
    get_loci_from_file(file: str) -> List
    search_candidates(output: str, cpus: int) -> None
    score_locus(locus: str, files: List) -> Tuple[int, bool, List]
    process_final_loci(final_loci: List, minimum_loci: int, output: str) -> None
    do_thresholding(
        target_loci: List, presence_threshold: float, maximum_loci: int, output: str
    ) -> List
    threshold_loci(target_loci: List, thresholds: List, output: str) -> List
    write_pa_table(pa_table: List, loci: List, output: str) -> None
    get_target_proteins(
        checkpoint: Checkpoint, output: str, seed: str, thresholds: List,
        cpus: int, random_seed_number: int
    ) -> None
'''
import os
import glob
import logging
from collections import Counter
import random
from typing import List, Tuple

from getphylo.ext import diamond
from getphylo.utils import io
from getphylo.utils.checkpoint import Checkpoint
from getphylo.utils.errors import(
    FileAlreadyExistsError,
    InsufficientLociError,
    NoCandidateLociError
)

def get_unique_hits_from_tsv(file: str) -> List:
    '''
    Finds unique hits from the provided dmnd result file.
        Arguments:
            file: path to tsv file
        Returns:
            unique_hits: list of hits that are unique in the genome
    '''
    hits = []
    unique_hits = []
    lines = io.read_tsv(file)
    for line in lines:
        hits.append(line[0])
    counter = Counter(hits)
    for hit in counter:
        if counter[hit] == 1:
            unique_hits.append(hit)
    return unique_hits

def get_seed_paths(seed: str, output: str) -> Tuple[str, str, str]:
    '''
    Take the path to a genbank file and return files with .fasta, .dmnd and .tsv extensions.
    Arguments:
        seed: path to the seed genome
        output: path to the output directory
    Returns:
        seed_fasta, seed_dmnd, seed_tsv
    '''
    seed = os.path.basename(seed)
    seed_fasta = io.change_extension(seed, "fasta")
    seed_fasta = os.path.join(output, 'fasta', seed_fasta)
    seed_dmnd = io.change_extension(seed, "dmnd")
    seed_dmnd = os.path.join(output, 'dmnd', seed_dmnd)
    seed_tsv = io.change_extension(seed, "tsv")
    seed_tsv = os.path.join(output, 'tsv', seed_tsv)
    return seed_fasta, seed_dmnd, seed_tsv

def get_singletons_from_seed(seed, output, thresholds, random_seed_number):
    '''
    Use diamond to identify singletons in the seed genome.
        Arguments:
            seed: path to the seed genbank file
            output: path to the output directory
            thresholds: list of thresholds from the parser
                [args.find, args.minlength, args.maxlength,
                args.presence, args.minloci, args.maxloci]
        Returns:
            candidate_loci:
                List of candidates selected from the seed genome
    '''
    io.make_folder(os.path.join(output, 'tsv'))
    logging.info("Identifying singletons in seed genome...")
    seed_fasta, seed_dmnd, seed_tsv = get_seed_paths(seed, output)
    seed_fasta_contents = io.read_file(seed_fasta)
    diamond.run_diamond_search(seed_fasta, seed_dmnd, seed_tsv)
    unique_loci = get_unique_hits_from_tsv(seed_tsv)
    logging.info("Found %s loci in the seed genome!", str(len(unique_loci)))
    if random_seed_number is None: #random, random if no random seed set
        random.shuffle(unique_loci)
    else: #use the random seed provided
        random.Random(random_seed_number).shuffle(unique_loci)
    loci = 0
    loci_fasta = []
    candidate_loci = []
    loci_to_find, loci_min_length, loci_max_length, _, _, _ = thresholds
    for locus in unique_loci:
        if loci_to_find < 0 or loci < loci_to_find:
            sequence = io.get_locus(seed_fasta_contents, locus)
            if loci_max_length > len(sequence) > loci_min_length:
                candidate_loci.append(locus)
                loci_fasta.append(">" + locus)
                loci_fasta.append(sequence)
                loci += 1
        else:
            break
    txt_path = os.path.join(output, 'tsv/candidate_loci.txt')
    fasta_path = os.path.join(output, 'tsv/candidate_loci.fasta')
    logging.info('%s singletons found in the seed genome!' % loci)
    io.write_to_file(txt_path, candidate_loci)
    io.write_to_file(fasta_path, loci_fasta)
    return candidate_loci

def get_loci_from_file(file: str) -> List:
    '''
    Gets a list of loci from the provided text file.
        Arguments:
            file: the path of the file to be searched
        Returns:
            loci: a list of loci extracted from the text file
    '''
    loci = io.read_file(file)
    loci = [locus.strip() for locus in loci]
    return loci

def search_candidates(output: str, cpus: int) -> None:
    '''
    Uses diamond blastP to search for the candidates in all other genomes.
        Arguments:
            output: path to the output folder
            cpus: the number of cpus avaliable
        Returns:
            None
    '''
    tsvs_folder = os.path.join(output, 'tsvs')
    io.make_folder(tsvs_folder)
    diamond_databases = glob.glob(os.path.join(output, 'dmnd/*.dmnd'))
    args_list = []
    for database in diamond_databases:
        tsv_name = io.change_extension(os.path.basename(database), 'tsv')
        tsv_name = os.path.join(tsvs_folder, tsv_name)
        candidate_loci_path = os.path.join(output, 'tsv/candidate_loci.fasta')
        args_item = [candidate_loci_path, database, tsv_name]
        args_list.append(args_item)
    io.run_in_parallel(diamond.run_diamond_search, args_list, cpus)

def score_locus(locus: str, files: List) -> Tuple[int, bool, List]:
    '''
    Scores the presence and uniqueness of a locus.
        Arguments:
            locus:
                the name of the locus being screened
            files:
                list of paths to blastP results
    '''
    pa_data = []
    presence_counter = 0
    unique_flag = True
    for file in files:
        counter = 0
        lines = io.read_file(file)
        for line in lines:
            if locus in line:
                counter += 1
        if counter > 0:
            presence_counter += 1
        if counter > 1:
            unique_flag = False
        pa_data.append(counter)
    return presence_counter, unique_flag, pa_data

def process_final_loci(final_loci: List, minimum_loci: int, output: str) -> None:
    '''
    Assess that the number of loci meets the user defined threshold to continue analysis
    and write the list of loci to file
        Arguments:
            output: path to the output directory
            final_loci: list of loci names to be used for the alignment
            minimum_loci:
                the minimum number of loci that need to be selected to continue the analysis
        Returns:
            None
    '''
    number_of_loci = len(final_loci)
    logging.info('%s loci selected for MLST!', number_of_loci)
    if number_of_loci < minimum_loci:
        logging.warning("The number of loci below defined threshold. Exiting...")
        raise InsufficientLociError('The number of loci selected are below the defined threshold.')
    filename = os.path.join(output, 'final_loci.txt')
    io.write_to_file(filename, final_loci)

def do_thresholding(
        target_loci: List, presence_threshold: float, maximum_loci: int, output: str
    ) -> List:
    '''
    Apply thresholding to finalise the loci to use for analysis.
    and write presence absence table for each loci.
    Arguments:
        target_loci: list of loci from the seed the genome to compare against other genomes
        presence_threshold:
            the percentage of genomes the loci needs to be present in to be selected for analysis
    Returns:
        None
    '''
    final_loci = []
    pa_table = []
    files = glob.glob(os.path.join(output, 'tsvs/*.tsv'))
    pa_table.append([os.path.splitext(os.path.basename(file))[0] for file in files])
    thresholding_data = ["locus;" + "presence;" + "unique"]
    if len(target_loci) < maximum_loci:
        maximum_loci = len(target_loci)
    for locus in target_loci:
        logging.debug(
            "final = %s, max = %s, targets = %s",
            len(final_loci), maximum_loci, len(target_loci))
        presence, unique, pa_data = score_locus(locus, files)
        pa_table.append(pa_data)
        number_of_loci = len(files)
        presence_percent = (presence / number_of_loci) * 100
        thresholding_string = str(locus) + ";" + str(presence_percent) + ";" + str(unique)
        thresholding_data.append(thresholding_string)
        if presence_percent >= presence_threshold and unique:
            final_loci.append(locus)
        if len(final_loci) >= maximum_loci:
            logging.info('The maximum number of loci was reached.')
            break
    if len(final_loci) <= maximum_loci:
        logging.warning('Number of loci selected is lower than the maximum defined.')
    filename = os.path.join(output, 'thresholding_data')
    io.write_to_file(filename, thresholding_data)
    write_pa_table(pa_table, target_loci, output)
    return final_loci


def threshold_loci(target_loci: List, thresholds: List, output: str) -> List:
    '''
    Score loci for singleton status and presence in dataset
        Arguments:
            target_loci: list of loci from the seed the genome to compare against other genomes
            thresholds: list of thresholds from the parser
                [args.find, args.minlength, args.maxlength,
                args.presence, args.minloci, args.maxloci]
            output: path to the output directory
        Returns:
            final_loci: list of loci names to be used in the downstream analysis
    '''
    _, _, _, presence_threshold, minimum_loci, maximum_loci = thresholds
    filename = os.path.join(output, 'final_loci.txt')
    if os.path.exists(filename):
        logging.error('%s already exists. Exiting!', filename)
        raise FileAlreadyExistsError(
            'File already exists. Please remove and restart from checkpoint.', filename
            )
    final_loci = do_thresholding(target_loci, presence_threshold, maximum_loci, output)
    process_final_loci(final_loci, minimum_loci, output)
    return final_loci

def write_pa_table(pa_table: List, loci: List, output: str) -> None:
    '''
    Transposes and writes the presence absence table to the specified folder
        Arguments:
            pa_table:
                list of lists, list is the presence absence data corresponding to a single locus
            loci: list of loci names
            output: path of the output directory
        Returns:
            None
    '''
    new_table = []
    header = ['strain']
    header.extend(loci)
    new_table.append(";".join(header))
    transposed_data = list(zip(*pa_table))
    for data in transposed_data:
        data_string = ';'.join([str(datum) for datum in data])
        new_table.append(data_string)
    filename = os.path.join(output, 'presence_absence_table.csv')
    io.write_to_file(filename, new_table)

def get_target_proteins(
        checkpoint: Checkpoint, output: str, seed: str, thresholds: List,
        cpus: int, random_seed_number: int
    ) -> None:
    '''
    The main routine for screen.py
        Arguments:
            checkpoint: the checkpoint provided by the user
            output: path of the output directory
            seed: the name of the file corresponding to the seed genome
            thresholds: list of arguments containing threholding information:
                [args.find, args.minlength, args.maxlength,
                args.presence, args.minloci, args.maxloci]
        Returns:
            None
    '''
    candidate_loci = None
    final_loci = None
    logging.debug('The output directory is: %s', output)
    if checkpoint < Checkpoint.SINGLETONS_IDENTIFIED:
        candidate_loci = get_singletons_from_seed(seed, output, thresholds, random_seed_number)
    logging.info("CHECKPOINT: SINGLETONS_IDENTIFIED")
    #candidate loci will not exist if restarted from a checkpoint
    if not candidate_loci:
        try:
            candidate_loci = get_loci_from_file(os.path.join(output, 'tsv/candidate_loci.txt'))
        except Exception:
            raise NoCandidateLociError(
                'Candidate loci could not be read from candidate_loci.txt.'
                'If restarting from a checkpoint ensure there is a final_loci.txt '
                'file in the specified output folder.'
            )
    #continue sequential analysis
    if checkpoint < Checkpoint.SINGLETONS_SEARCHED:
        logging.info("Screening candidate loci against other genomes...")
        search_candidates(output, cpus)
    logging.info("CHECKPOINT: SINGLETONS_SEARCHED")
    if checkpoint < Checkpoint.SINGLETONS_THRESHOLDED:
        logging.info("Thresholding candidate loci...")
        final_loci = threshold_loci(candidate_loci, thresholds, output)
    logging.info("CHECKPOINT: SINGLETONS_THRESHOLDED")
    return final_loci
