'''
Create alignments from a folder of fasta files.

Functions:
    get_locus_from_tsv(locus: str, fasta_name: str) -> Tuple[str, str]
    make_fasta_for_alignments(loci_list: List[str], output:str) -> None
    do_alignments(output:str) -> None
    get_locus_length(alignment: List[str]) -> int
    get_locus_alignment(alignment, taxon_name) -> str
    make_combined_alignment(gbks: List, output: str) -> None
    make_alignments(checkpoint: Checkpoint, output: str, loci: List, gbks: List) -> None
'''
import logging
import os
import glob
from typing import List, Tuple
from getphylo.ext import muscle
from getphylo.utils import io
from getphylo.utils.checkpoint import Checkpoint
from getphylo.utils.errors import FileAlreadyExistsError, BadLocusError

def get_locus_from_tsv(locus: str, fasta_name: str) -> Tuple[str, str]:
    '''
    Extracts specific loci for the alignment.
        Arguments:
            locus: locus to extract
            fasta_name: the name of the .fasta file
        Returns:
            sequence_name: the locus and taxon name
            sequence: the sequence for aligning
    '''
    tsv_name = fasta_name.replace("/fasta/", "/tsvs/")
    tsv_name = io.change_extension(tsv_name, "tsv")
    tsv = io.read_tsv(tsv_name)
    fasta_contents = io.read_file(fasta_name)
    for line in tsv:
        if locus in line[0]:
            sequence = io.get_locus(fasta_contents, line[1])
            organism, _ = os.path.splitext(os.path.basename(fasta_name))
            sequence_name = f'>{organism}_{line[1]}'
            return sequence_name, sequence
    raise BadLocusError("Locus %s not found in file: %s" % (locus, tsv_name))

def make_fasta_for_alignments(loci_list: List[str], output: str) -> None:
    '''
    Builds a .fasta file from sequences where there is a hit in the diamond search results.
    Arguments:
        loci_list: a list of locus tags to extract hits
        output: the path of the output directory
    Returns:
        None
    '''
    io.make_folder(os.path.join(output, 'unaligned_fasta'))
    files = glob.glob(os.path.join(output, 'fasta/*.fasta'))
    assert files
    logging.debug(loci_list)
    for locus in loci_list:
        write_lines = []
        for filename in files:
            try:
                logging.debug('locus: %s; filename: %s' % (locus, filename))
                write_lines.extend(get_locus_from_tsv(locus, filename))
                logging.debug('write_lines %s', write_lines)
            except BadLocusError as error:
                logging.warning(error)

        outfile = os.path.join(output, 'unaligned_fasta', locus + '.fasta')
        io.write_to_file(outfile, write_lines)

def do_alignments(output: str, cpus: int) -> None:
    '''
    Runs the pre-aligned fasta file through the muscle module.
        Arguments:
            output: the path of the outut directory
        Returns:
            None
    '''
    io.make_folder(os.path.join(output, 'aligned_fasta'))
    args_list = []
    for filename in glob.glob(os.path.join(output, 'unaligned_fasta/*.fasta')):
        outfile = os.path.join(output, 'aligned_fasta', os.path.basename(filename))
        args_list.append([filename, outfile])
    io.run_in_parallel(muscle.run_muscle, args_list, cpus) #add cpu arg

def get_locus_length(alignment: List[str]) -> int:
    '''
    Returns the length of the fasta locus by counting sequence lines.
        Arguments:
            alignment: a list containing each line of a fasta file
        Returns:
            alignment_length: an int reprisenting the alignment length
    '''
    alignment_length = 0
    if not alignment:
        raise ValueError('An alignment cannot be empty.')
    for line in alignment[1:]:
        if '>' not in line:
            alignment_length += len(line.strip())
        else:
            break
    return alignment_length

def get_locus_alignment(alignment, taxon_name) -> str:
    '''
    Extract an alignment for a loci for a taxa.
        Arguments:
            alignment: List of lines representing the alignment fasta
            taxon_name: the name of the taxon to extract lines for
        Returns:
            extracted_alignment_string:
                the alignment of the specific loci and taxa combined as a string
    '''
    seq_flag = False
    extracted_alignment = []
    for line in alignment:
        if line[0] == '>':
            seq_flag = False
        if '>' in line and taxon_name in line:
            seq_flag = True
        if seq_flag and taxon_name not in line:
            extracted_alignment.extend(line[:-1])
    extracted_alignment_string = ''.join(extracted_alignment)
    return extracted_alignment_string

def format_partition_data(partition_data: List) -> List:
    '''
    Takes the partition data [locus, locus_length] and reformats it as a partition file 
    (e.g. p1 = 1-50, p2 = 51-110)
        Arguments:
            partition_data: a list of [locus, locus_length] lists
        Returns:
            partition_lines: list of lines in the partition format
    '''
    model = 'WAG'
    partition_lines = []
    partition_start = 1
    for partition in partition_data:
        locus, length = partition[0], partition[1]
        locus = os.path.splitext(os.path.basename(locus))[0]
        partition_end = partition_start + length - 1
        partition_lines.append('%s, %s = %s-%s' % (model, locus, partition_start, partition_end))
        partition_start += length
    return partition_lines

def make_combined_alignment(gbks: List, output: str) -> None:
    '''
    Create a combined alignment from single locus alignments
        Arguments:
            gbks: a list of genbank files
            output: path  to output directory
        Returns:
            None
    '''
    combined_alignment_path = os.path.join(output, 'aligned_fasta/combined_alignment.fasta')
    partition_path = os.path.join(output, 'partition.txt')
    if os.path.exists(combined_alignment_path):
        logging.error(
            'ALERT: aligned_fasta/combined_alignment.fasta already exists. Exiting!'
            )
        raise FileAlreadyExistsError('%s alread exists.' % combined_alignment_path)
    combined_alignment = []
    partition_data = []
    taxa = glob.glob(gbks)
    assert taxa, gbks
    loci = glob.glob(os.path.join(output, 'aligned_fasta/*.fasta'))
    first_taxa = True
    for taxon in taxa:
        taxon_name = os.path.splitext(os.path.basename(taxon))[0]
        sequence_data = []
        for locus in loci:
            alignment = io.read_file(locus)
            locus_length = get_locus_length(alignment)
            if first_taxa is True:
                partition_data.append([locus, locus_length])
            locus_alignment = get_locus_alignment(alignment, taxon_name)
            if len(locus_alignment) == locus_length:
                sequence_data.append(locus_alignment)
            else:
                sequence_data.append('?' * locus_length)
        sequence_string = ''.join(sequence_data)
        assert sequence_data
        if sequence_string.count('?') == len(sequence_string):
            logging.error('[ALERT]: %s has no sequence data and has been removed.', taxon_name)
        else:
            combined_alignment.append(f'>{taxon_name}')
            combined_alignment.append(sequence_string)
        first_taxa = False
    partition_data = format_partition_data(partition_data)
    io.write_to_file(combined_alignment_path, combined_alignment)
    io.write_to_file(partition_path, partition_data)

def make_alignments(checkpoint: Checkpoint, output: str, loci: List, gbks: List, cpus: int) -> None:
    '''
    Main routine for align.
        Arguments:
            checkpoint: checkpoint defined by the user
            output: path to the output directory
            loci: list of loci to align
            gbks: list of the input genbank files
        Returns:
            None
    '''
    if checkpoint < Checkpoint.SINGLETONS_EXTRACTED:
        logging.info("Extracting sequences for alignment...")
        make_fasta_for_alignments(loci, output)
    logging.info("CHECKPOINT: SINGLETONS_EXTRACTED")
    if checkpoint < Checkpoint.SINGLETONS_ALIGNED:
        logging.info("Aligning sequences...")
        do_alignments(output, cpus)
    logging.info("CHECKPOINT: SINGLETONS_ALIGNED")
    if checkpoint < Checkpoint.ALIGNMENTS_COMBINED:
        logging.info("Making combined alingnment...")
        make_combined_alignment(gbks, output)
    logging.info("CHECKPOINT: ALIGNMENTS_COMBINED")
