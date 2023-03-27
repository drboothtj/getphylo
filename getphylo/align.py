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
    for line in tsv:
        if locus in line[0]:
            sequence = io.get_locus(fasta_name, line[1])
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

def do_alignments(output: str) -> None:
    '''
    Runs the pre-aligned fasta file through the muscle module.
        Arguments:
            output: the path of the outut directory
        Returns:
            None
    '''
    io.make_folder(os.path.join(output, 'aligned_fasta'))
    for file in glob.glob(os.path.join(output, 'unaligned_fasta/*.fasta')):
        outfile = os.path.join(output, 'aligned_fasta', os.path.basename(file))
        muscle.run_muscle(file, outfile)

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

def make_combined_alignment(gbks: List, output: str) -> None:
    '''Create a combined alignment from single locus alignments
        Arguments:
            gbks: a list of genbank files
            output: path  to output directory
        Returns:
            None
        '''
    combined_alignemnt_path = os.path.join(output, 'aligned_fasta/combined_alignment.fasta')
    if os.path.exists(combined_alignemnt_path):
        logging.error(
            'ALERT: aligned_fasta/combined_alignment.fasta already exists. Exiting!'
            )
        raise FileAlreadyExistsError('%s alread exists.' % combined_alignemnt_path)
    else:
        combined_alignment = []
        taxa = glob.glob(gbks)
        loci = glob.glob(os.path.join(output, 'aligned_fasta/*.fasta'))
        for taxon in taxa:
            taxon_name = os.path.splitext(taxon)[0]
            sequence_data = []
            for locus in loci:
                alignment = io.read_file(locus)
                locus_length = get_locus_length(alignment)
                locus_alignment = get_locus_alignment(alignment, taxon_name)
                if len(locus_alignment) == locus_length:
                    sequence_data.append(locus_alignment)
                else:
                    sequence_data.append('X' * locus_length)
            sequence_string = ''.join(sequence_data)
            if sequence_string.count('X') == len(sequence_string):
                logging.error('[ALERT]: %s has no sequence data and has been removed.', taxon_name)
                #make fatal and add option to ignore
            else:
                combined_alignment.append(f'>{taxon_name}')
                combined_alignment.append(sequence_string)
        io.write_to_file(combined_alignemnt_path, combined_alignment)
        #provide partition data!

def make_alignments(checkpoint: Checkpoint, output: str, loci: List, gbks: List) -> None:
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
        logging.info("Checkpoint: Extracting sequences for alignment...")
        make_fasta_for_alignments(loci, output)
    if checkpoint < Checkpoint.SINGLETONS_ALIGNED:
        logging.info("Checkpoint 6: Aligning sequences...")
        do_alignments(output)
    if checkpoint < Checkpoint.ALIGNMENTS_COMBINED:
        logging.info("Checkpoint 7: Combining alignments")
        make_combined_alignment(gbks, output)
