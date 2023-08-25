'''
Create alignments from a folder of fasta files.

Functions:
    get_locus_from_tsv(locus: str, fasta_name: str) -> Tuple[str, str]
    make_fasta_for_alignments(loci_list: List[str], output:str) -> None
    do_alignments(output:str) -> None
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

def read_fasta(filename: str) -> dict[str, str]:
    '''
    Reads a fasta file into a dictionary
        Arguments:
            filename: the path to the FASTA file to read
        Returns:
            a dictionary mapping sequence ID to sequence
    '''
    pairs = {}
    with open(filename, "r", encoding="utf-8") as fasta:
        current_seq: list[str] = []
        current_id: str = ""
        for line in fasta:
            line = line.strip()
            if not line:
                continue
            if line[0] == '>':
                if current_id and current_seq:
                    pairs[current_id] = "".join(current_seq)
                current_id = line[1:].replace(" ", "_")
                current_seq.clear()
            else:
                if not current_id:
                    raise ValueError("Sequence before identifier in fasta file")
                current_seq.append(line)
    if current_id and current_seq:
        pairs[current_id] = "".join(current_seq)
    if not pairs:
        raise ValueError("Fasta file contains no sequences")
    return pairs


def make_combined_alignment(gbks: str, output: str) -> None:
    '''
    Create a combined alignment from single locus alignments
        Arguments:
            gbks: a path of genbank files (e.g. 'inputs/*.gbk')
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
        raise FileAlreadyExistsError('%s already exists.' % combined_alignment_path)
    combined_alignment = []
    partition_data = []
    taxa = glob.glob(gbks)
    assert taxa, gbks

    # create a system to work out the taxa name from a fasta identifier
    taxa_by_locus_tag = {}
    all_taxa = set()
    for t in taxa:
        taxa_name = os.path.splitext(os.path.basename(t))[0]
        all_taxa.add(taxa_name)
        assert t.endswith("fasta"), "handle more file types later"
        try:
            content = read_fasta(t)
        except ValueError:
            print(t)
            raise

        for l in content:
            l = f"{taxa_name}_{l}"
            taxa_by_locus_tag[l] = taxa_name

    assert "GCF_008931305.1_ASM893130v1_genomic_NZ_CP042324_1_FQ762_RS02160" in taxa_by_locus_tag

    # every fasta header in aligned fastas is a unique combo of filename, record name, and locus tag
    loci = glob.glob(os.path.join(output, 'aligned_fasta/*.fasta'))
    # read all the fastas and get the maximum length of any of them
    alignments: list[dict[str, str]] = []
    max_seq_length = 0
    for locus in loci:
        # load the file
        this_alignment = read_fasta(locus)
        # we only care about the taxa each sequence came from, not the full identifier
        alignments.append({taxa_by_locus_tag[k]: v for k, v in this_alignment.items()})
        assert len(alignments[-1]) == len(this_alignment)  # and no taxa present multiple times in an alignment
        # each fasta better have the same length for all sequences
        seqs = list(this_alignment.values())
        length = len(seqs[0])  # should be the same for sequence in this alignment
        assert all(len(s) == length for s in seqs)
        # if it's longer than any other file's sequence length, update the max
        if length > max_seq_length:
            max_seq_length = length

    # find all taxas in any alignments
    taxas_in_any_alignment = set()
    for alignment in alignments:
        taxas_in_any_alignment.update(alignment)

    # for any taxa, if there's no member of it present in any alignment, warn the user
    missing_taxa = all_taxa - taxas_in_any_alignment
    for taxon_name in sorted(missing_taxa):
        logging.error('[ALERT]: %s was not present in any alignment and has been removed.', taxon_name)

    # build a combined fasta of all alignments, padding the end of shorter sequences with '?'
    # so that they're all the same length
    for taxon_name in taxas_in_any_alignment:
        for alignment in alignments:
            if taxon_name not in alignment:
                alignment[taxon_name] = len(list(alignment.values())[0]) * "?"

    # every taxa should now be present in every alignment
    assert all(len(alignment) == len(taxas_in_any_alignment) for alignment in alignments)

    # concatenate the alignments
    concatenated = {}
    for taxon_name in taxas_in_any_alignment:
        line = []
        for alignment in alignments:
            line.append(alignment[taxon_name])
        concatenated[taxon_name] = "".join(line)

    combined_alignment = [f">{taxon_name}\n{seq}" for taxon_name, seq in concatenated.items()]

    io.write_to_file(combined_alignment_path, combined_alignment)
    # TODO
    # partition_data = format_partition_data(partition_data)
#    io.write_to_file(partition_path, partition_data)

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
