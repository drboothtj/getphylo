'''
Build fasta and diamond databases from genbank files

Functions:
build_diamond_databases(output:str) -> None
extract_cdses(
    gbks: str, output: str, tag_label: str, ignore_bad_annotations: bool, ignore_bad_records: bool
    ) -> None
get_cds_from_genbank(
    filename: str, output: str, tag_label: str, ignore_bad_annotations: bool,
    ignore_bad_records: bool) -> None
extract_data(
    checkpoint: Checkpoint, output: str, gbks: str, tag_label: str,
    ignore_bad_annotations: bool, ignore_bad_records: bool
    ) -> None
'''
import logging
import os
import glob
from getphylo.ext import diamond
from getphylo.utils import io
from getphylo.utils.checkpoint import Checkpoint
from getphylo.utils.errors import BadAnnotationError, BadRecordError

def build_diamond_databases(output: str, cpus: int) -> None:
    '''
    Create diamond databases from from all the fasta files in .output/fasta/*.fasta
        Arguments:
            output: path to the output folder
            cpus: number of cpus available
        Returns:
            None
    '''
    dmnd_folder = os.path.join(output, 'dmnd')
    io.make_folder(dmnd_folder)
    fasta_files_path = os.path.join(output, 'fasta/*.fasta')
    args_list = []
    for filename in glob.glob(fasta_files_path):
        dmnd_database = os.path.basename(io.change_extension(filename, "dmnd"))
        dmnd_database = os.path.join(dmnd_folder, dmnd_database)
        args = [filename, dmnd_database]
        args_list.append(args)
    io.run_in_parallel(diamond.make_diamond_database, args_list, cpus)

def extract_cdses(
        gbks: str, output: str, tag_label: str,
        ignore_bad_annotations: bool, ignore_bad_records: bool, cpus: int
    ) -> None:
    '''
    Produce a fasta file from each genbank provided
        Arguments:
            gbks: search string for genbank files
            output: path to the output directory
            tag_args:
            cpus: number of cpus available
        Returns:
            None
    '''
    io.make_folder(os.path.join(output, 'fasta'))
    filenames = glob.glob(gbks)
    args_list = [[
        filename, output, tag_label, ignore_bad_annotations, ignore_bad_records
        ] for filename in filenames]
    io.run_in_parallel(get_cds_from_genbank, args_list, cpus)

def get_cds_from_genbank(
        filename: str, output: str, tag_label: str, ignore_bad_annotations: bool,
    ignore_bad_records: bool) -> None:
    '''
    Extract CDS translations from genbank files into ./fasta/*.fasta
        Arguments:
            filename: the name of the genbank file being read
            output: path to the output folder
            tag_label: the string defining the tag label (e.g. 'locus_tag')
            ignore_bad_annotations:
                bool flagging whether to ignore features with missing annotations
        Returns: None
    '''
    logging.debug('Extracting CDS annotations from %s', filename)
    lines = []
    seen = set()
    warning_flag = False
    try:
        records = io.get_records_from_genbank(filename)
        for record in records:
            for feature in record.features:
                try:
                    if feature.type == "CDS":
                        locus_tag = f'{record.id}_{feature.qualifiers.get(tag_label)[0]}'
                        if locus_tag in seen:
                            raise BadAnnotationError(f'{filename} contains duplicate: {locus_tag}')
                        seen.add(locus_tag)
                        lines.append(">" + locus_tag.replace(".", "_"))
                        translation = str(feature.qualifiers.get("translation")[0])
                        if translation == "":
                            raise BadAnnotationError(
                                f'{locus_tag} in {filename} contains an empty translation!'
                                )
                        else:
                            lines.append(translation)
                except TypeError:
                    if feature.qualifiers.get(tag_label) is None:
                        logging.warning(
                            'Missing %s in %s.', tag_label, record.id
                            )
                        if not ignore_bad_annotations:
                            raise BadAnnotationError(
                                f'Some features are missing the {tag_label} annotations.'
                                'Ensure the genbank file is correctly annotated or use'
                                '--ignore-bad-annotations flag.'
                            )
                    else:
                        warning_flag = True
            if warning_flag is True:
                logging.warning('%s has missing translations!', record.id)
    except ValueError as error:
        if not ignore_bad_records:
            raise BadRecordError(error)
    if not lines:
        if not ignore_bad_annotations:
            raise BadAnnotationError(f'No CDS Features in {filename}')
    filename = io.change_extension(os.path.basename(filename), "fasta")
    filename = os.path.join(output, 'fasta', filename)
    io.write_to_file(filename, lines)

def extract_data(
        checkpoint: Checkpoint, output: str, gbks: str, tag_label: str,
        ignore_bad_annotations: bool, ignore_bad_records: bool, cpus: int
    ) -> None:
    '''
    Called from main to build fasta and diamond databases from the provided genbankfiles
        Arguments:
            checkpoint: Checkpoint to begin the analysis
            output: path to the output folder
            gbks: search string for the genbank files
            tag_label: the string defining the tag label (e.g. 'locus_tag')
            ignore_bad_annotations:
                bool flagging whether to ignore features with missing annotations
        Returns: None
    '''
    if checkpoint < Checkpoint.FASTA_EXTRACTED:
        logging.info("Extracting CDS in fasta format...")
        extract_cdses(gbks, output, tag_label, ignore_bad_annotations, ignore_bad_records, cpus)
    logging.info("CHECKPOINT:FASTA_EXTRACTED")
    if checkpoint < Checkpoint.DIAMOND_BUILT:
        logging.info("Building diamond databases...")
        build_diamond_databases(output, cpus)
    logging.info("CHECKPOINT:DIAMOND_BUILT")
