'''
Build fasta and diamond databases from genbank files

Functions:
build_diamond_databases(output:str) -> None:
extract_cdses(
    gbks: str, output: str, tag_label: str, ignore_bad_annotations: bool, ignore_bad_records: bool
    ) -> None:
get_cds_from_genbank(
    filename: str, output: str, tag_label: str, ignore_bad_annotations: bool
    ) -> None:
extract_data(
    checkpoint: Checkpoint, output: str, gbks: str, tag_label: str,
    ignore_bad_annotations: bool, ignore_bad_records: bool
    ) -> None:
'''
import logging
import os
import glob
from getphylo.ext import diamond
from getphylo.utils import io
from getphylo.utils.checkpoint import Checkpoint
from getphylo.utils.errors import GetphyloError, BadAnnotationError, BadRecordError

def build_diamond_databases(output: str) -> None:
    '''Create diamond databases from from all the fasta files in .output/fasta/*.fasta
        Arguments:
            output: path to the output folder
        Returns:
            None'''
    logging.info("Creating diamond databases from extracted cdses...")
    dmnd_folder = os.path.join(output, 'dmnd')
    io.make_folder(dmnd_folder)
    fasta_files_path = os.path.join(output, 'fasta/*.fasta')
    for file in glob.glob(fasta_files_path):
        dmnd_database = os.path.basename(io.change_extension(file, "dmnd"))
        dmnd_database = os.path.join(dmnd_folder, dmnd_database)
        diamond.make_diamond_database(file, dmnd_database)

def extract_cdses(
        gbks: str, output: str, tag_label: str,
        ignore_bad_annotations: bool, ignore_bad_records: bool
    ) -> None:
    '''Produce a fasta file from each genbank provided
        Arguments:
            gbks: search string for genbank files
            output: path to the output directory
            tag_args:
        Returns:
            None
        '''
    io.make_folder(os.path.join(output, 'fasta'))
    files_valid = True
    for filename in glob.glob(gbks):
        logging.info('Extracting CDS annotations from %s', filename)
        try:
            get_cds_from_genbank(filename, output, tag_label, ignore_bad_annotations)
        except GetphyloError as error:
            logging.warning('%s raised error: %s', filename, error)
            files_valid = False

    if not files_valid:
        if ignore_bad_records:
            raise BadAnnotationError(
                'Some files were not parsed correctly, check the logging file for more information.'
                'If you want to skip these files,'
                'rerun the analysis with the --ignore-bad-records flag.'
                )

def get_cds_from_genbank(
        filename: str, output: str, tag_label: str, ignore_bad_annotations: bool
    ) -> None:
    '''Extract CDS translations from genbank files into ./fasta/*.fasta
        Arguments:
            filename: the name of the genbank file being read
            output: path to the output folder
            tag_label: the string defining the tag label (e.g. 'locus_tag')
            ignore_bad_annotations:
                bool flagging whether to ignore features with missing annotations
        Returns: None
    '''
    lines = []
    seen = set()

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
                        lines.append(str(feature.qualifiers.get("translation")[0]))
                except TypeError:
                    if feature.qualifiers.get(tag_label) is None:
                        logging.warning(
                            'Missing %s in %s.', tag_label, record.id #feature id?
                            )
                        if not ignore_bad_annotations:
                            raise BadAnnotationError(
                                f'Some features are missing the {tag_label} annotations.'
                                'Ensure the genbank file is correctly annotated or use'
                                '--ignore-bad-annotations flag.'
                            )
                    else:
                        logging.warning(
                            '%s has no translation!', feature.qualifiers.get(tag_label)[0]
                            #only show once!
                            )
    except ValueError as error:
        raise BadRecordError(error)
    if not lines:
        raise BadAnnotationError(f'No CDS Features in {filename}')
    filename = io.change_extension(os.path.basename(filename), "fasta")
    filename = os.path.join(output, 'fasta', filename)
    io.write_to_file(filename, lines)

def extract_data(
        checkpoint: Checkpoint, output: str, gbks: str, tag_label: str,
        ignore_bad_annotations: bool, ignore_bad_records: bool
    ) -> None:
    '''Called from main to build fasta and diamond databases from the provided genbankfiles
        Arguments:
            checkpoint: Checkpoint to begin the analysis
            output: path to the output folder
            gbks: search string for the genbank files
            tag_label: the string defining the tag label (e.g. 'locus_tag')
            ignore_bad_annotations:
                bool flagging whether to ignore features with missing annotations
        Returns: None'''
    if checkpoint < Checkpoint.FASTA_EXTRACTED:
        logging.info("Checkpoint: Extracting CDSs...")
        extract_cdses(gbks, output, tag_label, ignore_bad_annotations, ignore_bad_records)
    if checkpoint < Checkpoint.DIAMOND_BUILT:
        logging.info("Checkpoint: Building diamond databases...")
        build_diamond_databases(output)
