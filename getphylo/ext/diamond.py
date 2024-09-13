'''
Make DIAMOND databases and run searches using DIAMOND

Functions:
    make_diamond_database(filename: str, dmnd_database=None) -> None
    run_diamond_search(filename: str, dmnd_database=None, outname=None) -> None
'''
import logging
from getphylo.utils import io

def make_diamond_database(filename: str, dmnd_database=None, diamond_location='diamond') -> None:
    '''
    Create a DIAMOND database from a fasta file.
        Arguments:
            filename: path to the input file
            dmnd_database: path to output directory
        Returns:
            None
    '''
    if dmnd_database is None:
        database = filename.split('.')[0] + ".dmnd"
    else:
        database = dmnd_database
    command = [
        diamond_location, "makedb",
        "--db", database,
        "--in", filename,
    ]
    logging.debug(command)
    io.run_in_command_line(command)

def run_diamond_search(
    filename: str, dmnd_database=None, outname=None, diamond_location='diamond'
    ) -> None:
    '''
    Run BLASTP through DIAMOND.
        Arguments:
            filename: path to the input fasta file
            dmnd_database: name of a database file to read
            outname: name of the output file
        Returns:
            None
    '''
    if dmnd_database is None:
        database = filename.split('.')[0] + ".dmnd"
    else:
        database = dmnd_database
    if outname is None:
        output = dmnd_database.split('.')[0] + ".tsv"
    else:
        output = outname
    command = [
        diamond_location, "blastp",
        "--query", filename,
        "--db", database,
        "--out", output,
        "--outfmt", "6", "qseqid", "sseqid", "pident",
    ]
    logging.debug(command)
    io.run_in_command_line(command)
