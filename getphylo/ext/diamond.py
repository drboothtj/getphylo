'''
Make DIAMOND databases and run searches using DIAMOND

Functions:
    make_diamond_database(filename: str, dmnd_database=None) -> None
    run_diamond_search(filename: str, dmnd_database=None, outname=None) -> None
'''
import logging
from getphylo.utils import io

def make_diamond_database(filename: str, dmnd_database=None) -> None:
    '''
    Create a DIAMOND database from a fasta file.
        Arguments:
            filename: path to the input file
            dmnd_database: path to output directory
        Returns:
            None
    '''
    if dmnd_database is None:
        database = " --db " + filename.split('.')[0] + ".dmnd"
    else:
        database = " --db " + dmnd_database
    makedb = "diamond makedb"
    input_ = " --in " + filename + ""
    command = makedb + database + input_
    logging.debug(command)
    io.run_in_command_line(command)

def run_diamond_search(filename: str, dmnd_database=None, outname=None) -> None:
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
        database = " --db " + filename.split('.')[0] + ".dmnd"
    else:
        database = " --db " + dmnd_database
    if outname is None:
        output = " --out " + dmnd_database.split('.')[0] + ".tsv"
    else:
        output = " --out " + outname
    blastp = "diamond blastp"
    query = " --query " + filename
    outfmt = " --outfmt 6 qseqid sseqid pident"
    command = blastp + database + query + output + outfmt
    logging.debug(command)
    io.run_in_command_line(command)
