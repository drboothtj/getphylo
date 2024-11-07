'''
Make DIAMOND databases and run searches using DIAMOND

Functions:
    make_diamond_database(filename: str, dmnd_database=None) -> None
    run_diamond_search(
    filename: str, dmnd_database=None, outname=None, diamond_args=['diamond',None,None,None]
    ) -> None
'''
import logging
from getphylo.utils import io

def make_diamond_database(infile: str, dmnd_database=None, diamond_location='diamond') -> None:
    '''
    Create a DIAMOND database from a fasta file.
        Arguments:
            filename: path to the input file
            dmnd_database: path to output directory
        Returns:
            None
    '''
    if dmnd_database is None:
        database_name = infile.split('.')[0] + ".dmnd"
    else:
        database_name = dmnd_database
    command = [
        diamond_location, "makedb",
        "--db", database_name, 
        "--in", infile
        ]
    logging.debug(command)
    io.run_in_command_line(command)

def run_diamond_search(
    filename: str, dmnd_database=None, outname=None, diamond_args=['diamond',None,None,None]
    ) -> None:
    '''
    Run BLASTP through DIAMOND.
        Arguments:
            filename: path to the input fasta file
            dmnd_database: name of a database file to read
            outname: name of the output file
            diamond_args: list of arguments for diamond
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
        diamond_args[0], "blastp",
        "--db", database,
        "--query", filename,
        "--out", output,
        "--outfmt", "6", "qseqid", "sseqid", "pident"
        ]
    #there must be a nicer way of doing this but it works for now!
    if diamond_args[1] is not None:
        command.append("--id")
        command.append(str(diamond_args[1]))
    if diamond_args[2] is not None:
        command.append("--query-cover")
        command.append(str(diamond_args[2]))
    if diamond_args[3] is not None:
        command.append("--subject-cover")
        command.append(str(diamond_args[3]))
    logging.debug(command)
    io.run_in_command_line(command)
