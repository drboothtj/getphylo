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
        database = " --db " + filename.split('.')[0] + ".dmnd"
    else:
        database = " --db " + dmnd_database
    makedb = " makedb"
    input_ = " --in " + filename + ""
    command = diamond_location + makedb + database + input_
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
        database = " --db " + filename.split('.')[0] + ".dmnd"
    else:
        database = " --db " + dmnd_database
    if outname is None:
        output = " --out " + dmnd_database.split('.')[0] + ".tsv"
    else:
        output = " --out " + outname
    diamond = diamond_args[0]
    blastp = " blastp"
    query = " --query " + filename
    outfmt = " --outfmt 6 qseqid sseqid pident"
    command = diamond + blastp + database + query + output + outfmt
    #there must be a nicer way of doing this but it works for now!
    if diamond_args[1] is not None:
        identity = " --id " + str(diamond_args[1])
        command += identity
    if diamond_args[2] is not None:
        q_cover = " --query_cover " + str(diamond_args[2])
        command += q_cover
    if diamond_args[3] is not None:
        s_cover = " --subject_cover " + str(diamond_args[3])
        command += s_cover
    logging.debug(command)
    io.run_in_command_line(command)
