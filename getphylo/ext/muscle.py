'''
Runs MUSCLE on a provided fasta file.

Functions:
    run_muscle(filename, outname=None) -> None
'''
from getphylo.utils import io

def run_muscle(filename: str, outname=None) -> None:
    '''
    Run MUSCLE aligner on protein fasta file.
        Arguments:
            filename: path to unaligned sequences
            outname: path for the alignment
        Returns:
            None'''
    if outname is None:
        out = " -out aligned_" + filename
    else:
        out = " -out " + outname
    _in = " -in " + filename
    command = "muscle" + _in + out
    io.run_in_command_line(command)
