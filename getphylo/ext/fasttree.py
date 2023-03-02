'''
Run fasttree.

Functions:
    run_fasttree(filename, outfile=None) -> None

'''
from getphylo.utils import io

def run_fasttree(filename, outfile=None) -> None:
    '''
    Run fasttree on a protein alignment.
        Arguments:
            filename: path to the alignment
            outfile: path to the output file
        Returns:
            None
    '''
    if outfile is None:
        out = " -out " + io.change_extension(filename, "tree") + " "
    else:
        out = " -out " + outfile + " "
    command = "fasttree"
    io.run_in_command_line(command + out + filename)
