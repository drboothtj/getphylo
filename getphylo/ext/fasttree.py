'''
Run fasttree.

Functions:
    run_fasttree(filename, outfile=None, fasttree_location='fasttree') -> None

'''
from getphylo.utils import io

def run_fasttree(filename, outfile=None, fasttree_location='fasttree') -> None:
    '''
    Run fasttree on a protein alignment.
        Arguments:
            filename: path to the alignment
            outfile: path to the output file
        Returns:
            None
    '''

    if outfile is None:
        out = io.change_extension(filename, "tree")
    else:
        out = outfile
    command = [
        fasttree_location,
        "-out", out,
        filename
    ]
    io.run_in_command_line(command)
