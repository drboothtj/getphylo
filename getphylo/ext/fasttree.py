'''
Run fasttree.

Functions:


'''
from getphylo.utils import io

def run_fasttree(filename, outfile=None):
    '''Run fasttree on a protein alignment'''
    if outfile is None:
        out = " -out " + io.change_extension(filename, "tree") + " "
    else:
        out = " -out " + outfile + " "
    command = "fasttree"
    io.run_in_command_line(command + out + filename)
