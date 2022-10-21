'''
Runs MUSCLE on a provided fasta file.

Functions:

'''
from getphylo import console

def run_muscle(filename, outname=None):
    '''Run MUSCLE aligner on protein fasta file'''
    if outname is None:
        out = " -out aligned_" + filename
    else:
        out = " -out " + outname
    _in = " -in " + filename
    command = "muscle" + _in + out
    console.run_in_command_line(command)
