'''
Run iqtree.

Functions:
    run_iqtree(alignment_path, partition_path, outfile=None) -> None:

'''
from getphylo.utils import io

def run_iqtree(alignment_path, partition_path) -> None:
    '''
    Run fasttree on a protein alignment.
        Arguments:
            alignment_path: path to the alignment
            partition_path: path to the partition file
            outfile: path to the output file
        Returns:
            None
    '''
    command = 'iqtree'
    alignment = " ".join(['-s', alignment_path])
    partition = " ".join(['-spp', partition_path])
    model = '-m MFP'
    bootstraps = '-bb 1000'

    io.run_in_command_line(" ".join([command, alignment, partition, model, bootstraps]))
