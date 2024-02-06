'''
Run iqtree.

Functions:
    run_iqtree(alignment_path: str, partition_path: str, out_path: str) -> None

'''
from getphylo.utils import io

def run_iqtree(alignment_path: str, out_path: str, partition_path: str=None) -> None:
    '''
    Run fasttree on a protein alignment.
        Arguments:
            alignment_path: path to the alignment
            partition_path: path to the partition file
            out_path: path to the output file
        Returns:
            None
    '''
    command = 'iqtree'
    alignment = " ".join(['-s', alignment_path])
    out_path = " ".join(['-pre', out_path])
    model = '-m MFP'
    bootstraps = '-bb 1000'
    #do a run with or without partition file
    if partition_path is not None:
        partition = " ".join(['-spp', partition_path])
        io.run_in_command_line(" ".join([command, alignment, partition, out_path, model, bootstraps]))
    else:
        io.run_in_command_line(" ".join([command, alignment, out_path, model, bootstraps]))
