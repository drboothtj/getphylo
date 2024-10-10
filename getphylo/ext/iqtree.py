'''
Run iqtree.

Functions:
    run_iqtree(alignment_path: str, partition_path: str, out_path: str) -> None

'''
from getphylo.utils import io

def run_iqtree(
    alignment_path: str, out_path: str, partition_path: str=None, iqtree_location: str='iqtree'
    ) -> None:
    '''
    Run fasttree on a protein alignment.
        Arguments:
            alignment_path: path to the alignment
            partition_path: path to the partition file
            out_path: path to the output file
        Returns:
            None
    '''
    command = [
            iqtree_location,
            '-s', alignment_path,
            '-pre', out_path, 
            '-m', 'MFP',
            '-bb', '1000'
            ]
    if partition_path is not None:
        command.append('-spp')
        command.append(partition_path)
    io.run_in_command_line(command)
