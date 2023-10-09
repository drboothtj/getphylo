'''
Input and output operations for getphylo.

Functions:
    get_locus(file:str, locus:str) -> str
    count_files(directory: str) -> int
    change_extension(filename: str, new_extension: str) -> str
    get_records_from_genbank(filename: str) -> List
    make_folder(name: str) -> None
    read_file(filename: str) -> List[str]
    read_tsv(filename: str) -> List[str]
    run_in_command_line(command: str)
    run_in_parallel(function: Callable, args_list: Iterable[List], cpus: int) -> List
    write_to_file(filename: str, write_lines: List[str]) -> None
'''
import csv
import glob
import multiprocessing
import os
import subprocess
import logging
from typing import Callable, Iterable, List
from Bio import SeqIO

from getphylo.utils.errors import GetphyloError, FolderExistsError

def get_locus(fasta: List[str], locus: str) -> str:
    '''
    Returns a sequence from a fasta file with the provided locus name.
        Arguments:
            fasta: list of lines from in fasta format
            locus: locus to search for
        Returns:
            sequence: the sequence of the locus
    '''
    line_number = 0
    sequence = ''
    for line in fasta:
        line_number += 1
        if locus in line:
            sequence = fasta[line_number]
            break
    assert sequence
    return sequence

def count_files(directory: str) -> int:
    '''
    Counts the number of files in a directory.
        Arguments:
            directory: the path to the directory in which to count
        Returns:
            number_of_files: the number of files in that directory as an integer
    '''
    number_of_files = len(glob.glob(directory))
    return number_of_files

def change_extension(filename: str, new_extension: str) -> str:
    '''
    Takes a file name and changes the extension to the string provided.
        Arguments:
            filename: the name of the file with the old extension
            new_extension: the new extension to replace with
        Returns:
            new_filename: the filename with the new extension
    '''
    new_filename = os.path.splitext(filename)[0] + '.' + new_extension
    return new_filename

def get_records_from_genbank(filename: str) -> List:
    '''
    Use BioPython to get genbank records from a given file.
        Arguments:
            filename: the filename of the genbank file
        Returns:
            records: a list genbank records
    '''
    records = SeqIO.parse(filename, "genbank")
    return records

def make_folder(name: str) -> None:
    '''
    Attempts to make a folder with the given name but raises and exception if it already exists.
        Arguments:
            name: the name of the folder being created
        Returns:
            None
    '''
    if os.path.exists(name):
        raise FolderExistsError(
            f'ALERT: The directory {name} already exists.'
            )
    os.mkdir(name)

def read_file(filename: str) -> List[str]:
    '''
    Return a files contents as a list of lines
        Arguments:
            filename: The file to be read
        Returns:
            _file.readlines(): list of strings for each line of the file
    '''
    with open(filename) as _file:
        return _file.readlines()

def read_tsv(filename: str) -> List[str]:
    '''
    Read lines from a .tsv file.
        Arguments:
            filename: the path to the file to be read
        Returns:
            list containing the lines of the .tsv file
    '''
    contents = []
    with open(filename) as file:
        tsv_file = csv.reader(file, delimiter="\t")
        for line in tsv_file:
            contents.append(line)
    return contents

def run_in_command_line(command: str):
    '''
    Convert a string into a command and run in the terminal.
        Aruments:
            command: string containing the command for the terminal
        Returns:
            process: the process being run
    '''
    command = command.split(" ")
    logging.debug(command)
    with subprocess.Popen(command, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT) as process:
        process.communicate()
        if process.returncode != 0:
            raise RuntimeError('Failed to run: ' + str(command))
        return process

def run_in_parallel(function: Callable, args_list: Iterable[List], cpus: int) -> List:
    '''
    Run a given function on avaliable cpus. If only 1 cpu is available, run as normal.
        Arguments:
            function: the function to be called
            args_list: Iterable of lists containing the arguments for each call of the function
            cpus: the number of cpus available
        Returns:
            return_value: a list of return values for each call of the function
    '''
    if cpus <= 1:
        return_value = [function(*args) for args in args_list]
    else:
        for item in args_list:
            try:
                iter(item)
            except TypeError as error:
                raise GetphyloError from error
        try:        
            with multiprocessing.Pool(cpus) as pool:
                results = pool.starmap_async(function, args_list)
                return_value = results.get()
        except Exception:
            raise
    return return_value

def write_to_file(filename: str, write_lines: List[str]) -> None:
    '''
    Write a list line by line into a new file.
        Arguments:
            filename: path to the new file being written
            write_lines: list of strings to be written to the file
        Returns:
            None
    '''
    with open(filename, "a") as _file:
        for line in write_lines:
            _file.write(line + '\n')
        _file.close()
