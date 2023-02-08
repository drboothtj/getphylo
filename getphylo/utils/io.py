'''
Read and write files

Functions:
    count_files(directory: str) -> int
    change_extension(filename:str, new_extension:str) -> str
    get_records_from_genbank(filename:str) -> List[]
    make_folder(name: str) -> None
    read_tsv(database_file) -> database
    write_to_file(write_lines, write_name)
'''
import csv
import glob
import os
import subprocess
from typing import List
from Bio import SeqIO
from getphylo.utils import console

def count_files(directory: str) -> int:
    '''Counts the number of files in a directory.
        Arguments:
            directory: the path to the directory in which to count
        Returns:
            number_of_files: the number of files in that directory as an integer
    '''
    number_of_files = len(glob.glob(directory))
    return number_of_files

def change_extension(filename: str, new_extension: str) -> str:
    '''Takes a file name and changes the extension to the string provided.
        Arguments:
            filename: the name of the file with the old extension
            new_extension: the new extension to replace with
        Returns:
            new_filename: the filename with the new extension
    '''
    new_filename = os.path.splitext(filename)[0] + '.' + new_extension
    return new_filename

def get_records_from_genbank(filename: str) -> List:
    '''Use BioPython to get genbank records from a given file.
        Arguments: 
            filename: the filename of the genbank file
        Returns:
            records: a list genbank records'''
    records = SeqIO.parse(filename, "genbank")
    return records

def make_folder(name: str) -> None:
    '''Attempts to make a folder with the given name but raises and exception if it already exists.
        Arguments: 
            name: the name of the folder being created
        Returns:
            None'''
    if os.path.exists(name):
        raise Exception (
            f'ALERT: The directory {name} already exists.'
            )
    else:
        os.mkdir(name)

def read_file(filename: str) -> List[str]:
    '''Return a files contents as a list of lines
        Arguments: 
            filename: The file to be read
        Returns:
            _file.readlines(): list of strings for each line of the file
    '''
    _file = open(filename, "r")
    return _file.readlines()

def read_tsv(filename: str):
    '''read lines from a .tsv file'''
    contents = []
    with open(filename) as file:
        tsv_file = csv.reader(file, delimiter="\t")
        for line in tsv_file:
            contents.append(line)
    return contents

def run_in_command_line(command: str):
    '''Convert a string into a command and run in the terminal'''
    command = command.split(" ")
    process = subprocess.Popen(command, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
    process.communicate()
    return process

def write_to_file(filename: str, write_lines: list):
    '''Write a list line by line into a new file'''
    file = open(filename, "a")
    for line in write_lines:
        file.write(line + '\n')
    file.close()
