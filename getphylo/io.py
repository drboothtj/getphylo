'''
Read and write files

Functions:
    count_files(filename) -> number of files
    change_extension(filename, new_extension) -> new_filename
    get_records_from_genbank(filename) -> records
    read_tsv(database_file) -> database
    write_to_file(write_lines, write_name)
'''
import csv
import glob
import os
from Bio import SeqIO
from getphylo import console


def count_files(filename: str):
    '''count the number of files in a glob'''
    number_of_files = len(glob.glob(filename))
    return number_of_files

def change_extension(filename: str, new_extension: str):
    '''Takes a file name and changes the extension to the string provided'''
    new_filename = filename.split('.')[0] + '.' + new_extension
    return new_filename

def get_records_from_genbank(filename: str):
    '''use BioPython to get genbank records from a given file'''
    records = SeqIO.parse(filename, "genbank")
    return records

def make_folder(name):
    '''attempts to make a folder with the given name and returns an error if it already exists'''
    try:
        os.mkdir(name)
    except OSError as error:
        if error.errno == 17:
            console.print_to_system(
                f'ALERT: The directory {name} already exists. Exiting!'
                )
            exit()
        else:
            raise

def read_file(filename: str):
    '''read a file from the filename'''
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

def write_to_file(filename: str, write_lines: list):
    '''Write a list line by line into a new file'''
    file = open(filename, "a")
    for line in write_lines:
        file.write(line + '\n')
    file.close()
