'''
Build fasta and diamond databases from genbank files

Functions:
build_diamond_databases()
extract_cdses()
get_cds_from_genbank(filename)
extract_data(genbank_path, checkpoint)
'''
import os
import glob
from getphylo import console, diamond, io

def build_diamond_databases(output):
    '''create diamond databases from from all the fasta files in ./fasta/*.fasta'''
    console.print_to_system("Creating diamond databases from extracted cdses...")
    io.make_folder(f'{output}/dmnd')
    for file in glob.glob(f'{output}/fasta/*.fasta'):
        dmnd_database = io.change_extension(file, "dmnd").split('/')[2]
        dmnd_database = (f'{output}/dmnd/' + dmnd_database)
        diamond.make_diamond_database(file, dmnd_database)

def extract_cdses(gbks, output, tag_args):
    '''produce a fasta file from each genbank provided'''
    try:
        os.mkdir(f'{output}/fasta/')
    except OSError as error:
        if error.errno == 17:
            console.print_to_system(f'ALERT: The directory {output}/fasta already exists. Exiting!')
            exit()
        else:
            raise
    else:
        seen = set()
        for filename in glob.glob(gbks):
            console.print_to_system('Extracting CDS annotations from ' + filename)
            get_cds_from_genbank(filename, output, seen, tag_args)

def get_cds_from_genbank(filename, output, seen, tag_args):
    '''extract CDS translations from genbank files into ./fasta/*.fasta'''
    tag = tag_args[0]
    ignore = tag_args[1]
    lines = []
    records = io.get_records_from_genbank(filename)
    for record in records:
        for feature in record.features:
            try:
                if feature.type == "CDS":
                    locus_tag = f'{record.id}_{feature.qualifiers.get(tag)[0]}'
                    if locus_tag in seen:
                        raise ValueError(f'{filename} contains duplicate: {locus_tag}')
                    seen.add(locus_tag)
                    lines.append(">" + locus_tag.replace(".", "_"))
                    lines.append(str(feature.qualifiers.get("translation")[0]))
            except TypeError:
                if feature.qualifiers.get(tag) is None:
                    console.print_to_system(
                        f'ALERT: Missing {tag} in {record.id}.'
                        )
                    if not ignore:
                        exit()
                else:
                    console.print_to_system(
                        "ALERT: " + feature.qualifiers.get("locus_tag")[0] + 'has no translation!'
                        )
    if not lines:
        console.print_to_system("ALERT: No CDS Features in " + filename)
    filename = io.change_extension(filename, "fasta")
    filename = f'{output}/fasta/{filename}'
    io.write_to_file(filename, lines)

def extract_data(checkpoint, output, gbks, tag):
    '''called from main to build fasta and diamond databases from the provided genbankfiles'''
    if checkpoint < 1:
        console.print_to_system("CHECKPOINT 0: Extracting CDSs...")
        extract_cdses(gbks, output, tag)
    if checkpoint < 2:
        console.print_to_system("CHECKPOINT 1: Building diamond databases...")
        build_diamond_databases(output)
