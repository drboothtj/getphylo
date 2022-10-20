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
from getphylo import console, diamond, io, parser

def build_diamond_databases():
    '''create diamond databases from from all the fasta files in ./fasta/*.fasta'''
    console.print_to_system("Creating diamond databases from extracted cdses...")
    try:
        os.mkdir('dmnd')
    except OSError as error:
        if error.errno == 17:
            console.print_to_system(
                "ALERT: The directory './dmnd/' already exists. Exiting!"
                )
            exit()
        else:
            raise
    else:
        for file in glob.glob("fasta/*.fasta"):
            dmnd_database = io.change_extension(file, "dmnd").split('/')[1]
            dmnd_database = ("dmnd/" + dmnd_database)
            diamond.make_diamond_database(file, dmnd_database)

def extract_cdses():
    '''produce a fasta file from each genbank provided'''
    try:
        os.mkdir('fasta')
    except OSError as error:
        if error.errno == 17:
            console.print_to_system("ALERT: The directory './fasta/' already exists. Exiting!")
            exit()
        else:
            raise
    else:
        seen = set()
        for filename in glob.glob(parser.get_gbks()):
            console.print_to_system('Extracting CDS annotations from ' + filename)
            get_cds_from_genbank(filename, seen)

def get_cds_from_genbank(filename, seen):
    '''extract CDS translations from genbank files into ./fasta/*.fasta'''
    lines = []
    records = io.get_records_from_genbank(filename)
    for record in records:
        for feature in record.features:
            try:
                if feature.type == "CDS":
                    locus_tag = f'{record.id}_{feature.qualifiers.get("locus_tag")[0]}'
                    if locus_tag in seen:
                        raise ValueError(f'{filename} contains duplicate: {locus_tag}')
                    seen.add(locus_tag)
                    lines.append(">" + locus_tag.replace(".","_"))
                    lines.append(str(feature.qualifiers.get("translation")[0]))
            except TypeError:
                console.print_to_system(
                    "ALERT: " + feature.qualifiers.get("locus_tag")[0] + 'has no translation!'
                    )
    if not lines:
        console.print_to_system("ALERT: No CDS Features in" + filename)
    filename = io.change_extension(filename, "fasta")
    filename = "fasta/" + filename
    io.write_to_file(filename, lines)

def extract_data(checkpoint):
    '''called from main to build fasta and diamond databases from the provided genbankfiles'''
    if checkpoint < 1:
        console.print_to_system("CHECKPOINT 0: Extracting CDSs...")
        extract_cdses()
    if checkpoint < 2:
        console.print_to_system("CHECKPOINT 1: Building diamond databases...")
        build_diamond_databases()
