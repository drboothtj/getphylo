'''
Create alignments from a folder of fasta files.

Functions:

'''
import os
import glob
from getphylo import console, io, muscle, parser, screen

def get_fasta_for_alignments(loci_list):
    '''extracts sequences from fasta files for downstream alignment'''
    try:
        os.mkdir('unaligned_fasta')
    except OSError as error:
        if error.errno == 17:
            console.print_to_system(
                "ALERT: The directory './unaligned_fasta/' already exists. Exiting!"
                )
            exit()
        else:
            raise
    else:
        files = glob.glob("fasta/*.fasta")
        for locus in loci_list:
            write_lines = []
            for filename in files:
                tsv_name = "tsvs/" + io.change_extension(filename.split("/")[1], "tsv")
                tsv = io.read_tsv(tsv_name)
                for line in tsv:
                    if locus in line[0]:
                        sequence = screen.get_locus(filename, line[1]) #move to io?
                        sequence_name = f">{filename.split('/')[1].split('.')[0]}_{line[1]}"
                        write_lines.append(sequence_name)
                        write_lines.append(sequence)
                        break
            io.write_to_file("unaligned_fasta/" + locus + ".fasta", write_lines)


def do_alignments():
    '''run alignment file through muscle module'''
    try:
        os.mkdir('aligned_fasta')
    except OSError as error:
        if error.errno == 17:
            console.print_to_system(
                "ALERT: The directory './aligned_fasta/' already exists. Exiting!"
                )
            exit()
        else:
            raise
    else:
        for file in glob.glob("./unaligned_fasta/*.fasta"):
            outfile = "aligned_fasta/" + file.split("/")[2]
            muscle.run_muscle(file, outfile)

def get_locus_length(alignment):
    '''return the length of the fasta locus'''
    length = 0
    for line in alignment[1:]:
        if '>' not in line:
            length += len(line)-1
        else:
            break
    return length

def get_locus_alignment(alignment, taxon_name):
    '''extract individual loci alignments for the taxa'''
    seq_flag = False
    extracted_alignment = []
    for line in alignment:
        if line[0] == '>':
            seq_flag = False
        if '>' in line and taxon_name in line:
            seq_flag = True
        if seq_flag is True and taxon_name not in line:
            extracted_alignment.extend(line[:-1])
    return ''.join(extracted_alignment)

def make_combined_alignment():
    '''Create a combined alignment from single locus alignments'''
    if os.path.exists('aligned_fasta/combined_alignment.fasta'):
        console.print_to_system(
            'ALERT: aligned_fasta/combined_alignment.fasta already exists. Exiting!'
            )
        exit()
    else:
        combined_alignment = []
        taxa = glob.glob(parser.get_gbks())
        loci = glob.glob('aligned_fasta/*')
        for taxon in taxa:
            taxon_name = taxon.split('.')[0]
            combined_alignment.append(f'>{taxon_name}')
            sequence_data = []
            for locus in loci:
                alignment = io.read_file(locus)
                locus_length = get_locus_length(alignment)
                locus_alignment = get_locus_alignment(alignment, taxon_name)
                if len(locus_alignment) == locus_length:
                    sequence_data.append(locus_alignment)
                else:
                    sequence_data.append('N' * locus_length)
            sequence_string = ''.join(sequence_data)
            combined_alignment.append(sequence_string)
        io.write_to_file('aligned_fasta/combined_alignment.fasta', combined_alignment)
        #add check to recommend removing data that is mostly Ns

def make_alignments(checkpoint, loci):
    '''Main routine for align'''
    if checkpoint < 6:
        console.print_to_system("CHECKPOINT 5: Extracting sequences for alignment...")
        get_fasta_for_alignments(loci)
    if checkpoint < 7:
        console.print_to_system("CHECKPOINT 6: Aligning sequences...")
        do_alignments()
    if checkpoint < 8:
        console.print_to_system("CHECKPOINT 7: Combining alignments")
        make_combined_alignment()
