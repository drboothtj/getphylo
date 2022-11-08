'''
Create alignments from a folder of fasta files.

Functions:

'''
import os
import glob
from getphylo import console, io, muscle, screen #remove dependency on screen

def get_fasta_for_alignments(loci_list, output):
    '''extracts sequences from fasta files for downstream alignment'''
    io.make_folder(f'{output}/unaligned_fasta')
    files = glob.glob(f'{output}/fasta/*.fasta')
    for locus in loci_list:
        write_lines = []
        for filename in files:
            tsv_name = f'{output}/tsvs/{io.change_extension(filename.split("/")[2], "tsv")}'
            tsv = io.read_tsv(tsv_name)
            for line in tsv:
                if locus in line[0]:
                    sequence = screen.get_locus(filename, line[1]) 
                    sequence_name = f'>{filename.split("/")[2].split(".")[0]}_{line[1]}'
                    write_lines.append(sequence_name)
                    write_lines.append(sequence)
        io.write_to_file(f'{output}/unaligned_fasta/{locus}.fasta', write_lines)

def do_alignments(output):
    '''run alignment file through muscle module'''
    io.make_folder(f'{output}/aligned_fasta')
    for file in glob.glob(f'{output}/unaligned_fasta/*.fasta'):
        outfile = f'{output}/aligned_fasta/{file.split("/")[2]}'
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

def make_combined_alignment(gbks, output):
    '''Create a combined alignment from single locus alignments'''
    if os.path.exists('aligned_fasta/combined_alignment.fasta'):
        console.print_to_system(
            'ALERT: aligned_fasta/combined_alignment.fasta already exists. Exiting!'
            )
        exit()
    else:
        combined_alignment = []
        taxa = glob.glob(gbks)
        loci = glob.glob(f'{output}/aligned_fasta/*.fasta')
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
        io.write_to_file(f'{output}/aligned_fasta/combined_alignment.fasta', combined_alignment)
        #add check to recommend removing data that is mostly Ns

def make_alignments(checkpoint, output, loci, gbks):
    '''Main routine for align'''
    if checkpoint < 6:
        console.print_to_system("CHECKPOINT 5: Extracting sequences for alignment...")
        get_fasta_for_alignments(loci, output)
    if checkpoint < 7:
        console.print_to_system("CHECKPOINT 6: Aligning sequences...")
        do_alignments(output)
    if checkpoint < 8:
        console.print_to_system("CHECKPOINT 7: Combining alignments")
        make_combined_alignment(gbks, output)
