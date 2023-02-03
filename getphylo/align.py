'''
Create alignments from a folder of fasta files.

Functions:

'''
import os
import glob
from getphylo import console, io, muscle, screen #remove dependency on screen

def make_fasta_for_alignments(loci_list:[str], output:str) -> None:
    '''builds a .fasta file from sequences where there is a hit in the diamond search results
    
    Arguments:
        loci_list: a list of locus tags to extract hits
        output: the path of the output directory

    Returns:
        None
    '''
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

def do_alignments(output:str) -> None:
    '''runs alignment file through muscle module
    
    Arguments:
        output: the path of the outut directory

    Returns:
        None
    '''
    io.make_folder(f'{output}/aligned_fasta')
    for file in glob.glob(f'{output}/unaligned_fasta/*.fasta'):
        outfile = f'{output}/aligned_fasta/{file.split("/")[2]}'
        muscle.run_muscle(file, outfile)

def get_locus_length(alignment:[str]) -> int:
    '''returns the length of the fasta locus by counting sequence lines
    
    Arguments:
        alignment: a list containing each line of a fasta file

    Returns:
        alignment_length: an int of the alignment length
    '''
    for line in alignment[1:]:
        if '>' not in line:
            length = len(line)-1
        else:
            break
    return alignment_length

def get_locus_alignment(alignment, taxon_name):
    '''extract individual loci alignments for the taxa'''
    seq_flag = False
    extracted_alignment = []
    for line in alignment:
        if line[0] == '>':
            seq_flag = False
        if '>' in line and taxon_name in line:
            seq_flag = True
        if seq_flag and taxon_name not in line:
            extracted_alignment.extend(line[:-1])
    return ''.join(extracted_alignment)

def make_combined_alignment(gbks:[str], output:str) -> None:
    '''Create a combined alignment from single locus alignments
    
        Arguments:
        
        Returns: None
        '''
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
            sequence_data = []
            for locus in loci:
                alignment = io.read_file(locus)
                locus_length = get_locus_length(alignment)
                locus_alignment = get_locus_alignment(alignment, taxon_name)
                if len(locus_alignment) == locus_length:
                    sequence_data.append(locus_alignment)
                else:
                    sequence_data.append('X' * locus_length)
            sequence_string = ''.join(sequence_data)
            if sequence_string.count('X') == len(sequence_string):
                console.print_to_system(f'[ALERT]: {taxon_name} has no sequence data and has been removed.') #maybe convert to a fatal error.
            else:
                combined_alignment.append(f'>{taxon_name}')
                combined_alignment.append(sequence_string)
        io.write_to_file(f'{output}/aligned_fasta/combined_alignment.fasta', combined_alignment)
        #provide partition data!

def make_alignments(checkpoint, output, loci, gbks):
    '''Main routine for align'''
    if checkpoint < 6:
        console.print_to_system("CHECKPOINT 5: Extracting sequences for alignment...")
        make_fasta_for_alignments(loci, output)
    if checkpoint < 7:
        console.print_to_system("CHECKPOINT 6: Aligning sequences...")
        do_alignments(output)
    if checkpoint < 8:
        console.print_to_system("CHECKPOINT 7: Combining alignments")
        make_combined_alignment(gbks, output)
