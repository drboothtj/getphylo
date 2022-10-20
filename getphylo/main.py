'''
Main routine for getphylo

Functions:
main()
set_seed()

'''
import glob
from random import shuffle
from getphylo import align, console, extract, parser, screen, trees

def set_seed():
    '''get the seed genome name'''
    seed = parser.get_seed()
    if seed is None:
        gbks = glob.glob(parser.get_gbks())
        shuffle(gbks)
        seed = gbks[0]
    return seed

def main():
    '''main routine for getphylo'''
    console.print_to_system("Running getphylo version 0.0.1")
    seed = set_seed()
    console.print_to_system("The seed genome is " + seed)
    checkpoint = parser.get_checkpoint()
    if checkpoint < 2:
        extract.extract_data(checkpoint)
    if checkpoint < 5:
        final_loci = screen.get_target_proteins(checkpoint, seed)
    if checkpoint > 4:
        final_loci = screen.get_loci_from_file('final_loci.txt')
    if checkpoint < 8:
        align.make_alignments(checkpoint, final_loci)
    if checkpoint <9: 
        trees.make_trees()
    console.print_to_system("Analysis complete. Thank you for using getphylo!")

    #if checkpoint < Y align.main(checkpoint)
    #if checkpoint < Y tree.main(checkpoint)
    #add custom error types for files existing etc.
    #add logging
    #allow screening multiple genomes for target loci (loop screening)
    #structure as single diamond database
    #generate bio utils library
    #add reading and writing details at each cp
    #add file/directory check as io function!
if __name__ == "__main__":
    main()
