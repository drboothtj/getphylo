'''
Main routine for getphylo

Functions:
main()
set_seed()

'''
import glob
from random import shuffle
from getphylo import console, extract, parser, screen

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
        extract.extract_data()
    if checkpoint < 3:
        screen.get_target_proteins(seed)

    #if checkpoint < Y align.main(checkpoint)
    #if checkpoint < Y tree.main(checkpoint)
    #add custom error types for files existing
    #add logging
if __name__ == "__main__":
    main()
