'''
Main routine for getphylo

Functions:
main()
set_seed()

'''
import glob
from getphylo import align, console, extract, io, parser, screen, trees

def set_seed(seed, checkpoint, gbks):
    '''get the seed genome name'''
    if seed is None:
        if checkpoint > 0:
            console.print_to_system(
                'ALERT: A checkpoint has been set! Please ensure the seed is defined. Exiting!'
                )
            exit()
        else:
            gbks = glob.glob(gbks)
            seed = gbks[0]
            console.print_to_system(
                f'No seed defined. Using first file in glob ({seed}) as seed.'
                )
    return seed

def main():
    '''main routine for getphylo'''
    console.print_to_system("Running getphylo version 0.0.1")
    #store arguments as args and get checkpoint and output directory, then get args as required
    args = parser.parse_args()
    checkpoint = args.checkpoint
    output = args.output
    try:
        io.make_folder(output)
    except:
        console.print_to_system(
            f'ALERT: {output} already exists. Continuing analysis in that directory.'
            )
    #set seed
    gbks = args.gbks
    seed = args.seed
    seed = set_seed(seed, checkpoint, gbks)
    console.print_to_system("The seed genome is " + seed)
    #checkpoint 0 and 1 extracting the CDSs
    if checkpoint < 2:
        extract.extract_data(checkpoint, output, gbks)
    #checkpoint 2, 3 and 4 screening the target proteins
    thresholds = [args.find, args.minlength, args.maxlength, args.presence, args.minloci]
    if checkpoint < 5:
        final_loci = screen.get_target_proteins(checkpoint, output, seed, thresholds)
    if checkpoint > 4:
        final_loci = screen.get_loci_from_file(f'{output}/final_loci.txt')
    #checkpoint 5, 6 and 7 extracting sequences and alignment
    if checkpoint < 8:
        align.make_alignments(checkpoint, output, final_loci, gbks)
    #checkpoint 8 building trees
    if checkpoint < 9:
        trees.make_trees(output)
    console.print_to_system("Analysis complete. Thank you for using getphylo!")

    #add logging

if __name__ == "__main__":
    main()
