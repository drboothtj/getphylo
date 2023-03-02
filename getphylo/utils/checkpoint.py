'''
Checkpoint class used by getphylo.
'''
from enum import IntEnum

class Checkpoint(IntEnum):
    '''Checkpoint class as IntEnum'''
    START = 0
    FASTA_EXTRACTED = 100
    DIAMOND_BUILT = 200
    SINGLETONS_IDENTIFIED = 300
    SINGLETONS_SEARCHED = 400
    SINGLETONS_THRESHOLDED = 500
    SINGLETONS_EXTRACTED = 600
    SINGLETONS_ALIGNED = 700
    ALIGNMENTS_COMBINED = 800
    TREES_BUILT = 900
    DONE = 1000
