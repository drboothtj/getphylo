'''
Unique errors for getphylo.
'''

class GetphyloError(Exception):
    '''General class of errors unique to getphylo'''
    pass

class BadInputError(GetphyloError):
    '''Called when user provides bad input'''
    pass

class BadSeedError(GetphyloError):
    '''Called when a seed cannot be correctly set'''
    pass

class NoFinalLociError(GetphyloError):
    '''Called when final_loci is empty and cannot be read from final_loci.txt'''
    pass

class NoCandidateLociError(GetphyloError):
    '''Called when candidate_loci is empty and cannot be read from final_loci.txt'''
    def __init__(self, path):
        self.path = path
        super().__init__(
            f'Candidate loci could not be read from {self.path}. '
            'If restarting from a checkpoint ensure there is a final_loci.txt '
            'file in the specified output folder.'
            )

class BadAnnotationError(GetphyloError):
    '''
        Called when a genbank files is poorly annotated
        (e.g. duplicate locus tags or missing annotations)
    '''
    pass

class BadRecordError(GetphyloError):
    '''Called when BioPython cannot read records due to misformatting'''
    pass

class FolderExistsError(GetphyloError):
    '''Called by getphylo.utils.io.make_folder when a folder exists.'''
    pass

class FileAlreadyExistsError(GetphyloError):
    '''Called by getphylo.screen when a attempting to write a file and that file already exists.'''
    pass

class InsufficientLociError(GetphyloError):
    '''Called in screen if the number of loci are below the threshold defined by the user'''
    pass

class BadLocusError(GetphyloError):
    '''Called in align when a locus is not present.'''
    pass

class BadExecutableError(GetphyloError):
    '''Called when a non-existant executable path is provided'''
    def __init__(self, command):
        self.command = command
        super().__init__(
        f'getphylo could not find the executable "{self.command[0]}", ' +
        'please ensure the correct paths to all executables are provided'
        )

class BadMethodError(GetphyloError):
    '''
    Called if a phylogentic tool is defined that is not 'fasttree' or 'iqtree'
    Note: It shouldn't be feasable for the user.
    '''
    def __init__(self, method):
        self.method = method
        super().__init__(
        f'You somehow selected {self.method} as the method.'
        'It should not be possible for you to generate this error - please report!'
        )
