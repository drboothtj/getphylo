'''
Unique errors for getphylo.
'''

class GetphyloError(Exception):
    '''General class of errors unique to getphylo'''
    pass

class BadSeedError(GetphyloError):
    '''Called when a seed cannot be correctly set'''
    pass

class NoFinalLociError(GetphyloError):
    '''Called when final_loci is empty and cannot be read from final_loci.txt'''
    pass

class NoCandidateLociError(GetphyloError):
    '''Called when candidate_loci is empty and cannot be read from final_loci.txt'''
    pass

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
