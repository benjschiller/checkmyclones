'''
Created on Aug 11, 2011

@author: ben
'''
from cogent.align.algorithm import sw_align
from cogent import Alignment
from cogent.core.sequence import SequenceI
#from cogent.parse.record import FileFormatError

def align_clone_to_ref(clone, reference):
    '''Aligns two pycogent Sequence objects with Smith-Waterman algorithm'''
    if not isinstance(clone, SequenceI):
        raise ValueError('clone must be a cogent.SequenceI object')
    if not isinstance(reference, SequenceI):
        raise ValueError('reference must be a cogent.SequenceI object')
    aligned = sw_align(clone, reference)
    ref_len = len(reference)
    clone_matched = str(aligned[0].parseOutGaps()[1])
    ref_matched = str(aligned[1].parseOutGaps()[1])
    if ref_matched == '': raise AlignmentError('No alignment')
    
    aln = CloneAlignment(aligned)
    # find where we are in the reference
    aln.first_ref_pos = str(reference).index(ref_matched)
    aln.first_clone_pos = str(clone).index(clone_matched)
    aln.last_ref_pos = aln.first_ref_pos + len(ref_matched)
    aln.last_clone_pos = aln.first_clone_pos + len(clone_matched)
    match_len = len(ref_matched)
    aln.reference_len = ref_len
    aln.is_truncated = not match_len == ref_len
    aln.has_gaps = aligned[1].isGapped()
    aln.has_mismatches = not aligned[1].canMatch(aligned[0])
    
    aln.Seqs[0].Name = clone.Name
    aln.Seqs[1].Name = reference.Name
    return aln

class CloneAlignment(Alignment):
    '''
    Class for representing two pycogent Sequence objects aligned
        with Smith-Waterman algorithm
    has some methods for judging alignment

    the constructor follows the same format as Alignment and SequenceCollection
    (i.e. takes the data argument directly)
    **The order of sequences matters**
    the first is assumed to be the sequencing result
    the second is assumed to be the reference sequence
    '''

    def __init__(self, data, *args, **kwargs):
        '''
        Constructor
        '''
        super(CloneAlignment, self).__init__(data, *args, **kwargs)
        self.Clone = self.Seqs[0]
        self.Reference = self.Seqs[1]

    def is_match(self):
        """
        Returns true only if clone matched the full reference
        without mismatches or gaps
        """
        return (not self.is_truncated) and (not self.has_gaps) and \
               (not self.has_mismatches) 
            
    def aligned_fragment(self):
        """
        returns the fully aligned clone fragment matching the reference sequence
        raises GapError if the alignment contains a gap
        raises TruncationError if the alignment is truncated
        """
        if self.has_gaps:
            raise GapError('The alignment contains gaps')
        elif self.is_truncated:
            raise TruncationError('The alignment was truncated at pos %s' % 
                                  self.last_good_pos)
        else:
            return self[self.first_good_pos:self.last_good_pos]
    
    def aligned_reference(self):
        return self.Seqs[0]
    
    def aligned_clone(self):
        return self.Seqs[1]
            
class AlignmentError(ValueError):
    def __init__(self, msg):
        super(AlignmentError, self).__init__(msg)

class GapError(AlignmentError): 
    def __init__(self, msg):
        super(GapError, self).__init__(msg)
    
class TruncationError(AlignmentError): 
    def __init__(self, msg):
        super(TruncationError, self).__init__(msg)