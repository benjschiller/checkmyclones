'''
Created on Aug 11, 2011

@author: ben
'''
from cogent.align.algorithm import sw_align
from cogent import Alignment
from cogent.core.sequence import SequenceI
#from cogent.parse.record import FileFormatError

class CloneToRefAlignment(Alignment):
    '''
    Aligns two pycogent Sequence objects with Smith-Waterman algorithm
    has some methods for judging alignment

    the first is assumed to be the sequencing result
    the second is assumed to be the reference sequence
    '''

    def __init__(self, clone_seq, ref_seq):
        '''
        Constructor
        '''
        if not isinstance(clone_seq, SequenceI):
            raise ValueError('clone_seq must be a cogent.SequenceI object')
        if not isinstance(ref_seq, SequenceI):
            raise ValueError('ref_seq must be a cogent.SequenceI object')
        aligned = sw_align(clone_seq, ref_seq)
        ref_len = len(ref_seq)
        # find where we are in the reference
        clone_matched = str(aligned[0].parseOutGaps()[1])
        ref_matched = str(aligned[1].parseOutGaps()[1])
        if ref_matched == '': raise AlignmentError('No alignment')
        match_len = len(ref_matched)
        self.first_ref_pos = str(ref_seq).index(ref_matched)
        self.first_clone_pos = str(clone_seq).index(clone_matched)
        self.last_ref_pos = self.first_ref_pos + len(ref_matched)
        self.last_clone_pos = self.first_clone_pos + len(clone_matched)
        
        self.is_truncated = not match_len == ref_len
        self.has_gaps = aligned[1].isGapped()
        self.has_mismatches = not aligned[1].canMatch(aligned[0])
        
        super(CloneToRefAlignment, self).__init__(aligned)
        self.Seqs[0].Name = clone_seq.Name
        self.Seqs[1].Name = ref_seq.Name


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