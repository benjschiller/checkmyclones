'''
Created on Aug 11, 2011

@author: ben
'''
from cogent.align.algorithm import nw_align
from cogent import Alignment

class Comparison(object):
    '''
    Compares two pycogent Sequence objects
    the first is assumed to be the reference sequence
    and the second is assumed to be the sequencing data
    '''

    def __init__(self, ref, seq):
        '''
        Constructor
        '''
        aligned = nw_align(ref, seq)
        ref_aligned, seq_aligned = aligned
        alignment = Alignment(aligned)
        self.ref = ref_aligned
        self.seq = seq_aligned
        