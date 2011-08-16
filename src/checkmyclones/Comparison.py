'''
Created on Aug 11, 2011

@author: ben
'''
from cogent.align.algorithm import nw_align
from cogent import Alignment
import operator

class Comparison(object):
    '''
    Compares two pycogent Sequence objects
    the first is assumed to be the sequencing result
    the second is assumed to be the reference sequence
    '''

    def __init__(self, clone_seq, ref_seq):
        '''
        Constructor
        '''
        aligned = nw_align(clone_seq, ref_seq)
        ref_len = len(ref_seq)
        clone_len = len(clone_seq)
        ref_aligned, clone_aligned = aligned
        aln = Alignment(aligned)
        self.aln = aln
#        aln_len = len(aln)
        aln_eq = map(operator.eq, *aligned)
        seq_to_aln_map = aln.Seqs[0].getGappedSeq().gapMaps()[0]
        first_good_pos = seq_to_aln_map[0]
        self.ref_aligned = ref_aligned
        self.clone_aligned = clone_aligned
        self.is_truncated = False
        self.has_mismatches = False
        self.aligned_frag = None
        if sum(aln_eq) == ref_len:
            # we matched the whole reference, but there might be insertions
            last_good_pos = seq_to_aln_map[ref_len - 1]
            self.aln_frag = aln[first_good_pos:last_good_pos]
            if last_good_pos - first_good_pos + 1 == ref_len:
                # perfect match
                self.has_gaps = False
                self.longest_aln = self.aln_frag
                return
            else:
                self.has_gaps = True
                spans = self.aln_frag.Seqs[0].map.nongap().spans
                longest_span = sorted(enumerate(map(len, spans)),
                                      key=operator.itemgetter(1))[-1]
                start, end = longest_span.Start, longest_span.End
                self.longest_aln = self.aln_frag[start:end]
                # Let the user see the gaps for now
                return
        else:
            # we did not match the whole reference sequence
            aln_frag = aln[first_good_pos:]
            aln_frag_eq = map(operator.eq, *aln_frag)
            uneq_count = aln_frag_eq.count(False)
            ref_frag_aln = aln_frag.Seqs[0]
            spans = ref_frag_aln.map.nongap().spans
            longest_span = sorted(enumerate(map(len, spans)),
                                  key=operator.itemgetter(1))[-1]
            start, end = longest_span.Start, longest_span.End
            self.longest_aln = ref_frag_aln[start:end]
            gap_count = len(ref_frag_aln.gaps())
            self.has_gaps = gap_count > 0
            # maybe we are truncated?
            if ref_len > clone_len - first_good_pos:
                self.is_truncated = True
                if uneq_count == 0: return
                else:
                    # there are mismatches or gaps
                    self.has_mismatches = not gap_count == uneq_count
                    return
            else:
                # we are not truncated but there are mismatches
                self.has_mismatches = True
                return
                
        