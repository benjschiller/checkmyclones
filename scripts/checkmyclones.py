#!/usr/bin/env python
from glob import glob
import os.path
import cogent
from cogent import LoadSeqs, DNA
import argparse
from twobitreader import TwoBitFile
import scripter
from scripter import Usage, get_logger, debug
from pkg_resources import get_distribution
import platform
import signal
from logging import WARNING
from pprint import pprint
from itertools import groupby
import multiprocessing
from clonechecker.align import align_clone_to_ref, AlignmentError
from clonechecker.filetools import load_seqs, read_bed_file, find_2bit_file
from operator import itemgetter, attrgetter
try:
    from cPickle import dumps, loads
except ImportError:
    from pickle import dumps, loads
__version__ = get_distribution('checkmyclones').version
VERSION = __version__

def load_all_seqs(glob_path, moltype=DNA, recursive=True):
    """
    loads all available sequences that match a path
    wildcards are allowed
    all sequences are assumed to be DNA
    
    if the path is a directory, then we will search recursively
    """
    seqs = []
    if type(glob_path) is list:
        for p in glob_path:
            seqs.extend(load_all_seqs(p, moltype=DNA, recursive=True))
        return seqs
    for foo in glob(glob_path):
        print 'Trying to import sequences from %s' % foo
        if os.path.isfile(foo):
            try:
                myseqs = load_seqs(foo, moltype=moltype)
            except:
                print 'ERROR: Loading %s failed' % foo
                continue
            for seq in myseqs:
                seq.Name = '%s (%s)' % (seq.Name, foo)
            seqs.extend(myseqs)
        elif recursive and os.path.isdir(foo):
            for boo in glob(os.path.join(foo,'*')):
                seqs.extend(load_all_seqs(boo, moltype=moltype, recursive=True))
    return seqs

def real_name(name):
    """
    strip coordinates from a name so that we can compare to original 
    """
    end = name.find('(')
    if not end == -1: name = name[0:end-1]
    return name

def main():
    """
    runs the main checkmyclones script
    """
    e = scripter.Environment(doc=__doc__, version=VERSION, handle_files=False)
    parser = e.argument_parser
    parser.add_argument('--path-to-gbdb', default='/gbdb',
help='Location of "gdbdb" or 2bit files. If gbdb is not in /gbdb or C:\gbdb, specify the path here')
    ggroup = parser.add_mutually_exclusive_group()
    ggroup.add_argument('--genome',
               help='Use 2bit file foo as reference genome (Looks also for {path-to-gbdb}/foo/foo.2bit))')
    ggroup.add_argument('--hg18', const='hg18', action='store_const',
               help='Shortcut for --genome hg18') 
    ggroup.add_argument('--hg19', const='hg19', action='store_const',
               help='Shortcut for --genome hg19') 
    ggroup.add_argument('--mm9', const='mm9', action='store_const',
               help='Shortcut for --genome mm9')
    parser.add_argument('--reverse-orientation', action='store_true',
                        help='Check only the reverse orientation')
    parser.add_argument('--both-orientations', action='store_true',
                        help='Check forward and reverse orientations')
    parser.add_argument('--clones', nargs='+',
                        help='list of files that contain clone sequences')
    parser.add_argument('--references', nargs='*',
                        help='list of files that contain reference sequences')
    parser.add_argument('--bed-reference',
                        help='Use the regions listed in the bed file as reference sequences')
    parser.add_argument('--only-use-references', nargs='*',
                        help='Use the only regions with the following names')
    parser.set_defaults(**{'genome': 'hg19', 'logging_level': WARNING})
    args = parser.parse_args()
    context = vars(args)
    scripter.LOGGER.setLevel(context['logging_level'])
    clones = load_all_seqs(context['clones'], recursive=context['recursive'])
    ref_seqs = []
    if len(clones) == 0:
        raise Usage('Could not find any clone sequences')
    if context['references'] is None and context['bed_reference'] is None:
        raise Usage('No reference sequences specified')
    else:
        if context['bed_reference'] is not None:
            genome = find_2bit_file(context['genome'], context['path_to_gbdb'])
            print 'Fetching sequences from %s using %s' % (context['bed_reference'],
                                                           genome)
            ref_seqs.extend(read_bed_file(context['bed_reference'],
                                          genome=genome))
        if context['references'] is not None:
            ref_seqs.extend(load_all_seqs(context['references'],
                                          recursive=context['recursive']))
        specified_references = context['only_use_references']
        if specified_references is not None:
            good_name = lambda ref: real_name(ref.Name) in specified_references
            ref_seqs = filter(good_name, ref_seqs)
        if len(ref_seqs) == 0:
            raise Usage('Could not find any reference sequences')
    signal.signal(signal.SIGCHLD, signal.SIG_DFL)
    debug('multiprocessing enabled')
    p = multiprocessing.Pool(processes=context['num_cpus'])
    debug('Initialized pool of %d workers', context['num_cpus'])
    results = []
    forward = not context['reverse_orientation']
    rc = context['reverse_orientation'] or context['both_orientations'] or False
    for ref in ref_seqs:
        print 'Loaded reference %s' % ref.Name
    for clone in clones:
        p.apply_async(announce_first, (clone,), context)
        for ref in ref_seqs:
            if forward:
                r = p.apply_async(compare_clone_to_ref,
                                  (clone, ref), context)
                results.append(r)
            if rc:
                r = p.apply_async(compare_clone_to_ref,
                                  (clone.rc(), ref), context)
                results.append(r)
    p.close()
    p.join()
    result_values = []
    for r in results:
        current_pickle = r.get()
        current_result = loads(current_pickle)
        if current_result is None: continue
        else: result_values.append(current_result)
    all_matches = []
    for clone_name, group in groupby(result_values, key=itemgetter(0)):
        alns = map(itemgetter(1), list(group))
        is_matched = lambda aln: not aln.is_truncated and not aln.has_gaps
        matches = filter(is_matched, alns)
        if len(matches) > 0:
            all_matches.extend(alns)
            continue
        elif len(alns) == 0:
            print 'No match for %s' % clone_name
            continue
        else:
            print_good_alns(alns)
    print_matched_alns(all_matches)
    return

def print_good_alns(alns):
    """
    takes a list of CloneAlignments and prints any alignment of length
    at least the reference length
    
    if there are no such alignments, print_good_alns will print the longest
    available alignment
    """
    alns_sorted = sorted(alns, key=len, reverse=True)
    best_aln = alns_sorted[0]
    max_length = len(best_aln)
    ref_len = best_aln.reference_len
    clone_name = best_aln.Clone.Name
    print 'Partially matched %s to %s with length %s / %s' % (
            clone_name, best_aln.Reference.Name, max_length, ref_len)
    fasta_print(best_aln.Reference, name=best_aln.Reference.Name)
    fasta_print(best_aln.Clone, name=clone_name)
    if max_length > ref_len: return
    for aln in alns_sorted[1:]:
        if max_length > aln.reference_len:
            print 'Also partially matched %s to %s with length %s / %s' % (
                     clone_name, aln.Reference.Name, len(aln), ref_len)
            fasta_print(aln.Reference, name=aln.Reference.Name)
            fasta_print(aln.Clone, name=clone_name)
            print aln
        else: break
    return

def print_matched_alns(matched_alns):
    """
    print matched alignments to stdout
    """
    clones = []
    references = []
    print '='*60
    for ref, group in groupby(matched_alns, key=attrgetter('Reference')):
        ref_name = ref.Name
        print 'Matches to %s:' % ref_name
        alns = list(group)
        for aln in alns:
            if aln.has_mismatches:
                msg = '\t%s with mismatches (%s / %s)'
                print msg % (aln.Clone.Name, len(aln),
                             len(aln.degap().Seqs[1]))
                mismatches = [(i, x) for i, x in
                              enumerate(aln.iterPositions()) if not x[0]==x[1]]
                try:
                    coords = ref_name[ref_name.find('(') + 1:-1]
                    chr = coords[0:coords.find(':')]
                    start = coords[coords.find(':')+1:coords.find('-')]
                    for m in mismatches:
                        msg = '\t\t%s %d (1-based position %d) %s->%s'
                        print msg % (chr, int(start) + m[0] + 1, m[0] + 1,
                                     m[1][1], m[1][0])
                except ValueError, TypeError:
                    for m in mismatches:
                        msg = '\t\t1-based position %d %s->%s'
                        print msg % (m[0] + 1, m[1][0],
                                     m[1][1])
            else:
                print '%s perfectly' % aln.Clone.Name
        fasta_print(ref, name=ref_name)
        for aln in alns:
            fasta_print(aln.Clone, name=aln.Clone.Name)
        print '='*60

def fasta_print(a, name='', width=60):
    """print an Aligned sequence object in fasta format"""
    print '>%s' % name
    j = 0
    line = []
    for x in a:
        if j == width:
            print ''.join(line)
            line = []
            j = 0
        line.append(x)
        j += 1
    if len(line) > 0: print ''.join(line)

def announce_first(clone, logging_level=20, **kwargs):
    """
    Announce to stderr that we are about to start comparing a clone to
    all available references
    """
    logger = get_logger(logging_level)
    logger.info('Comparing %s to references', clone.Name)

def compare_clone_to_ref(clone, ref, logging_level=20, **kwargs):
    """
    process-safe comparison of a clone Sequence to reference Sequence
    aligns using clonechecker.align.align_clone_to_ref
    returns the pickled tuple (clone, Alignment)
            or None (if AlignmentError is raised)
    """
    logger = get_logger(logging_level)
    try:
        aln = align_clone_to_ref(clone, ref)
    except AlignmentError: return
    if aln.is_match():
        logger.info('match found %s, %s', clone.Name, ref.Name)
    if aln.is_truncated:
        logger.debug('alignment is truncated (%s, %s)', clone.Name, ref.Name)
        if not aln.has_gaps and not aln.has_mismatches:
            logger.info('truncated match found %s, %s', clone.Name, ref.Name)
    if aln.has_gaps:
        logger.debug('alignment has gaps (%s, %s)', clone.Name, ref.Name)
    if aln.has_mismatches:
        logger.debug('alignment has mismatches (%s, %s)', clone.Name, ref.Name)
        if not aln.is_truncated and not aln.has_gaps:
            logger.info('mutated match found %s, %s', clone.Name, ref.Name)
    return dumps((aln.Clone.Name, aln))

if __name__=='__main__': main()