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
from itertools import groupby
import clonechecker
import clonechecker.AlignClone
import multiprocessing
from clonechecker.AlignClone import CloneToRefAlignment
from clonechecker.filetools import load_one_seq, read_bed_file, find_2bit_file
from operator import itemgetter
__version__ = get_distribution('checkmyclones').version
VERSION = __version__

def load_all_seqs(glob_path, moltype=DNA, recursive=True):
    seqs = []
    if type(glob_path) is list:
        for p in glob_path:
            seqs.extend(load_all_seqs(p, moltype=DNA, recursive=True))
        return seqs
    for foo in glob(glob_path):
        print 'Trying to import sequence %s' % glob_path
        if os.path.isfile(foo):
            seq =load_one_seq(foo, moltype=moltype)
            seq.Name = '%s (%s)' % (seq.Name, foo)
            seqs.append(seq)
        elif recursive and os.path.isdir(foo):
            for boo in glob(os.path.join(foo,'*')):
                seqs.extend(load_all_seqs(boo, moltype=moltype, recursive=True))
    return seqs

def real_name(name):
    end = name.find('(')
    if not end == -1: name = name[0:end-1]
    return name

def main():
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
    parser.add_argument('--clones', nargs='+',
                        help='list of files that contain clone sequences')
    parser.add_argument('--references', nargs='*',
                        help='list of files that contain reference sequences')
    parser.add_argument('--bed-reference',
                        help='Use the regions listed in the bed file as reference sequences')
    parser.add_argument('--only-use-references', nargs='*',
                        help='Use the only regions with the following names')
    parser.set_defaults(**{'genome': 'hg19'})
    args = parser.parse_args()
    context = vars(args)
    scripter.LOGGER.setLevel(context['logging_level'])
    clones = load_all_seqs(context['clones'], recursive=context['recursive'])
    ref_seqs = []
    if len(clones) == 0:
        raise Usage('Could not find any clone sequences')
    if context['references'] is None and context['bed_reference'] is None:
        raise Usage('No reference sequeneces specified')
    else:
        if context['bed_reference'] is not None:
            genome = find_2bit_file(context['genome'], context['path_to_gbdb'])
            debug('Using genome %s', genome)
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
    for ref in ref_seqs:
        print 'Loaded reference %s' % ref.Name
    for clone in clones:
        p.apply_async(announce_first, (clone,), context)
        for ref in ref_seqs:
            r = p.apply_async(compare_clone_to_ref, (clone, ref), context)
            results.append(r)
    p.close()
    p.join()
    result_values = []
    for r in results:
        result_values.append(r.get())
    for clone_name, g in groupby(result_values, key=itemgetter(0)):
        group = list(g)
        matches = map(itemgetter(5), group)
        if any(matches):
            for length in group:
                cx, rx, lx, lx2, aln, good = length
                if good:
                    print 'Matched %s to %s with mismatches (%s / %s)' % (
                                      clone_name, rx, lx, lx2)
                    print aln
            continue
        lengths = sorted(group, key=itemgetter(2), reverse=True)
        if lengths[0][2] == 0:
            print 'No match for %s' % clone_name
        else:
            max_length = lengths[0][2]
            print 'Partially matched %s to %s with length %s / %s' % (
                              clone_name, lengths[0][1], max_length,
                              lengths[0][3])
            print lengths[0][4]
            for length in lengths[1:]:
                cx, rx, lx, lx2, aln, good = length
                if lx == max_length:
                    print 'Also partially matched %s to %s with length %s / %s' % (
                                      clone_name, rx, lx, lx2)
                    print aln
                else: break
    return

def announce_first(clone, logging_level=20, **kwargs):
    logger = get_logger(logging_level)
    logger.info('Comparing %s to references', clone.Name)

def compare_clone_to_ref(clone, ref, logging_level=20, **kwargs):
    logger = get_logger(logging_level)
    aln = CloneToRefAlignment(clone, ref)
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
    return (clone.Name, ref.Name, len(aln),
            len(aln.degap().Seqs[1]), str(aln), not aln.is_truncated and not aln.has_gaps)

if __name__=='__main__': main()