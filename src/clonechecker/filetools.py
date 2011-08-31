'''
Created on Aug 20, 2011

@author: ben
'''
from cogent.parse.record import FileFormatError
from cogent import LoadSeqs, DNA
import os.path
from twobitreader import TwoBitFile
from clonechecker.TabFile import BedFile
import platform
from errno import ENOENT


def load_one_seq(foo, moltype=DNA):
    """
    load one sequence from path foo if it is a FASTQ or plain-text file
    """   
    try:
        seqs = LoadSeqs(foo, aligned=False, format='fasta', moltype=moltype)
    except FileFormatError:
        fh = open(foo, 'rU')
        s = fh.read().strip()
        name = os.path.basename(foo)
        seqs = LoadSeqs(data=[(name, s)], aligned=False, moltype=moltype)
    return seqs.Seqs[0]

def load_seqs(foo, moltype=DNA):
    """
    load one or more sequences from path foo if it is a FASTQ or plain-text file
    
    returns a list of Seq objects
    """   
    try:
        seqs = LoadSeqs(foo, aligned=False, format='fasta', moltype=moltype)
    except FileFormatError:
        fh = open(foo, 'rU')
        s = fh.read().strip()
        name = os.path.basename(foo)
        seqs = LoadSeqs(data=[(name, s)], aligned=False, moltype=moltype)
    return [seq for seq in seqs.Seqs]

def write_to_fasta(file_handle, sequence, name=''):
    """write a record to a fasta file without closing it"""
    # write first line
    file_handle.write('>{!s}\n'.format(name))
    # write sequence
    num_lines = (len(sequence)+ 59 ) / 60
    for j in xrange(num_lines):
        start = 0 + j*60
        end = 60 + j*60
        file_handle.write('{!s}\n'.format(sequence[start:end]))
    return

def read_bed_file(foo, moltype=DNA, genome=None, write_fasta=True):
    """read in all regions from a BED file"""
    t = TwoBitFile(genome)
    bf = BedFile(foo)
    seqs = []
    if write_fasta: fasta_output = open(os.path.basename(foo) + '.fa', 'w')
    for row in bf:
        chr, start, end = row.chrom(), row.chrom_start(), row.chrom_end()
        region_seq = t[chr][long(start):long(end)]
        region_name = row.name() + ' (%s:%s-%s)' % (chr, start, end)
        if write_fasta: write_to_fasta(fasta_output, region_seq, name=region_name)
        s = LoadSeqs(data=[(region_name, region_seq)], aligned=False,
                     moltype=moltype)
        seqs.append(s.Seqs[0])
    return seqs

def find_2bit_file(ref_genome, path_to_gbdb=None):
    """find a .2bit file with genome named ref_genome"""
    if ref_genome is None: raise ValueError("No reference genome specified")
    if os.path.exists(ref_genome): return ref_genome
    fname = "{!s}{!s}2bit".format(ref_genome, os.extsep)
    if os.path.exists(fname): return fname
    if path_to_gbdb is not None:
        # check in its own dir
        specified_path = os.path.join(path_to_gbdb, ref_genome, fname)
        if os.path.exists(specified_path): return specified_path
        # check in the main dir
        specified_path = os.path.join(path_to_gbdb, fname)
        if os.path.exists(specified_path): return specified_path
    # try absolute paths
    new_path = os.path.join("gbdb", ref_genome, fname)
    if platform.system() == 'Windows':
        home_drive = os.environ['HOMEDRIVE']
        windows_path = "{!s}{!s}{!s}".format(home_drive, os.sep, new_path)
        if os.path.exists(windows_path): return windows_path
    else: # assume *nix
        unix_path = os.sep + new_path
        if os.path.exists(unix_path): return unix_path
    # give up
    raise IOError(ENOENT, os.strerror(ENOENT), fname)