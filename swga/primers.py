# -*- coding: utf-8 -*-
"""primers.py

Contains the table definitions for the Primer/Set database, as well as functions
that help manage the database. Also contains functions for counting (k)mers in a
genome.

"""
from __future__ import with_statement, division
import subprocess
import os
import re
import struct
import swga
import swga.core
import swga.locate
import swga.resources as resources
from swga.database import Primer
from pyfaidx import Fasta


def count_kmers(k, genome_fp, cwd, threshold=1):
    assert isinstance(threshold, int)
    dsk = resources.get_dsk()
    genome = genome_fp.split(os.sep).pop()
    out = '%s-%dmers' % (genome, k)
    outfile = os.path.join(cwd, out + '.solid_kmers_binary')
    if os.path.isfile(outfile):
        swga.message("Binary kmer file already exists at %s; parsing..."
                     % outfile)
    else:
        cmdstr = ("{dsk} {genome_fp} {k} -o {out} -t {threshold}"
                  .format(**locals()))
        swga.message("In {cwd}:\n> {cmdstr}".format(**locals()))
        try:
            subprocess.check_call(cmdstr, shell=True, cwd=cwd)
        except:
            if os.path.isfile(outfile):
                os.remove(outfile)
            raise
    primers = dict((kmer, freq)
                   for kmer, freq in parse_kmer_binary(outfile)
                   if freq >= threshold)
    return primers
        

def parse_kmer_binary(fp):
    # Adapted from `dsk/parse_results.py`
    with open(fp, 'rb') as f:
        try:
            kmer_nbits = struct.unpack('i', f.read(4))[0]
            k = struct.unpack('i', f.read(4))[0]
        except struct.error:
            if os.path.isfile(fp):
                os.remove(fp)
            raise
        try:
            while True:
                kmer_binary = struct.unpack('B' * (kmer_nbits // 8),
                                            f.read(kmer_nbits // 8))
                freq = struct.unpack('I', f.read(4))[0]
                kmer = ""
                for i in xrange(k):
                    kmer = "ACTG"[(kmer_binary[i//4] >> (2 * (i%4))) % 4] + kmer
                yield kmer, freq
        except struct.error:
            pass

        
def read_primer_list(lines, fg_genome_fp, bg_genome_fp):
    '''
    Reads in a list of primers, one per line, and returns the corresponding
    records from the primer database. If the primer doesn't exist in the db,
    tries to create it manually. If the primer doesn't appear in the fg genome,
    it skips it with a warning.
    '''
    seqs = [re.split(r'[ \t]+', line.strip('\n'))[0] for line in lines]
    primers = list(Primer.select().where(Primer.seq << seqs).execute())
    if len(primers) < len(seqs):
        primer_seqs = [p.seq for p in primers]
        missing = [_ for _ in seqs if _ not in primer_seqs]
        for seq in missing:
            try:
                primers.append(create_primer_from_seq(seq, fg_genome_fp,
                                                      bg_genome_fp)) 
            except swga.core.SWGAError:                
                swga.warn(seq + " does not bind to the foreground genome; skipping.")
    return primers


def create_primer_from_seq(seq, fg_genome_fp, bg_genome_fp):
    """
    Instead of using DSK, count the number of times a primer occurs in the
    foreground and genome background manually. 
    For user-entered primer sets that contain primers that don't exist in the
    database already. 
    """
    
    fg = Fasta(fg_genome_fp)
    bg = Fasta(bg_genome_fp)
    fg_freq = swga.locate.count(seq, fg)
    if fg_freq == 0:
        raise swga.core.SWGAError("Primer not found in foreground genome")
    bg_freq = swga.locate.count(seq, bg) if fg_freq > 0 else 0
    ratio = fg_freq/float(bg_freq) if bg_freq > 0 else 0
    return Primer(seq=seq, fg_freq=fg_freq, bg_freq=bg_freq,
                  ratio=ratio).execute() 
