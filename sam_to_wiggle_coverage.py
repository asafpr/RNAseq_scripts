#!/usr/bin/env python

"""
Read a sam file and generae a wiggle file with the coverage on both strands
"""

import sys
import optparse
import pysam
from collections import defaultdict
import csv
import logging

def process_command_line(argv):
    """
    Return a 2-tuple: (settings object, args list).
    `argv` is a list of arguments, or `None` for ``sys.argv[1:]``.
    """
    if argv is None:
        argv = sys.argv[1:]

    # initialize the parser object:
    parser = optparse.OptionParser(
        formatter=optparse.TitledHelpFormatter(width=78),
        add_help_option=None)


    parser.add_option(
        '-s', '--samfile',
        help='Input sam or bam file.')
    parser.add_option(
        '-t', '--title',
        help='Title of wig file.')
    parser.add_option(
        '-d', '--description',
        help='Description of wiggle file.')
    parser.add_option(
        '-r', '--reverse', action='store_true', default=False,
        help='Reverse strand of mapping.')
    parser.add_option(
        '-f', '--first_pos', action='store_true', default=False,
        help='Plot the sum of reads that start in each position instead of '
        'number of reads covering the position.')
    parser.add_option(      # customized description; put --help last
        '-h', '--help', action='help',
        help='Show this help message and exit.')

    settings, args = parser.parse_args(argv)

    # check number of arguments, verify values, etc.:
    if args:
        parser.error('program takes no command-line arguments; '
                     '"%s" ignored.' % (args,))

    # further process settings & args if necessary

    return settings, args



def get_paired_pos(read, rev=False):
    """
    Given a read which is the positive of the pairs return a strand and
    start and end positions
    Arguments:
    - `read`: A read from pysam
    """
    strand = '+'
    if rev!=read.is_read2:
        strand = '-'
    fpos = read.pos
    tpos = read.tlen + fpos
    return strand, fpos, tpos

def get_single_pos(read, rev=False):
    """
    Given a read which is the positive of the pairs return a strand and
    start and end positions
    Arguments:
    - `read`: A read from pysam
    """
    strand = '+'
    if rev!=read.is_reverse:
        strand = '-'
    fpos = read.pos
    tpos = read.qlen + fpos
    return strand, fpos, tpos




def generate_wig(samfile, rev=False, first_pos=False):
    """
    Go over the samfile and return two histograms (for + and - strands) of
    coverage
    
    Arguments:
    - `samfile`: A pysam object
    - `rev`: reverse the strand of the read
    - `first_pos`: Count only the first position of each read
    """
    # Build the structure of the dictionary chromosome->strand->list of 0
    coverage = {}
    for i, rfg in enumerate(samfile.references):
        rlen = samfile.lengths[i]
        coverage[rfg] = {'-':[0] * rlen, '+':[0] * rlen}
    for read in samfile.fetch():
        if read.is_paired:
            if read.is_reverse or read.is_unmapped or\
                read.mate_is_unmapped or\
                read.is_reverse==read.mate_is_reverse or\
                not read.is_proper_pair:
                continue
        else: # single end
            if  read.is_unmapped:
                continue
        # Take only the forward mate
        try:
            chrname = samfile.getrname(read.tid)
        except ValueError:
            logging.warn("Read has no valid chr name %s"%(str(read)))
            continue
        # Get the positions of the fragment
        if read.is_paired:
            strand, fpos, tpos = get_paired_pos(read, rev=rev)
        else:
            strand, fpos, tpos = get_single_pos(read, rev=rev)
        rrange = range(fpos, tpos)
        if first_pos:
            if strand == '+':
                rrange = [fpos]
            else:
                rrange = [tpos]
        for i in rrange:
            try:
                coverage[chrname][strand][i] += 1
            except IndexError:
                logging.warn("IndexError: trying to set index %d on chr %s, bu length is only %d"%(i, chrname, len(coverage[chrname][strand])))
    return coverage


def print_wiggle(coverage, title, description, outf):
    """
    Print the coverage into an open wiggle file
    Arguments:
    - `coverage`: returned from generate_wig
    - `title`: Title of wig track
    - `description`: Description of track
    - `outf`: Open file to write to
    """
    for chr_name in coverage:
        outf.write('track type=wiggle_0 name=%s_PLUS description="%s PLUS" visibility=full color=0,0,255\n'%(
            title, description))
        outf.write("fixedStep chrom=%s start=1 step=1\n"%chr_name)
        for c in coverage[chr_name]['+']:
            outf.write("%d\n"%c)
        

        outf.write('track type=wiggle_0 name=%s_MINUS description="%s MINUS" visibility=full color=255,0,0\n'%(
            title, description))
        outf.write("fixedStep chrom=%s start=1 step=1\n"%chr_name)
        for c in coverage[chr_name]['-']:
            outf.write("%d\n"%c)
    

def main(argv=None):
    settings, args = process_command_line(argv)
    samfile = pysam.Samfile(settings.samfile)
    coverage = generate_wig(
        samfile, rev=settings.reverse, first_pos=settings.first_pos)
    print_wiggle(coverage, settings.title, settings.description, sys.stdout)        
        
    # application code here, like:
    # run(settings, args)
    return 0        # success

if __name__ == '__main__':
    status = main()
    sys.exit(status)
