#!/usr/bin/env python

"""
This script is intendent to count the number of sequencing fragments that
overlap each gene (or feature) in the genome. This was written with a bacterial
gene expression pattern in mind meaning there are long transcription units
that contain more than one gene. Since some genes could be small they might be
contained inside operons and thus not be represented in the reads (the ends
of the fragment).
Count the fragment for a gene only if the overlap is higher than a predefined
variable (l). If the fragment is mapped to different locations count each of
these locations.
Input is a gtf (gff) file, read only the lines that are "exon" (column 2),
identifier is quoted after gene_id in column 8.

22/10/14: Added single-end support
07/12/14: Added -r option for reverse reads (useful mainly in single-end)
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
        '-g', '--gtf',
        help='GTF file containing the features.')
    parser.add_option(
        '-s', '--samfile',
        help='Input sam or bam file.')
    parser.add_option(
        '-r', '--reverse', action='store_true', default=False,
        help='Count the genes on the reverse strand of the mapping.')
    parser.add_option(
        '-f', '--feature', default='exon',
        help='Name of features to count on the GTF file (column 2).')
    parser.add_option(
        '-i', '--identifier', default='gene_id',
        help='Name of identifier to print (in column 8 of the GTF file).')
    parser.add_option(
        '-o', '--overlap', type='int', default=5,
        help='Minimal required overlap between the fragment and the feature.')
    parser.add_option(
        '-t', '--strict', default=False, action='store_true',
        help="Assign the read to a feature only if it covers only one feature.")
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

def read_gtf(gtf_file, feature, identifier):
    """
    Read a GTF file and return a list in the length of the genome in which each
    position contains a list of features that overlap this position
    Arguments:
    - `gtf_file`: An open gtf_file
    - `feature`: the name of the feature to index
    - `identifier`: The identifier to use from column 8
    """
    # First initialize a dictionary
    pos_feat = defaultdict(lambda: defaultdict(set))
    # Get all the names of the features
    all_features = set(['~~intergenic', '~~antisense'])
    for line in csv.reader(gtf_file, delimiter='\t'):
        if line[2] != feature:
            continue
        ids_dict = {}
        for id_pair in line[8].strip().split(';'):
            try:
                k, v = id_pair.strip().split(' ')
            except ValueError:
                pass
            ids_dict[k] = v.replace('"','')
        fid = ids_dict[identifier]
        all_features.add(fid)
        # Change to 0-based coordinates and add this feature to all the
        # positions it convers
        for i in range(int(line[3])-1, int(line[4])):
            # Concat the strand tot he name of the chromosome
            pos_feat[line[0]+line[6]][i].add(fid)
    # Change the dictionary to list
    pos_feat_list = {}
    for chrom, data in pos_feat.items():
        maxpos = max(data.keys())
        list_of_sets = []
        for k in range(maxpos+1):
            list_of_sets.append(list(data[k]))
        pos_feat_list[chrom] = list_of_sets
    return pos_feat_list, all_features


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




def count_features(features_lists, samfile, overlap, rev=False, strict=False):
    """
    Go over the samfile and for each pair of reads find the features that
    overlap the fragment with at least 'overlap' nucleotides. Add 1 to the count
    of these features
    
    Arguments:
    - `features_lists`: The list of features returned from the read_gtf function
    - `samfile`: A pysam object
    - `overlap`: The minimal overlap between the feature and read
    - `rev`: reverse the strand of the read
    """
    fcounts = defaultdict(int)
    counter = 0
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
        counter += 1
        if counter%100000==0:
            sys.stderr.write("Processed %i fragments\n"%counter)
        try:
            chrname = samfile.getrname(read.tid)
        except ValueError:
            sys.stderr.write(str(read)+"\n")
        # Get the positions of the fragment
        if read.is_paired:
            strand, fpos, tpos = get_paired_pos(read, rev=rev)
        else:
            strand, fpos, tpos = get_single_pos(read, rev=rev)
        
        # Count the number of times a feature intersects with the fragmen
        rcounts = defaultdict(int)
        for fset in features_lists[chrname+strand][fpos:tpos]:

            for el in fset:
                rcounts[el] += 1
        # Go over the list of features, if the number of counts is above the
        # Threshold add 1 to the count of this feature
        is_counted = False
        if strict and len(rcounts)>1:
            fcounts['~~multiple'] += 1
            continue
        for feature, counts in rcounts.items():
            if counts >= overlap:
                fcounts[feature] += 1
                is_counted = True
        if not is_counted:
            # Test if antisense
            rev_str = '-'
            if strand == '-':
                rev_str = '+'
            rev_counts = defaultdict(int)
            for fset in features_lists[chrname+rev_str][fpos:tpos]:
                for el in fset:
                    rev_counts[el] += 1
            is_antis = False
            for feature, counts in rev_counts.items():
                if counts >= overlap:
                    is_antis = True
                    break
            if is_antis:
                fcounts['~~antisense'] += 1
            else:
                fcounts['~~intergenic'] += 1
    return fcounts


def main(argv=None):
    settings, args = process_command_line(argv)
    features_lists, all_features = read_gtf(
        open(settings.gtf), settings.feature, settings.identifier)
    samfile = pysam.Samfile(settings.samfile)
    fcounts = count_features(
        features_lists, samfile, settings.overlap, rev=settings.reverse,
        strict=settings.strict)
    for k in sorted(list(all_features)):
        print("%s\t%i"%(k, fcounts[k]))
        
    # application code here, like:
    # run(settings, args)
    return 0        # success

if __name__ == '__main__':
    status = main()
    sys.exit(status)
