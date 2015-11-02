#!/usr/bin/env python

"""
Read a list of fastq (or fastq.gz) files of paired-end or single-end
experiments, split them according to their indices and/or barcodes given in
files in the format:
name[tab]i7_sequence[tab]i5_sequence[tab]internal_barcode 
if either i7 or i5 weren't read or not important leave these cells empty.
For instance, if only internal barcodes were used the table should look like:
exp1_rep1[tab][tab][tab]ACGATGA
exp1_rep2[tab][tab][tab]AATCCAG
...
if only i7 was read:
exp1_rep1[tab]ACATG
exp1_rep2[tab]GATGA
...
etc.


Added the option for one mismatch (including N) using -n

The internal barcode only (last column of table) can be in different length
"""
if __name__ == '__main__' and __package__ is None:
    from os import sys, path
    sys.path.append(path.dirname(path.dirname(path.abspath(__file__))))
    
import sys
import optparse
import gzip
import logging
import FastSeqIO
import subprocess

if sys.version.startswith("3"):
    import io
    io_method = io.BytesIO
else:
    import cStringIO
    io_method = cStringIO.StringIO


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

    # define options here:
    parser.add_option(
        '-1', '--files1',
        help='A comma separated list of files of the first pair (containing the'
        ' barcode.')
    parser.add_option(
        '-2', '--files2',
        help='The pairs of the input files int he same order.')
    parser.add_option(
        '-b', '--barcodes_file',
        help='Barcodes table in the format:'
        'name\ti7_sequence\ti5_sequence\tbarcode or one/some of them.')
    parser.add_option(
        '-m', '--mismatch', action='store_true', default=False,
        help='Allow one mismatch.')
    parser.add_option(
        '-o', '--output_dir', default='.',
        help='Directory to write files to.')
    parser.add_option(
        '-u', '--unknown', action='store_true', default=False,
        help='For unknown barcodes make a file for each barcode instead of'
        ' putting them in one big file. Beware that there is a limit on the'
        ' number of open files, when this limit reaches no more files will'
        ' be opened.')
    parser.add_option(
        '-s', '--strict', default=False, action='store_true',
        help="Don't write reads with no barcode match, overrides --unknown")
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

def read_barcodes(barfile, mismatch=False):
    """
    Read the barcodes file and return a dictionary with the concatenated
    sequence of all barcodes (i7, i5, internal) as key and lib name as value.
    When two mismatched barcodes interfere remove them from the dictionary.
    If a primary barcode (without mismatches) interfere with a mismatched
    barcode remove the later.
    Arguments:
    - `barfile`: The file containing the table
    - `mismatch`: Allow one mismath
    Return:
    - `barcodes`: A dictionary from seq to library name
    - `cols_lens`: i7, i5 and barcode sequence lengths
    """
    barcodes = {}
    # Mismatched indices to remove - they fit more than one barcode
    toremove = set()
    # initial barcodes - they wont be removed
    initials = set()
    cols_lens = [set(), set(), set()]
    with open(barfile) as barin:
        for line in barin:
            if line.startswith('#'):
                continue
            bseqs = line.strip().split('\t')
            name = bseqs.pop(0)
            for i, c in enumerate(bseqs):
                cols_lens[i].add(len(c))
            seq = ''.join(bseqs)
            initials.add(seq)
            name = name.replace(' ', '_').replace('/', '_').replace(':', '_')
            barcodes[seq] = name
            if mismatch:
                for i in range(len(seq)):
                    for c in ('A', 'C', 'T', 'G', 'N'):
                        if seq[i] != c:
                            s2 = seq[:i] + c
                            if i+1 < len(seq):
                                s2 += seq[i+1:]
                            if s2 in barcodes:
                                logging.warn(
                                    "Mismatched barcode %s matches libraries %s"
                                    "and %s\n"%(s2, barcodes[s2],  name))
                                toremove.add(s2)
                            else:
                                barcodes[s2] = name
    for s in toremove-initials:
        del barcodes[s]
    return barcodes, cols_lens


def open_file(fname):
    """
    Open a fastq or fastq.gz file and return a fastIO object
    Arguments:
    - `fname`: file name to open
    """
    try:
        if fname.endswith('.gz'):
            p = subprocess.Popen(["zcat", fname], stdout = subprocess.PIPE)
            fin = p.stdout
#            fin = io_method(p.communicate()[0])
            # fin = gzip.open(fname)
        else:
            fin = open(fname)
    except IOError:
        logging.error("Can't open file %s for reading."%fname)
    f1_seqs = FastSeqIO.fastIO(fin)
    return f1_seqs


def get_possible_bcseqs(read, has_i7, has_i5, has_ib, blens):
    """
    Return a set of possible barcode sequences. If blens has only one
    value (including 0) one possible barcode will be returned, if it has
    multiple values, more than one will be returned
    Arguments:
    - `read`: A Read object
    - `has_i7`: i7 sequences were read from the table
    - `has_i5`: i5 sequences were read from the table
    - `has_ib`: internal barcodes were read from the table
    - `blens`: Possible lengths of internal barcodes
    """
    bcseq = ''
    if has_i7 and has_i5:
        try:
            bc1, bc2 = read.name.split(':')[-1].split('+')
            bcseq = bc1 + bc2
        except ValueError: # i7 and/or i5 wasn't read
            logging.error("i7 and i5 were given at the index table"
                          " but are not read by the machine")
            return None
    elif has_i7:
        bc1 = read.name.split(':')[-1]
        if '+' in bc1: # Both read but i5 is not neede
            bc1 = bc1.split('+')[0]
        bcseq = bc1
    elif has_i5: # i7 is read but doesn't matter
        bc2 = read.name.split(':')[-1].split('+')[1]
        bcseq = bc2
    len_i75 = len(bcseq)
    all_possible_seqs = set()
    if has_ib:
        for bl in blens:
            all_possible_seqs.add(bcseq + read.seq[:bl])
    else:
        all_possible_seqs.add(bcseq)
    return all_possible_seqs, len_i75

def split_reads(
    files1, files2, barcodes_file, mismatch=False, outdir=".", unknown=False,
    strict=False, bufsize=10000000):
    """
    The main function (can be called from other scripts)
    Arguments:
    - `files1`: A list of fastq files of read 1
    - `files2`: A list of fastq files of the mate reads
                if single end set as None
    - `barcodes_file`: A file with the indices and barcodes table
    - `mismatch`: Allow one mismatch in the barcode sequence
    - `outdir`: Where to write the files to
    - `unknown`: If encountered a barcode that is not recognized write
                 the sequence to a file with the barcode sequence as file name.
                 in default write all of the sequences to the same file.
    - `strict`: disregard reads with no matching index
    - `bufsize`: Size of buffer (in bytes), it's higher than normal since
                 there is a lot of memory and writing is expensive
    """
    barcodes, cols_lens = read_barcodes(barcodes_file, mismatch=mismatch)
    f2open = True
    if not files2:
        f2open = False
        files2 = files1 # For simplicity
    barfiles = {}
    for bcname in set(barcodes.values()):
        try:
            barfiles[bcname] = [open(
                "%s/%s_1.fastq"%(outdir,bcname), 'w', bufsize)]
            if f2open:
                barfiles[bcname].append(
                    open("%s/%s_2.fastq"%(outdir,bcname), 'w', bufsize))
        except IOError:
            logging.critical("Can't open fastq files for writing.")
            raise
    if not strict:
        if unknown:
            unseq_files = {}
        else:
            unm_files = [
                open('%s/nobarcode_reads_1.fastq'%outdir, 'w', bufsize)]
            if f2open:
                unm_files.append(
                    open('%s/nobarcode_reads_2.fastq'%outdir, 'w', bufsize))
    # Different lengths of barcodes
    blens = cols_lens[2]
    # Which barcodes to look for: i7, i5 and internal (ib)
    has_i7 = any(cols_lens[0])
    has_i5 = any(cols_lens[1])
    has_ib = any(cols_lens[2])
    for f1, f2 in zip(files1, files2):
        f1_seqs = open_file(f1)
        if f2open:
            f2_seqs = open_file(f2)
        for s1 in f1_seqs:
            if f2open:
                try:
                    s2 = f2_seqs.next()
                except StopIteration:
                    logging.error(
                        "Uneven number of reads in R1 and R2 in file %s"%f1)
                    break
            all_possible_bcseqs, i75_len = get_possible_bcseqs(
                s1, has_i7, has_i5, has_ib, blens)
            # Test if the barcode sequences match one (or more) from the table
            bcmatch = None
            for bcs in all_possible_bcseqs:
                if bcs in barcodes:
                    if bcmatch:
                        logging.warn("read %s mapped to more than one barcodes"
                                     " removing from library"%s1.name)
                        bcmatch = None
                        break
                    else:
                        bcmatch = bcs
            if bcmatch:
                s1[len(bcmatch)-i75_len:].write(barfiles[barcodes[bcmatch]][0])
                if f2open:
                    s2.write(barfiles[barcodes[bcmatch]][1])
            elif not strict:
                if unknown:
                    # Grab the longest barcode and write to file
                    wrt_bcseq = max(all_possible_bcseqs, key=len)
                    if wrt_bcseq not in unseq_files:
                        try:
                            unseq_files[wrt_bcseq] = [open(
                                "%s/%s_1.fastq"%(outdir,wrt_bcseq),'w',bufsize)]
                        except IOError:
                            # too many open files probably, ignore and go on
                            logging.error("Too many open files, won't open")
                            unseq_files[wrt_bcseq] = None
                        else:
                            if f2open:
                                try:
                                    unseq_files[wrt_bcseq].append(open(
                                        "%s/%s_2.fastq"%(outdir,wrt_bcseq),
                                        'w', bufsize))
                                except IOError:
                                    logging.error(
                                        "Too many open files, won't open")
                                    unseq_files[wrt_bcseq] = None
                    if unseq_files[wrt_bcseq]:
                        s1[len(wrt_bcseq)-i75_len:].write(
                            unseq_files[wrt_bcseq][0])
                        if f2open:
                            s2.write(unseq_files[wrt_bcseq][1])
                else:
                    s1.write(unm_files[0])
                    if f2open:
                        s2.write(unm_files[1])
    # Close all the output files
    for sl in barfiles.values():
        for sf in sl:
            sf.close()
    if not strict:
        if unknown:
            for sl in unseq_files.values():
                for sf in sl:
                    sf.close()
        else:
            for sf in unm_files:
                sf.close()
    
def main(argv=None):
    settings, args = process_command_line(argv)
    if not settings.files2:
        f2list = None
    else:
        f2list = settings.files2.split(',')
    split_reads(
    settings.files1.split(','), f2list, settings.barcodes_file,
    settings.mismatch, settings.output_dir, settings.unknown, settings.strict)
    return 0        # success

if __name__ == '__main__':
    status = main()
    sys.exit(status)
