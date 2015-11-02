#!/usr/bin/env python

"""
Module docstring.
"""

import sys
import optparse
from subprocess import call, Popen, PIPE
from tempfile import NamedTemporaryFile
import time
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

    # define options here:
    parser.add_option(
        '-1', '--p5_reads',
        help="Name of fastq(.gz) file containing the 5' reads (R1)")
    parser.add_option(
        '-2', '--p3_reads',
        help="Name of fastq(.gz) file containing the 3' reads (R2)")
    parser.add_option(
        '-f', '--fasta_genome',
        help='An indexed fasta file.')
    parser.add_option(
        '-p', '--output_prefix',
        help='Output files prefix.')
    parser.add_option(
        '-g', '--gtf',
        help='GTF file to be used for counting.')
    parser.add_option(
        '-n', '--mismatches', type='int', default=2,
        help='Number of mismatches for bwa aln command.')
    parser.add_option(
        '-b', '--bwa_cmd', default='/usr/bin/bwa',
        help='bwa executable.')
    parser.add_option(
        '-S', '--samtools_cmd', default='samtools',
        help='Samtools executable.')
    parser.add_option(
        '-a', '--params_aln', default='-t 8 -k 1 -R 200 -l 20',
        help='Additional parameters for aln function of bwa.')
    parser.add_option(
        '-s', '--sampe_params', default='-a 1500 -P',
        help='Additional parameters for sampe function of bwa.')
    parser.add_option(
        '--samse_params', default=' ',
        help='Additional parameters for samse function of bwa.')
    parser.add_option(
        '-o', '--output_dir', default=".",
        help='Output directory name.')
    parser.add_option(
        '-c', '--count_cmd',
        default='/home/users/assafp/lib/python/RNAseq_pipeline/count_PE_fragments.py',
        help='Executable of count_PE_fragments.py.')
    parser.add_option(
        '-u', '--params_count', default=None,
        help='Additional parameters to be transferred to count_PE_fragments.py')
        
    
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


def run_bwa(bwa_cmd, fname1, fname2, output_dir, output_prefix, mismatches,
            fasta_genome, params_aln, params_sampe, params_samse, samtools_cmd):
    """
    Run bwa on paired or single end fastq files. write a sorted  bam and a bai
    file
    Arguments:
    - `bwa_cmd`: executable of bwa
    - `fname1`: The R1 fastq file
    - `fname2`: The R2 fastq file
    - `output_dir`: Where to write the files to
    - `output_prefix`: Name of bam file (without the extension
    - `mismatches`: Allowed mismatches
    - `fasta_genome`: genome fasta file. Must be indexed!
    - `params_aln`: extra parametrs for aln command
    - `params_sampe`: extra paarmeters for sampe execution
    - `params_samse`: extra parameters for samse execution
    - `samtools_cmd`: samtools command
    """
    # Run aln of bwa on both files
    sai1 = NamedTemporaryFile(dir=output_dir)
    call([bwa_cmd, 'aln', '-n', str(mismatches), params_aln,
          fasta_genome, fname1], stdout=sai1)
    if fname2:
        sai2 = NamedTemporaryFile(dir=output_dir)
        next_cmd = [bwa_cmd, 'aln', '-n', str(mismatches),
                    params_aln, fasta_genome, fname2]
        logging.info("Executing %s", ' '.join(next_cmd))
        call(next_cmd, stdout=sai2)
    # Run sampe on both sai files and convert to bam on the fly

    samtobam = [samtools_cmd, 'view', '-Sb', '-']
    bamsort = [samtools_cmd, 'sort', '-',
               "%s/%s"%(output_dir, output_prefix)]
    if fname2:
        next_cmd = ' '.join([bwa_cmd, 'sampe', params_sampe,
                       fasta_genome, sai1.name, sai2.name,
                       fname1, fname2] + ['|'] + samtobam + ['|'] + bamsort)
        logging.info("Executing %s"%next_cmd)
        call(next_cmd, shell=True)
    else: # Single end
        next_cmd = ' '.join(
            [bwa_cmd, 'samse', fasta_genome, sai1.name, fname1] + \
            ['|'] + samtobam + ['|'] + bamsort)
        logging.info("Executing %s"%next_cmd)
        call(next_cmd, shell=True)
    index_cmd = [samtools_cmd, 'index', "%s/%s.bam"%(output_dir, output_prefix)]
    logging.info("Indexing bam file %s"%' '.join(index_cmd))
    call(index_cmd)


def main(argv=None):
    settings, args = process_command_line(argv)
    run_bwa(settings.bwa_cmd, settings.p5_reads, settings.p3_reads,
            settings.output_dir,
            settings.output_prefix, settings.mismatches, settings.fasta_genome,
            settings.params_aln, settings.sampe_params, settings.samse_params,
            settings.samtools_cmd)
    # Count the reads per gene
    count_out = open("%s/%s.counts"%(settings.output_dir, settings.output_prefix), 'w')
    count_array = [settings.count_cmd, '-g', settings.gtf,
          '-s', "%s/%s.bam"%(settings.output_dir, settings.output_prefix)]
    if settings.params_count:
        count_array.append(settings.params_count)
          
    call(count_array, stdout=count_out)
          
    # application code here, like:
    # run(settings, args)
    return 0        # success

if __name__ == '__main__':
    status = main()
    sys.exit(status)
