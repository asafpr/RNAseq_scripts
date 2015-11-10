#!/usr/bin/env python

"""
This script takes a nextseq run which its bcl files are under /cs/nextseq/ and
automatically process it with:
1. Split the run using indices (i7, i5) or/and barcodes in the 5' region
2. Adapter removal and low quality 3' end removal with cutadapt
3. fastqc files, copied to the www directory (accessible from home)
4. Mapping to the genome with bwa
5. Count the number of reads per gene
6. Generate wig file with all the tracks, send an email with a link to upload
   the wig file to microbial UCSC genome browser
An email will be sent after the run is over with a summary and a link to the
fastq html files and wig files.

The files will be written to ~asafp5/lab_nextseq_complete_data/data/[experimentalist name]/[experiment name]/[Date]/
It will have the structure:
/fastq_files - fastq.gz files after barcode splitting and adapter trimming
               etc. reads with no barcode will not be saved (they can be
               recovered with bcl2fastq)
/fastqc_results - html and zip files of fastqc. They will have links from
                  www directory under ~asafp5/www/lab_nextseq_data/
/bwa_mapping - bam, bam.bai, wig.gz and counts files including a table that
               summarizes all the counts for all the libraries in the run
               resulting from bwa mapping. The wig file will be linked to
               ~asafp5/www/lab_nextseq_data

Input:
* In order to run this script and keep the data organized you should first
prepare a barcode index file in the following format:
library name[tab]i7_sequence[tab]i5_sequence[tab]internal_barcode
If you used only internal barcodes (they are the 5' ends of the reads) you
should leave columns 2 and 3 empty. If you read indices i7 and i5 but the
sequence of i7 doesn't matter you should only give a sequence in column 3.
The sequences should be the sequences read by the machine so sometimes they are
the reverse complement of the oligos used in the experiment.
Using Jonathan Livny's protocol there's an extra T after the barcode, it should
be included as well in the sequence. 
* The next thing you need to know is the exact date the run took place, you
should enter it in the format YYMMDD. For instance a run from July 15th 2014
will be 140715.
* Adapter sequences - If your cDNA is short and the read is long you will
probably run into the adapter at the end of the read. You should give the
sequence of the read adapter (which is the reverse-complement of the adapter)
with -a, if the sequencing is paired-end insert the 5' adapter with -A. If you
used standard Illumina adapters you will see if they are read in the fastqc
results.
* Reverse complement - Using Jonathan Livny's protocol the reads are the reverse
complement of the actual RNA, use -r to treat them as such so the counts per
gene and wig file will be correct, the barcode sequences, as mentioned above,
should be rev-comped manually.
* To keep things organized the name of the experimentalist and the experiment
should be the first two parameters.

How to run:
process_nextseq_run experimentalist experiment_name date genome_fsata genome_gff
                    -i index_file -a adapter_3' -A adapter_5' [-r]
                    -e yourmail@host,sysadminmail@host...
genome_fasta -  a fasta file with the genome sequence indexed by bwa. The genome
             of E. coli is under referece_genomes/genome.fa
genome_gff - a file with gene positions, the one of E. coli is in:
           reference_genomes/refseq_ver2_genes_chr.gff
email addresses - it will send emails with a summary including links to the
                  fastqc html, the wig file to upload to microbial UCSC browser
                  and the location of the files on the disk.

Special cases:
Libraries from other sources:
* Libraries which are already split to fastq files should be put in the
fastq_files directory, the options: --skip_bcl2fastq and --skip_split should be
used. Make sure the files are fastq and not fastq.gz, if they are paired they
should be _1.fastq and _2.fastq
* Libraries which are not split should come in files Undetermined*_R1*.fastq.gz
and Undetermined*_R2*.fastq.gz. They should be placed under the base directory
(under the date) and --skip_bcl2fastq should be used.
* Remapping: If you want to remap the libraries only use --remap, you can use
different parameters like --allowed_mismatches (number of allowed mismatches)
or --force_single to use only the first read if you see that the pair is useless
(full of poly-G or poly-N)

"""
if __name__ == '__main__' and __package__ is None:
    from os import sys, path
    sys.path.append(path.dirname(path.dirname(path.abspath(__file__))))
sys.path.append("/cs/icore/asafp5/lib/python") # For pysam

import os
try:
    os.environ["PYTHONPATH"] += ":/usr/local/icore-hm/x86_64.debian64_5775/python2.7/lib/python2.7/site-packages:/cs/icore/asafp5/lib/python"
except KeyError:
    os.environ["PYTHONPATH"] = ":/usr/local/icore-hm/x86_64.debian64_5775/python2.7/lib/python2.7/site-packages:/cs/icore/asafp5/lib/python"
#print os.environ["PYTHONPATH"]
import sys
import argparse
import logging
import time
import glob
from subprocess import call
from tempfile import NamedTemporaryFile, SpooledTemporaryFile
import pysam
import gzip
import csv
import stat
import shutil

from smtplib import SMTPAuthenticationError
import index_splitter
import map_and_count
import executable_defaults as defs
import count_PE_fragments
import sam_to_wiggle_coverage
import sendmail

def process_command_line(argv):
    """
    Return a settings object.
    `argv` is a list of arguments, or `None` for ``sys.argv[1:]``.
    """
    if argv is None:
        argv = sys.argv[1:]

    # initialize the parser object, replace the description
    parser = argparse.ArgumentParser(
        description='Process a nextseq run',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument(
        'benchman',
        help='Experimentalist name, will be the first subdirectory name.')
    parser.add_argument(
        'experiment_name',
        help='Will be the second subdirectory under [benchman].')
    parser.add_argument(
        'date',
        help='The date of the nextseq run. Format: YYMMDD.'
        ' Will be used as the third'
        ' subdirectory under experiment name and also to find your run under'
        ' /cs/nextseq/. If more than one run have the same date you will be '
        ' prompt to choose the right one (look at basespace to see the flowcell'
        ' id and select the right one).')
    parser.add_argument(
        'genome_fasta',
        help='Name of genome fasta file. The file must be indexed using'
        ' bwa index command prior to this run.')
    parser.add_argument(
        'genes_gff',
        help='Name of gff file to count the reads per gene.')
    parser.add_argument(
        '-r', '--reverse_complement', default=False,
        action='store_true',
        help='Treat the reads as reverse complement only when counting'
        ' number of reads per gene and generating wig files. The resulting BAM'
        ' files will be the original ones.')
    parser.add_argument(
        '-i', '--index_table',
        help='A tab-delimited file containing library name in first column, '
        'i7 index sequence in the second column, i5 sequence in the third '
        "and internal barcode sequence (at the 5' end). optional."
        " If i7, i5 or internal barcode are not used leave the column blank")
    parser.add_argument(
        '-e', '--email_address',
        help='Send a short summary to these email addresses when the run ends.'
        ' The addresses should be comma separated.')
    parser.add_argument(
        '-a', '--adapter_seq', default=defs.ADAPTER_1,
        help="Sequence of adapter to remove from 3' end of read."
        " If multiple sequences are given (comma separated) use all of them.")
    parser.add_argument(
        '-A', '--adapter_pair_seq', default=defs.ADAPTER_2,
        help='Adapter sequence of the paired read.')
    parser.add_argument(
        '-m', '--minimal_length', type=int, default=defs.MINLEN,
        help='Minimal lenth of read to keep after cutadapt. default: %d'%defs.MINLEN)
    parser.add_argument(
        '-q', '--quality_cutoff', type=int, default=defs.QUAL_CUTOFF,
        help='Used by cutadapt and BWA if no adapter is given to remove '
        "nucleotides with low quality from the 3' end. default: %d"%defs.QUAL_CUTOFF)
    parser.add_argument(
        '--choose', type=int, default=None,
        help='If there are multiple libraries choose this one instead of prompting.')
    parser.add_argument(
        '--dont_delete', default=False, action='store_true',
        help="Don't delete intermediate fastq files.")
    parser.add_argument(
        '--force_dir', default=None,
        help='Extract run data from this directory.')
    parser.add_argument(
        '--skip_bcl2fastq', action='store_true', default=False,
        help='Start from Undetermined*.fastq.gz files.')
    parser.add_argument(
        '--skip_split', action='store_true', default=False,
        help='Skip the index splitting part, start from fastq files.'
        ' you should add --skip_bcl2fastq as well.')
    parser.add_argument(
        '--skip_cutadapt', action='store_true', default=False,
        help='Skip cutadapt part of the analysis.')
    parser.add_argument(
        '--keep_nobarcode', action='store_true', default=False,
        help='Keep the reads with no matching barcode.')
    parser.add_argument(
        '--remap', action='store_true', default=False,
        help='Only remap splitted sequences. Use the cutadapt files.')
    parser.add_argument(
        '--skip_mapping', action='store_true', default=False,
        help='Skip the mapping part, useful with --remap to recount and'
        ' generate wig files.')
    parser.add_argument(
        '--force_single', action='store_true', default=False,
        help='When mapping use only the first mate, ignore the second mate.')
    parser.add_argument(
        '--allowed_mismatches', type=int, default=defs.MISMATCHES,
        help="Allowed mismatches for BWA mapping, default=%s"%defs.MISMATCHES)
    parser.add_argument(
        '--plot_first_position', action='store_true', default=False,
        help='Plot the number of reads starting in each position instead'
        ' of the coverage.')
    parser.add_argument(
        '--dbname', default='eschColi_K12',
        help='The name of the UCSC microbial DB name, ecoli is the default.'
        'Salmonella LT2 for instance is salmTyph_LT2')
    parser.add_argument(
        '--gene_overlap', type=int, default=defs.OVERLAP,
        help='Number of nucleotides that should overlap between the gene and'
        ' the read to add one to the number of reads per gene.')
    parser.add_argument(
        '--bwa_exec', default=defs.BWA,
        help='bwa command, default: %s.'%defs.BWA)
    parser.add_argument(
        '--bwa_aln_params', default=defs.PARAMS_ALN,
        help='Parameters of bwa aln.')
    parser.add_argument(
        '--bwa_samse_params', default=defs.PARAMS_SAMSE,
        help='Parameters for bwa samse (used for single-end only).')
    parser.add_argument(
        '--bwa_sampe_params', default=defs.PARAMS_SAMPE,
        help='Parameters for bwa sampe (used for paired-end only).')
    parser.add_argument(
        '--fastqc_exec', default=defs.FASTQC,
        help='fastqc executable, default: %s.'%defs.FASTQC)
    parser.add_argument(
        '--samtools_cmd', default=defs.SAMTOOLS,
        help='samtools executable.')
    parser.add_argument(
        '--cutadapt_cmd', default=defs.CUTADAPT_CMD,
        help='cutadapt executable.')
    parser.add_argument(
        '--cutadapt_params', default=defs.CUTADAPT_ADD,
        help='Additional parameters for cutadapt.')
    settings = parser.parse_args(argv)

    return settings

# Create directories with subdirectories and chmod them
def supermakedirs(path, mode):
    if not path or os.path.exists(path):
        return []
    (head, tail) = os.path.split(path)
    res = supermakedirs(head, mode)
    try:
        os.mkdir(path)
        os.chmod(path, mode)
    except OSError, ose:
        if ose.errno == 17: # File exists
            logging.warn(
                "Base directory %s already exists, data will be overwritten"%\
                path)
        elif ose.errno == 13: # Permission denied
            logging.critical(
                "Base directory can't be initialized, permission denied: %s"%\
                path)
            raise
    res += [path]
    return res

def make_outdir(
    settings, wwwdir, homedir=defs.DATA,
    subdirs=(defs.FASTQ_DIR, defs.QCRESULTS_DIR, defs.BWA_DIR)):
    """
    Makes the directory if it doesn't exist
    Arguments:
    - `settings`: A argparse object with the parametrs
    - `subdirs`: A list of subdirectories to make
    Returns:
    - `basedir`: The main output directory
    """
    basedir = "%s/%s/%s/%s"%(
        homedir,settings.benchman, settings.experiment_name, settings.date)
    try:
        supermakedirs(basedir, 0777)
    except OSError, ose:
        if ose.errno == 17: # File exists
            logging.warn(
                "Base directory %s already exists, data will be overwritten"%\
                basedir)
        elif ose.errno == 13: # Permission denied
            logging.critical(
                "Base directory can't be initialized, permission denied: %s"%\
                basedir)
            raise
    for subd in subdirs:
        try:
            os.mkdir("%s/%s"%(basedir, subd))
            os.chmod("%s/%s"%(basedir, subd), 0777)
        except OSError:
            pass # A warning was already issued when base dir was recreated

    try:
        os.makedirs(wwwdir)
        os.chmod(
            wwwdir,
            stat.S_IXUSR+stat.S_IRUSR+stat.S_IWUSR+stat.S_IXOTH+stat.S_IROTH)
                                 
    except OSError, ose:
        if ose.errno == 17: # File exists
            logging.warn(
                "www directory %s already exists, data will be overwritten"%\
                wwwdir)
        elif ose.errno == 13: # Permission denied
            logging.error(
                "www directory can't be initialized, permission denied: %s"%\
                wwwdir)

    return basedir


def run_fastqc(files, outdir, wwwdir, fastqc_cmd):
    """
    Run fastqc on all the raw files, write them to outdir and 
    copy to www so they can be viewed from outside.
    Arguments:
    - `files`: A list of fastq(.gz) files
    - `outdir`: Where to write the files to
    - `wwwdir`: Where to link the files to
    - `fastqc_cmd`: fastqc executable
    """
    fc_cmd = [fastqc_cmd, '-o', outdir, '-t', '32']
    fc_cmd.extend(files)
    logging.info("Running: %s"%(' '.join(fc_cmd)))
    call(fc_cmd)
    for ext in defs.FASTQC_FILES:
        for cf in glob.glob("%s/*.%s"%(outdir, ext)):
            try:
                shutil.copy(cf, wwwdir)
            except OSError:
                logging.error("Can't copy file %s"%cf)
            
        
    
    
    
def concat_gzipped(flist, fout, bufsize=10000000):
    """
    Copy the gzipped files in the list to the output file.
    """
    outf = open(fout, 'w', bufsize)
    for fr1 in flist:
        fin = gzip.open(fr1)
        for chunk in iter(lambda: fin.read(bufsize), ''):
            outf.write(chunk)
        fin.close()
    outf.close()

def run_bcl2fatq(
    basedir, date, force_dir=None, choose=None, csbase=defs.BASE,
    command=defs.BCL2FASTQ, defpar=defs.BCLPARAMS):
    """
    Run bcl2fastq to the appropriate dir, the files (Undetermined...fastq.gz)
    will eventually be erased but the rest will remain
    Arguments:
    - `basedir`: output dir home
    - `date`: Date of run
    - `force_dir`: Take data from this dir
    - `choose`: Choose one library if there are multiple
    - `csbase`: base directory of nextseq data (/cs/nextseq/ by default)
    - `command`: The bcl2fastq executable full path
    - `defpar`: default parameters for bcl2fastq
    """
    if force_dir:
        data_dir = force_dir
    else:
        # Detect the right run
        fnames = glob.glob("%s/%s*"%(csbase, date)) +\
                 glob.glob("%s/NGS/%s*"%(csbase, date))
        fnum = 0
        if len(fnames) > 1 and (choose is None):
            # prompt the user which one to use
            print("There are more than one nextseq run for the date specified.")
            print("Please select the number of the run you want to process:")
            for i, f in enumerate(fnames):
                print("%d - %s"%(i, f))
            try:
                fnum = int(input("> "))
            except ValueError:
                logging.critical("A number should be selected")
                raise
        elif len(fnames) > 1:
            fnum = choose
        data_dir = fnames[fnum]
    cmd_exec = [command, '-R', data_dir, '-o', basedir,
                '--with-failed-reads', '--interop-dir', "%s/InterOp"%basedir]
    cmd_exec.extend(defpar.split())
    logging.info("Running: %s"%' '.join(cmd_exec))
    spl = SpooledTemporaryFile(dir=defs.TMP_DIR)
    call(cmd_exec, stdout = spl)
    spl.seek(0)
    for line in spl:
        logging.info(line.strip())
    spl.close()
    logging.info("bcl2fastq done")
            
        
def run_cutadapt(fname1, is_paired, adapter_seq, adapter_pair_seq, minlen,
                 qualval, cutadapt_cmd, cutadapt_add, dont_delete):
    """
    Run cutadapt on the input file, trimming low quality 3' ends and adapters
    Arguments:
    - `fname1`: Name of R1 file, assume ends with _1.fastq
    - `is_paired`: The run is paired, assume file name of R2 is same as R1
                   but with _2.fastq
    - `adapter_seq`: optional (set to None otherwise), remove this adapter
    - `adapter_pair_seq`: adapter of paired sequence
    - `minlen`: Minimal length of read to report
    - `qualval`: Quality minimal value
    - `cutadapt_cmd`: executable of cutadapt
    - `cutadapt_add`: Additional parameters for cutadapt
    """
#    os.environ["PYTHONPATH"] += ":/usr/local/icore-hm/x86_64.debian64_5775/python2.7/lib/python2.7/site-packages"
    cutadapt_cmd = [cutadapt_cmd, '-m', str(minlen),
                    '-q', str(qualval)]
    if is_paired:
        # First pass
        fname2 = fname1.rsplit("_", 1)[0] + "_2.fastq"
        tmp_outf1 = NamedTemporaryFile(dir=defs.TMP_DIR, suffix='.fastq')
        tmp_outf2 = NamedTemporaryFile(dir=defs.TMP_DIR, suffix='.fastq')
        add_params = ['-o', tmp_outf1.name, '-p', tmp_outf2.name]
        if adapter_seq:
            for ads in adapter_seq.split(','):
                add_params.extend(['-a', ads])
        add_params.extend(cutadapt_add.split())
        logging.info("First cutadapt pass over the pair %s"%fname1)
        logging.info("Running %s"%' '.join(
            cutadapt_cmd + add_params + [fname1, fname2]))
        spl = SpooledTemporaryFile(dir=defs.TMP_DIR)
        call(cutadapt_cmd + add_params + [fname1, fname2], stdout=spl)
        spl.seek(0)
        for line in spl:
            logging.info(line.strip())
        spl.close()
        outf1_name = fname1.rsplit("_", 1)[0] + "_cutadapt_1.fastq.gz"
        outf2_name = fname1.rsplit("_", 1)[0] + "_cutadapt_2.fastq.gz"
        add_params = ['-o', outf2_name, '-p', outf1_name]
        if adapter_pair_seq:
            add_params.extend(['-a', adapter_pair_seq])
        logging.info("Second cutadapt pass over the pair %s"%fname2)
        logging.info("Running: %s"%' '.join(
            cutadapt_cmd + add_params + [tmp_outf2.name, tmp_outf1.name]))
        spl = SpooledTemporaryFile(dir=defs.TMP_DIR)
        call(cutadapt_cmd + add_params + [tmp_outf2.name, tmp_outf1.name],
             stdout=spl)
        spl.seek(0)
        for line in spl:
            logging.info(line.strip())
        spl.close()
        # Remove the input files and close the temp files
        tmp_outf1.close()
        tmp_outf2.close()
        if not dent_delete:
            os.unlink(fname1)
            os.unlink(fname2)
    else:
        # One pass
        outf1_name = fname1.rsplit("_", 1)[0] + "_cutadapt_1.fastq.gz"
        add_params = ['-o', outf1_name]
        if adapter_seq:
            for ads in adapter_seq.split(','):
                add_params.extend(['-a', ads])
        logging.info("Running cutadapt on %s"%fname1)
        logging.info("Running: %s"%(
            cutadapt_cmd + add_params + [fname1]))
        spl = SpooledTemporaryFile(dir=defs.TMP_DIR)
        call(cutadapt_cmd + add_params + [fname1], stdout=spl)
        spl.seek(0)
        for line in spl:
            logging.info(line.strip())
        spl.close()
        if not dont_delete:
            os.unlink(fname1)
    logging.info("Finished cutadapt on file %s" %fname1)


def generate_summary(basedir, wwwdir, dbname):
    """
    Generate an html file that will be sent through email to the user.
    It will contain links to the fastqc summaries and wig file.
    Arguments:
    - `basedir`: Dir of this run
    - `wwwdir`: Where the files are available through internet
    Returns:
    - `content`: An email with all the details
    """
    p1, wdir = wwwdir.split('/www/', 1)
    username = p1.rsplit('/', 1)[1]
    wigname = glob.glob("%s/*.wig.gz"%wwwdir)[0].rsplit('/', 1)[1]
    content = """
    <p>Your nextseq auto process is done for library %s.
    FastQC results are here:<br>
    """%(basedir.rsplit('/', 1)[1])
    for hfile in glob.glob("%s/*.html"%wwwdir):
        content += "http://www.cs.huji.ac.il/~%s/%s/%s<br>\n"%(
            username, wdir, hfile.rsplit('/',1)[1])
    content += """
    You can find your bam and counts files in: %s/%s <br>
    You can upload this wiggle file to UCSC genome browser using this link:<br>
    http://microbes.ucsc.edu/cgi-bin/hgTracks?db=%s&hgt.customText=http://www.cs.huji.ac.il/~%s/%s/%s</p>
    """%(basedir, defs.BWA_DIR, dbname, username, wdir, wigname)
    return content
    



def main(argv=None):
    settings = process_command_line(argv)
    # Open a log file unique for this run
    run_id = "%s_%s_%s"%(
        settings.benchman, settings.experiment_name, settings.date)
    logging.basicConfig(
        filename="%s/%s_%s.log"%(
        defs.LOGS, run_id, time.strftime("%y-%m-%d-%H-%M")), filemode='w',
        level=logging.DEBUG, format='%(asctime)s %(levelname)s: %(message)s')
    logging.info(' '.join(sys.argv))
    wwwdir = "%s/%s"%(defs.WWWDIR, run_id)
    # Create the output directory and subdirectories
    basedir = make_outdir(settings, wwwdir)
    if settings.remap:
        settings.skip_bcl2fastq = True
        settings.skip_split = True
        settings.skip_cutadapt = True
    if not settings.skip_bcl2fastq:
        # Run bcl2fastq
        run_bcl2fatq(
            basedir, settings.date, force_dir=settings.force_dir,
            choose=settings.choose, csbase=defs.BASE,
            command=defs.BCL2FASTQ, defpar=defs.BCLPARAMS)

    if not settings.skip_split:
        # Get the Undetermined files
        undet_R1 = sorted(glob.glob("%s/%s*_R1*.fastq.gz"%(basedir,defs.UNDET)))
        undet_R2 = sorted(glob.glob("%s/%s*_R2*.fastq.gz"%(basedir,defs.UNDET)))
        is_paired = len(undet_R2) > 0
        # Run fastqc on all the raw files

        run_fastqc(
            undet_R1+undet_R2, "%s/%s"%(basedir, defs.QCRESULTS_DIR), wwwdir,
            settings.fastqc_exec)
        # Split by indices
        if settings.index_table:
            logging.info("splitting files: %s %s using index file %s"%(
                undet_R1, undet_R2, settings.index_table))
            index_splitter.split_reads(
                undet_R1, undet_R2, settings.index_table, mismatch=True,
                outdir="%s/%s"%(basedir, defs.FASTQ_DIR), unknown=False,
                strict=(not settings.keep_nobarcode))
            logging.info("Done index_splitter")
            # Copy the barcode table to the basedir
            try:
                shutil.copy(settings.index_table, basedir)
            except OSError:
                pass
        else:
            # Move all the files to fastq dir and concatenate them to one file
            concat_gzipped(
                undet_R1,
                "%s/%s/%s_1.fastq"%(basedir, defs.FASTQ_DIR, defs.SINGLE_NAME))
            if undet_R2:
                concat_gzipped(
                    undet_R2,
                    "%s/%s/%s_2.fastq"%(
                    basedir, defs.FASTQ_DIR, defs.SINGLE_NAME))
        # Remove the Undetermined files
        if not settings.dont_delete:
            for fname in undet_R1 + undet_R2:
                os.unlink(fname)
                logging.info("Removing raw file %s"%fname)
    else:
        is_paired = len(
            glob.glob("%s/%s/*_2.fastq"%(basedir, defs.FASTQ_DIR))) > 0
        
    if not settings.skip_cutadapt:
        # Get the list of files in the fastq directory
        input_files_1 = glob.glob("%s/%s/*_1.fastq"%(basedir, defs.FASTQ_DIR))
        # Remove adapter from 3' end and low quality bases
        for fname1 in input_files_1:
            run_cutadapt(
                fname1, is_paired, settings.adapter_seq,
                settings.adapter_pair_seq, settings.minimal_length,
                settings.quality_cutoff, settings.cutadapt_cmd,
                settings.cutadapt_params, settings.dont_delete)

    # After cutadapt the list should be different with cutadapt in the middle
    input_files_1 = glob.glob("%s/%s/*_1.fastq.gz"%(basedir, defs.FASTQ_DIR))
    # Run fastqc on all the fastq files after cutadapt
    input_files_2 = glob.glob("%s/%s/*_2.fastq.gz"%(basedir, defs.FASTQ_DIR))
    run_fastqc(
        input_files_1+input_files_2, "%s/%s"%(basedir, defs.QCRESULTS_DIR),
        wwwdir, settings.fastqc_exec)

    # Run bwa on all files
    if settings.force_single:
        is_paired = False
    if not settings.skip_mapping:
        for fname1 in input_files_1:
            if is_paired:
                fname2 = fname1.rsplit("_", 1)[0] + "_2.fastq.gz"
            else:
                fname2 = None
            map_and_count.run_bwa(
                settings.bwa_exec, fname1, fname2,
                "%s/%s"%(basedir, defs.BWA_DIR),
                fname1.rsplit('/', 1)[1].rsplit("_", 1)[0] + "_bwa",
                settings.allowed_mismatches, settings.genome_fasta,
                settings.bwa_aln_params, settings.bwa_sampe_params,
                settings.bwa_samse_params, settings.samtools_cmd)
        
    # count and prepare wig files
    features, feat_list = count_PE_fragments.read_gtf(
        open(settings.genes_gff), defs.FEAT_IDEF, defs.GENE_IDEF)
    feat_list = sorted(list(feat_list))
    bamfiles = glob.glob("%s/%s/*.bam"%(basedir, defs.BWA_DIR))
    all_counts = {}
    wigout = gzip.open("%s/%s/%s"%(basedir, defs.BWA_DIR, defs.WIGFILE), 'w')
    for bf in bamfiles:
        samfile = pysam.Samfile(bf)
        lib_name = bf.rsplit("/", 1)[1].rsplit(".", 1)[0]
        logging.info("Counting reads per gene for %s"%lib_name)
        all_counts[lib_name] = count_PE_fragments.count_features(
            features, samfile, settings.gene_overlap,
            rev=settings.reverse_complement)
        logging.info("Generate wiggle track for %s"%lib_name)
        coverage = sam_to_wiggle_coverage.generate_wig(
            samfile, rev=settings.reverse_complement,
            first_pos=settings.plot_first_position)
        sam_to_wiggle_coverage.print_wiggle(
            coverage, "%s_%s"%(lib_name, settings.date),
            "%s_%s"%(lib_name, settings.date),wigout)
        # Write to file
        with open(bf.rsplit(".", 1)[0] + ".counts", 'w') as counts_file:
            for g in feat_list:
                counts_file.write("%s\t%d\n"%(g, all_counts[lib_name][g]))
    wigout.close()
    # Make a link to the www dir so it can be uploaded to UCSC genome browser
    dest_www_file = "%s/%s.wig.gz"%(wwwdir, run_id)
    try:
        shutil.copy("%s/%s/%s"%(basedir, defs.BWA_DIR, defs.WIGFILE), "%s/"%wwwdir)
    except IOError, ose:
        logging.error("Can't copy file to www directory %s"%dest_www_file)
    # Write a table with all the counts in the experiment
    lib_names = sorted(all_counts.keys())
    with open("%s/%s/%s"%(basedir, defs.BWA_DIR, defs.COUNTSFILE), 'w') as acounts:
        writer = csv.writer(acounts, delimiter="\t")
        header = [''] + lib_names
        writer.writerow(header)
        for g in feat_list:
            writer.writerow([g] + [all_counts[l][g] for l in lib_names])
    if settings.email_address:
        body = generate_summary(basedir, wwwdir, settings.dbname)
        logging.info(
            "Sending email notification to %s"%(settings.email_address))
        try:
            sendmail.send_html(
                defs.SENDER, defs.PASSWORD, settings.email_address.split(','),
                defs.SUBJECT + " For %s"%run_id, body)
        except SMTPAuthenticationError, err:
            logging.error("Can't send the email: %s"%str(err))
        

    return 0        # success

if __name__ == '__main__':
    status = main()
    sys.exit(status)
