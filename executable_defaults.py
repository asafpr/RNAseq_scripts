# Default path of executables needed to run process_nextseq_run.py
BWA = "/usr/bin/bwa"
FASTQC = "/cs/icore/asafp5/software/FastQC/fastqc"
FASTQC_FILES = ("html", "zip")
BCL2FASTQ = "/usr/local/bcl2fastq/2.15.0.4/bin/bcl2fastq"
SAMTOOLS = "samtools"
BCLPARAMS = "-p 8 -d 6 -r 4 -w 4"
DATA = "/mnt/lustre/hms-01/fs01/asafp5/lab_nextseq_complete_data/data"#"/cs/icore/asafp5/lab_nextseq_complete_data/data"
LOGS = "/mnt/lustre/hms-01/fs01/asafp5/lab_nextseq_complete_data/logs"#"/cs/icore/asafp5/lab_nextseq_complete_data/logs"
BASE = "/cs/nextseq/"
UNDET = "" #"Undetermined"
FASTQ_DIR = "fastq_files"
QCRESULTS_DIR = "fastqc_results"
BWA_DIR = "bwa_mapping"
SINGLE_NAME = "all_run_reads"
#CUADAPT_PATH = "/usr/local/icore-hm/x86_64.debian64_5775/python2.7/lib/python2.7/site-packages/"
CUTADAPT_CMD = "/cs/icore/asafp5/lib/python/cutadapt"
MINLEN = 21
QUAL_CUTOFF = 15
MISMATCHES = 2
PARAMS_ALN = '-t 32 -k 1 -R 200 -l 20'
PARAMS_SAMPE = '-a 1500 -P'
PARAMS_SAMSE = ' '
OVERLAP = 10
FEAT_IDEF = "exon"
GENE_IDEF = "gene_id"
WIGFILE = "all_tracks.wig.gz"
COUNTSFILE = "all_counts_table.txt"
WWWDIR = "/cs/icore/asafp5/www/lab_nextseq_files"
SENDER = "nextseqautorun@gmail.com"
PASSWORD = "nextseqisthebest"
SUBJECT = "Hurray! your nextseq data is ready!"
ADAPTER_1 = "AGATCGGAAGAGC"
ADAPTER_2 = "AGATCGGAAGAGC"
LUSTER_DIR = "/mnt/lustre/hms-01/fs01/asafp5/lab_nextseq_complete_data/data"
TMP_DIR = "/mnt/lustre/hms-01/fs01/asafp5/lab_nextseq_complete_data/tmp"