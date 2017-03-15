"""
Read and write fastq files faster than SeqIO of Biopython
"""
class Read(object):
    """
    represents one read
    """
    
    def __init__(self, name, seq, qual):
        """
        Initialize a read with a name, sequence string and quality string
        Arguments:
        - `name`: The read name (line 0)
        - `seq`: The sequence (line 1)
        - `qual`: Quality line (line 3)
        """
        self.name = name
        self.seq = seq
        self.qual = qual

    def __getitem__(self, key):
        """
        Return a portion of the read
        Arguments:
        - `key`: The indices
        """
        return Read(self.name, self.seq[key], self.qual[key])

    def write(self, outf):
        """
        Print the read to a fastq file
        Arguments:
        - `outf`: An open file handle for writing
        """
        outf.write("%s\n%s\n+\n%s\n"%(self.name, self.seq, self.qual))


class fastIO(object):
    """
    Read a fastq file
    """
    
    def __init__(self, fin):
        """
        Use an open file to generate sequences
        Arguments:
        - `fin`: An open file handle
        """
        self._fin = fin

    def __iter__(self, ):
        """
        """
        return self

    def next(self, ):
        """
        Return the next sequence
        """
        name = self._fin.next().rstrip()
        seq = self._fin.next().rstrip()
        _ = self._fin.next()
        qual = self._fin.next().rstrip()
        return Read(name, seq, qual)


                   
def seqwrite(rlist, outf):
    """
    Write the list of reads to a file
    Arguments:
    - `outf`: An open file
    - `rlist`: A list of reads
    """
    for r in rlist:
        r.write(outf)

