from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import GC
import sys

fasta_sequences = SeqIO.parse(open(sys.argv[1),'fasta')
for fasta in fasta_sequences:
    name = fasta.description
    print name
    print GC(fasta.seq)
    print len(fasta.seq)