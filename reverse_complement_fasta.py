#######################################################################################################
#                                                                                                     #
## AUTHOR: Melinda M Ashcroft                                                                         #
## AFFILIATION: Beatson Lab | SCMB - University of Queensland                                         #
#                                                                                                     #
## PURPOSE: Reverse complement a fasta sequence							      #
## Date: 19 December 2014                                                                            #
#                                                                                                     #

#######################################################################################################

# Import modules
from Bio.Seq import Seq
from Bio.Alphabet import generic_nucleotide
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import sys

with open(sys.argv[1], 'r') as fin: # Open the infile (fasta)
    with open(sys.argv[2], 'w') as fout: # Create and open the outfile (fasta)
	# To open a file format other than fasta, change "fasta"
	# I.e. "gb" or "embl"
	# If your file is multifasta, this will reverse complement all sequences 
        for seq_record in SeqIO.parse(fin, "fasta"): 
            forward_seq = seq_record.seq 			# obtain sequence
            reverse_seq = forward_seq.reverse_complement()	# reverse complement
            reverse_rec = SeqRecord(reverse_seq, seq_record.id)	# create new record
            SeqIO.write(reverse_rec, fout, "fasta")		# write out reverse complemented sequence
