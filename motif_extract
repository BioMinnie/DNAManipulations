## AUTHOR: Melinda Ashcroft 01-09-14
## AFFILIATIONS: Beatson lab | SCMB - UQ
## NOTE: Thank you Areej for your help

## Background: This script takes in a fasta file, searches for a motif, converting
## taking into account IUPAC code for ambiguous nucleotides, finds the motif, start
## and end positions and the distance between motifs

## Warning: This script is still a work in progress
####################################################################################

# Import libraries
from Bio import SeqIO
from Bio.Seq import Seq
import re
import csv
import sys
import numpy as np

# Specify dictionary of IUPAC code for nucleotides
iupacdict = {'A':'A',
'C':'C',
'G':'G',
'T':'T',
'M':'[AC]',
'R':'[AG]',
'W':'[AT]',
'S':'[CG]',
'Y':'[CT]',
'K':'[GT]',
'V':'[ACG]',
'H':'[ACT]',
'D':'[AGT]',
'B':'[CGT]',
'X':'[ACGT]',
'N':'[ACGT]'}

# Motif to use for search
motif = str(sys.argv[4])

# Create empty lists
distance = []
sequence = []
seqstart = []
seqend = []
merged = []

with open(sys.argv[1], 'r') as fin: # open infile
     
    # Make the mregex a for loop using the motif user input and the iupacdict
    mregex = ''
    for n in motif:
        mregex += iupacdict.get(n)


    tosearch = re.compile(str(mregex), re.I)  # Compile the final regex, make sure to include re.Ignorecase!
    for record in SeqIO.parse(fin, 'fasta'):  # Parse the input fasta file
        cur = str(record.seq)                 # Convert the sequence to a string
        for match in tosearch.finditer(cur):  # Iterate through the string with the final regex
            if match:                         # If a match is found
                sequence.append(str(match.group())) # Save the matched motif to the sequence variable
                seqstart.append(str(match.start())) # Save the matched start position to the seqstart variable
                seqend.append(str(match.end()))     # Save the matched end position to the seqend variable
            else:                             # If no match is found raise an exception
                raise Exception('motif not found')
    
    # Add the first value zero to the distance list as there is no distance until 
    # we find the distance between the first and second motifs
    distance.append(0)
    # Enumerate through seqend, append distance list
    for i, item in enumerate(seqend):
        if (i < len(seqend)-1): # Tells enumerate where to stop iterating
            distance.append(int(seqend[i]) - int(seqstart[i+1])) # Saves the distance between motifs
            # The distance section is a work in progress. Currently saves in a non-standard format and as a negative value
            # Due to non-standard format, cannot merge with sequence, seqstart and seqend lists
            # Thus requiring writing to a second file and merging the two files using awk

# Merge sequence, seqstart and seqend lists
for row in zip(sequence, seqstart, seqend):
    merged.append(' '.join(row))
 
# Create a function to write output
def write_to_csv(fout, items_to_write):
    csvwrite = open(fout, 'wb')
    writer = csv.writer(csvwrite, quotechar=' ', quoting=csv.QUOTE_ALL) # Doesn't work if quoting is set to anything else
    for item in items_to_write: # Iterate through list of items to write        
        writer.writerow([item]) # Write each item of each list        

# Write output to files
write_to_csv(sys.argv[2], merged)     # first file
write_to_csv(sys.argv[3], distance)   # second file

"""
Have to merge sys.argv[2] and sys.argc[3] files using awk
# Awk command below:
paste <(awk '{print $1, $2, $3}' file2) <(awk '{print $1}' file3) > merged_file.txt
"""
