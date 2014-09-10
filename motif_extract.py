## Melinda Ashcroft 01-09-14
## Beatson lab
## Special thanks to Areej for her help

# Import modules
from Bio import SeqIO
from Bio.Seq import Seq
import re
import csv
import sys

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
example_motif = str(sys.argv[3])

# Create empty lists
distance = []
sequence = []
seqstart = []
seqend = []
merged = []

with open(sys.argv[1], 'r') as fin1: # Open infile
        
    # Make the mregex a for loop using the example_motif user input and the iupacdict
    mregex = ''
    for n in example_motif:
        mregex += iupacdict.get(n)


    tosearch = re.compile(str(mregex), re.I)  # Compile the final regex, make sure to include re.Ignorecase!
    for record in SeqIO.parse(fin1, 'fasta'): # Parse the input fasta file
        cur = str(record.seq)                 # Convert the sequence to a string
        for match in tosearch.finditer(cur):  # Iterate through the string with the final regex
            if match:                         # If a match is found:
                sequence.append(str(match.group())) # Save the matched sequence to the sequence variable
                seqstart.append(str(match.start())) # Save the matched start position to the seqstart variable
                seqend.append(str(match.end()))     # Save the matched end position to the seqend variable
            else:                             # If no match is found, raise an exception
                raise Exception('name not found')

    # Manually add the distance between the final motif and the first motif (bacterial genome is circular)
    distance.append(str(652))

    # Enumerate through seqend, append distance list
    for i, item in enumerate(seqend):
        if (i < len(seqend)-1):
            distance.append(str(int(seqend[i]) - int(seqstart[i+1])))

    # Merge sequence, seqstart and seqend and distance lists
    for row in zip(sequence, seqstart, seqend, distance):
        merged.append('\t'.join(row))
       
# Create a function to write output
def write_to_csv(output_file, items_to_write):
    csvwrite = open(output_file, 'w')
    writer = csv.writer(csvwrite, delimiter="\t", quotechar=' ', quoting=csv.QUOTE_ALL) # Doesn't work if quoting is set to anything else
    headers = ('motif', 'start', 'end', 'distance')
    writer.writerow(headers)
    for item in items_to_write: # Iterate through list of items to write        
        writer.writerow([item]) # Write each item of each list        

write_to_csv(sys.argv[2], merged) # Write to output_file  
