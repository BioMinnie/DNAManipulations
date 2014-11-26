## AUTHOR: Melinda Ashcroft 09-09-14
## AFFILIATIONS: Beatson lab | SCMB - UQ

## Used to extract the motif, start and end coordinates
## and distance between motifs from the first input file
## when iterating through a second input file of coordinates

# Import modules
import csv
import sys

# Create empty lists
motif = []
dist = []
mge_start = []
mge_end = []
motif_start = []
motif_end = []
new_motif_start = []
new_motif_end = []
mge = []
new_mge = []
merged = []

# Open both input files
# file1 is the list of motifs, coordinates and distances in whole genome
# file2 is the list of coordinates (ie. MGEs)

with open(sys.argv[1], 'r') as fin1, \
open(sys.argv[2], 'r') as fin2:
        
    # Specify action for first input file
    motifs_file = csv.reader(fin1, delimiter='\t', quotechar='', quoting=csv.QUOTE_NONE)
    motifs_headers = fin1.readline() # Reads and ignores the header line 
    for line in motifs_file:         # Iterates through the list of motifs
        motif.append(line[0])        # Motif is in the first column (python counts from zero)
        motif_start.append(line[1])  # Coordinate is in the second column
        motif_end.append(line[2])    # Coordinate is in the third column
        dist.append(str(line[3]))    # Distance between motifs is in the fourth column
   
    # Specify action for second input file
    mobile = csv.reader(fin2, delimiter='\t', quotechar='', quoting=csv.QUOTE_NONE)
    for row in mobile:                  # Iterates through the list of MGEs
        mge_start.append(row[1])        # Stores start position of MGE
        mge_end.append(row[2])          # Stores end position of MGE
        mge.append(row[4])              # Stores name of MGE
   
    # Extract coordinates and store in lists
    for i in range(len(mge_start)):     # Iterate through smallest list of coordinates first
        for j in range(len(motif_start)):
            # If start of motif is within mge and end of motif is within mge:
            if (int(motif_start[j]) >= int(mge_start[i]) and int(motif_end[j]) <= int(mge_end[i])):
                new_motif_start.append(str(motif_start[j]))      # Stores start position of motif
                new_motif_end.append(str(int(motif_start[j])+3)) # Stores end position of motif
                new_mge.append(mge[i])                           # Stores mge associated with each motif
   
    # Merge motif, new_motif_start, new_motif_end and new_mge lists
    for row in zip(motif, new_motif_start, new_motif_end, new_mge):
        merged.append('\t'.join(row))
    
# Create a function to write output
def write_to_csv(fout, items_to_write):
    csvwrite = open(fout, 'w')  # Open outfile
    # Must include an escapechar if using csv.QUOTE_NONE
    writer = csv.writer(csvwrite, escapechar='', quotechar="'", quoting=csv.QUOTE_NONE)
    for item in items_to_write: # Iterate through list of items to write        
        writer.writerow([item]) # Write each item of list or array

#Write output
write_to_csv(sys.argv[3], merged)
