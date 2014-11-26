## AUTHOR: Melinda Ashcroft 01-09-14
## AFFILIATIONS: Beatson lab | SCMB - UQ
## Used to alter coordinates in a csv file as python is zero based
## Used to convert string of distances to integers

# Import libraries
import csv
import sys

offset1 = 1 # As python is zero based - update coordinates

with open(sys.argv[1],'r') as fin1: # Open infile
    with open(sys.argv[2], 'w') as fout: # Open outfile
        writer = csv.writer(fout, delimiter='\t', quotechar='', quoting=csv.QUOTE_NONE)
        # Change headers to whatever headings you require for each column
        headers = ('motif', 'updated_start', 'end', 'updated_distance')
        writer.writerow(headers)        # Write headers
        file_headers = fin1.readline()  # Reads and ignores the header line 
        lines = fin1.readlines()        # Iterate through the lines
        for i in range(0, len(lines)):  # For each line in infile
            motif, start, end, distance = lines[i].strip().split("\t")
            updated_distance = abs(int(distance)) # Convert string to absolute integer
            updated_start = int(start)+offset1    # Convert string to integer and add offset
            write = motif, updated_start, end, updated_distance # Save variables to write to a single variable
            writer.writerow(write)                # Write the write variable to the output file
