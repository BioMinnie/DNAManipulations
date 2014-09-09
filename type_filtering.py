## AUTHOR: Melinda Ashcroft 09-09-14
## AFFILIATIONS: Beatson lab | SCMB - UQ

## Background: This script takes in a tab delimited text file, extracts the
## columns, searches column 3 for two strings and only writes (to a csv file)
## column 1 if those strings are not found.

# Import modules
import csv
import sys

# Open both input file and output file
with open(sys.argv[1], 'r') as fin:
    with open(sys.argv[2], 'w') as fout:
        
        # Specify the csv writer features and create a header row
        writer = csv.writer(fout, delimiter='\t', quotechar='', quoting=csv.QUOTE_NONE)
        headers = ('file_names',)
        writer.writerow(headers) # write headers to outfile
        
        # Specify action for the input file
        for line in fin:
            line = line.strip()        # split each line
            columns = line.split("\t") # split each tab and save to columns variable
            if columns[2] not in ('131', '-'): # if column 3 (python zero based) does not contain these strings
                writer.writerow([columns[0]])  # write output to outfile
