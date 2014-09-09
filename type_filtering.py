## AUTHOR: Melinda Ashcroft 09-09-14
## AFFILIATIONS: Beatson lab | SCMB - UQ

## Background: This script takes in a tab delimited text file, extracts the
## columns, searches column 3 for two strings and only writes (to a csv file)
## column 1 if those strings are not found.

# Import modules
import csv
import sys

# Create empty lists
file_name = []
typing = []

# Open both input files and output file
with open(sys.argv[1], 'r') as fin:
    with open(sys.argv[2], 'w') as fout:
        
        # Specify the csv writer features and create a header row
        writer = csv.writer(fout, delimiter='\t', quotechar='', quoting=csv.QUOTE_NONE)
        headers = ('file_names',)
        writer.writerow(headers)
        
        # Specify action for first input file
        for line in fin:
            line = line.strip()
            columns = line.split("\t")
            if columns[2] not in ('131', '-'):
                writer.writerow([columns[0]])