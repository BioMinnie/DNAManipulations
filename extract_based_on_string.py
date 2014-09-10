## AUTHOR: Melinda Ashcroft 10-09-14
## AFFILIATIONS: Beatson lab | SCMB - UQ

# Import modules
import csv
import sys

# Open both input file and output file
with open(sys.argv[1], 'r') as fin1:
    with open(sys.argv[2], 'w') as fout:
        
        # Specify the csv writer features and create a header row
        writer = csv.writer(fout, escapechar='', quotechar="'", quoting=csv.QUOTE_NONE)
        
        # Specify action for first input file
        for line in fin1.readlines():
            line = line.strip()         # split each line on the linebreak
            columns = line.split("\t")  # split each line on the tab delimiter
            # This line is interchangeable: Use < for less than 1000bp or >= for greater than or equal to 1000bp
            if int(columns[3]) >= 1000: 
                writer.writerow([line])
