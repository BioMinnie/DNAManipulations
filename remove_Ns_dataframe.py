import pandas as pd
import sys

# read in file as a 1xN dataframe (delimiter doesn't matter)
# change filename/path of input file below
df=pd.read_csv(sys.argv[1], header=None)


nIndex=set() # create a new set to hold indexes containing 'N'

# loop through each row (genome)
for index, row in df.iterrows():
    print("Inspecting genome " + str(index+1) + "...")
    # loop through each position in the genome
    for pos in range(0, len(row[0])):
        # record indexes containing 'N'
        if row[0][pos] == 'N':
            nIndex.add(pos)

print("Inspection complete")

# we now have a set of all unique indexes that contain an 'n'
# however python strings are immutable, so we can't just delete them
# instead let's loop through the df again, this time reconstructing genomes
# and store those in a new df, which we can then save

new_rows=[]
for index, row in df.iterrows():
    print("Constructing revised genome " + str(index+1) + "...")
    new_row=""
    for pos in range(0, len(row[0])):
        if pos not in nIndex:
            new_row=new_row+row[0][pos]
    new_rows.append(new_row)

print("Reconstruction complete")

new_df=pd.DataFrame(columns=df.columns)
new_df=new_df.append(new_rows)
#Change filename/path of output file below
new_df.to_csv(path_or_buf=sys.argv[2])

print("Reconstruction saved")