import glob
from Bio import SeqIO
file = glob.glob('*.gb')
#print file
for f in file:
    outf = '.'.join(f.split('.')[0:-1])+'.fa'
    SeqIO.convert(f, "genbank", outf, "fasta")
    #print outf
print 'finish converting gb to fasta'

	