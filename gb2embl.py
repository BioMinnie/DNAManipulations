import glob
from Bio import SeqIO
file = glob.glob('gb/*.gb')
#print file
for f in file:
    outf = '.'.join(f.split('.')[0:-1])+'.embl'
    SeqIO.convert(f, "genbank", outf, "embl")
    #print outf
print 'finish converting gb to embl'