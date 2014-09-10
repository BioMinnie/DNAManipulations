## ORIGINAL AUTHOR: Benjamin Tovar 21/04/2012
## MODIFIED BY: Melinda Ashcroft 10/09/2014
## MELINDA AFFILIATIONS: Beatson lab | SCMB - UQ

# Import library
require('ape')

# Import R file containing list of functions
source(file.choose()) # Requires this file: "motifOccurrence.R" - Author Benjamin Tovar

# Current genbank acccession is for strain EC958
genome <- read.GenBank("HG941718.1")           # Read in genbank file from NCBI, given an accession number
genomeDNA <- as.character.DNAbin(genome)       # Convert data to a character vector
genomeDNA <- genomeDNA[[1]]                    # Converts character vector to string of bases
motifEcoli <- c("g", "a", "g", "a", "c", "c")  # Motif to search for: GAGACC - this can be changed to any non-fuzzy motif

# Run coordMotif function in motifOccurence.R on DNA to extract the coordinates of the motif
genomeDNAcoord <- coordMotif(genomeDNA,motifEcoli)
# run computeDistance function in motifOccurence.R on coordinates to extract the average distance between motifs
genomeDNAmotifDistance <- computeDistance(genomeDNAcoord)
genomeDNAmotifDistance # prints the average distance between motifs

# Compute the number of occurrences of the motif among the input DNA
occur <- (length(genomeDNAcoord)-1)
occur # prints the number of occurrences
