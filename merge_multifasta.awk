awk '
    /^>/ { 
        # print the first header
        if (c++ == 0) {print; print ""} 
        next
    } 
    /^$/ {next} 
    {printf "%s", $0} 
    END {print ""}
' /Users/BioMinnie/Desktop/working/kSNP3/Laptop/input/E_coli_ST101/NDM26.contigs.fasta > /Users/BioMinnie/Desktop/working/NDM_ST101/seq/NDM26_merged.fasta