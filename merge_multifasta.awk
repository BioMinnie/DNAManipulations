awk '
    /^>/ { 
        # print the first header
        if (c++ == 0) {print; print ""} 
        next
    } 
    /^$/ {next} 
    {printf "%s", $0} 
    END {print ""}
' input_file.fasta > output_file.fasta
