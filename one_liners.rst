Bash/Awk one liners - useful for data manipulations
===================================================

Remove every other line:

awk 'NR % 2 == 0' file > newfile

Remove every 10th line:

awk 'NR % 10 != 0' file > newfile
