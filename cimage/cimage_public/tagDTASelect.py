#!/usr/bin/env python
#
# tag each peptide line in DTASelect output file with its IPI name at beginning

import sys
from sys import argv

if len(argv) != 2:
    print 'Usage: %s <DTASelect-filter.txt>'%argv[0]
    print 'tag each peptide line in DTASelect output file with its IPI name at beginning'
    sys.exit(-1)

## list to save ipi and peptide lines in case multiple protein with the same peptides
ipi_lines = []
pep_lines = []
tag=''
last_line_is_peptide = False
## control header and tail printing
print_on = True

for line in open(argv[1]):
    line = line.rstrip()
    words = line.split()
    if len(words) <= 3:
        print
        continue
    # none peptide entry line
    if words[3].find('%') != -1 or words[0]=='Proteins':
        if words[0]=='Proteins':
            print_on = True
        else:
            print_on = False

        if last_line_is_peptide:
            # print out saved lines
            for ipi in ipi_lines:
                print 'cimageipi-'+ipi
                # find tag
                ipi_words = ipi.split()
                i = ipi_words[0].find("|")
                if i>0:
                    tag = ipi_words[0][0:i]
                else:
                    tag = ipi_words[0]
                # print out tagged lines
                for pep in pep_lines:
                    print 'cimagepep-'+tag, pep[1:]
            # emtpy lists
            ipi_lines = []
            pep_lines = []

        if not print_on:
            ipi_lines.append(line)
            last_line_is_peptide = False
    else:
        if not print_on:
            pep_lines.append(line)
            last_line_is_peptide = True

    if print_on:
        print line