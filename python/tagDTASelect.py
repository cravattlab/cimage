#!/usr/bin/env python
#
# tag each peptide line in DTASelect output file with its IPI name at beginning

import sys
from sys import argv

if len(argv) != 2:
    print 'Usage: %s <DTASelect-filter.txt>'%argv[0]
    print 'tag each peptide line in DTASelect output file with its IPI name at beginning'
    sys.exit(-1)

# get current IPI tag
tag=''
for line in open(argv[1]):
    line = line.rstrip()
    words = line.split()
    # none peptide entry line
    if len(words) <= 5:
        tag=''
        print line
        continue
    if words[0].find('IPI') != -1:
        i = words[0].find("|")
        if i>0:
            tag = words[0][0:i]
        else:
            tag = words[0]
        print line
    else:
        print tag, line



