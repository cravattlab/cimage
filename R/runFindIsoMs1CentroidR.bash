#!/bin/bash

# a wrapper script to run findIsoMs1Centroid.R each individual mzXML file to avoid overflowing memory

if [ $# -lt 1 ]; then
    echo Usage: $0 input_mzXML_file [fast]
    exit -1;
fi

R --vanilla --args $@ < /home/chuwang/svnrepos/R/findIsoMs1Centroid.R > $1.Rout

exit 0
