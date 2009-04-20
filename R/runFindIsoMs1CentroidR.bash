#!/bin/bash

# a wrapper script to run findIsoMs1Centroid.R each individual mzXML file to avoid overflowing memory

if [ $# -lt 1 ]; then
    echo Usage: $0 input_mzXML_file
    exit -1;
fi

for p in $@
  do
  echo $p
  cat ~/svnrepos/R/findIsoMs1Centroid.R | sed "s/bash\.input\.file/$p/g" > tmp.R
  R CMD BATCH --vanilla tmp.R
  mv tmp.Rout $p.Rout
done

exit 0
