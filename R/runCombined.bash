#!/bin/bash

if [ $# -lt 2 ]; then
    echo Usage: $0 text_file dir1 ...
    exit -1
fi

txt=$1
shift
allpt=""
outname="combined"
for p in $@
do
    outname=$outname\_$p
    pt=$(find $p -name $txt)
    allpt="$allpt $pt"
done
dirs=$(echo $allpt | sed "s/$txt//g")

R --vanilla --args $txt $dirs < /home/chuwang/svnrepos/R/combined.R > combined.Rout

mv combined.txt $outname.txt

/home/chuwang/svnrepos/perl/textTableCombinedToHtml.pl $outname.txt

