#!/bin/bash

if [ $# -lt 2 ]; then
    echo Usage: $0 [by_protein] text_file dir1 ...
    exit -1
fi

by_protein=""
if [ $1 == "by_protein" ]; then
    by_protein=$1
    shift
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

if [ "$by_protein" == "by_protein" ]; then
    R --vanilla --args $txt $dirs < /home/chuwang/svnrepos/R/combined_by_protein.R > $outname.by_protein.Rout
else
    R --vanilla --args $txt $dirs < /home/chuwang/svnrepos/R/combined.R > $outname.Rout
fi

mv combined.txt $outname.txt
mv combined.png $outname.png
mv combined_vennDiagram.png $outname.vennDiagram.png

/home/chuwang/svnrepos/perl/textTableCombinedToHtml.pl $outname.txt $outname $allpt

