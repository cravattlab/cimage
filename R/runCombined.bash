#!/bin/bash

if [ $# -lt 1 ]; then
    echo Usage: $0 [by_protein] [exclude_singleton] text_file dir1 ...
    exit -1
fi

by_protein=""
if [ $1 == "by_protein" ]; then
    by_protein=$1
    shift
fi
exclude_singleton=""
if [ $1 == "exclude_singleton" ]; then
    exclude_singleton=$1
    shift
fi

tmptxt=$(echo $1 | cut -c1-7)
if [ $tmptxt == "output_rt" ]; then
    txt=$1
    shift
else
    txt="output_rt_10_sn_2.5.to_excel.txt"
fi
allpt=""
outname="combined"
for p in $@
do
    outname=$outname\_$p
    pt=$(find $p -name $txt)
    allpt="$allpt $pt"
done
dirs=$(echo $allpt | sed "s/$txt//g")

nchar=$(echo $outname | wc -c)
if [ "$nchar" -gt 20 ]; then
    outname="combine_all"
fi

if [ "$by_protein" == "by_protein" ]; then
    R --vanilla --args $exclude_singleton $txt $dirs < $CIMAGE_PATH/R/combined_by_protein.R > $outname.by_protein.Rout
else
#    echo  R --vanilla --args $txt $dirs
    R --vanilla --args $txt $dirs < $CIMAGE_PATH/R/combined.R > $outname.Rout
fi

mv combined.txt $outname.txt
mv combined.png $outname.png
if [ -s combined_vennDiagram.png ]; then
    mv combined_vennDiagram.png $outname.vennDiagram.png
fi

cwd=$(pwd)

if [ "$by_protein" == "by_protein" ]; then
    $CIMAGE_PATH/perl/textTableCombinedToHtml_by_protein.pl $outname $cwd $allpt
else
    $CIMAGE_PATH/perl/textTableCombinedToHtml.pl $outname $cwd $allpt
fi

