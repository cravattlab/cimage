#!/bin/bash

if [ $# -lt 1 ]; then
    echo Usage: $0 [by_protein] [exclude_singleton] [descending]  text_file dir1 ...
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

descending=""
if [ $1 == "descending" ]; then
    descending=$1
    shift
fi

# filter by my_list.txt
mylist=$(find ./ -name "my_list.txt" 2> /dev/null | wc -l)
if [ "$mylist" -eq 1 ];
then
    echo User provides a customized list in my_list.txt -- filtering...
    listpath=$(find $PWD -name "my_list.txt")
else
    listpath="mylist_none"
fi

tmptxt=$(echo $1 | cut -c1-9)
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
##if [ "$nchar" -gt 20 ]; then
##    outname="combine_all"
##fi

if [ "$by_protein" == "by_protein" ]; then
    echo "R --vanilla --args $exclude_singleton $descending $listpath $txt $dirs < $CIMAGE_PATH/R/combined_by_protein.R > $outname.by_protein.Rout"
    R --vanilla --args $exclude_singleton $descending $listpath $txt $dirs < $CIMAGE_PATH/R/combined_by_protein.R > $outname.by_protein.Rout
else
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

