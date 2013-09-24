#!/bin/bash

if [ $# -lt 1 ]; then
    echo Usage: $0 data_dir
    exit -1
fi


cd $1

hst=$(hostname)

wd=$(pwd)
echo current working dir: $wd

#nhtml=$(\ls ./*.html | wc -l);
#if [ "$nhtml" -eq 0 ]; then
#    echo no html results found! Exit...
#    exit -1;
#fi

rhost="www.scripps.edu"
rdir="/web/docs/cravatt/cimage/tempul/"

out="temp"
#echo compress txt, html and png files into $out.zip
#find ./ -name "*.txt" -o -name "*.html" -o -name "*.png" | zip $out -@
echo $hst:$wd > $out.txt
if ! grep public_html $out.txt; then
    echo ERROR: Your data is not accessible from public_html folder! Exit ...
    exit -1
fi

find -L ./ -name "combine*.html" |sort >> $out.txt
chmod a+rw $out.txt

echo transfer $out.txt to server $rhost:$rdir\(use your TSRI password\)
rsync -avzr $out.txt $rhost:$rdir

echo transfer complete! Add it to CIMAGE database at URL below!
echo http://www.scripps.edu/cravatt/cimage/

exit 0
