#!/bin/bash

if [ $# -lt 1 ]; then
    echo Usage: $0 data_dir
    exit -1
fi


cd $1

wd=$(pwd)
echo current working dir: $wd

nhtml=$(\ls ./combine*.html | wc -l);
if [ "$nhtml" -eq 0 ]; then
    echo no combined html results found! Exit...
    exit -1;
fi

rhost="137.131.5.161"
ruser="samba"
rdir="/srv/www/htdocs/cimage/tempul/"

out="temp"
echo compress txt, html and png files into $out.zip
find ./ -name "*.txt" -o -name "*.html" -o -name "*.png" | zip $out -@

chmod a+rw $out.zip

echo transfer $out.zip to server \(use password \"cravatt\"\) and this may take some time
rsync -ar --progress --remove-source-files $out.zip $ruser@$rhost:$rdir

echo transfer complete! Add it to CIMAGE database at URL below!
echo http://bfclabcomp3.scripps.edu/cimage/

exit 0
