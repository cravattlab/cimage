#!/bin/bash

if [ $# -lt 1 ]; then
    echo Usage: $0 set_1 set_2 ...
    exit -1
fi

mzxml=$@

echo -n > tmp.ipi_name
echo -n > tmp.key_scan
echo -n > tmp.seq_mass
# tagged and fwd
echo Parsing DTASelect files
for p1 in $mzxml
do
    echo $p1
    for p2 in $(ls DTASelect-filter_$p1\_*.txt);
    do
	/home/chuwang/svnrepos/python/tagDTASelect.py $p2 > $p2.tagged;
	#cat $p2.tagged | grep "^IPI" | sed 's/\ \*//g' > $p2.tagged.fwd;
	cat $p2.tagged | grep $p1 | sed 's/IPI://g' | awk '{print $1}' > $p2.tmp.ipi
	cat $p2.tagged | grep $p1 | sed 's/IPI://g' | awk '{print $2}' > $p2.tmp.FileName
	cat $p2.tagged | grep $p1 | sed 's/IPI://g' | awk '{print $3}' > $p2.tmp.xcorr
	cat $p2.tagged | grep $p1 | sed 's/IPI://g' | awk '{print $NF}' > $p2.tmp.peptide
	cat $p2.tmp.peptide | sed 's/\*//g' | cut -f2 -d \.  > $p2.tmp.sequence
	for p3 in $(cat $p2.tmp.sequence | sort | uniq )
	do
	    echo -n "$p3 "
	    /home/chuwang/svnrepos/python/peptideCalcMass.py $p3 mono;
	done > $p2.tmp.uniq.mass
	for p3 in $(cat $p2.tmp.sequence)
	do
	    grep "^$p3 " $p2.tmp.uniq.mass | awk '{print $2}'
	done > $p2.tmp.mass
	##cat $p2.tagged | grep $p1 | sed 's/IPI://g' | awk '{print $6}' > $p2.tmp.mass
	cat $p2.tmp.FileName | awk -F "." '{print $1}' | awk -F "_" '{print $NF}' > $p2.tmp.segment
	cat $p2.tmp.FileName | awk -F "." '{print $2}'  > $p2.tmp.scan
	cat $p2.tmp.FileName | awk -F "." '{print $NF}' > $p2.tmp.charge
	cat $p2.tmp.FileName | awk -v run=$p1 '{print run}' > $p2.tmp.run
	paste -d":" $p2.tmp.ipi $p2.tmp.peptide $p2.tmp.charge $p2.tmp.segment > $p2.tmp.key
	paste -d " " $p2.tmp.run $p2.tmp.scan $p2.tmp.mass $p2.tmp.xcorr $p2.tmp.key >> tmp.key_scan
	rm -rf $p2.tmp.*
	cat $p2 | grep Gene_Symbol | sed 's/IPI://g' | awk '{print $1} '| awk -F"|" '{print $1}'> tmp.ipi
	cat $p2 | grep Gene_Symbol | sed 's/IPI://g' | cut -f3 -d"=" | cut -c1-50 | sed -e s/^\-/_/g > tmp.name
	paste tmp.ipi tmp.name >> tmp.ipi_name
    done
done

echo Creating input files for xcms
## create ipi number to protein name map
echo "name" > ipi_name.table
cat tmp.ipi_name | tr -d '\r' | sort | uniq | sed -e s/\'//g | sed -e s/\"//g | sed -e s/\;/\ /g >> ipi_name.table

## table with scans from different dataset
scanfiles=""
keys=$(cat tmp.key_scan | awk '{print $NF}' | sort | uniq )
for p1 in $mzxml
do
    echo $p1 > $p1.tmp.scan
    for key in $keys
    do
	match=$(cat tmp.key_scan | grep -F $key | grep "^$p1 " | wc -l)
	#echo $key $match
	if [ "$match" != "0" ]; then
	    cat tmp.key_scan | grep -F $key | grep "^$p1 " | sort -k 4 -rn | head -1 | awk '{print $2}' >> $p1.tmp.scan
	else
	    echo "0" >> $p1.tmp.scan
	fi
    done

    scanfiles="$scanfiles $p1.tmp.scan"
done

echo "key mass" > tmp.seq_mass
for key in $keys
do
    cat tmp.key_scan | grep -F $key | head -1 | awk '{print $NF, $3}' >> tmp.seq_mass
done
paste -d " " tmp.seq_mass $scanfiles > cross_scan.table

## table with all ms2 scans
echo "key run scan" > all_scan.table
cat tmp.key_scan | awk '{print $NF, $1, $2}' >> all_scan.table

rm -rf tmp.ipi tmp.name tmp.ipi_name tmp.key_scan tmp.seq_mass  *.tmp.scan

echo Running xcms to extract chromatographic peaks
R --vanilla --args $mzxml < /home/chuwang/svnrepos/R/findMs1AcrossSetsFromDTASelect.R > findMs1AcrossSetsFromDTASelect.Rout

echo Generating graphs
mkdir -p output/PNG
cd output
for p in $(\ls *.ps | sed 's/\.ps//g')
do
    convert -rotate 90 $p.ps $p.png
    mv $p*.png ./PNG/
    ps2pdf $p.ps
done
cd ..

mkdir -p output/TEXT
cp ipi_name.table output/TEXT
cp cross_scan.table output/TEXT
cp all_scan.table output/TEXT

echo Finished
exit 0