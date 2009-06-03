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
	cat $p2.tagged | grep $p1 | sed 's/IPI://g' | awk '{print $2}' > $p2.tmp.scan
	cat $p2.tagged | grep $p1 | sed 's/IPI://g' | awk '{print $NF}' > $p2.tmp.peptide
	cat $p2.tmp.peptide | sed 's/\*//g' | cut -f2 -d \.  > $p2.tmp.sequence
	for p3 in $(cat $p2.tmp.sequence)
	do
	    /home/chuwang/svnrepos/python/peptideCalcMass.py $p3 mono;
	done > $p2.tmp.mass
	cat $p2.tmp.scan | awk -F "." '{print $1}' | awk -F "_" '{print $NF}' > $p2.tmp.segment
	cat $p2.tmp.scan | awk -F "." '{print $NF}' > $p2.tmp.charge
	paste -d"_" $p2.tmp.peptide $p2.tmp.ipi $p2.tmp.segment $p2.tmp.charge | sed 's/\*/z/g' | sed 's/\-/z/g'> $p2.tmp.key
	paste -d " " $p2.tmp.ipi $p2.tmp.scan $p2.tmp.peptide $p2.tmp.mass $p2.tmp.key >> tmp.key_scan
	rm -rf $p2.tmp.*
	cat $p2 | grep Gene_Symbol | sed 's/IPI://g' | awk '{print $1} '| awk -F"|" '{print $1}'> tmp.ipi
	cat $p2 | grep Gene_Symbol | sed 's/IPI://g' | cut -f3 -d"=" | cut -c1-50 | sed -e s/^\-/_/g > tmp.name
	paste tmp.ipi tmp.name >> tmp.ipi_name
    done
done
echo ipi.name > ipi_name.table
cat tmp.ipi_name | tr -d '\r' | sort | uniq | sed -e s/\'//g | sed -e s/\"//g | sed -e s/\;/\ /g >> ipi_name.table

scanfiles=""
keys=$(cat tmp.key_scan | awk '{print $NF}' | sort | uniq )
for run in $mzxml
do
    echo $run > $run.tmp.scan
    for key in $keys
    do
	match=$(cat tmp.key_scan | grep $key | grep "$run\_" | wc -l)
	#echo $key $match
	if [ "$match" != "0" ]; then
	    cat tmp.key_scan | grep $key | grep "$run\_" | head -1 | awk '{print $2}' >> $run.tmp.scan
	else
	    echo none >> $run.tmp.scan
	fi
    done

    scanfiles="$scanfiles $run.tmp.scan"
done

echo "peptide mass ipi" > tmp.seq_mass
for key in $keys
do
    cat tmp.key_scan | grep $key | head -1 | awk '{print $3, $4, $1 }' >> tmp.seq_mass
done

paste -d " " tmp.seq_mass $scanfiles > cross_scan.table

rm -rf tmp.ipi tmp.name tmp.ipi_name tmp.seq_mass tmp.key_scan *.tmp.scan

echo running xcms to extract chromatographic peaks
R --vanilla --args ipi_name.table cross_scan.table  < /home/chuwang/svnrepos/R/findMs1AcrossSetsFromDTASelect.R > $infile.findMs1AcrossSetsFromDTASelect.Rout

echo generating graphs
mkdir -p ipi_name.table_output/PNG
cd ipi_name.table_output
for p in $(\ls *.ps | sed 's/\.ps//g')
do
    convert $p.ps $p.png
    mv $p*.png ./PNG/
    ps2pdf $p.ps
done
cd ..

mkdir -p ipi_name.table_output/TEXT
cp ipi_name.table ipi_name.table_output/TEXT
cp cross_scan.table ipi_name.table_output/TEXT