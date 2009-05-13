#!/bin/bash

if [ $# -lt 3 ]; then
    echo Usage: $0 input_spreadsheet.txt set_1 set_2 ...
    exit -1
fi

infile=$(echo $1 | sed 's/\.txt$//g')
shift
mzxml=$@

cat $infile.txt | grep -v "Description" | awk 'BEGIN{FS="\t"}{print $3}' > $infile.peptide_tagged
cat $infile.txt | grep -v "Description" | awk 'BEGIN{FS="\t"}{print $3}' | sed 's/\*//g' | cut -f2 -d \.  > $infile.peptide

# tagged and fwd
for p1 in $mzxml
do
    for p2 in $(ls DTASelect-filter_$p1\_*.txt);
    do
	/home/chuwang/svnrepos/python/tagDTASelect.py $p2 > $p2.tagged;
	cat $p2.tagged | grep "^IPI:IPI" | sed 's/\ \*//g' > $p2.tagged.fwd;
    done
done

scanfiles=""
# pull out scan number
for f in $mzxml
do
    for p in $(cat $infile.peptide_tagged )
    do
	pp=$(cat DTASelect-filter_$f\_*.fwd | awk -v frag=$p '(frag==$NF){print $1, $2, ";"}' );
	echo $p $pp;
    done  > DTASelect-filter_$f.match
    echo $f > $infile\_$f.scan
    cat DTASelect-filter_$f.match  | awk '{if (NF==1) print "none";else print $3}' >> $infile\_$f.scan
    scanfiles="$scanfiles $infile"_$f.scan
done

# mass
echo "mass" > $infile.peptide_mass
for p in $(cat $infile.peptide);
do
    /home/chuwang/svnrepos/python/peptideCalcMass.py $p mono;
done >> $infile.peptide_mass

echo "peptide" > $infile.peptide_seq
cat $infile.peptide >> $infile.peptide_seq

paste $infile.peptide_seq $infile.peptide_mass $scanfiles > $infile.txt.seq_mass_scan

R --vanilla --args $infile.txt $infile.txt.seq_mass_scan < /home/chuwang/svnrepos/R/findMs1AcrossSets.R > findMs1AcrossSets.Rout

cd output
ps2pdf *.ps
cd ..
