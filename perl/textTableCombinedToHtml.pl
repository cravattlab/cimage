#!/usr/bin/perl -w

die "Usage: $0 text_table\n" if (@ARGV <1 );

$intable=$ARGV[0];
##$title=$ARGV[1];

open(INFILE,$intable) || die "cannot open $intable: $!\n";
@txt = <INFILE>;
close(INFILE);
$nrow = @txt;
@header = split(/\t/,$txt[0]);
$ncol = @header;


$outhtml=$intable;
$outhtml =~ s/\.txt$/\.html/g;

open(OUTFILE, ">$outhtml") || die "cannot write to $outhtml: $!\n";
##print OUTFILE "Content-type: text/html\n\n";
print OUTFILE <<ENDOFHEADER;
<HTML><HEAD>
<STYLE TYPE="text/css">
    * {
        font-family: arial, helvetica, myriad;
        font-size: 13px;
    }
</STYLE>
ENDOFHEADER
print OUTFILE "<TABLE border=2 frame=\"border\" rules=\"groups\" summary=$intable>\n";
##print OUTFILE "<CAPTION> <b>$title</b> </CAPTION>\n";
print OUTFILE "<COLGROUP span=1>\n";
print OUTFILE "<COLGROUP span=4>\n";
print OUTFILE "<COLGROUP span=1>\n";
print OUTFILE "<COLGROUP span=3>\n";
print OUTFILE "<COLGROUP span=3>\n";
print OUTFILE "<COLGROUP span=3>\n";
print OUTFILE "<COLGROUP span=1>\n";
print OUTFILE "<COLGROUP span=1>\n";
print OUTFILE "<TBODY>\n";
print OUTFILE "<TR>\n";
for ($i=0; $i<$ncol; ++$i) {
    print OUTFILE "<TH align=\"center\" bgcolor=\"HoneyDew\">$header[$i]\n";
}
print OUTFILE "<TBODY>\n";
@bgcolormap=@header;
for ($i=0; $i<$ncol; ++$i) {
    $bgcolormap[$i]=""
}
for ($i=6; $i<9; ++$i) {
    $bgcolormap[$i]="bgcolor=\"#DCDCDC\""
}
for ($i=1; $i<$nrow; ++$i) {
    print OUTFILE "<TR>";
    @line = split('\t',$txt[$i]);
    if ($line[0] =~ /\d+/ ) {
	$bold1="<b>";
	$bold2="</b>";
	$anchor1 = 1;
    } else {
	$bold1="";
	$bold2="";
	$anchor1 = 0;
    }
    for ($j=0; $j<$ncol; ++$j) {
	$_ = $line[$j];
	if ( /^=HYPERLINK/ ) {
	    /^=HYPERLINK\("(.*)","(\d+.\d+)"\)$/;
	    print OUTFILE "<TD><A HREF=\"$1\">$2</A>";
	} elsif ( /^([a-zA-Z]+)$/ && $anchor1 ) {
	    print OUTFILE "<TD align=\"center\" $bgcolormap[$j]> <A NAME=\"$1\"></A> $bold1 $1 $bold2";
	} else {
	    print OUTFILE "<TD align=\"center\" $bgcolormap[$j]> $bold1 $line[$j] $bold2";
	}
    }
    print OUTFILE "\n";
}
print OUTFILE '</TABLE>'."\n";
