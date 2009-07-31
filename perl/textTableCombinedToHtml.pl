#!/usr/bin/perl -w

die "Usage: $0 text_table [ratio.png] [run_dirs] \n" if (@ARGV <1 );

$intable=$ARGV[0];
##$title=$ARGV[1];

open(INFILE,$intable) || die "cannot open $intable: $!\n";
@txt = <INFILE>;
close(INFILE);
$nrow = @txt;
@header = split(/\t/,$txt[0]);
$ncol = @header;
@bgcolormap=@header;
$nset=0;
for ($i=0; $i<$ncol; ++$i) {
    if ($header[$i] =~ /^mr\./) {
	$bgcolormap[$i]="bgcolor=\"#DCDCDC\"";
	$nset++;
    } else {
	$bgcolormap[$i]=""
    }
}


$outhtml=$intable;
$outhtml =~ s/\.txt$/\.html/g;

open(OUTFILE, ">$outhtml") || die "cannot write to $outhtml: $!\n";
##print OUTFILE "Content-type: text/html\n\n";
print OUTFILE <<ENDOFHEADER;
<HTML><HEAD>
<STYLE TYPE="text/css">
*{
    font-family: arial, helvetica, myriad;
    font-size: 13px;
}
table#sample TD {
    text-align: center;
}
</STYLE>
</HEAD>
<BODY>
ENDOFHEADER

if(@ARGV>1) {
    print OUTFILE "<A HREF=\"$ARGV[1]\">ratio plot</A><BR><BR>\n";
}
if(@ARGV>2) {
    for ( $i=2; $i<@ARGV; $i++ ) {
	$j = $i-1;
	print OUTFILE "run $j: <A HREF=\"$ARGV[$i]\">$ARGV[$i]</A><BR><BR>";
  }
}
print OUTFILE "<TABLE id=\"sample\" border=2 frame=\"border\" rules=\"groups\" summary=$intable>\n";
##print OUTFILE "<CAPTION> <b>$title</b> </CAPTION>\n";
print OUTFILE "<COLGROUP span=1>\n";
print OUTFILE "<COLGROUP span=4>\n";
print OUTFILE "<COLGROUP span=1>\n";
print OUTFILE "<COLGROUP span=$nset>\n";
print OUTFILE "<COLGROUP span=$nset>\n";
print OUTFILE "<COLGROUP span=3>\n";
print OUTFILE "<COLGROUP span=1>\n";
print OUTFILE "<COLGROUP span=1>\n";
print OUTFILE "<TBODY>\n";
print OUTFILE "<TR>\n";
for ($i=0; $i<$ncol; ++$i) {
    #print OUTFILE "<TH align=\"center\" bgcolor=\"HoneyDew\">$header[$i]\n";
    print OUTFILE "<TH bgcolor=\"HoneyDew\">$header[$i]\n";
}
print OUTFILE "<TBODY>\n";
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
	    print OUTFILE "<TD $bgcolormap[$j]> <A NAME=\"$1\"></A> $bold1 $1 $bold2";
	} else {
	    print OUTFILE "<TD $bgcolormap[$j]> $bold1 $line[$j] $bold2";
	}
    }
    print OUTFILE "\n";
}
print OUTFILE '</TABLE>'."\n";
print OUTFILE '</BODY>'."\n";
print OUTFILE '</HTML>'."\n";
