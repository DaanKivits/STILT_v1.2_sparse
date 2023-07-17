#!/bin/csh -f

#  cvs ID: $Id: arl-chk.csh,v 1.6 2008/01/15 19:13:25 trn Exp $

if ($?TMPTRN) then
 set tmpdir=$TMPTRN
else
 set tmpdir=/home/awoude/STILT/WRF-ARL-converter/NWP/arl-format/arl/tmp
endif
if ($?NWPDIR) then
 set NWPDIR=$NWPDIR
else
 set NWPDIR=/home/awoude/STILT/WRF-ARL-converter/NWP/
endif
if ($?SXTEEN) then
 set SXTEEN=$SXTEEN
else
 set SXTEEN=F
endif

switch ($#argv)
case 0:
	echo "$0 aborting because supplied $#argv arguments"
	echo "It needs 1 or more:"
	echo "1 (or more) - ARL input files"
	echo "Also optional are environment variables:"
	echo "TMPTRN (default used otherwise: $tmpdir)"
	echo "NWPDIR (default used otherwise: $NWPDIR)"
	echo "SXTEEN (T/F for 16-bit; default used otherwise: $SXTEEN)"
	exit (1)
default:
	set infiles = ( $* )
endsw

mkdir $tmpdir
mkdir $tmpdir/chk-arl
cd $tmpdir/chk-arl

@ count=0
foreach file ( $infiles )

 echo "$file" >! arl-chk.in
 echo "$SXTEEN" >>  arl-chk.in
 echo "4" >>  arl-chk.in
 echo "-1 1 1 -1 1 1" >>  arl-chk.in
 $NWPDIR/arl-format/chk_data.x < arl-chk.in > ! arl-chk.txt

 set arg2="append"
 if ($count == 0) set arg2="overwrite"
 $NWPDIR/arl-format/split_chk_data.csh arl-chk.txt $arg2
 @ count++

end

echo "source('$NWPDIR/arl-format/read.split.chkdata.r')" >! arl-plot.in
echo "tmp.list <- read.all.split.chkdata('./')" >>! arl-plot.in
echo "plot.all.split.chkdata(tmp.list)" >>! arl-plot.in
set Rcmd=`which R`
echo $Rcmd | grep 'not found'
set found=$status
if ($found == 0) then
   if ($?R_HOME) then
      set Rcmd="${R_HOME}/bin/R"
   else
      echo "R not found in path, will not be run"
      set Rcmd=echo
   endif
endif
if ("a$Rcmd" != "aecho") set Rcmd="$Rcmd --no-save"
echo "Rcmd=$Rcmd"
$Rcmd < arl-plot.in >&! arl-plot.out
