#!/bin/csh -f
# File id: $Id: split_chk_data.csh,v 1.8 2013/09/27 15:12:01 trn Exp $

# Purpose: see help message below

set append=0
switch ($#argv)
case 2:
   if ($2 == "append") set append=1
   if ($2 == "overwrite") set append=0
case 1:
   set fname=$1
   breaksw
default: 
   echo "$0 help message:"
   echo "This script splits output from chk_data by variable, to allow ingest into"
   echo "Splus or IDL for further analysis"
   echo "Output will be written to cwd, in files of the form chk_data_xxxx.txt"
   echo "NOTE: If these files exist, they will be overwritten unless arg 2 specifies otherwise"
   echo "This script needs 1-2 argument:"
   echo "1 - chk_data output file name"
   echo "2 (optional) - either 'append' or 'overwrite'"
   exit 1
endsw

set skip_line=`grep -n INDX $fname | head -1 | sed -e 's/:.*//'`
set xorno=`tail -n +$skip_line $fname | grep -v INDX | cut -c11,11 | head -1`
if ("a$xorno" == "aX") then
   set vars=( `tail -n +$skip_line $fname | grep -v INDX | cut -c17-20 | sort | uniq` )
else
   set vars=( `tail -n +$skip_line $fname | grep -v INDX | cut -c15-18 | sort | uniq` )
endif

foreach var ( $vars )
  set outfile=chk_data_${var}.txt
  if ($append) then
    set outstr=", appending to"
  else
    if (-e $outfile) /bin/rm -f $outfile    
    set outstr=", creating"
  endif
  echo "Processing var=$var $outstr outfile=$outfile"
  if ("a$xorno" == "aX") then
     tail -n +$skip_line $fname | grep "${var}" | grep ": min" | cut -c1-14,26-28,44-88 >>! $outfile
  else
     tail -n +$skip_line $fname | grep "${var}" | grep ": min" | cut -c1-12,24-26,42-86 >>! $outfile
  endif
end
