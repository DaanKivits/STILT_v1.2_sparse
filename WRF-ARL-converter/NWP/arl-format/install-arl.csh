#!/bin/csh -f

# Add additional compiler types as needed
# Issue from arl-format directory (NWP/arl-format), w/o args for help message

set help=yes
set supptypes = ( pgf90 xlf ifort ivyifort )
if ($#argv == 1) then
   set comptype=$1
   foreach supptype ( $supptypes )
    if ("a$comptype" == "a$supptype") set help=no
   end
endif

if ("a$help" == "ayes") then
   echo "$0 help message:"
   echo "Need to provide one argument: compiler type, one of:"
   echo " $supptypes"
   echo "If needed, add additional code for new machines"
   echo "Use environment variable STILTSRC to specify stilt/src pathname (default: ../p1229/stilt/src)"
   echo "    - only needed for interp_terrain.x, not needed for grib/ncdf_to_arl converters"
   echo ""
   exit 1
endif

# ifort, as on Cartesius Linux x86_64, Intel compiler

if ($?STILTSRC == 0) setenv STILTSRC "/home/awoude/STILT/STILT_Model/merged_stilt_hysplit"
if ("a$comptype" == "aifort") then
  set compname="ifort"
  set compcmd="ifort"
  set cflags0="-FR -convert big_endian -I${STILTSRC}"
  set cflags1="-FR -convert big_endian"
  set cflags2="-convert big_endian"
  set cflags3="-FR -O -convert big_endian -I."
  set confarg="F77=ifort"
endif
       
#preferentially use nf-config over nc-config, allow specification via NETCDFDIR environement variable:
set configs = ( nf-config nc-config )
set iconfig = 0
set notfound = 1
while ($notfound && $iconfig < $#configs)
 @ iconfig++
 set configcmd=$configs[$iconfig]
 if ($?NETCDFDIR) then
  if ("a$NETCDFDIR" != "aIGNORE") then
   set configcmd="${NETCDFDIR}/bin/${configcmd}"
  endif
 endif
 set configcmd=`which ${configcmd}`
 set notfound=$status
end
set libnet=""
set flgnet=""
if ($notfound) then
  echo "nc-config not found, 'which' returned: $configcmd"
  echo "Trying to use NETCDF environment variable instead, may not work for netcdf4"
  if ($?NETCDF == 0) then
   echo "Missing environment variable NETCDF, trying `which ncdump`"
   set nethome=`which ncdump | sed -e 's/\/bin\/ncdump//'`
   set testnh=`echo $nethome | cut -c 1,1`
   if ("a$testnh" == "a/") then
      setenv NETCDF $nethome
      echo "Found NETCDF=$NETCDF"
   else
      echo "NETCDF not found, linking will likely fail for wrfnc2arl"
   endif
  endif
  set libnet="-L${NETCDF}/lib  -lnetcdf"
  set flgnet="-I${NETCDF}/include"
else
  set libnet=`$configcmd --flibs`
  set flgnet=`$configcmd --fflags`
endif

echo "Using netcdf flags:"
echo " flgnet=$flgnet"
echo " libnet=$libnet"
set topmake="make -k MODULEDIR=$STILTSRC"
echo "Using topmake=$topmake"

if ("a$comptype" == "apgf90") then
# Linux machine with pgf90
 set compname="pgf90"
 set compcmd="pgf90"
 set cflags0="-O -Mfree -byteswapio -I. -I${STILTSRC}"
 set cflags1="-O -Mfree -byteswapio -I."
 set cflags2="-O -byteswapio"
 set cflags3="-O -Mfree -byteswapio -I./ -Mlfs"
 set confarg="F77=pgf90 CC=pgcc"
endif

# IBM (such as AER's humboldt)
if ("a$comptype" == "axlf") then
 set compname="xlf"
 set compcmd="xlf90_r"
 set cflags0="-qstrict -d -qarch=auto -w -qspill=20000  -qmaxmem=32767 -I. -I${STILTSRC}"
 set cflags1="-qstrict -d -qarch=auto -w -qspill=20000  -qmaxmem=32767 -I."
 set cflags2="-qstrict -d -qarch=auto -w -qspill=20000  -qmaxmem=32767 -I. -qfixed"
 set cflags3="-qstrict -d -qarch=auto -w -qspill=20000  -qmaxmem=32767 -I. -qcheck"
 set confarg="F77=xlf CC=xlc"
endif

# ifort, as on SGI machines at Pleiades or rotor
if ("a$comptype" == "aifort") then
 set compname="ifort"
 set compcmd="ifort"
 set cflags0="-free -w -I. -assume byterecl -I${STILTSRC}"
 set cflags1="-free -w -I. -assume byterecl"
 set cflags2="-w -I. -assume byterecl"
 set cflags3="-free -O -w -I. -assume byterecl"
 set confarg="F77=ifort"
endif

# Ivy bridge (Pleiades) SGI machines at Pleiades
if ("a$comptype" == "aivyifort") then
 set compname="ifort"
 set compcmd="ifort"
 set cflags0="-free -O -axAVX -xSSE4.1 -w -I. -assume byterecl -I${STILTSRC}"
 set cflags1="-free -O -axAVX -xSSE4.1 -w -I. -assume byterecl"
 set cflags2="-w -I. -assume byterecl"
 set cflags3="-free -O -axAVX -xSSE4.1 -w -I. -assume byterecl"
 set confarg="F77=ifort"
endif

$topmake FC="$compcmd" FFLAGS="$cflags0"
cd data2arl/libfcsubs/
./configure $confarg
make
make test
ln -s lib${compname}fcsubs.a libfcsubs.a
cd ../libhysplit/
make FC="$compcmd" CFLAGS="$cflags1"
cd ../libcmapf/
make FC="$compcmd" FFLAGS="$cflags2"
cd ..
mkdir library
cd library/
ln -s ../libfcsubs/libfcsubs.a .
ln -s ../libcmapf/libcmapf.a .
ln -s ../libhysplit/libhysplit.a .
cd ../wrfgrib2arl
make FC="$compcmd" CFLAGS="$cflags3 -I../libhysplit"
cd ../wrfnc2arl
make FC="$compcmd" FFLAGS="$cflags3 -I../libhysplit $flgnet" LIB="-L../library -lhysplit $libnet"
cd ../geos2arl
make FC="$compcmd" FFLAGS="$cflags3 -I../libhysplit $flgnet -I../wrfnc2arl" LIB="-L../library -lhysplit $libnet"
