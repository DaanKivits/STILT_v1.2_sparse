#------------------------------------------------------------------
# Makefile for the grib2arl data converter 
#------------------------------------------------------------------

SHELL = /bin/sh

PGM=grib2arl

SUN5 = -O -free
AIXC = -O3 
AIX8 = -O3 -qarch=auto -qmaxmem=8192 -qstrict
AIX6 = -O -qarch=com -qmaxmem=8192
AIX5 = -O -qarch=com -qmaxmem=8192
AIX4 = -O -qarch=pwr2 -qmaxmem=8192
AIX3 = -O -qcharlen=3072 -Pv -Wp,-gv
DEC3 = -O -assume byterecl
SGI5 = -O -bytereclen -freeform
SGI3 = -O -old_rl
ABSF = -O1 -f free 
PGF9 = -O -Mfree -byteswapio
INTL = -O -FR -assume byterecl -convert big_endian -I.
GF95 = -O2 -fendian=big -ffree-form 
GFOR = -O2 -fconvert=big-endian -frecord-marker=4 -ffree-form 
CFOR = -fconvert=big-endian -frecord-marker=4 -ffree-form -Wall -fbounds-check -fbacktrace -finit-local-zero

CFLAGS = $(GFOR)
LFLAGS = -vr
FC = gfortran

$(PGM) : $(PGM).f 
	$(FC) -o $@ $(CFLAGS) $? -L../library -lhysplit 

clean :
	rm -f $(PGM)  
