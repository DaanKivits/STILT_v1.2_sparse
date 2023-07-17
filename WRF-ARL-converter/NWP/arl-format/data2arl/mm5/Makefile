#---------------------------------------------------------
# Master Makefile for the ../data2arl/mm5 directory
# Last Revised: 08 Sep 2003
#---------------------------------------------------------

SHELL = /bin/sh

EXE = mm5toarl
LIB = ../library

SUN5 = -O  -c -I. -free
AIX8 = -O3 -c -I. -qarch=auto -qmaxmem=8192 -qstrict
AIX6 = -O  -c -I. -qarch=com  -qmaxmem=8192
AIX5 = -O  -c -I. -qarch=com  -qmaxmem=8192
AIX4 = -O  -c -I. -qarch=pwr2 -qmaxmem=8192
AIX3 = -O  -c -I. -qcharlen=3072 -Pv -Wp,-gv
DEC3 = -O  -c -I. -assume byterecl
SGI5 = -O  -c -I. -bytereclen -freeform
SGI3 = -O  -c -I. -old_rl
ABSF = -O1 -c -I. -f free
PGF9 = -O  -c -I. -Mfree -byteswapio

CFLAGS = $(AIX5)
LFLAGS = -vr
FC = xlf90  

.f.o:
	rm -f $@
	$(FC) $(CFLAGS) $*.f

MODULES = module_consts.o    \
	  module_file.o

FOBJS = mm5toarl.o    \
	convert.o     \
	crs2dot.o     \
	rdmm5v3.o     \
	setupv3.o     \
	setgrid.o     \
	setvarb.o     \
	wrtarl.o 

$(EXE) : $(MODULES) $(FOBJS)
	$(FC) -o $@ $(FOBJS) $(MODULES) -L$(LIB) -lhysplit -lcmapf

clean:
	rm -f $(EXE)
	rm -f *.o
	rm -f *.mod
