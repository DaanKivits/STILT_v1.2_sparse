# $RCSfile: Makefile,v $
#--------------------------------------------------------
# Version: $Revision: 1.4 $ Last Revised: $Date: 2013/07/31 19:45:14 $
#--------------------------------------------------------

MODULEDIR = ../p1229/stilt/multi

SHELL = /bin/sh

SUN5 = -O -free -I./ -M$(MODULEDIR)
AIX8 = -O3 -qarch=auto -qmaxmem=8192 -I./ -qstrict -I$(MODULEDIR)
AIX6 = -O -bmaxdata:2000000000 -bmaxstack:256000000 -qarch=com -qmaxmem=16384 -I./ -I$(MODULEDIR)
AIX5 = -O -qarch=com  -qmaxmem=8192 -I./ -I$(MODULEDIR)
AIX4 = -O -qarch=pwr2 -I./ -I$(MODULEDIR)
AIX3 = -O -qcharlen=3072 -Wv -Wp,-gv -I$(MODULEDIR)
DEC3 = -O -assume byterecl -I./ -I$(MODULEDIR)
SGI5 = -O -bytereclen -I./ -freeform -I$(MODULEDIR)
SGI3 = -O -old_rl -I./ -I$(MODULEDIR)
ABSF = -O1 -f free -I./ -I$(MODULEDIR)
PGF9 = -O -Mfree -byteswapio -I.
HUM = -qsource -O2 -qmaxmem=-1 -q64 -C #AER-humboldt for xlf compilers; -C does bounds checking

GFOR = -O2 -fconvert=big-endian -frecord-marker=4 -fbounds-check -ffree-form -I.
IFOR = -FR -O -convert big_endian -I./ -I$(MODULEDIR)

FFLAGS = $(IFOR)
FC       = ifort

EXECS = chk_data.x interp_terrain.x interp_receptors.x

OBJECTS = unpack_subs.o interp_subs.o $(MODULEDIR)/constants_module.o $(MODULEDIR)/misc_definitions_module.o \
	$(MODULEDIR)/mprintf.o $(MODULEDIR)/module_map_utils.o

all : $(EXECS)

chk_data.x : chk_data.o unpack_subs.o
	$(FC) $(FFLAGS) -o $@ $^

interp_terrain.x : interp_terrain.o $(OBJECTS)
	$(FC) $(FFLAGS) -o $@ $^

interp_receptors.x : interp_receptors.o $(OBJECTS)
	$(FC) $(FFLAGS) -o $@ $^

clean :
	rm -f *.o *.mod
	rm -f $(EXECS)

interp_terrain.o: unpack_subs.o interp_subs.o 
interp_receptors.o: unpack_subs.o interp_subs.o 
chk_data.o: unpack_subs.o
