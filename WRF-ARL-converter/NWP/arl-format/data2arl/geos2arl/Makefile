# $RCSfile: Makefile,v $
#--------------------------------------------------------
# Version: $Revision: 1.3 $ Last Revised: $Date: 2016/04/12 14:54:44 $
#--------------------------------------------------------

SHELL = /bin/sh
RM = /bin/rm -f

LIB = -L../library -lhysplit -L$(NETCDF)/lib  -lnetcdf 

SUN5 = -O -free -I./ -I../libhysplit 
AIX8 = -O3 -qarch=auto -qmaxmem=8192 -I./ -I../libhysplit -qstrict
AIX6 = -O -bmaxdata:2000000000 -bmaxstack:256000000 -qarch=com -qmaxmem=16384 -I./ -I../libhysplit
AIX5 = -O -qarch=com  -qmaxmem=8192 -I./ -I../libhysplit  
AIX4 = -O -qarch=pwr2 -I./ -I../libhysplit 
AIX3 = -O -qcharlen=3072 -Wv -Wp,-gv 
DEC3 = -O -assume byterecl -I./ -I../libhysplit 
SGI5 = -O -bytereclen -I./ -I../libhysplit -freeform 
SGI3 = -O -old_rl -I./ -I../libhysplit
ABSF = -O1 -f free -I./ -I../libhysplit
PGF9 = -O -Mfree -byteswapio -I./ -I../libhysplit

FFLAGS = $(AIX5) -I$(NETCDF)/include
FC = xlf90

all: geos2arl

WRFNC_OBJECTS = ../wrfnc2arl/datcnv.o ../wrfnc2arl/datini.o ../wrfnc2arl/getdim.o ../wrfnc2arl/read_header.o ../wrfnc2arl/setvar.o
OBJECTS = get_fname.o make_index.o set_lev.o set_map.o process_onetime.o 

$(WRFNC_OBJECTS) :
	(cd ../wrfnc2arl; $(MAKE) FC="$(FC)" FFLAGS="$(FFLAGS)" $(subst ../wrfnc2arl/,,$(WRFNC_OBJECTS)))

geos2arl : geos2arl.f $(OBJECTS) $(WRFNC_OBJECTS)
	$(FC) $(FFLAGS) -o $@ $@.f $(OBJECTS) $(WRFNC_OBJECTS) $(LIB)

clean :
	$(RM) *.o *.mod
	$(RM) geos2arl

#geos2arl.o: $(WRFNC_OBJECTS) process_onetime.o 
get_fname.o:
make_index.o:
process_onetime.o: $(WRFNC_OBJECTS) get_fname.o make_index.o set_lev.o set_map.o
set_map.o: $(WRFNC_OBJECTS)
set_lev.o: $(WRFNC_OBJECTS)
