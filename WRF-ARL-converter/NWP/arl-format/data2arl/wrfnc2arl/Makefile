# $RCSfile: Makefile,v $
#--------------------------------------------------------
# Version: $Revision: 1.2 $ Last Revised: $Date: 2016/04/12 14:54:44 $
#--------------------------------------------------------

SHELL = /usr/bin/sh
#NETCDF = /usr/local/netcdf-3.6.3-gfortran
NETCDF = /sw/arch/RedHatEnterpriseServer7/EB_production/2019/software/netCDF/4.6.1-intel-2018b
#NETCDF = /hpc/sw/netcdf-fortran-4.4.2-intel-seq/bin
LIB = -L../library -lhysplit -L$(NETCDF)/lib64  -lnetcdf -lnetcdff

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

GFOR = -O2 -fconvert=big-endian -frecord-marker=4 -fbounds-check -ffree-form -I.
IFOR = -FR -O -convert big_endian -I./ -I$(MODULEDIR)


#FFLAGS = $(GFOR) -I$(NETCDF)/include
FFLAGS = $(IFOR)/include -I../libhysplit
#FC = xlf90
FC = ifort

OBJECTS = copy_arr.o datcnv.o datini.o getdim.o gettime.o makndx.o read_header.o setmap.o setvar.o \
 hybrid_vertical_coords.o xtract.o

wrfnc2arl : wrfnc2arl.f $(OBJECTS)
	$(FC) $(FFLAGS) -o $@ $@.f $(OBJECTS) $(LIB)

clean :
	rm -f *.o *.mod
	rm -f wrfnc2arl

makndx.o: setmap.o 
wrfnc2arl.o: read_header.o xtract.o 
xtract.o: getdim.o gettime.o makndx.o setvar.o 
gettime.o: getdim.o setvar.o 
makndx.o: setmap.o 
wrfnc2arl.o: read_header.o xtract.o
hybrid_vertical_coords.o: setvar.o
xtract.o: copy_arr.o datcnv.o datini.o getdim.o gettime.o makndx.o setvar.o hybrid_vertical_coords.o
