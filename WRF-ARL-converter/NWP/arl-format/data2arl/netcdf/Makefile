#------------------------------------------------------------------
# Master Makefile for the ../data2arl/netcdf directory
# Last Revised: 08 Sep 2003
#------------------------------------------------------------------
SHELL = /bin/sh

NETCDF = /sw/arch/RedHatEnterpriseServer7/EB_production/2019/software/netCDF/4.6.1-intel-2018b
#NETCDF = /hpc/sw/netcdf-fortran-4.4.2-intel-seq/bin
#LIB = -L../library -lhysplit -L$(NETCDF)/lib64  -lnetcdf -lnetcdff
#LIB = -L../library -L$(NETCDF)/lib64  -lnetcdf -lnetcdff
LIB = -I../library -L$(NETCDF)/lib64  -lnetcdf -lnetcdff

SUN5 = -O -free -I./
AIX8 = -O3 -qarch=auto -qmaxmem=8192 -I./ -qstrict
AIX6 = -O -bmaxdata:2000000000 -bmaxstack:256000000 -qarch=com -qmaxmem=16384 -I./
AIX5 = -O -qarch=com  -qmaxmem=8192 -I./
AIX4 = -O -qarch=pwr2 -I./
AIX3 = -O -qcharlen=3072 -Wv -Wp,-gv
DEC3 = -O -assume byterecl -I./
SGI5 = -O -bytereclen -I./ -freeform
SGI3 = -O -old_rl -I./
ABSF = -O1 -f free -I./
PGF9 = -O -Mfree -byteswapio -I./
IFOR = -FR -O -convert big_endian -I./ -I$(MODULEDIR)

CFLAGS = $(IFOR)/include -I../libhysplit
FC = ifort

CODE = ncscan.f \
       nc2dat.f \
       ncjulh.f \
       gaussl.f \
       intp2d.f \
       setxy.f

all : cdf2arl

$(LIB) : $(CODE)
	for member in $?; do make member MODULE=$$member;done

member : $(MODULE)
	$(FC) -c $(CFLAGS) $(FPATH)$?                                 
	$(AR) $(LFLAGS) $(LIB) $(?F:.f=.o)
	rm -f $(?F:.f=.o)

cdf2arl : cdf2arl.f $(LIB) 
#cdf2arl : cdf2arl.f
	$(FC) -o $@ $(CFLAGS) $? -L./ -lnetcdf -L$(DIR) -lfcsubs -lhysplit -lcmapf -L$(LIB)
#	$(FC) -o $@ $(CFLAGS) $? -L./ -lnetcdf -L$(DIR) -lfcsubs -L$(LIB)

clean :
	rm -f *.o
	rm -f cdf2arl
#	rm -f $(LIB)
