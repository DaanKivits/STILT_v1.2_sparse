#------------------------------------------------------------------
# Makefile for the ../data2arl/arw2arl data converter 
# Last Revised: 28 Jan 2007 
#               09 Aug 2007 - support for gfortran
#               16 May 2013 - use netcdff for version 4
#               25 Aug 2014 - netCDF libaries previously exported
#------------------------------------------------------------------

SHELL = /bin/sh

PGM=arw2arl
EXE=../../exec
LIB=libarw.a

# exported in calling script
# NETLIB = /usr/local/lib
# NETINC = /usr/local/include

# library name for netCDF3=netcdf and for netCDF4=netcdff
  LINKS = -L$(NETLIB) -lnetcdf  -L../../library -lhysplit
# LINKS = -L$(NETLIB) -lnetcdff -L../../library -lhysplit

# standard compilation flags
SUN5 = -O -free
AIX8 = -O3 -qarch=auto -qmaxmem=8192 -qstrict
DEC3 = -O -assume byterecl
SGI5 = -O -bytereclen -freeform
ABSF = -O1 -f free 
PGF9 = -O -Mfree -byteswapio
INTL = -O -FR -assume byterecl -convert big_endian -I.
GF95 = -O2 -fendian=big -ffree-form 
GFOR = -O2 -fconvert=big-endian -frecord-marker=4 -ffree-form 

LFLAGS = -vr

CFLAGS = $(GFOR)
FC = gfortran
#CFLAGS = $(INTL)
#FC = ifort 
#CFLAGS = $(PGF9)
#FC = pgf90 

CODE = cfgrec.f \
       cfgvar.f \
       get1dv.f \
       get2dv.f \
       get3dv.f \
       get4dv.f \
       getdim.f \
       setmap.f \
       setvar.f 

all : $(PGM) 

$(LIB) : $(CODE)
	for member in $?; do make member MODULE=$$member;done

member : $(MODULE)
	$(FC) -c $(CFLAGS) -I$(NETINC) $?                                 
	$(AR) $(LFLAGS) $(LIB) $(?F:.f=.o)
	ranlib $(LIB)
	rm -f $(?F:.f=.o)

$(PGM) : $(PGM).f $(LIB) 
	$(FC) -o $(EXE)/$@ $(CFLAGS) -I$(NETINC) $? $(LINKS)  
	rm -f $(LIB)

clean :
	rm -f $(EXE)/$(PGM)  
	rm -f $(LIB)
