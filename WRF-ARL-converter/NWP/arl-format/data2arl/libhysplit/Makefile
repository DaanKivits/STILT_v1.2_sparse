#-----------------------------------------------------------------
# Master Makefile
# Last Revised: 17 Nov 2003
#-----------------------------------------------------------------

SHELL = /bin/sh

DIR=.
LIB=$(DIR)/libhysplit.a

SUN5 = -O -free -I. 
AIX8 = -O3 -qarch=auto -qmaxmem=8192 -I. -qstrict
AIX6 = -O -bmaxdata:2000000000 -bmaxstack:256000000 -qarch=com -qmaxmem=8192 -I.
AIX5 = -O -qarch=com -qmaxmem=8192 -I.
AIX4 = -O -qarch=pwr2 -qmaxmem=8192 -I.  
AIX3 = -O -qcharlen=3072 -Pv -Wp,-gv 
DEC3 = -O -assume byterecl -I. 
SGI5 = -O -bytereclen -I. -freeform 
SGI3 = -O -old_rl -I.
ABSF = -O1 -f free -I. 
PGF9 = -O -Mfree -byteswapio -I.

GFOR = -O2 -fconvert=big-endian -frecord-marker=4 -fbounds-check -ffree-form -I.
IFOR = -FR -O -convert big_endian -I./ -I$(MODULEDIR)

CFLAGS = $(GFOR)
LFLAGS = -vr           
FC = gfortran

CODE =          \
       pakinp.f \
       pakout.f \
       pakout_16.f \
       pakset.f \
       pakset_16.f \
       pakrec.f \
       pakrec_16.f \
       pakini.f \
       pakndx.f \
       pakndx_16.f \
       tm2day.f \
       tm2jul.f \
       tm2min.f \
       tminit.f \
       tmplus.f \
       w3code.f 

$(LIB) : $(CODE)
	for member in $?; do make member MODULE=$$member;done

member : $(MODULE)
	$(FC) -c $(CFLAGS) $?                                 
	$(AR) $(LFLAGS) $(LIB) $(?F:.f=.o)
	rm -f $(?F:.f=.o)

clean :
	rm -f *.o
	rm -f $(LIB)
