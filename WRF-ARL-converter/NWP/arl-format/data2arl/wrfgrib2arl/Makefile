# $RCSfile: Makefile,v $
#--------------------------------------------------------
# Version: $Revision: 1.2 $ Last Revised: $Date: 2016/04/12 14:54:44 $
#--------------------------------------------------------

SHELL = /bin/sh

LIB = ../library 

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

CFLAGS = $(AIX5)
FC = xlf90

wrfgrib2arl : wrfgrib2arl.f delta_time.o
	$(FC) $(CFLAGS) -o $@ $@.f delta_time.o -L$(LIB) -lfcsubs -lhysplit -lcmapf 

delta_time.o : delta_time.f
	$(FC) -c $(CFLAGS) $^

clean :
	rm -f *.o
	rm -f wrfgrib2arl
