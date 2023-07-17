#--------------------------------------------------------
# Last Revised: 17 Nov 2003
#--------------------------------------------------------

SHELL = /bin/sh

LIB = ../library 

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

CFLAGS = $(AIX5)
FC = xlf90

eta40arl : eta40arl.f 
	$(FC) $(CFLAGS) -o $@ $@.f -L$(LIB) -lfcsubs -lhysplit -lcmapf 

eta12arl : eta12arl.f 
	$(FC) $(CFLAGS) -o $@ $@.f -L$(LIB) -lfcsubs -lhysplit -lcmapf 

eta04arl : eta04arl.f 
	$(FC) $(CFLAGS) -o $@ $@.f -L$(LIB) -lfcsubs -lhysplit -lcmapf 

clean :
	rm -f *.o
	rm -f eta40arl
	rm -f eta12arl
	rm -f eta04arl
