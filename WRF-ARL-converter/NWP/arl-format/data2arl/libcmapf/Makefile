SHELL = /bin/sh

F_OBJECTS = stlmbr stcm1p eqvlat stcm2p \
	cc2gxy cc2gll cg2cxy cg2cll \
	cgszll cgszxy ccrvxy ccrvll \
	cpolxy cpolll cll2xy cxy2ll \
	cnllxy cnxyll cspanf

LIB = libcmapf.a

library: $(LIB)

$(LIB): $(F_OBJECTS:=.o) 
	$(AR) $(ARFLAGS) $(LIB) $?

test: cmaptst
	cmaptst >fortest.txt

cmaptst: $(LIB) cmaptst.f
	$(FC) -O -o cmaptst cmaptst.f $(LIB) 

clean:
	rm  *.o
	rm *.a
	rm cmaptst
	rm fortest.txt
