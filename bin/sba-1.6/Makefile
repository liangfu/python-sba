#
# Makefile for Sparse Bundle Adjustment library & demo program
#
CC=gcc
#ARCHFLAGS=-march=pentium4 # YOU MIGHT WANT TO UNCOMMENT THIS FOR P4
ARCHFLAGS=-isysroot /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX.sdk -arch x86_64
CFLAGS=$(ARCHFLAGS) -O3 -funroll-loops -Wall #-Wstrict-aliasing #-g -pg
OBJS=sba_levmar.o sba_levmar_wrap.o sba_lapack.o sba_crsm.o sba_chkjac.o
SRCS=sba_levmar.c sba_levmar_wrap.c sba_lapack.c sba_crsm.c sba_chkjac.c
AR=ar
RANLIB=ranlib
MAKE=make

all: libsba.so dem

libsba.a: $(OBJS)
	$(AR) crv libsba.a $(OBJS)
	$(RANLIB) libsba.a

libsba.so: 
	$(CC) -fPIC -c sba.h sba_levmar.c sba_levmar_wrap.c sba_lapack.c sba_crsm.c sba_chkjac.c $(CFLAGS)
	$(CC) -shared -Wall -o libsba.so $(OBJS) -llapack -lblas $(CFLAGS)

install: 
	@echo To install libsba.so, do the following as root:
	@echo cp libsba.so /usr/local/lib/libsba.so.1.6.1
	@echo chmod a-x /usr/local/lib/libsba.so.1.6.1
	@echo rm /usr/local/lib/libsba.so
	@echo ln -s /usr/local/lib/libsba.so.1.6.1 /usr/local/lib/libsba.so
	@echo cp sba.h /usr/local/include/sba.h
	@echo ldconfig

sba_levmar.o: sba.h sba_chkjac.h compiler.h
sba_levmar_wrap.o: sba.h
sba_lapack.o: sba.h compiler.h
sba_crsm.o: sba.h
sba_chkjac.o: sba.h sba_chkjac.h compiler.h

dem:
	cd demo; $(MAKE)

clean:
	@rm -f $(OBJS)
	cd demo; $(MAKE) clean
	cd matlab; $(MAKE) clean

realclean cleanall: clean
	@rm -f libsba.a

depend:
	makedepend -f Makefile $(SRCS)

# DO NOT DELETE THIS LINE -- make depend depends on it.
