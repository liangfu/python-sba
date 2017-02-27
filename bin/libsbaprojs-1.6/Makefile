#
# Makefile for Sparse Bundle Adjustment projections and Jacobians
# modified from Lourakis' sba demo eucsbademo.c
#
CC=gcc
ARCHFLAGS=-isysroot /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX.sdk -arch x86_64 -I/usr/local/include -L/usr/local/lib
CFLAGS=$(ARCHFLAGS) -I.. -O3 -Wall #-g -pg
OBJS=sbaprojs.o imgproj.o 
SRCS=sbaprojs.c imgproj.c


LAPACKLIBS=-llapack -lblas #-lf2c # On systems with a FORTRAN (not f2c'ed) version of LAPACK, -lf2c is
                                 # not necessary; on others -lf2c is equivalent to -lF77 -lI77

#LAPACKLIBS=-L/usr/local/atlas/lib -llapack -lcblas -lf77blas -latlas -lf2c # This works with the ATLAS updated lapack and Linux_P4SSE2
                                                                            # from http://www.netlib.org/atlas/archives/linux/

#LAPACKLIBS=-llapack -lgoto -lpthread -lf2c # This works with GotoBLAS
                                            # from http://www.tacc.utexas.edu/resources/software/

#LAPACKLIBS=-L/opt/intel/mkl/8.0.1/lib/32/ -lmkl_lapack -lmkl_ia32 -lguide -lf2c # This works with MKL 8.0.1 from
                                            # http://www.intel.com/cd/software/products/asmo-na/eng/perflib/mkl/index.htm

LIBS=-lsba $(LAPACKLIBS) -lm
LDFLAGS=$(ARCHFLAGS) -L..

libsbaprojs.so: 
	$(CC) -fPIC -c sbaprojs.h sbaprojs.c imgproj.c $(CFLAGS)
	$(CC) -shared -Wall $(LDFLAGS) $(OBJS) -o libsbaprojs.so $(LIBS)
# do i need to include -lsba and have libsba.so somewhere?!
# answer is yes - they must be in the default Ubuntu 12.04 places
# /usr/local/lib/libsba.so
# /usr/local/include/sba.h

install:
	@chmod a-x libsbaprojs.so
	@cp libsbaprojs.so /usr/local/lib/libsbaprojs.so.1.6.0
	@rm -f /usr/local/lib/libsbaprojs.so
	@ln -s /usr/local/lib/libsbaprojs.so.1.6.0 /usr/local/lib/libsbaprojs.so
	@cp sbaprojs.h /usr/local/include
	@ldconfig
# this must be done as root (sudo)

clean:
	@rm -f $(OBJS)

realclean cleanall: clean
	@rm -f libsbaprojs.so

depend:
	makedepend -f Makefile $(SRCS)

# DO NOT DELETE THIS LINE -- make depend depends on it.
