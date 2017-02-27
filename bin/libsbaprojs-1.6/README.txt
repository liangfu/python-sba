========================================================
libsbaprojs - Projections and Jacobians for use with sba
========================================================

This is a C library (shared object) of the projections and
Jacobians provided by Lourakis in sba/demo/eucsbademo.c. The 
functions have the correct return type and argument list for
use with Lourakis' sparse bundle adjustment routines.  This
shared library can be wrapped with Python (example provided
here... though the real wrapping happens in python-sba in
the projections.py submodule). 

The most recent version can be obtained from bitbucket via::

    hg clone ssh://hg@bitbucket.org/devangel77b/libsbaprojs

Installation
============
Prerequisites to build this library you must have libsba installed (see
python-sba; this requires compling and installing Lourakis' sba library 
as libsba.so, installing it in /usr/local/lib, and copying the header sba.h 
to /usr/local/include.)

To actually build the library::
    make libsbaprojs.so

To copy and link the library in the default places for Ubuntu 12.04::
    sudo make install 

The last line requires administrator permissions and will copy the
library to /usr/local/lib/libsbaprojs.so.1.6.0; create symbolic link 
/usr/local/lib/libsbaprojs.so; set the permissions and copy sbaprojs.h
to /usr/local/include. 

Known limitations
=================
This has only been tested on Ubuntu 12.04 with mostly vanilla installs
of everything, straight from apt.  I have never tested it with any other
flavors of C, LAPACK, etc etc etc and do not intend to. 

Acknowledgements
================

The original sba C library was written by Manolis Lourakis and is 
described in Lourakis, Manolis I A and Antonis A Argyros (2004), "The design 
and implementation of a generic sparse bundle adjustment software package 
based on the Levenberg-Marquardt algorithm", FOURTH_ICS TR-340. The projection 
functions provided here are taken straight from Lourakis' demo program 
eucsbademo.c.  


If using this package in research work, we would appreciate you citing it: D The
riault, N Fuller, B Jackson, E Bluhm, D Evangelista, Z Wu, M Betke, and T Hedric
k (in review). A method for accurate multi-camera field videography. J exp Bio
l. The BibTeX entry is::

    @article{Theriault:2014,
      author = {Theriault, D and Fuller, N and Jackson, B and Bluhm, E and Evang
elista, D and Wu, Z and Betke, M and Hedrick, T},
      title = {A method for accurate multi-camera field videography}
      journal = {J exp Biol},
      year = {2014},
      volume = {217},
      pages = {1843--1848}
    }
