===============================================================
Python wrapper for Lourakis' sparse bundle adjustment C library 
===============================================================

Enjoy! The most recent version can be obtained from bitbucket via::

    hg clone ssh://hg@bitbucket.org/devangel77b/python-sba

Typical usage
=============

The main way to use this is as follows::

    import sba

    cameras = sba.Cameras.fromTxt('cams.txt')
    points = sba.Points.fromTxt('pts.txt',cameras.ncameras)
    newcams, newpts, info = sba.SparseBundleAdjust(cameras,points)

If you wish to alter the default and autodetected options, you can
create an Options object and change it, and then pass it to sba::

    options = sba.Options.fromInputs(points,cameras)
    # can also update options.XXX to appropriate values
    newcams,newpts,info = sba.SparseBundleAdjust(cameras,points,options)

Hopefully this is cleaner than the original way to call it in C. 
    
Contributors
============

The original sba C library was written by Manolis Lourakis and is 
described in Lourakis, Manolis I A and Antonis A Argyros (2004), "The design 
and implementation of a generic sparse bundle adjustment software package 
based on the Levenberg-Marquardt algorithm", FOURTH_ICS TR-340.


Thanks also to
==============

Manolis Lourakis and Antonis Argyros, Ty Hedrick, Evan Bluhm, my mom and the academy
