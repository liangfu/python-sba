#!/usr/bin/env python

import argparse
import logging
import numpy as np
import sba

__version__ = "0.0"



def sbaDriver(camname,ptsname,intrname=None,camoname=None,ptsoname=None):
    logging.debug("sbaDriver() {0} {1}".format(camname,ptsname))

    if intrname is None:
        cameras = sba.Cameras.fromTxt(camname) # Load cameras from file
    else:
        # Intrinsic file for all identical cameras not yet implemented
        logging.warning("Fixed intrinsic parameters not yet implemented")
        cameras = sba.Cameras.fromTxtWithIntr(camname,intrname)
    points = sba.Points.fromTxt(ptsname,cameras.ncameras) # Load points

    options = sba.Options.fromInput(cameras,points)
    #options.camera = sba.OPTS_CAMS_NODIST # only use this if passing in 5+3+3
    options.nccalib=sba.OPTS_FIX2_INTR # fix all intrinsics
    # If you wish to fix the intrinsics do so here by setting options
    options.ncdist = sba.OPTS_FIX5_DIST # fix all distortion coeffs

    newcameras,newpoints,info = sba.SparseBundleAdjust(cameras,points,options)
    info.printResults()

    if camoname:
        newcameras.toTxt(camoname)
    if ptsoname:
        newpoints.toTxt(ptsoname)

    # other ways to save output:
    #  - all can be pickled 
    #  - save options to yaml (and load from)?






if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Sparse bundle adjustment using Lourakis' sba library")

    parser.add_argument("camname",
                        help="fx px py AR s r2 r4 t1 t2 r6 q0 qi qj qk tx ty tz")
    parser.add_argument("ptsname",
                        help="X Y Z Nvisible Cam0 x0 y0 Cam1 x1 y1")

    parser.add_argument("camo",nargs="?",default=None,help="Output filename for cameras txt")
    parser.add_argument("ptso",nargs="?",default=None,help="Output filename for 3D points txt")

    parser.add_argument("--intrname",help="Input intrinsics, not implemented")
    parser.add_argument("-v","--verbose",action="store_true",
                        help="Toggle verbosity")
    parser.add_argument("--version",action="version",
                        version = "%(prog)s {0}".format(__version__))
    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    sbaDriver(args.camname,args.ptsname,camoname=args.camo,ptsoname=args.ptso)



    
