#!/usr/bin/env/python
"""
drivers.py
Dennis Evangelista, 2013
"""


import sba.routines
import sba.options
import sba.info
import logging
import numpy as np
import ctypes
import time
import copy

ERROR = -1

def SparseBundleAdjust(cameras,points,
                       options=None,
                       projection=None, jacobian=None):
    """
    Take Cameras object, Points object and optional Options object
    and performs sparse bundle adjustment using Lourakis' sba C library

    Returns updated newcameras, newpoints, and an Information object
    giving results from the sparse bundle adjustment. 

    Additional options can be specified with an optional Options input.

    Also, Python callable functions can be specified for projection and
    Jacobian; these must conform to the return type and argtypes expected
    by sba in C - see sba.projections for more information.
    """

    logging.debug("sba.SparseBundleAdjust using {0}".format(points))
    logging.debug("sba.SparseBundleAdjust using {0}".format(cameras))

    if options==None:
        options = sba.options.Options.fromInput(cameras,points)

    info = sba.info.Information.fromInput(cameras,points,options)
    logging.debug("sba.SparseBundleAdjust created Information() object")
    newinfo = sba.info.CInfoType()

    routine = sba.routines.whichRoutine(options)
    info.method = routine.__str__
    logging.debug("sba.SparseBundleAdjust using routine {0}".format(routine.__str__))

    if projection==None:
        projection = sba.projections.whichProjection(options)
        logging.debug("sba.SparseBundleAdjust using projection {0}".format(projection.__str__))

    if jacobian==None:
        jacobian = sba.projections.whichJacobian(options)
        logging.debug("sba.SparseBundleAdjust using Jacobian {0}".format(jacobian.__str__))


    # start timer here
    logging.debug("sba.SparseBundleAdjust starting timer and calling C")
    info.t0 = time.clock()

    # dispatcher
    if options.expert:
        if options.motstruct==sba.options.OPTS_MOTSTRUCT:
            cameras.zeroLocalRotationEstimates()
            P = copy.deepcopy(np.hstack((cameras.Aravel(),points.Bravel())))
            newX = points.XtoC() # Have to do this to get back numpy array
            print("n {0}".format(points.n))
            print("m {0}".format(cameras.ncameras))
            print("x {0}".format(newX[0:5]))
            print("p {0}".format(P[0:5]))
            print("all x {0}".format(newX))
            print("length x {0}".format(len(newX)))
            n = routine(ctypes.c_int(points.n), #const int
                        ctypes.c_int(0), #const int
                        ctypes.c_int(cameras.ncameras), # const int
                        ctypes.c_int(0), #const int
                        points.vmaskToC(),
                        P.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                        ctypes.c_int(options.cnp), # const int
                        ctypes.c_int(options.pnp), # const int
                        newX.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                        points.covxToC(),
                        ctypes.c_int(options.mnp), # const int
                        sba.FuncProjection(projection), # pointer
                        sba.FjacProjJacobian(jacobian), # pointer
                        options.globsToC(), # void*?
                        options.maxIterToC(), # const
                        options.verboseToC(), # const
                        options.optsToC(), # const
                        newinfo) # pointer, i want it back
        elif options.motstruct==sba.options.OPTS_MOT:
            cameras.zeroLocalRotationEstimates()
            options.ptparams = points.Bravel()
            P = copy.deepcopy(cameras.Aravel())
            newX = points.XtoC()
            n = routine(ctypes.c_int(points.n),
                        ctypes.c_int(cameras.ncameras),
                        ctypes.c_int(0),
                        points.vmaskToC(),
                        P.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                        ctypes.c_int(options.cnp),
                        newX.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                        points.covxToC(),
                        ctypes.c_int(options.mnp),
                        sba.FuncProjection(projection),
                        sba.FjacProjJacobian(jacobian),
                        options.globsToC(),
                        options.maxIterToC(),
                        options.verboseToC(),
                        options.optsToC(),
                        newinfo)
        elif options.motstruct==sba.options.OPTS_STRUCT:
            P = copy.deepcopy(points.Bravel())
            options.camparams = cameras.Aravel()
            newX = points.XtoC()
            n = routine(ctypes.c_int(points.n),
                        ctypes.c_int(0),
                        ctypes.c_int(cameras.ncameras),
                        points.vmaskToC(),
                        P.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
#                        ctypes.c_int(options.cnp),
                        ctypes.c_int(options.pnp),
                        newX.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                        points.covxToC(),
                        ctypes.c_int(options.mnp),
                        sba.FuncProjection(projection),
                        sba.FjacProjJacobian(jacobian),
                        options.globsToC(),
                        options.maxIterToC(),
                        options.verboseToC(),
                        options.optsToC(),
                        newinfo)
        else:
            raise sba.errors.SbaError("Unrecognized sba option")
    else:
        logging.warning("Simple mode untested!")
        if options.motstruct==sba.options.OPTS_MOTSTRUCT:
            cameras.zeroLocalRotationEstimates()
            P = np.hstack((cameras.Aravel(),points.Bravel())).copy()
            newX = points.XtoC() # Have to do this to get back numpy array
            n = routine(ctypes.c_int(points.n), #const int
                        ctypes.c_int(0), #const int
                        ctypes.c_int(cameras.ncameras), # const int
                        ctypes.c_int(0), #const int
                        points.vmaskToC(),
                        P.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                        ctypes.c_int(options.cnp), # const int
                        ctypes.c_int(options.pnp), # const int
                        newX.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                        points.covxToC(),
                        ctypes.c_int(options.mnp), # const int
                        sba.FuncProjectionSimple(projection), # pointer
                        sba.FjacProjJacobianSimple(jacobian), # pointer
                        options.globsToC(), # void*?
                        options.maxIterToC(), # const
                        options.verboseToC(), # const
                        options.optsToC(), # const
                        newinfo) # pointer, i want it back
        elif options.motstruct==sba.options.OPTS_MOT:
            cameras.zeroLocalRotationEstimates()
            P = copy.deepcopy(cameras.Aravel())
            options.ptsparams = points.Bravel()
            newX = points.XtoC()
            n = routine(ctypes.c_int(points.n),
                        ctypes.c_int(cameras.ncameras),
                        ctypes.c_int(0),
                        points.vmaskToC(),
                        P.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                        ctypes.c_int(options.cnp),
                        ctypes.c_int(options.pnp),
                        newX.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                        points.covxToC(),
                        ctypes.c_int(options.mnp),
                        sba.FuncProjectionSimple(projection),
                        sba.FjacProjJacobianSimple(jacobian),
                        options.globsToC(),
                        options.maxIterToC(),
                        options.verboseToC(),
                        options.optsToC(),
                        newinfo)
        elif options.motstruct==sba.options.OPTS_STRUCT:
            P = copy.deepcopy(points.Bravel())
            options.camparams = cameras.Aravel()
            newX = points.XtoC()
            n = routine(ctypes.c_int(points.n),
                        ctypes.c_int(0),
                        ctypes.c_int(cameras.ncameras),
                        points.vmaskToC(),
                        P.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                        ctypes.c_int(options.cnp),
                        ctypes.c_int(options.pnp),
                        newX.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                        points.covxToC(),
                        ctypes.c_int(options.mnp),
                        sba.FuncProjectionSimple(projection),
                        sba.FjacProjJacobianSimple(jacobian),
                        options.globsToC(),
                        options.maxIterToC(),
                        options.verboseToC(),
                        options.optsToC(),
                        newinfo)
        else:
            raise sba.errors.SbaError("Unrecognized sba option")
        pass
    
    # stop timer here
    info.tf = time.clock()
    logging.debug("sba.SparseBundleAdjust returned from C and stopped timer")
    
    if n == ERROR:
        raise sba.errors.SbaError("sparseBundleAdjust() C call returned error {0}".format(n))





    # unpack some stuff
    logging.debug("sba.SparseBundleAdjust unpacking results")
    newcameras = copy.deepcopy(cameras)
    newpoints = copy.deepcopy(points)
    info.fromC(newinfo)
    if options.motstruct==sba.options.OPTS_MOTSTRUCT:
        newcameras.camarray = P[0:options.ncameras*options.cnp].reshape(options.ncameras,options.cnp) # unpack A from P first part
        newcameras.combineLocalAndInitialRotation(options.rot0params) # combine local rotation estimates with initial ones
        newpoints.B = P[options.ncameras*options.cnp:].reshape(newpoints.n,options.pnp) # unpack B from P second part

    elif options.motstruct==sba.options.OPTS_MOT:
        newcameras.camarray = P.reshape(options.ncameras,options.cnp) # unload A camera params only
        newcameras.combineLocalAndInitialRotation(options.rot0params) # combine local rotation estimates with initial ones

    elif options.motstruct==sba.options.OPTS_STRUCT:
        newpoints.B = P.reshape(options.n,options.pnp) # unload B structure params only

    else:
        raise sba.errors.SbaError("Unrecognized option for howto: {0}".format(options.motstruct))
    # refined motion and structure are now in motstruct
    
    logging.debug("sba.SparseBundleAdjust returning results")
    return newcameras,newpoints,info












