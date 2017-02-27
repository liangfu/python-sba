#!/usr/bin/env python
"""
routines.py
Python wrapper for sba C library by Manolis Lourakis, ver 1.6
Wrapper by Dennis Evangelista, 2013

This part wraps the actual routines, eg sba_motstr_levmar_x
"""

import ctypes
import logging
import numpy as np
from sba.crsm import Crsm
from sba.projections import FuncProjection, FjacProjJacobian, FuncProjectionSimple, FjacProjJacobianSimple
from sba.options import OPTS_SIZE
from sba.info import INFO_SIZE
import sba.options
import sys


if sys.platform == 'linux2':
    # This line actually loads the sba shared object
    _libsba = ctypes.CDLL("/usr/local/lib/libsba.so")
    # On Ubuntu 12.04 this is /usr/local/lib/libsba.so, which is a symbolic link
    # to /usr/local/lib/libsba.so.1.6.4 currently

elif sys.platform == 'darwin':
    # On Mac you may wish to use .dylib extension instead; be sure to rename
    # the .so in the Makefile or when installing and use the correct Mac 
    # locations here instead.
    try:
        _libsba = ctypes.CDLL("libsba.dylib")
    except Exception:
        _libsba = ctypes.CDLL("libsba.so")
    # THIS MAY NEED DEBUGGING

elif sys.platform == 'win32' or sys.platform == 'win64':
    # on Windows,
    _libsba = ctypes.CDLL("libsba.dll")
    # THIS MAY NEED DEBUGGING

# "expert" drivers from sba.h (BU modified version 1.6.1)
# function prototype for sba_motstr_levmar_x
_libsba.sba_motstr_levmar_x.argtypes = [ctypes.c_int, # n num of points
                                        ctypes.c_int, # ncon number not to be modified
                                        ctypes.c_int, # m num of images
                                        ctypes.c_int, # mcon number not to be modified
                                        ctypes.c_char_p, # visibility mask
                                        ctypes.POINTER(ctypes.c_double), # initial parameter vector
                                        ctypes.c_int, # cnp number of params for ONE camera
                                        ctypes.c_int, # pnp number of params for ONE point eg 3
                                        ctypes.POINTER(ctypes.c_double), # x measurements vector
                                        ctypes.POINTER(ctypes.c_double), # covx measurements covariance
                                        ctypes.c_int, # mnp number of params for EACH measurement eg 2
                                        FuncProjection, # functional relation describing measurements?
                                        FjacProjJacobian, # evaluate sparse Jacobian dX/dp
                                        ctypes.c_void_p, # adata, additional data to func, fjac
                                        ctypes.c_int, # itmax max iterations
                                        ctypes.c_int, # verbosity
                                        ctypes.c_double * OPTS_SIZE, # opts options
                                        ctypes.c_double * INFO_SIZE] # info 
_libsba.sba_motstr_levmar_x.restype = ctypes.c_int 
# Python handle
_MotionStructureExpert = _libsba.sba_motstr_levmar_x
_MotionStructureExpert.__doc__ = """
function prototype for sba routine from libsba
long int sba_motstr_levmar_x(long int n, long int ncon, long int m, long int mcon,
char *vmask, double *p, long int cnp, long int pnp, double *x, double *covx,
long int mnp, (func - a projection function), (fjac - a Jacobian function), 
void *adata, long int itmax, long int verbose, double opts[OPTS_SIZE], 
double info[INFO_SIZE]

This does sba adjustment with motion and structure in expert mode. 
"""
_MotionStructureExpert.__str__ = "MotionStructureExpert"


# function prototype for sba_mot_levmar_x
_libsba.sba_mot_levmar_x.argtypes = [ctypes.c_int, # n num of points
                                     ctypes.c_int, # m num of images
                                     ctypes.c_int, # mcon number not to be modified
                                     ctypes.c_char_p, # visibility mask
                                     ctypes.POINTER(ctypes.c_double), # initial parameter vector
                                     ctypes.c_int, # cnp number of params for ONE camera
                                     #ctypes.c_int, # pnp number of params for ONE point eg 3
                                     ctypes.POINTER(ctypes.c_double), # x measurements vector
                                     ctypes.POINTER(ctypes.c_double), # covx measurements covariance
                                     ctypes.c_int, # mnp number of params for EACH measurement eg 2
                                     FuncProjection, # functional relation describing measurements?
                                     FjacProjJacobian, # evaluate sparse Jacobian dX/dp
                                     ctypes.c_void_p, # adata, additional data to func, fjac
                                     ctypes.c_int, # itmax max iterations
                                     ctypes.c_int, # verbosity
                                     ctypes.c_double * OPTS_SIZE, # opts options
                                     ctypes.c_double * INFO_SIZE] # info 
_libsba.sba_mot_levmar_x.restype = ctypes.c_int 
# Python handle
_MotionExpert = _libsba.sba_mot_levmar_x
_MotionExpert.__doc__ = """
function prototype for sba routine from libsba
long int sba_mot_levmar_x(long int n, long int m, long int mcon,
char *vmask, double *p, long int cnp, long int pnp, double *x, double *covx,
long int mnp, (func - a projection function), (fjac - a Jacobian function), 
void *adata, long int itmax, long int verbose, double opts[OPTS_SIZE], 
double info[INFO_SIZE]

This does sba adjustment with motion in expert mode. 
"""
_MotionExpert.__str__ = "MotionExpert"

# function prototype for sba_str_levmar_x
_libsba.sba_str_levmar_x.argtypes = [ctypes.c_int, # n num of points
                                     ctypes.c_int, # ncon number not to be modified
                                     ctypes.c_int, # m num of images
                                     ctypes.c_char_p, # visibility mask
                                     ctypes.POINTER(ctypes.c_double), # initial parameter vector
                                     #ctypes.c_int, # cnp number of params for ONE camera
                                     ctypes.c_int, # pnp number of params for ONE point eg 3
                                     ctypes.POINTER(ctypes.c_double), # x measurements vector
                                     ctypes.POINTER(ctypes.c_double), # covx measurements covariance
                                     ctypes.c_int, # mnp number of params for EACH measurement eg 2
                                     FuncProjection, # functional relation describing measurements?
                                     FjacProjJacobian, # evaluate sparse Jacobian dX/dp
                                     ctypes.c_void_p, # adata, additional data to func, fjac
                                     ctypes.c_int, # itmax max iterations
                                     ctypes.c_int, # verbosity
                                     ctypes.c_double * OPTS_SIZE, # opts options
                                     ctypes.c_double * INFO_SIZE] # info 
_libsba.sba_str_levmar_x.restype = ctypes.c_int 
# Python handle
_StructureExpert = _libsba.sba_str_levmar_x
_StructureExpert.__doc__ = """
function prototype for sba routine from libsba
long int sba_motstr_levmar_x(long int n, long int ncon, long int m
char *vmask, double *p, long int cnp, long int pnp, double *x, double *covx,
long int mnp, (func - a projection function), (fjac - a Jacobian function), 
void *adata, long int itmax, long int verbose, double opts[OPTS_SIZE], 
double info[INFO_SIZE]

This does sba adjustment with structure in expert mode. 
"""
_StructureExpert.__str__ = "StructureExpert"




# Simple drivers
_libsba.sba_motstr_levmar.argtypes = [ctypes.c_int, # n
                                      ctypes.c_int, # ncon
                                      ctypes.c_int, # m
                                      ctypes.c_int, # mcon
                                      ctypes.c_char_p, # *vmask
                                      ctypes.POINTER(ctypes.c_double), # *p
                                      ctypes.c_int, # cnp
                                      ctypes.c_int, # pnp
                                      ctypes.POINTER(ctypes.c_double), # *x
                                      ctypes.POINTER(ctypes.c_double), # *covx
                                      ctypes.c_int, # mnp
                                      FuncProjectionSimple,
                                      FjacProjJacobianSimple,
                                      ctypes.c_void_p, # *adata
                                      ctypes.c_int, # itmax
                                      ctypes.c_int, # verbose
                                      ctypes.c_double * OPTS_SIZE, # double *opts
                                      ctypes.c_double * INFO_SIZE] # double *info
_libsba.sba_motstr_levmar.restype = ctypes.c_int
_MotionStructureSimple = _libsba.sba_motstr_levmar
_MotionStructureSimple.__doc__ = """
function prototype for sba routine from libsba
long int sba_motstr_levmar(long int n, long int ncon, long int m, long int mcon,
char *vmask, double *p, long int cnp, long int pnp, double *x, double *covx,
long int mnp, (func - a projection function), (fjac - a Jacobian function), 
void *adata, long int itmax, long int verbose, double opts[OPTS_SIZE], 
double info[INFO_SIZE]

This does sba adjustment with motion and structure in simple mode. 
"""
_MotionStructureSimple.__str__ = "MotionStructureSimple"

_libsba.sba_mot_levmar.argtypes = [ctypes.c_int, # n
                                   ctypes.c_int, # m
                                   ctypes.c_int, # mcon
                                   ctypes.c_char_p, # *vmask
                                   ctypes.POINTER(ctypes.c_double), # *p
                                   ctypes.c_int, # cnp
                                   ctypes.c_int, # pnp
                                   ctypes.POINTER(ctypes.c_double), # *x
                                   ctypes.POINTER(ctypes.c_double), # *covx
                                   ctypes.c_int, # mnp
                                   FuncProjectionSimple,
                                   FjacProjJacobianSimple,
                                   ctypes.c_void_p, # *adata
                                   ctypes.c_int, # itmax
                                   ctypes.c_int, # verbose
                                   ctypes.c_double * OPTS_SIZE, # double *opts
                                   ctypes.c_double * INFO_SIZE] # double *info
_libsba.sba_mot_levmar.restype = ctypes.c_int
_MotionSimple = _libsba.sba_mot_levmar
_MotionSimple.__doc__ = """
function prototype for sba routine from libsba
long int sba_mot_levmar(long int n, long int m, long int mcon,
char *vmask, double *p, long int cnp, long int pnp, double *x, double *covx,
long int mnp, (func - a projection function), (fjac - a Jacobian function), 
void *adata, long int itmax, long int verbose, double opts[OPTS_SIZE], 
double info[INFO_SIZE]

This does sba adjustment with motion and structure in simple mode. 
"""
_MotionSimple.__str__ = "MotionSimple"

_libsba.sba_str_levmar.argtypes = [ctypes.c_int, # n
                                   ctypes.c_int, # ncon
                                   ctypes.c_int, # m
                                   ctypes.c_char_p, # *vmask
                                   ctypes.POINTER(ctypes.c_double), # *p
                                   ctypes.c_int, # cnp
                                   ctypes.c_int, # pnp
                                   ctypes.POINTER(ctypes.c_double), # *x
                                   ctypes.POINTER(ctypes.c_double), # *covx
                                   ctypes.c_int, # mnp
                                   FuncProjectionSimple,
                                   FjacProjJacobianSimple,
                                   ctypes.c_void_p, # *adata
                                   ctypes.c_int, # itmax
                                   ctypes.c_int, # verbose
                                   ctypes.c_double * OPTS_SIZE, # double *opts
                                   ctypes.c_double * INFO_SIZE] # double *info
_libsba.sba_str_levmar.restype = ctypes.c_int
_StructureSimple = _libsba.sba_str_levmar
_StructureSimple.__doc__ = """
function prototype for sba routine from libsba
long int sba_motstr_levmar(long int n, long int ncon, long int m,
char *vmask, double *p, long int cnp, long int pnp, double *x, double *covx,
long int mnp, (func - a projection function), (fjac - a Jacobian function), 
void *adata, long int itmax, long int verbose, double opts[OPTS_SIZE], 
double info[INFO_SIZE]

This does sba adjustment with motion and structure in simple mode. 
"""
_StructureSimple.__str__ = "StructureSimple"



"""
Not implemented:
CRSM routines 
Linear algebra helper routines
"""

def whichRoutine(options):
    """
    Dispatcher for which sba routine to use based on options
    
    input:
    options - Options object specifying sba options

    returns:
    a Python pointer to a C function (ctypes.funcptr) for sba projection.
    e.g. sba_levmar_motstr_x, sba_levmar_str_x, etc. 

    Raises SbaError if options are not recognized
    """
    #logging.debug("routines.whichRoutine() called with {0}".format(options))
    if options.expert:
        if options.motstruct==sba.options.OPTS_MOTSTRUCT:
            return _MotionStructureExpert
        elif options.motstruct==sba.options.OPTS_MOT:
            return _MotionExpert
        elif options.motstruct==sba.options.OPTS_STRUCT:
            return _StructureExpert
        else:
            raise SbaError("Unrecognized sba option")
    else:
        if options.motstruct==sba.options.OPTS_MOTSTRUCT:
            return _MotionStructureSimple
        elif options.motstruct==sba.options.OPTS_MOT:
            return _MotionSimple
        elif options.motstruct==sba.options.OPTS_STRUCT:
            return _StructureSimple
        else:
            raise SbaError("Unrecognized sba option")
