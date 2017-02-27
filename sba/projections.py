#!/usr/bin/env python
"""
projections.py
Python wrapper for projections and jacobians from eucsbademo by Lourakis
Wrapper by Dennis Evangelista, 2013
"""

import ctypes
import logging
import numpy as np
import sba.options
from sba.crsm import Crsm
#from sba.options import Globs_
from sba.errors import SbaError
import sys


if sys.platform == 'linux2':
    # This line actually loads the sba shared object from /usr/local/lib
    _libsbaprojs = ctypes.CDLL("/usr/local/lib/libsbaprojs.so")
    # On Ubuntu 12.04, this should be in /usr/local/lib/libsbaprojs.so as a link
    # to /usr/local/lib/libsbaprojs.so.1.6.4

elif sys.platform == 'darwin':
    # On Mac, you may wish to use .dylib rather than .so and set the paths
    # accordingly. 
    _libsbaprojs = ctypes.CDLL("libsbaprojs.dylib")
    # THIS MAY NEED DEBUGGING

elif sys.platform == 'win32' or sys.platform == 'win64':
    # On Windows, 
    _libsbaprojs = ctypes.CDLL("libsbaprojs.dll")
    # THIS MAY NEED DEBUGGING


# This is used for going from Python into C
FuncProjection = ctypes.CFUNCTYPE(None, # return type void
                                  ctypes.POINTER(ctypes.c_double), #p
                                  ctypes.POINTER(Crsm), #idxij
                                  ctypes.POINTER(ctypes.c_int), #rcidxs
                                  ctypes.POINTER(ctypes.c_int), #rcsubs
                                  ctypes.POINTER(ctypes.c_double), #hx
                                  ctypes.c_void_p) #adata
FuncProjection.__doc__ = """
Function prototype for callable projection func in expert mode sba:
void (*func) (double *p, struct sba_crsm *idxij, int *rcidxs, 
int *rcsubs, double *hx, void *adata)

Use to wrap a Python function when sending a projection function from 
Python into C.
"""

# This is used for going from Python into C
FjacProjJacobian = ctypes.CFUNCTYPE(None, # return type void
                                    ctypes.POINTER(ctypes.c_double), #p
                                    ctypes.POINTER(Crsm), #idxij
                                    ctypes.POINTER(ctypes.c_int), #rcidxs
                                    ctypes.POINTER(ctypes.c_int), #rcsubs
                                    ctypes.POINTER(ctypes.c_double), #jac
                                    ctypes.c_void_p) #adata
FjacProjJacobian.__doc__ = """
Function prototype for callable Jacobian fjac in expert mode sba:
void (*fjac) (double *p, struct sba_crsm *idxij, int *rcidxs, 
int *rcsubs, double *hx, void *adata)

Use to wrap a Python function when sending a Jacobian function from 
Python into C.
"""

FuncProjectionSimple = ctypes.CFUNCTYPE(None, # return type void
                                        ctypes.c_int, # j
                                        ctypes.c_int, # i
                                        ctypes.POINTER(ctypes.c_double), # *bi
                                        ctypes.POINTER(ctypes.c_double), # *xij
                                        ctypes.c_void_p) # adata
FuncProjectionSimple.__doc__ = """
Function prototype for callable projection func in simple mode sba:
void (*proj) (int j, int i, double *bi, double *xij, void *adata)

Use to wrap a Python function when sending a projection function from 
Python into C.
"""

FjacProjJacobianSimple = ctypes.CFUNCTYPE(None, # return type void
                                          ctypes.c_int, # j
                                          ctypes.c_int, # i
                                          ctypes.POINTER(ctypes.c_double), # *bi
                                          ctypes.POINTER(ctypes.c_double), # *xij
                                          ctypes.c_void_p) # adata
FuncProjectionSimple.__doc__ = """
Function prototype for callable Jacobian fjac in simple mode sba:
void (*projac) (int j, int i, double *bi, double *xij, void *adata)

Use to wrap a Python function when sending a Jacobian function from 
Python into C.
"""






# projection functions and their Jacobians
# full sparse bundle adjustment
_libsbaprojs.img_projsRTS_x.argtypes = [ctypes.POINTER(ctypes.c_double),
                                        ctypes.POINTER(Crsm),
                                        ctypes.POINTER(ctypes.c_int),
                                        ctypes.POINTER(ctypes.c_int),
                                        ctypes.POINTER(ctypes.c_double),
                                        ctypes.c_void_p]
_libsbaprojs.img_projsRTS_x.restype = None
_libsbaprojs.img_projsRTS_x.__str__ = "MotionsStructureExpert_NoCameras"
# Python handle
_MotionStructureExpert_NoCameras = _libsbaprojs.img_projsRTS_x
_MotionStructureExpert_NoCameras.__doc__="""
Function prototype for callable projection func in expert mode sba:
void (*func) (double *p, struct sba_crsm *idxij, int *rcidxs, 
int *rcsubs, double *hx, void *adata)

Does projection from P to X. 
"""

_libsbaprojs.img_projsRTS_jac_x.argtypes = [ctypes.POINTER(ctypes.c_double),
                                            ctypes.POINTER(Crsm),
                                            ctypes.POINTER(ctypes.c_int),
                                            ctypes.POINTER(ctypes.c_int),
                                            ctypes.POINTER(ctypes.c_double),
                                            ctypes.c_void_p]
_libsbaprojs.img_projsRTS_jac_x.restype = None
_libsbaprojs.img_projsRTS_jac_x.__str__ = "MotionsStructureJacobianExpert_NoCameras"
# Python handle
_MotionStructureJacobianExpert_NoCameras = _libsbaprojs.img_projsRTS_jac_x
_MotionStructureJacobianExpert_NoCameras.__doc__="""
Function prototype for callable Jacobian fjac in simple mode sba:
void (*fjac) (int j, int i, double *bi, double *xij, void *adata)

Does projection from P to X. 
"""

# sparse bundle adjustment with camera parameters only
_libsbaprojs.img_projsRT_x.argtypes = [ctypes.POINTER(ctypes.c_double),
                                       ctypes.POINTER(Crsm),
                                       ctypes.POINTER(ctypes.c_int),
                                       ctypes.POINTER(ctypes.c_int),
                                       ctypes.POINTER(ctypes.c_double),
                                       ctypes.c_void_p]
_libsbaprojs.img_projsRT_x.restype = None
_libsbaprojs.img_projsRT_x.__str__ = "MotionExpert_NoCameras"
# Python handle
_MotionExpert_NoCameras = _libsbaprojs.img_projsRT_x
_MotionExpert_NoCameras.__doc__="""
Function prototype for callable projection func in expert mode sba:
void (*func) (double *p, struct sba_crsm *idxij, int *rcidxs, 
int *rcsubs, double *hx, void *adata)

Does projection from P to X. 
"""

_libsbaprojs.img_projsRT_jac_x.argtypes = [ctypes.POINTER(ctypes.c_double),
                                           ctypes.POINTER(Crsm),
                                           ctypes.POINTER(ctypes.c_int),
                                           ctypes.POINTER(ctypes.c_int),
                                           ctypes.POINTER(ctypes.c_double),
                                           ctypes.c_void_p]
_libsbaprojs.img_projsRT_jac_x.restype = None
_libsbaprojs.img_projsRT_jac_x.__str__ = "MotionJacobianExpert_NoCameras"
# Python handle
_MotionJacobianExpert_NoCameras = _libsbaprojs.img_projsRT_jac_x
_MotionJacobianExpert_NoCameras.__doc__="""
Function prototype for callable Jacobian fjac in simple mode sba:
void (*fjac) (int j, int i, double *bi, double *xij, void *adata)

Does projection from P to X. 
"""

# sparse bundle adjustment with structure parameters only
_libsbaprojs.img_projsS_x.argtypes = [ctypes.POINTER(ctypes.c_double),
                                      ctypes.POINTER(Crsm),
                                      ctypes.POINTER(ctypes.c_int),
                                      ctypes.POINTER(ctypes.c_int),
                                      ctypes.POINTER(ctypes.c_double),
                                      ctypes.c_void_p]
_libsbaprojs.img_projsS_x.restype = None
_libsbaprojs.img_projsS_x.__str__ = "StructureExpert_NoCameras"
# Python handle
_StructureExpert_NoCameras = _libsbaprojs.img_projsS_x
_StructureExpert_NoCameras.__doc__="""
Function prototype for callable projection func in expert mode sba:
void (*func) (double *p, struct sba_crsm *idxij, int *rcidxs, 
int *rcsubs, double *hx, void *adata)

Does projection from P to X. 
"""

_libsbaprojs.img_projsS_jac_x.argtypes = [ctypes.POINTER(ctypes.c_double),
                                          ctypes.POINTER(Crsm),
                                          ctypes.POINTER(ctypes.c_int),
                                          ctypes.POINTER(ctypes.c_int),
                                          ctypes.POINTER(ctypes.c_double),
                                          ctypes.c_void_p]
_libsbaprojs.img_projsS_jac_x.restype = None
_libsbaprojs.img_projsS_jac_x.__str__ = "StructureJacobianExpert_NoCameras"
# Python handle
_StructureJacobianExpert_NoCameras = _libsbaprojs.img_projsS_jac_x
_StructureJacobianExpert_NoCameras.__doc__="""
Function prototype for callable Jacobian fjac in simple mode sba:
void (*fjac) (int j, int i, double *bi, double *xij, void *adata)

Does projection from P to X. 
"""

# Cameras with distortion
_libsbaprojs.img_projsKDRTS_x.argtypes = [ctypes.POINTER(ctypes.c_double),
                                          ctypes.POINTER(Crsm),
                                          ctypes.POINTER(ctypes.c_int),
                                          ctypes.POINTER(ctypes.c_int),
                                          ctypes.POINTER(ctypes.c_double),
                                          ctypes.c_void_p]
_libsbaprojs.img_projsKDRTS_x.restype = None
_libsbaprojs.img_projsKDRTS_x.__str__ = "MotionsStructureExpert_Cameras"
# Python handle
_MotionStructureExpert_Cameras = _libsbaprojs.img_projsKDRTS_x
_MotionStructureExpert_Cameras.__doc__="""
Function prototype for callable projection func in expert mode sba:
void (*func) (double *p, struct sba_crsm *idxij, int *rcidxs, 
int *rcsubs, double *hx, void *adata)

Does projection from P to X. 
"""

_libsbaprojs.img_projsKDRTS_jac_x.argtypes = [ctypes.POINTER(ctypes.c_double),
                                              ctypes.POINTER(Crsm),
                                              ctypes.POINTER(ctypes.c_int),
                                              ctypes.POINTER(ctypes.c_int),
                                              ctypes.POINTER(ctypes.c_double),
                                              ctypes.c_void_p]
_libsbaprojs.img_projsKDRTS_jac_x.restype = None
_libsbaprojs.img_projsKDRTS_jac_x.__str__ = "MotionsStructureJacobianExpert_Cameras"
# Python handle
_MotionStructureJacobianExpert_Cameras = _libsbaprojs.img_projsKDRTS_jac_x
_MotionStructureJacobianExpert_Cameras.__doc__="""
Function prototype for callable projection fjac in simple mode sba:
void (*fjac) (int j, int i, double *bi, double *xij, void *adata)

Does projection from P to X. 
"""

_libsbaprojs.img_projsKDRT_x.argtypes = [ctypes.POINTER(ctypes.c_double),
                                         ctypes.POINTER(Crsm),
                                         ctypes.POINTER(ctypes.c_int),
                                         ctypes.POINTER(ctypes.c_int),
                                         ctypes.POINTER(ctypes.c_double),
                                         ctypes.c_void_p]
_libsbaprojs.img_projsKDRT_x.restype = None
_libsbaprojs.img_projsKDRT_x.__str__ = "MotionExpert_Cameras"
# Python handle
_MotionExpert_Cameras = _libsbaprojs.img_projsKDRT_x
_MotionExpert_Cameras.__doc__="""
Function prototype for callable projection func in expert mode sba:
void (*func) (double *p, struct sba_crsm *idxij, int *rcidxs, 
int *rcsubs, double *hx, void *adata)

Does projection from P to X. 
"""

_libsbaprojs.img_projsKDRT_jac_x.argtypes = [ctypes.POINTER(ctypes.c_double),
                                             ctypes.POINTER(Crsm),
                                             ctypes.POINTER(ctypes.c_int),
                                             ctypes.POINTER(ctypes.c_int),
                                             ctypes.POINTER(ctypes.c_double),
                                             ctypes.c_void_p]
_libsbaprojs.img_projsKDRT_jac_x.restype = None
_libsbaprojs.img_projsKDRT_jac_x.__str__ = "MotionJacobianExpert_Cameras"
# Python handle
_MotionJacobianExpert_Cameras = _libsbaprojs.img_projsKDRT_jac_x
_MotionJacobianExpert_Cameras.__doc__="""
Function prototype for callable Jacobian fjac in simple mode sba:
void (*fjac) (int j, int i, double *bi, double *xij, void *adata)

Does projection from P to X. 
"""

_libsbaprojs.img_projsKDS_x.argtypes = [ctypes.POINTER(ctypes.c_double),
                                        ctypes.POINTER(Crsm),
                                        ctypes.POINTER(ctypes.c_int),
                                        ctypes.POINTER(ctypes.c_int),
                                        ctypes.POINTER(ctypes.c_double),
                                        ctypes.c_void_p]
_libsbaprojs.img_projsKDS_x.restype = None
_libsbaprojs.img_projsKDS_x.__str__ = "StructureExpert_Cameras"
# Python handle
_StructureExpert_NoCameras = _libsbaprojs.img_projsKDS_x
_StructureExpert_NoCameras.__doc__="""
Function prototype for callable projection func in expert mode sba:
void (*func) (double *p, struct sba_crsm *idxij, int *rcidxs, 
int *rcsubs, double *hx, void *adata)

Does projection from P to X. 
"""

_libsbaprojs.img_projsKDS_jac_x.argtypes = [ctypes.POINTER(ctypes.c_double),
                                            ctypes.POINTER(Crsm),
                                            ctypes.POINTER(ctypes.c_int),
                                            ctypes.POINTER(ctypes.c_int),
                                            ctypes.POINTER(ctypes.c_double),
                                            ctypes.c_void_p]
_libsbaprojs.img_projsKDS_jac_x.restype = None
_libsbaprojs.img_projsKDS_jac_x.__str__ = "StructureJacobianExpert_Cameras"
# Python handle
_StructureJacobianExpert_Cameras = _libsbaprojs.img_projsKDS_jac_x
_StructureJacobianExpert_Cameras.__doc__="""
Function prototype for callable projection fjac in simple mode sba:
void (*fjac) (int j, int i, double *bi, double *xij, void *adata)

Does projection from P to X. 
"""

# Ideal cameras with no distortion
_libsbaprojs.img_projsKRTS_x.argtypes = [ctypes.POINTER(ctypes.c_double),
                                         ctypes.POINTER(Crsm),
                                         ctypes.POINTER(ctypes.c_int),
                                         ctypes.POINTER(ctypes.c_int),
                                         ctypes.POINTER(ctypes.c_double),
                                         ctypes.c_void_p]
_libsbaprojs.img_projsKRTS_x.restype = None
_libsbaprojs.img_projsKRTS_x.__str__ = "MotionStructureExpert_NoDistortion"
# Python handle
_MotionStructureExpert_NoDistortion = _libsbaprojs.img_projsKRTS_x
_MotionStructureExpert_NoDistortion.__doc__="""
Function prototype for callable projection func in expert mode sba:
void (*func) (double *p, struct sba_crsm *idxij, int *rcidxs, 
int *rcsubs, double *hx, void *adata)

Does projection from P to X. 
"""

_libsbaprojs.img_projsKRTS_jac_x.argtypes = [ctypes.POINTER(ctypes.c_double),
                                             ctypes.POINTER(Crsm),
                                             ctypes.POINTER(ctypes.c_int),
                                             ctypes.POINTER(ctypes.c_int),
                                             ctypes.POINTER(ctypes.c_double),
                                             ctypes.c_void_p]
_libsbaprojs.img_projsKRTS_jac_x.restype = None
_libsbaprojs.img_projsKRTS_jac_x.__str__ = "MotionStructureJacobianExpert_NoDistortion"
# Python handle
_MotionStructureJacobianExpert_NoDistortion = _libsbaprojs.img_projsKRTS_jac_x
_MotionStructureJacobianExpert_NoDistortion.__doc__="""
Function prototype for callable projection fjac in simple mode sba:
void (*fjac) (int j, int i, double *bi, double *xij, void *adata)

Does projection from P to X. 
"""

_libsbaprojs.img_projsKRT_x.argtypes = [ctypes.POINTER(ctypes.c_double),
                                        ctypes.POINTER(Crsm),
                                        ctypes.POINTER(ctypes.c_int),
                                        ctypes.POINTER(ctypes.c_int),
                                        ctypes.POINTER(ctypes.c_double),
                                        ctypes.c_void_p]
_libsbaprojs.img_projsKRT_x.restype = None
_libsbaprojs.img_projsKRT_x.__str__ = "MotionExpert_NoDistortion"
# Python handle
_MotionExpert_NoDistortion = _libsbaprojs.img_projsKRT_x
_MotionExpert_NoDistortion.__doc__="""
Function prototype for callable projection func in expert mode sba:
void (*func) (double *p, struct sba_crsm *idxij, int *rcidxs, 
int *rcsubs, double *hx, void *adata)

Does projection from P to X. 
"""

_libsbaprojs.img_projsKRT_jac_x.argtypes = [ctypes.POINTER(ctypes.c_double),
                                            ctypes.POINTER(Crsm),
                                            ctypes.POINTER(ctypes.c_int),
                                            ctypes.POINTER(ctypes.c_int),
                                            ctypes.POINTER(ctypes.c_double),
                                            ctypes.c_void_p]
_libsbaprojs.img_projsKRT_jac_x.restype = None
_libsbaprojs.img_projsKRT_jac_x.__str__ = "MotionJacobianExpert_NoDistortion"
# Python handle
_MotionJacobianExpert_NoDistortion = _libsbaprojs.img_projsKRT_jac_x
_MotionJacobianExpert_NoDistortion.__doc__="""
Function prototype for callable projection fjac in simple mode sba:
void (*fjac) (int j, int i, double *bi, double *xij, void *adata)

Does projection from P to X. 
"""

_libsbaprojs.img_projsKS_x.argtypes = [ctypes.POINTER(ctypes.c_double),
                                       ctypes.POINTER(Crsm),
                                       ctypes.POINTER(ctypes.c_int),
                                       ctypes.POINTER(ctypes.c_int),
                                       ctypes.POINTER(ctypes.c_double),
                                       ctypes.c_void_p]
_libsbaprojs.img_projsKS_x.restype = None
_libsbaprojs.img_projsKS_x.__str__ = "StructureExpert_NoDistortion"
# Python handle
_StructureExpert_NoDistortion = _libsbaprojs.img_projsKRT_x
_StructureExpert_NoDistortion.__doc__="""
Function prototype for callable projection func in expert mode sba:
void (*func) (double *p, struct sba_crsm *idxij, int *rcidxs, 
int *rcsubs, double *hx, void *adata)

Does projection from P to X. 
"""

_libsbaprojs.img_projsKS_jac_x.argtypes = [ctypes.POINTER(ctypes.c_double),
                                           ctypes.POINTER(Crsm),
                                           ctypes.POINTER(ctypes.c_int),
                                           ctypes.POINTER(ctypes.c_int),
                                           ctypes.POINTER(ctypes.c_double),
                                           ctypes.c_void_p]
_libsbaprojs.img_projsKS_jac_x.restype = None
_libsbaprojs.img_projsKS_jac_x.__str__ = "StructureJacobianExpert_NoDistortion"
# Python handle
_StructureJacobianExpert_NoDistortion = _libsbaprojs.img_projsKS_jac_x
_StructureJacobianExpert_NoDistortion.__doc__="""
Function prototype for callable projection fjac in simple mode sba:
void (*fjac) (int j, int i, double *bi, double *xij, void *adata)

Does projection from P to X. 
"""







# Simple mode
_libsbaprojs.img_projRTS.argtypes = [ctypes.c_int, # j
                                     ctypes.c_int, # i
                                     ctypes.POINTER(ctypes.c_double), # *bi
                                     ctypes.POINTER(ctypes.c_double), # *xij
                                     ctypes.c_void_p]
_libsbaprojs.img_projRTS.restype = None
_libsbaprojs.img_projRTS.__str__ = "MotionsStructureSimple_NoCameras"
# Python handle
_MotionStructureSimple_NoCameras = _libsbaprojs.img_projRTS
_MotionStructureSimple_NoCameras.__doc__="""
Function prototype for callable projection func in simple mode sba:
void (*proj) (int j, int i, double *bi, double *xij, void *adata)

Use to wrap a Python function when sending a projection function from 
Python into C.
"""
_libsbaprojs.img_projRT.argtypes = [ctypes.c_int, # j
                                    ctypes.c_int, # i
                                    ctypes.POINTER(ctypes.c_double), # *bi
                                    ctypes.POINTER(ctypes.c_double), # *xij
                                    ctypes.c_void_p]
_libsbaprojs.img_projRT.restype = None
_libsbaprojs.img_projRT.__str__ = "MotionSimple_NoCameras"
# Python handle
_MotionSimple_NoCameras = _libsbaprojs.img_projRT
_MotionSimple_NoCameras.__doc__="""
Function prototype for callable projection func in simple mode sba:
void (*proj) (int j, int i, double *bi, double *xij, void *adata)

Use to wrap a Python function when sending a projection function from 
Python into C.
"""
_libsbaprojs.img_projS.argtypes = [ctypes.c_int, # j
                                   ctypes.c_int, # i
                                   ctypes.POINTER(ctypes.c_double), # *bi
                                   ctypes.POINTER(ctypes.c_double), # *xij
                                   ctypes.c_void_p]
_libsbaprojs.img_projS.restype = None
_libsbaprojs.img_projS.__str__ = "StructureSimple_NoCameras"
# Python handle
_StructureSimple_NoCameras = _libsbaprojs.img_projS
_StructureSimple_NoCameras.__doc__="""
Function prototype for callable projection func in simple mode sba:
void (*proj) (int j, int i, double *bi, double *xij, void *adata)

Use to wrap a Python function when sending a projection function from 
Python into C.
"""

_libsbaprojs.img_projKDRTS.argtypes = [ctypes.c_int, # j
                                       ctypes.c_int, # i
                                       ctypes.POINTER(ctypes.c_double), # *bi
                                       ctypes.POINTER(ctypes.c_double), # *xij
                                       ctypes.c_void_p]
_libsbaprojs.img_projKDRTS.restype = None
_libsbaprojs.img_projKDRTS.__str__ = "MotionsStructureSimple_Cameras"
# Python handle
_MotionStructureSimple_Cameras = _libsbaprojs.img_projKDRTS
_MotionStructureSimple_Cameras.__doc__="""
Function prototype for callable projection func in simple mode sba:
void (*proj) (int j, int i, double *bi, double *xij, void *adata)

Use to wrap a Python function when sending a projection function from 
Python into C.
"""
_libsbaprojs.img_projKDRT.argtypes = [ctypes.c_int, # j
                                      ctypes.c_int, # i
                                      ctypes.POINTER(ctypes.c_double), # *bi
                                      ctypes.POINTER(ctypes.c_double), # *xij
                                      ctypes.c_void_p]
_libsbaprojs.img_projKDRT.restype = None
_libsbaprojs.img_projKDRT.__str__ = "MotionSimple_Cameras"
# Python handle
_MotionSimple_Cameras = _libsbaprojs.img_projKDRT
_MotionSimple_Cameras.__doc__="""
Function prototype for callable projection func in simple mode sba:
void (*proj) (int j, int i, double *bi, double *xij, void *adata)

Use to wrap a Python function when sending a projection function from 
Python into C.
"""
_libsbaprojs.img_projKDS.argtypes = [ctypes.c_int, # j
                                     ctypes.c_int, # i
                                     ctypes.POINTER(ctypes.c_double), # *bi
                                     ctypes.POINTER(ctypes.c_double), # *xij
                                     ctypes.c_void_p]
_libsbaprojs.img_projKDS.restype = None
_libsbaprojs.img_projKDS.__str__ = "StructureSimple_Cameras"
# Python handle
_StructureSimple_Cameras = _libsbaprojs.img_projKDS
_StructureSimple_Cameras.__doc__="""
Function prototype for callable projection func in simple mode sba:
void (*proj) (int j, int i, double *bi, double *xij, void *adata)

Use to wrap a Python function when sending a projection function from 
Python into C.
"""

_libsbaprojs.img_projKRTS.argtypes = [ctypes.c_int, # j
                                      ctypes.c_int, # i
                                      ctypes.POINTER(ctypes.c_double), # *bi
                                      ctypes.POINTER(ctypes.c_double), # *xij
                                      ctypes.c_void_p]
_libsbaprojs.img_projKRTS.restype = None
_libsbaprojs.img_projKRTS.__str__ = "MotionsStructureSimple_NoDistortion"
# Python handle
_MotionStructureSimple_NoDistortion = _libsbaprojs.img_projKRTS
_MotionStructureSimple_NoDistortion.__doc__="""
Function prototype for callable projection func in simple mode sba:
void (*proj) (int j, int i, double *bi, double *xij, void *adata)

Use to wrap a Python function when sending a projection function from 
Python into C.
"""
_libsbaprojs.img_projKRT.argtypes = [ctypes.c_int, # j
                                     ctypes.c_int, # i
                                     ctypes.POINTER(ctypes.c_double), # *bi
                                     ctypes.POINTER(ctypes.c_double), # *xij
                                     ctypes.c_void_p]
_libsbaprojs.img_projKRT.restype = None
_libsbaprojs.img_projKRT.__str__ = "MotionSimple_NoDistortion"
# Python handle
_MotionSimple_NoDistortion = _libsbaprojs.img_projKRT
_MotionSimple_NoDistortion.__doc__="""
Function prototype for callable projection func in simple mode sba:
void (*proj) (int j, int i, double *bi, double *xij, void *adata)

Use to wrap a Python function when sending a projection function from 
Python into C.
"""
_libsbaprojs.img_projKS.argtypes = [ctypes.c_int, # j
                                    ctypes.c_int, # i
                                    ctypes.POINTER(ctypes.c_double), # *bi
                                    ctypes.POINTER(ctypes.c_double), # *xij
                                    ctypes.c_void_p]
_libsbaprojs.img_projKS.restype = None
_libsbaprojs.img_projKS.__str__ = "StructureSimple_NoDistortion"
# Python handle
_StructureSimple_NoDistortion = _libsbaprojs.img_projKS
_StructureSimple_NoDistortion.__doc__="""
Function prototype for callable projection func in simple mode sba:
void (*proj) (int j, int i, double *bi, double *xij, void *adata)

Use to wrap a Python function when sending a projection function from 
Python into C.
"""

_libsbaprojs.img_projRTS.argtypes = [ctypes.c_int, # j
                                     ctypes.c_int, # i
                                     ctypes.POINTER(ctypes.c_double), # *bi
                                     ctypes.POINTER(ctypes.c_double), # *xij
                                     ctypes.c_void_p]
_libsbaprojs.img_projRTS.restype = None
_libsbaprojs.img_projRTS.__str__ = "MotionsStructureSimple_NoCameras"
# Python handle
_MotionStructureSimple_NoCameras = _libsbaprojs.img_projRTS
_MotionStructureSimple_NoCameras.__doc__="""
Function prototype for callable projection func in simple mode sba:
void (*proj) (int j, int i, double *bi, double *xij, void *adata)

Use to wrap a Python function when sending a projection function from 
Python into C.
"""
_libsbaprojs.img_projRT.argtypes = [ctypes.c_int, # j
                                    ctypes.c_int, # i
                                    ctypes.POINTER(ctypes.c_double), # *bi
                                    ctypes.POINTER(ctypes.c_double), # *xij
                                    ctypes.c_void_p]
_libsbaprojs.img_projRT.restype = None
_libsbaprojs.img_projRT.__str__ = "MotionSimple_NoCameras"
# Python handle
_MotionSimple_NoCameras = _libsbaprojs.img_projRT
_MotionSimple_NoCameras.__doc__="""
Function prototype for callable projection func in simple mode sba:
void (*proj) (int j, int i, double *bi, double *xij, void *adata)

Use to wrap a Python function when sending a projection function from 
Python into C.
"""
_libsbaprojs.img_projS.argtypes = [ctypes.c_int, # j
                                   ctypes.c_int, # i
                                   ctypes.POINTER(ctypes.c_double), # *bi
                                   ctypes.POINTER(ctypes.c_double), # *xij
                                   ctypes.c_void_p]
_libsbaprojs.img_projS.restype = None
_libsbaprojs.img_projS.__str__ = "StructureSimple_NoCameras"
# Python handle
_StructureSimple_NoCameras = _libsbaprojs.img_projS
_StructureSimple_NoCameras.__doc__="""
Function prototype for callable projection func in simple mode sba:
void (*proj) (int j, int i, double *bi, double *xij, void *adata)

Use to wrap a Python function when sending a projection function from 
Python into C.
"""

_libsbaprojs.img_projKDRTS_jac.argtypes = [ctypes.c_int, # j
                                           ctypes.c_int, # i
                                           ctypes.POINTER(ctypes.c_double), # *bi
                                           ctypes.POINTER(ctypes.c_double), # *xij
                                           ctypes.c_void_p]
_libsbaprojs.img_projKDRTS_jac.restype = None
_libsbaprojs.img_projKDRTS_jac.__str__ = "MotionsStructureJacobianSimple_Cameras"
# Python handle
_MotionStructureJacobianSimple_Cameras = _libsbaprojs.img_projKDRTS_jac
_MotionStructureJacobianSimple_Cameras.__doc__="""
Function prototype for callable projection fjac in simple mode sba:
void (*fjac) (int j, int i, double *bi, double *xij, void *adata)

Use to wrap a Python function when sending a Jacobian from 
Python into C.
"""
_libsbaprojs.img_projKDRT_jac.argtypes = [ctypes.c_int, # j
                                          ctypes.c_int, # i
                                          ctypes.POINTER(ctypes.c_double), # *bi
                                          ctypes.POINTER(ctypes.c_double), # *xij
                                          ctypes.c_void_p]
_libsbaprojs.img_projKDRT_jac.restype = None
_libsbaprojs.img_projKDRT_jac.__str__ = "MotionJacobianSimple_Cameras"
# Python handle
_MotionJacobianSimple_Cameras = _libsbaprojs.img_projKDRT_jac
_MotionJacobianSimple_Cameras.__doc__="""
Function prototype for callable projection fjac in simple mode sba:
void (*fjac) (int j, int i, double *bi, double *xij, void *adata)

Use to wrap a Python function when sending a projection function from 
Python into C.
"""
_libsbaprojs.img_projKDS_jac.argtypes = [ctypes.c_int, # j
                                         ctypes.c_int, # i
                                         ctypes.POINTER(ctypes.c_double), # *bi
                                         ctypes.POINTER(ctypes.c_double), # *xij
                                         ctypes.c_void_p]
_libsbaprojs.img_projKDS_jac.restype = None
_libsbaprojs.img_projKDS_jac.__str__ = "StructureJacobianSimple_Cameras"
# Python handle
_StructureJacobianSimple_Cameras = _libsbaprojs.img_projKDS_jac
_StructureJacobianSimple_Cameras.__doc__="""
Function prototype for callable projection fjac in simple mode sba:
void (*fjac) (int j, int i, double *bi, double *xij, void *adata)

Use to wrap a Python function when sending a projection function from 
Python into C.
"""

_libsbaprojs.img_projKRTS_jac.argtypes = [ctypes.c_int, # j
                                          ctypes.c_int, # i
                                          ctypes.POINTER(ctypes.c_double), # *bi
                                          ctypes.POINTER(ctypes.c_double), # *xij
                                          ctypes.c_void_p]
_libsbaprojs.img_projKRTS_jac.restype = None
_libsbaprojs.img_projKRTS_jac.__str__ = "MotionsStructureJacobianSimple_NoDistortion"
# Python handle
_MotionStructureJacobianSimple_NoDistortion = _libsbaprojs.img_projKRTS_jac
_MotionStructureJacobianSimple_NoDistortion.__doc__="""
Function prototype for callable projection fjac in simple mode sba:
void (*fjac) (int j, int i, double *bi, double *xij, void *adata)

Use to wrap a Python function when sending a projection function from 
Python into C.
"""
_libsbaprojs.img_projKRT_jac.argtypes = [ctypes.c_int, # j
                                         ctypes.c_int, # i
                                         ctypes.POINTER(ctypes.c_double), # *bi
                                         ctypes.POINTER(ctypes.c_double), # *xij
                                         ctypes.c_void_p]
_libsbaprojs.img_projKRT_jac.restype = None
_libsbaprojs.img_projKRT_jac.__str__ = "MotionJacobianSimple_NoDistortion"
# Python handle
_MotionJacobianSimple_NoDistortion = _libsbaprojs.img_projKRT_jac
_MotionJacobianSimple_NoDistortion.__doc__="""
Function prototype for callable projection fjac in simple mode sba:
void (*fjac) (int j, int i, double *bi, double *xij, void *adata)

Use to wrap a Python function when sending a projection function from 
Python into C.
"""
_libsbaprojs.img_projKS_jac.argtypes = [ctypes.c_int, # j
                                        ctypes.c_int, # i
                                        ctypes.POINTER(ctypes.c_double), # *bi
                                        ctypes.POINTER(ctypes.c_double), # *xij
                                        ctypes.c_void_p]
_libsbaprojs.img_projKS_jac.restype = None
_libsbaprojs.img_projKS_jac.__str__ = "StructureJacobianSimple_NoDistortion"
# Python handle
_StructureJacobianSimple_NoDistortion = _libsbaprojs.img_projKS_jac
_StructureJacobianSimple_NoDistortion.__doc__="""
Function prototype for callable projection fjac in simple mode sba:
void (*fjac) (int j, int i, double *bi, double *xij, void *adata)

Use to wrap a Python function when sending a projection function from 
Python into C.
"""






def whichProjection(options):
    """
    Dispatcher for providing the correct projection function from C

    input:
    options - Options object specifying sba options

    returns:
    a Python pointer to a C function (ctypes.funcptr) for sba projection.
    This can be passed to an sba driver, eg sba_levmar_motstr_x

    Raises SbaError if options are not recognized
    """
    #logging.debug("projections.whichProjection() called with {0}".format(options))
    if options.expert:
        if options.motstruct==sba.options.OPTS_MOTSTRUCT:
            if options.camera==sba.options.OPTS_CAMS:
                #return _libsbaprojs.img_projsKDRTS_x
                return _MotionStructureExpert_Cameras
            elif options.camera==sba.options.OPTS_CAMS_NODIST:
                #return _libsbaprojs.img_projsKRTS_x
                return _MotionStructureExpert_NoDistortion
            elif options.camera==sba.options.OPTS_NO_CAMS:
                return _MotionStructureExpert_NoCameras
            else:
                raise SbaError("Unrecognized projection option")
        elif options.motstruct==sba.options.OPTS_MOT:
            if options.camera==sba.options.OPTS_CAMS:
                return _MotionExpert_Cameras
            elif options.camera==sba.options.OPTS_CAMS_NODIST:
                return _MotionExpert_NoDistortion
            elif options.camera==sba.options.OPTS_NO_CAMS:
                return _MotionExpert_NoCameras
            else:
                raise SbaError("Unrecognized projection option")
        elif options.motstruct==sba.options.OPTS_STRUCT:
            if options.camera==sba.options.OPTS_CAMS:
                return _StructureExpert_Cameras
            elif options.camera==sba.options.OPTS_CAMS_NODIST:
                return _StructureExpert_NoDistortion
            elif options.camera==sba.options.OPTS_NO_CAMS:
                return _StructureExpert_NoCameras
            else:
                raise SbaError("Unrecognized projection option")
    else: # simple mode
        if options.motstruct==sba.options.OPTS_MOTSTRUCT:
            if options.camera==sba.options.OPTS_CAMS:
                return _MotionStructureSimple_Cameras
            elif options.camera==sba.options.OPTS_CAMS_NODIST:
                return _MotionStructureSimple_NoDistortion
            elif options.camera==sba.options.OPTS_NO_CAMS:
                return _MotionStructureSimple_NoCameras
            else:
                raise SbaError("Unrecognized projection option")
        elif options.motstruct==sba.options.OPTS_MOT:
            if options.camera==sba.options.OPTS_CAMS:
                return _MotionSimple_Cameras
            elif options.camera==sba.options.OPTS_CAMS_NODIST:
                return _MotionSimple_NoDistortion
            elif options.camera==sba.options.OPTS_NO_CAMS:
                return _MotionSimple_NoCameras
            else:
                raise SbaError("Unrecognized projection option")
        elif options.motstruct==sba.options.OPTS_STRUCT:
            if options.camera==sba.options.OPTS_CAMS:
                return _StructureSimple_Cameras
            elif options.camera==sba.options.OPTS_CAMS_NODIST:
                return _StructureSimple_NoDistortion
            elif options.camera==sba.options.OPTS_NO_CAMS:
                return _StructureSimple_NoCameras
            else:
                raise SbaError("Unrecognized projection option")

def whichJacobian(options):
    """
    Dispatcher for providing the correct Jacobian function

    input:
    options - Options object specifying sba options

    returns:
    a Python pointer to a C function (ctypes.funcptr) for sba Jacobian 
    This can be passed to an sba driver, eg sba_levmar_motstr_x

    Raises SbaError if options are not recognized
    """
    #logging.debug("projections.whichJacobian() called with {0}".format(options))
    if options.expert:
        if options.analjac:
            if options.motstruct==sba.options.OPTS_MOTSTRUCT:
                if options.camera==sba.options.OPTS_CAMS:
                    return _MotionStructureJacobianExpert_Cameras
                elif options.camera==sba.options.OPTS_CAMS_NODIST:
                    return _MotionStructureJacobianExpert_NoDistortion
                elif options.camera==sba.options.OPTS_NO_CAMS:
                    return _MotionStructureJacobianExpert_NoCameras
                else:
                    raise SbaError("Unrecognized Jacobian option")
            elif options.motstruct==sba.options.OPTS_MOT:
                if options.camera==sba.options.OPTS_CAMS:
                    return _MotionJacobianExpert_Cameras
                elif options.camera==sba.options.OPTS_CAMS_NODIST:
                    return _MotionJacobianExpert_NoDistortion
                elif options.camera==sba.options.OPTS_NO_CAMS:
                    return _MotionJacobianExpert_NoCameras
                else:
                    raise SbaError("Unrecognized Jacobian option")
            elif options.motstruct==sba.options.OPTS_STRUCT:
                if options.camera==sba.options.OPTS_CAMS:
                    return _StructureJacobianExpert_Cameras
                elif options.camera==sba.options.OPTS_CAMS_NODIST:
                    return _StructureJacobianExpert_NoDistortion
                elif options.camera==sba.options.OPTS_NO_CAMS:
                    return _StructureJacobianExpert_NoCameras
                else:
                    raise SbaError("Unrecognized Jacobian option")
        else:
            return None
    else: # simple mode
        if options.analjac:
            if options.motstruct==sba.options.OPTS_MOTSTRUCT:
                if options.camera==sba.options.OPTS_CAMS:
                    return _MotionStructureJacobianSimple_Cameras
                elif options.camera==sba.options.OPTS_CAMS_NODIST:
                    return _MotionStructureJacobianSimple_NoDistortion
                elif options.camera==sba.options.OPTS_NO_CAMS:
                    return _MotionStructureJacobianSimple_NoCameras
                else:
                    raise SbaError("Unrecognized Jacobian option")
            elif options.motstruct==sba.options.OPTS_MOT:
                if options.camera==sba.options.OPTS_CAMS:
                    return _MotionJacobianSimple_Cameras
                elif options.camera==sba.options.OPTS_CAMS_NODIST:
                    return _MotionJacobianSimple_NoDistortion
                elif options.camera==sba.options.OPTS_NO_CAMS:
                    return _MotionJacobianSimple_NoCameras
                else:
                    raise SbaError("Unrecognized Jacobian option")
            elif options.motstruct==sba.options.OPTS_STRUCT:
                if options.camera==sba.options.OPTS_CAMS:
                    return _StructureJacobianSimple_Cameras
                elif options.camera==sba.options.OPTS_CAMS_NODIST:
                    return _StructureJacobianSimple_NoDistortion
                elif options.camera==sba.options.OPTS_NO_CAMS:
                    return _StructureJacobianSimple_NoCameras
                else:
                    raise SbaError("Unrecognized Jacobian option")
            else:
                return None




