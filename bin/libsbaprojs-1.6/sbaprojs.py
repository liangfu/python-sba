#!/usr/bin/env python
"""
sbaprojs.py
Python wrapper for projections and jacobians from eucsbademo by Lourakis
Wrapper by Dennis Evangelista, 2013
"""

__version__ = "1.6.0"

import ctypes
import logging
import numpy as np

# This line actually loads the sba shared object
_libsbaprojs = ctypes.CDLL("libsbaprojs.so")
# hmmm where will this file go when it is actually installed... 

# what a horrible way to pass parameters to the projection routines?!
# define this in sba.py? 
class Globs_(ctypes.Structure):
    _fields_ = [("rot0params", ctypes.POINTER(ctypes.c_double)),
                ("intrcalib",ctypes.POINTER(ctypes.c_double)),
                ("nccalib",ctypes.c_int),
                ("ncdist",ctypes.c_int),
                ("cnp",ctypes.c_int),
                ("pnp",ctypes.c_int),
                ("mnp",ctypes.c_int),
                ("ptparams",ctypes.POINTER(ctypes.c_double)),
                ("camparams",ctypes.POINTER(ctypes.c_double))]

# This should be defined in sba.py?                 
class Crsm(ctypes.Structure):
    """Sparse matrix represntation using compressed row storage
    Python wrapper for libsba struct sba_crsm.
    nr, nc - rows, columns for the sparse matrix
    nnz - number of nonzero array elements
    val - storage for nonzero array elements, size nnz
    colidx - column indices for nonzero elements, size nnz
    rowptr - locations in val that start a row, size nr+1"""
    _fields_ = [("nr", ctypes.c_int),
                ("nc", ctypes.c_int),
                ("nnz", ctypes.c_int),
                ("val", ctypes.POINTER(ctypes.c_int)),
                ("colidx", ctypes.POINTER(ctypes.c_int)),
                ("rowptr", ctypes.POINTER(ctypes.c_int))]



# This is used for going from Python into C
FuncProjection = ctypes.CFUNCTYPE(ctypes.POINTER(ctypes.c_double), #p
                                  ctypes.POINTER(Crsm), #idxij
                                  ctypes.POINTER(ctypes.c_int), #rcidxs
                                  ctypes.POINTER(ctypes.c_int), #rcsubs
                                  ctypes.POINTER(ctypes.c_double), #hx
                                  ctypes.c_void_p) #adata
"""Function prototype for callable projection func in expert mode sba"""

# This is used for going from Python into C
FjacProjJacobian = ctypes.CFUNCTYPE(ctypes.POINTER(ctypes.c_double), #p
                                    ctypes.POINTER(Crsm), #idxij
                                    ctypes.POINTER(ctypes.c_int), #rcidxs
                                    ctypes.POINTER(ctypes.c_int), #rcsubs
                                    ctypes.POINTER(ctypes.c_double), #jac
                                    ctypes.c_void_p) #adata
"""Function prototype for callable Jacobian fjac in expert mode sba"""







# projection functions and their Jacobians
# full sparse bundle adjustment
_libsbaprojs.img_projsRTS_x.argtypes = [ctypes.POINTER(ctypes.c_double),
                                        ctypes.POINTER(Crsm),
                                        ctypes.POINTER(ctypes.c_int),
                                        ctypes.POINTER(ctypes.c_int),
                                        ctypes.POINTER(ctypes.c_double),
                                        ctypes.c_void_p]
_libsbaprojs.img_projsRTS_x.restype = None
# provide Python-friendly handle for using this
fullProjection = _libsbaprojs.img_projsRTS_x


_libsbaprojs.img_projsRTS_jac_x.argtypes = [ctypes.POINTER(ctypes.c_double),
                                            ctypes.POINTER(Crsm),
                                            ctypes.POINTER(ctypes.c_int),
                                            ctypes.POINTER(ctypes.c_int),
                                            ctypes.POINTER(ctypes.c_double),
                                            ctypes.c_void_p]
_libsbaprojs.img_projsRTS_jac_x.restype = None
# provide Python-friendly handle for using this
fullJacobian = _libsbaprojs.img_projsRTS_jac_x

# sparse bundle adjustment with camera parameters only
_libsbaprojs.img_projsRT_x.argtypes = [ctypes.POINTER(ctypes.c_double),
                                       ctypes.POINTER(Crsm),
                                       ctypes.POINTER(ctypes.c_int),
                                       ctypes.POINTER(ctypes.c_int),
                                       ctypes.POINTER(ctypes.c_double),
                                       ctypes.c_void_p]
_libsbaprojs.img_projsRT_x.restype = None
# provide Python-friendly handle for using this
justCamProjection = _libsbaprojs.img_projsRT_x

_libsbaprojs.img_projsRT_jac_x.argtypes = [ctypes.POINTER(ctypes.c_double),
                                           ctypes.POINTER(Crsm),
                                           ctypes.POINTER(ctypes.c_int),
                                           ctypes.POINTER(ctypes.c_int),
                                           ctypes.POINTER(ctypes.c_double),
                                           ctypes.c_void_p]
_libsbaprojs.img_projsRT_jac_x.restype = None
# provide Python-friendly handle for using this
justCamJacobian = _libsbaprojs.img_projsRT_jac_x

# sparse bundle adjustment with structure parameters only
_libsbaprojs.img_projsS_x.argtypes = [ctypes.POINTER(ctypes.c_double),
                                      ctypes.POINTER(Crsm),
                                      ctypes.POINTER(ctypes.c_int),
                                      ctypes.POINTER(ctypes.c_int),
                                      ctypes.POINTER(ctypes.c_double),
                                      ctypes.c_void_p]
_libsbaprojs.img_projsS_x.restype = None
# provide Python-friendly handle for using this
justStructureProjection = _libsbaprojs.img_projsS_x

_libsbaprojs.img_projsS_jac_x.argtypes = [ctypes.POINTER(ctypes.c_double),
                                          ctypes.POINTER(Crsm),
                                          ctypes.POINTER(ctypes.c_int),
                                          ctypes.POINTER(ctypes.c_int),
                                          ctypes.POINTER(ctypes.c_double),
                                          ctypes.c_void_p]
_libsbaprojs.img_projsS_jac_x.restype = None
# provide Python-friendly handle for using this
justStructureProjection = _libsbaprojs.img_projsS_jac_x








# Cameras with distortion
_libsbaprojs.img_projsKDRTS_x.argtypes = [ctypes.POINTER(ctypes.c_double),
                                          ctypes.POINTER(Crsm),
                                          ctypes.POINTER(ctypes.c_int),
                                          ctypes.POINTER(ctypes.c_int),
                                          ctypes.POINTER(ctypes.c_double),
                                          ctypes.c_void_p]
_libsbaprojs.img_projsKDRTS_x.restype = None
# provide Python-friendly handle for using this
fullCameraProjection = _libsbaprojs.img_projsKDRTS_x

_libsbaprojs.img_projsKDRTS_jac_x.argtypes = [ctypes.POINTER(ctypes.c_double),
                                              ctypes.POINTER(Crsm),
                                              ctypes.POINTER(ctypes.c_int),
                                              ctypes.POINTER(ctypes.c_int),
                                              ctypes.POINTER(ctypes.c_double),
                                              ctypes.c_void_p]
_libsbaprojs.img_projsKDRTS_jac_x.restype = None
# provide Python-friendly handle for using this
fullCameraJacobian = _libsbaprojs.img_projsKDRTS_jac_x

_libsbaprojs.img_projsKDRT_x.argtypes = [ctypes.POINTER(ctypes.c_double),
                                         ctypes.POINTER(Crsm),
                                         ctypes.POINTER(ctypes.c_int),
                                         ctypes.POINTER(ctypes.c_int),
                                         ctypes.POINTER(ctypes.c_double),
                                         ctypes.c_void_p]
_libsbaprojs.img_projsKDRT_x.restype = None
# provide Python-friendly handle for using this
justStructureCameraProjection = _libsbaprojs.img_projsKDRT_x

_libsbaprojs.img_projsKDRT_jac_x.argtypes = [ctypes.POINTER(ctypes.c_double),
                                             ctypes.POINTER(Crsm),
                                             ctypes.POINTER(ctypes.c_int),
                                             ctypes.POINTER(ctypes.c_int),
                                             ctypes.POINTER(ctypes.c_double),
                                             ctypes.c_void_p]
_libsbaprojs.img_projsKDRT_jac_x.restype = None
# provide Python-friendly handle for using this
justStructureCameraJacobian = _libsbaprojs.img_projsKDRT_jac_x

_libsbaprojs.img_projsKDS_x.argtypes = [ctypes.POINTER(ctypes.c_double),
                                        ctypes.POINTER(Crsm),
                                        ctypes.POINTER(ctypes.c_int),
                                        ctypes.POINTER(ctypes.c_int),
                                        ctypes.POINTER(ctypes.c_double),
                                        ctypes.c_void_p]
_libsbaprojs.img_projsKDS_x.restype = None
# provide Python-friendly handle for using this
justMotionCameraProjection = _libsbaprojs.img_projsKDS_x

_libsbaprojs.img_projsKDS_jac_x.argtypes = [ctypes.POINTER(ctypes.c_double),
                                            ctypes.POINTER(Crsm),
                                            ctypes.POINTER(ctypes.c_int),
                                            ctypes.POINTER(ctypes.c_int),
                                            ctypes.POINTER(ctypes.c_double),
                                            ctypes.c_void_p]
_libsbaprojs.img_projsKDS_jac_x.restype = None
# provide Python-friendly handle for using this
justMotionCameraJacobian = _libsbaprojs.img_projsKDS_jac_x






# Ideal cameras with now distortion
_libsbaprojs.img_projsKRTS_x.argtypes = [ctypes.POINTER(ctypes.c_double),
                                         ctypes.POINTER(Crsm),
                                         ctypes.POINTER(ctypes.c_int),
                                         ctypes.POINTER(ctypes.c_int),
                                         ctypes.POINTER(ctypes.c_double),
                                         ctypes.c_void_p]
_libsbaprojs.img_projsKRTS_x.restype = None
# provide Python-friendly handle for using this
fullNoDistortionProjection = _libsbaprojs.img_projsKRTS_x

_libsbaprojs.img_projsKRTS_jac_x.argtypes = [ctypes.POINTER(ctypes.c_double),
                                             ctypes.POINTER(Crsm),
                                             ctypes.POINTER(ctypes.c_int),
                                             ctypes.POINTER(ctypes.c_int),
                                             ctypes.POINTER(ctypes.c_double),
                                             ctypes.c_void_p]
_libsbaprojs.img_projsKRTS_jac_x.restype = None
# provide Python-friendly handle for using this
fullNoDistortionJacobian = _libsbaprojs.img_projsKRTS_jac_x

_libsbaprojs.img_projsKRT_x.argtypes = [ctypes.POINTER(ctypes.c_double),
                                        ctypes.POINTER(Crsm),
                                        ctypes.POINTER(ctypes.c_int),
                                        ctypes.POINTER(ctypes.c_int),
                                        ctypes.POINTER(ctypes.c_double),
                                        ctypes.c_void_p]
_libsbaprojs.img_projsKRT_x.restype = None
# provide Python-friendly handle for using this
justStructureNoDistortionProjection = _libsbaprojs.img_projsKRT_x

_libsbaprojs.img_projsKRT_jac_x.argtypes = [ctypes.POINTER(ctypes.c_double),
                                            ctypes.POINTER(Crsm),
                                            ctypes.POINTER(ctypes.c_int),
                                            ctypes.POINTER(ctypes.c_int),
                                            ctypes.POINTER(ctypes.c_double),
                                            ctypes.c_void_p]
_libsbaprojs.img_projsKRT_jac_x.restype = None
# provide Python-friendly handle for using this
justStructureNoDistortionJacobian = _libsbaprojs.img_projsKRT_jac_x

_libsbaprojs.img_projsKS_x.argtypes = [ctypes.POINTER(ctypes.c_double),
                                       ctypes.POINTER(Crsm),
                                       ctypes.POINTER(ctypes.c_int),
                                       ctypes.POINTER(ctypes.c_int),
                                       ctypes.POINTER(ctypes.c_double),
                                       ctypes.c_void_p]
_libsbaprojs.img_projsKS_x.restype = None
# provide Python-friendly handle for using this
justMotionNoDistortionProjection = _libsbaprojs.img_projsKS_x

_libsbaprojs.img_projsKS_jac_x.argtypes = [ctypes.POINTER(ctypes.c_double),
                                           ctypes.POINTER(Crsm),
                                           ctypes.POINTER(ctypes.c_int),
                                           ctypes.POINTER(ctypes.c_int),
                                           ctypes.POINTER(ctypes.c_double),
                                           ctypes.c_void_p]
_libsbaprojs.img_projsKS_jac_x.restype = None
# provide Python-friendly handle for using this
justMotionNoDistortionJacobian = _libsbaprojs.img_projsKS_jac_x

