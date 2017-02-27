#!/usr/bin/env python
"""
crsm.py
Python wrapper for compressed row storage in Lourakis sba code
"""

import ctypes

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

