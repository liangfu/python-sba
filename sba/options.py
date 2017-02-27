#!/usr/bin/env python
"""
options.py
Dennis Evangelista, 2013

When C calling routines ask for it, Options is used to help figure out:

1. Which projection and jacobian to call?
       use magic decoder function in sba.projections

2. What options array to pass
       use Options.optsArray() , which calls Options.opts.toC()

3. What maxiter, verbose to pass, converted to C...

4. What Globs to pass, converted to a C structure
"""

import ctypes
import numpy as np
import logging







# These are stopping criteria for Levenberg-Marquardt algorithm
OPTS_SIZE = 5 # Size of LM options array of double
_INIT_MU = 1.0E-3 # mu, LM parameter
_STOP_THRESH = 1.0E-12 # stopping threshold
_MAXITER = 1000 # max iterations
_MAXITER2 = 1500 # unused? 

COptsType = ctypes.c_double*OPTS_SIZE

class LMOptions(object):
    """
    Python friendly implementation of sba Levenberg-Marquardt options

    attributes:
    mu - default 1e-3, LM parameter
    J^Teps_inf - default 1e-12, LM parameter
    DP - default 1e-12, LM parameter
    eps - default 1e-12, LM parameter
    eps/stuff - default 0.0, LM parameter
    """
    def __init__(self):
        self.__mu = _INIT_MU # LM parameter
        self.__JTeps = _STOP_THRESH # LM parameter
        self.__DP = _STOP_THRESH # LM parameter
        self.__eps = _STOP_THRESH # LM parameter
        self.__epsrr = 0.0 # LM parameter
        self.__optsarray = np.array([self.__mu, self.__JTeps, self.__DP,
                                     self.__eps, self.__epsrr],dtype=np.double)
        
        super(LMOptions, self).__init__()

    def _getmu(self):
        return self.__mu
    def _setmu(self, value):
        self.__mu = value
        self.__optsarray[0] = value
    mu = property(_getmu, _setmu,doc="mu parameter for Levenberg-Marquardt")

    def _getJTeps(self):
        return self.__JTeps
    def _setJTeps(self, value):
        self.__JTeps = value
        self.__optsarray[1] = value
    JTeps = property(_getJTeps, _setJTeps,doc="J^Teps_inf parameter for LM")

    def _getDP(self):
        return self.__DP
    def _setDP(self, value):
        self.__DP = value
        self.__optsarray[2] = value
    DP = property(_getDP, _setDP,doc="DP parameter for Levenberg-Marquardt")

    def _geteps(self):
        return self.__eps
    def _seteps(self, value):
        self.__eps = value
        self.__optsarray[3] = value
    eps = property(_geteps, _seteps,doc="epsilon stopping param for LM")

    def _getepsrr(self):
        return self.__epsrr
    def _setepsrr(self, value):
        self.__epsrr = value
        self.__optsarray[4] = value
    epsrr = property(_getepsrr, _setepsrr,doc="eps/stuff stopping param for LM")

    def _getoptsarray(self):
        return self.__optsarray
    def _setoptsarray(self,newoptsarray):
        self.__optsarray = np.array(newoptsarray,dtype=np.double)
        self.__mu = self.__optsarray[0]
        self.__JTeps = self.__optsarray[1]
        self.__DP = self.__optsarray[2]
        self.__eps = self.__optsarray[3]
        self.__epsrr = self.__optsarray[4]
    optsarray = property(_getoptsarray, _setoptsarray,doc="5 element options array for LM")

    def toC(self):
        """Returns a C-passable (double ``*``) pointer to the options array"""
        return COptsType(self.__optsarray[0],self.__optsarray[1],
                         self.__optsarray[2],self.__optsarray[3],
                         self.__optsarray[4])


    









class Globs(object):
    """
    Globs for passing additional options to projection and Jacobian 
    functions in sba.

    attributes:
    rot0params, default None, specifies initial rotation parameters?
    intrcalib, default None, specifies intrinsics if all cams held identical
    nccalib, default OPTS_FIX4_INTR, how many intrinsics to hold fixed
    ncdist, default OPTS_FIX3_DIST, how many distortion coeffs to hold fixed
    ptparams, default None, used when adjusting for cam params only
    camparams, default None, used when adjusting for structure only

    dimensions:
    cnp, default 3+3+5+5, number of parameters for each camera
    pnp, default 2, number of parameters for each image point (uv)
    mnp, default 3, number of parameters for each world point (xyz)
    """
    def __init__(self):
        """Create Globs object with default parameters"""
        self.__rot0params = None # initial rotation parameters 
        self.intrcalib = None # calibration when all cams fixed and identical
        self.nccalib = OPTS_FIX4_INTR # hold this many intrinsics fixed
        self.ncdist = OPTS_FIX3_DIST # hold this many dist coeffs fixed
        self.ptparams = None # needed when adjusting for cam params only
        self.camparams = None # needed when adjusting for structure only

        self._cnp = 5+5+3+3 # must set this, number of parameters per camera
        self.pnp = 2 # numbers for an image point
        self.mnp = 3 # numbers for a world point

        super(Globs,self).__init__()

    def _getrot0params(self):
        return self.__rot0params
    def _setrot0params(self,value):
        if value is None:
            self.__rot0params = None
            self.__Crot0params = None
        else:
            self.__rot0params = np.array(value,dtype=np.double)
            self.__Crot0params = self.__rot0params.ravel()
    rot0params = property(_getrot0params,_setrot0params,doc="initial rotation parameters")

    @property
    def Crot0params(self):
        return self.__Crot0params.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        
    def _getcnp(self):
        return self._cnp
    def _setcnp(self,val):
        self._cnp = val
        if self._cnp != 16:
            self.ncdist = -9999; 
    cnp = property(_getcnp,_setcnp,doc="Number of parameters per camera")
    
    @property 
    def hasFixedCal(self):
        return self.nccalib == OPTS_FIX5_INTR 

    @property
    def hasDistortion(self):
        return self.cnp == 5+5+3+3

    def toC(self):
        print self.ncdist
        """Returns Globs as a C Structure globs``_`` for passing to sba C routines"""
        if type(self.intrcalib) is np.ndarray:
            Cintrcalib = self.intrcalib.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        else:
            Cintrcalib = None
        if type(self.ptparams) is np.ndarray:
            Cptparams = self.ptparams.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        else:
            Cptparams = None
        if type(self.camparams) is np.ndarray:
            Ccamparams = self.camparams.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        else:
            Ccamparams = None
        return ctypes.pointer(_Globs_(self.Crot0params, # * double
                                      Cintrcalib, # * double
                                      ctypes.c_int(self.nccalib), # int
                                      ctypes.c_int(self.ncdist),  # int
                                      ctypes.c_int(self.cnp), # int
                                      ctypes.c_int(self.pnp), # int
                                      ctypes.c_int(self.mnp), # int
                                      Cptparams,  # * double
                                      Ccamparams)) # * double





class _Globs_(ctypes.Structure):
    """ctypes Structure for passing Globs to sba only"""
    _fields_ = [("rot0params", ctypes.POINTER(ctypes.c_double)),
                ("intrcalib",ctypes.POINTER(ctypes.c_double)),
                ("nccalib",ctypes.c_int),
                ("ncdist",ctypes.c_int),
                ("cnp",ctypes.c_int),
                ("pnp",ctypes.c_int),
                ("mnp",ctypes.c_int),
                ("ptparams",ctypes.POINTER(ctypes.c_double)),
                ("camparams",ctypes.POINTER(ctypes.c_double))]

        














# SBA options

# 1. Use expert mode (default) or simple?
# 2. Allow sba to alter camera?
OPTS_CAMS = 0 # camera intrinsics and distortion (default)
OPTS_CAMS_NODIST = 1 # camera intrinsics without distortion
OPTS_NO_CAMS = 2 # don't optimize cameras

# 3. Allow sba to alter motion, structure, or both
OPTS_MOTSTRUCT = 0 # motion and structure, default
OPTS_MOT = 1 # just motion
OPTS_STRUCT = 2 # just structure

# 4. Use analytic Jacobian (default) or not? 

# 5. Verbose or not? 
# What to print? 
OPTS_PRINT_NONE = -1
OPTS_PRINT_MOTSTRUCT = 0
OPTS_PRINT_MOT = 1
OPTS_PRINT_STRUCT = 2

# 6. Which camera intrinsics to hold fixed?
OPTS_FIX0_INTR = 0 # all free 
OPTS_FIX1_INTR = 1 # skew is fixed
OPTS_FIX2_INTR = 2 # skew and aspect ratio are fixed
OPTS_FIX3_INTR = 3 # not meaningful
OPTS_FIX4_INTR = 4 # skew, AR, and principal point fixed
OPTS_FIX5_INTR = 5 # all fixed
OPTS_FREE_INTR = OPTS_FIX0_INTR # synonym for all free
OPTS_FIX_SKEW = OPTS_FIX1_INTR # synonym 
OPTS_FIX_SKEW_AR = OPTS_FIX2_INTR # synonym
OPTS_FIX_SKEW_AR_PP = OPTS_FIX4_INTR # synonym
OPTS_FIXALL_INTR = OPTS_FIX5_INTR # synonym

# 7. Which distortion coeffs to hold fixed? 
OPTS_FIX0_DIST = 0 # all free
OPTS_FIX1_DIST = 1 # R6 fixed
OPTS_FIX2_DIST = 2 # R6, T2 fixed
OPTS_FIX3_DIST = 3 # R6, T2, T1 fixed
OPTS_FIX4_DIST = 4 # R6, T2, T1, R4 fixed
OPTS_FIX5_DIST = 5 # all fixed
OPTS_FREE_DIST = OPTS_FIX0_DIST # synonym for all free
OPTS_FIX_R6 = OPTS_FIX1_DIST # synonym
OPTS_FIX_TANG = OPTS_FIX3_DIST # synonym
OPTS_FIXALL_DIST = OPTS_FIX5_DIST #synonym

class Options(LMOptions,Globs):
    """
    Sparse bundle adjustment options settings

    attributes:
    expert - Use expert mode (default True)
    camera - Alter intrinsics and distortion? (default OPTS_CAMS = both)
    motstruct - Do sba for motions and structure? (default OPTS_MOTSTRUCT = both)
    analjac - Use analytic Jacobian (default True)
    verbose - (default True)
    printwhat - (default OPTS_PRINT_MOTSTRUCT = both)
    maxiter - Max iterations (default 100)

    Also a subclass of LMOptions, which has the Levenberg-Marquardt options
    and Globs, which has additional options passed to the projection and
    Jacobian functions. See those for additional options provided:

    attributes inherited from Globs:
    rot0params, default None, specifies initial rotation parameters?
    intrcalib, default None, specifies intrinsics if all cams held identical
    nccalib, default OPTS_FIX4_INTR, how many intrinsics to hold fixed
    ncdist, default OPTS_FIX3_DIST, how many distortion coeffs to hold fixed
    ptparams, default None, used when adjusting for cam params only
    camparams, default None, used when adjusting for structure only
    cnp, default 3+3+5+5, number of parameters for each camera
    pnp, default 2, number of parameters for each image point (uv)
    mnp, default 3, number of parameters for each world point (xyz)
 
    attributes inherited from LMOptions:
    mu - default 1e-3, LM parameter
    J^Teps_inf - default 1e-12, LM parameter
    DP - default 1e-12, LM parameter
    eps - default 1e-12, LM parameter
    eps/stuff - default 0.0, LM parameter

    pro tip: the base C routines seem to get confused if you call 
    OPTS_CAMS_NODIST plus OPTS_FIX5_INTR or OPTS_FIX5_DIST... (for example)
    For our work at UNC we usually supply 5+5+3+3 camera parameters,
    if you're doing that, use OPTS_CAMS (the cameras do have dist coeffs),
    so that it knows how to properly interpret what you want fixed. 
    """
    def __init__(self): 
        self.expert = True # Use expert mode
        self.camera = OPTS_CAMS # sba can alter camera intrinsics+dist
        self.motstruct = OPTS_MOTSTRUCT # sba for motion and structure
        self.analjac = True # use analytic jacobian
        self.verbose = True # verbose output
        self.printwhat = OPTS_PRINT_MOTSTRUCT
        self.maxiter = _MAXITER2 # max iterations 100
        super(Options, self).__init__()
        #self.opts = LMOptions() # numerical options for Levenberg-Marquardt
        #self.globs = Globs()

    def maxIterToC(self):
        return ctypes.c_int(self.maxiter)
    def verboseToC(self):
        if self.verbose:
            return ctypes.c_int(1)
        else:
            return ctypes.c_int(0) 
    def globsToC(self):
        return Globs.toC(self) # Why can't I do this using super? 
    def optsToC(self):
        return LMOptions.toC(self) # Why can't I do this using super? 
    def toC(self):
        return (self.maxIterToC(),self.verboseToC(),self.globsToC(),
                self.optsToC())

    def __str__(self):
        """Human readable representation of sba.Options()"""
        return "<sba.Options() object, expert {0}, motstruct {1}>".format(self.expert,self.motstruct)

    @classmethod
    def fromInput(cls,cameras,points):
        """
        Alternate constructor (a class method) to 
        create default Options object with values compatible
        with input objects cameras and points
        """
        newoptions = cls()
        newoptions.ncameras = cameras.ncameras
        if cameras.cnp==5+5+3+3:
            newoptions.camera = OPTS_CAMS
        elif cameras.cnp==5+3+3:
            newoptions.camera = OPTS_CAMS_NODIST
        elif cameras.cnp==3+3:
            newoptions.camera = OPTS_NOCAMS
        newoptions.cnp = cameras.cnp
        newoptions.rot0params = cameras.rot0params
        newoptions.mnp = points.mnp
        newoptions.pnp = points.pnp
        return newoptions
    # fromInput = classmethod(_fromInput)

    # pretty print method? 












