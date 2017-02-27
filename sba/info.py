#!/usr/bin/env python
"""
info.py
Dennis Evangelista, 2013
Python implementation of Lourakis' sba code and parts of demo
This is for handling output info 
"""

import ctypes
import numpy as np
import sba.options

# Provides readable reason for algorithm termination
INFO_REASONS = {1:"stopped by small gradient J^T e", 
                2:"stopped by small dp",
                3:"stopped by itmax",
                4:"stopped by small relative reduction in e_2",
                5:"too many attempts to increase damping, restart with increased mu",
                6:"stopped by small e_2",
                7:"error, stopped by invalid NaN or Inf func values",
                0:"no reason, sometimes it just happens"}

INFO_SIZE = 10 # Size of info array double [INFO SIZE]

CInfoType = ctypes.c_double*INFO_SIZE

class Information(object):
    """
    Python friendly class for representing information output by sba

    attributes:
    eps0 - LM parameter, initial value
    epsF - LM parameter, final value
    JTeps - J^Teps_inf LM parameter
    DP - LM parameter
    mu_max - LM parameter
    iterations - number of iterations
    reason - for termination. see why() method 
    fevals - number of projection function evaluations
    jevals - number of Jacobian evaluations
    nsystems - number of linear systems solved
    
    infoarray - access to 10 element numpy array 

    additional attributes (copied in during SparseBundleAdjust)
    nprojs - number of projections (3D point mapped to a 2D point in 1 image)
    nvars - number of free variables adjusted
    method - method used
    analjac - analytic (True) or approximate (False)
    havecovx - with (True) or without (False)
    havedistortion, ncdist
    fixedcal, nccalib
    t0 - start time of last sparse bundle adjustment
    tf - end time of last sparse bundle adjustment
    """
    def __init__(self,infoarray=None):
        if np.all(infoarray): # constructor called with input
            self.__eps0 = infoarray[0]
            self.__epsF = infoarray[1]
            self.__JTeps = infoarray[2]
            self.__DP = infoarray[3]
            self.__mu_max = infoarray[4]
            self.__iterations = infoarray[5]
            self.__reason = infoarray[6]
            self.__fevals = infoarray[7]
            self.__jevals = infoarray[8]
            self.__nsystems = infoarray[9]
            self.__infoarray = np.array(infoarray,dtype=np.double)
        else: # constructor called with no input
            self.__eps0 = 0.
            self.__epsF = 0.
            self.__JTeps = 0.
            self.__DP = 0.
            self.__mu_max = 0.
            self.__iterations = 0
            self.__reason = 0
            self.__fevals = 0
            self.__jevals = 0
            self.__nsystems = 0
            self.__infoarray = np.array([self.__eps0,self.__epsF,self.__JTeps,
                                         self.__DP,self.__mu_max,
                                         self.__iterations,
                                         self.__reason,self.__fevals,
                                         self.__jevals,
                                         self.__nsystems],
                                        dtype=np.double)

            # additional attributes
            self.__n=0
            self.__ncameras = 0
            self.__nprojs = 0 # number of projections of 3D->2D in an image
            self.__nvars = 0 # number of free variables adjusted
            self.__method = "not defined" # e.g. MotionStructureExpert_Cameras
            self.__expert = False 
            self.__analjac = False # have analytic Jacobian
            self.__havecovx = False # have input covx for image poitns
            self.__havedist = False # using distortion?
            self.__ncdist = 0 # if so, how many distortion coeffs fixed? 
            self.__fixedcal = False # using fixed cal for all cameras? 
            self.__nccalib = 0 # if not, how many intrinsics are fixed?  
            self.__t0 = 0. # start time of last sparse bundle adjustment
            self.__tf = 0. # end time of last sparse bundle adjustment
            
            
    def why(self):
        """What is the sound of one hand clapping?"""
        return INFO_REASONS[int(self.__reason)]

    def toC(self):
        """Provides C format info array as double *"""
        return CInfoType(self.__infoarray[0],self.__infoarray[1],
                         self.__infoarray[2],self.__infoarray[3],
                         self.__infoarray[4],self.__infoarray[5],
                         self.__infoarray[6],self.__infoarray[7],
                         self.__infoarray[8],self.__infoarray[9])

    def fromC(self,cinfo):
        temp = np.array([cinfo[0],cinfo[1],cinfo[2],cinfo[3],cinfo[4],
                         cinfo[5],cinfo[6],cinfo[7],cinfo[8],cinfo[9]],
                        dtype=np.double)
        self.infoarray=temp

    def _geteps0(self):
        return self.__eps0
    def _seteps0(self,value):
        self.__eps0 = value
        self.__infoarray[0] = self.__eps0
    eps0 = property(_geteps0,_seteps0,doc="eps0 initial error")

    def _getepsF(self):
        return self.__epsF
    def _setepsF(self,value):
        self.__epsF = value
        self.__infoarray[1] = self.__epsF
    epsF = property(_getepsF,_setepsF,doc="epsF final error")

    def _getJTeps(self):
        return self.__JTeps
    def _setJTeps(self,value):
        self.__JTeps = value
        self.__infoarray[2] = self.__JTeps
    JTeps = property(_getJTeps,_setJTeps,doc="J^Teps_inf LM parameter")

    def _getDP(self):
        return self.__DP
    def _setDP(self,value):
        self.__DP = value
        self.__infoarray[3] = self.__DP
    DP = property(_getDP,_setDP,doc="DP LM parameter")

    def _getmu_max(self):
        return self.__mu_max
    def _setmu_max(self,value):
        self.__mu_max = value
        self.__infoarray[4] = self.__mu_max
    mu_max = property(_getmu_max,_setmu_max,doc="mu_max LM parameter")

    def _getiterations(self):
        return int(self.__iterations)
    def _setiterations(self,value):
        self.__iterations = value
        self.__infoarray[5] = self.__iterations
    iterations = property(_getiterations,_setiterations,doc="number of iterations")

    def _getreason(self):
        return int(self.__reason)
    def _setreason(self,value):
        self.__reason = value
        self.__infoarray[6] = self.__reason
    reason = property(_getreason,_setreason,doc="reason for stopping, see also INFO_REASONS")

    def _getfevals(self):
        return int(self.__fevals)
    def _setfevals(self,value):
        self.__fevals = value
        self.__infoarray[7] = self.__fevals
    fevals = property(_getfevals,_setfevals,doc="number of projection function evals")

    def _getjevals(self):
        return int(self.__jevals)
    def _setjevals(self,value):
        self.__jevals = value
        self.__infoarray[8] = self.__jevals
    jevals = property(_getjevals,_setjevals,doc="number of Jacobian function evals")

    def _getnsystems(self):
        return self.__nsystems
    def _setnsystems(self,value):
        self.__nsystems = value
        self.__infoarray[9] = self.__nsystems
    nsystems = property(_getnsystems,_setnsystems,doc="number of linear systems solved")

    def _getinfoarray(self):
        return self.__infoarray
    def _setinfoarray(self,newinfoarray):
        self.__infoarray = np.array(newinfoarray,dtype=np.double)
        self.__eps0 = self.__infoarray[0]
        self.__epsF = self.__infoarray[1]
        self.__JTeps = self.__infoarray[2]
        self.__DP = self.__infoarray[3]
        self.__mu_max = self.__infoarray[4]
        self.__iterations = self.__infoarray[5]
        self.__reason = self.__infoarray[6]
        self.__fevals = self.__infoarray[7]
        self.__jevals = self.__infoarray[8]
        self.__nsystems = self.__infoarray[9]
    infoarray = property(_getinfoarray,_setinfoarray,doc="10 element info array")

    def __str__(self):
        """Human readable representation of sba.Information()"""
        return "<sba.Information() object, last stopped for {0}>".format(self.why()) 

    def _getn(self):
        return self.__n
    def _setn(self,value):
        self.__n = value
    n = property(_getn,_setn,doc="Number of 3D points")

    def _getncameras(self):
        return self.__ncameras
    def _setncameras(self,value):
        self.__ncameras = value
    ncameras = property(_getncameras,_setncameras,doc="Number of cameras")

    def _getnprojs(self):
        return self.__nprojs
    def _setnprojs(self,value):
        self.__nprojs = value
    nprojs = property(_getnprojs,_setnprojs,
                      doc="""
    number of projections of 3D-2D, over all cams. In SparseBundleAdjust, this is computed
    as the sum of points.vmask""")

    def _getnvars(self):
        return self.__nvars
    def _setnvars(self,value):
        self.__nvars = value
    nvars = property(_getnvars,_setnvars,
                     doc="""
    number of free variables adjusted.  In SparseBundleAdjust, this is computed as:
    ncameras*cnp+npts3d*pnp is using motion and structure
    ncameras*cnp if using motion only
    npts3d*pnp if using structure only""")

    def _getmethod(self):
        return self.__method
    def _setmethod(self,value):
        self.__method = value
    method = property(_getmethod,_setmethod,
                      doc="""e.g. MotionStructureExpert_Cameras.  In SparseBundleAdjust, 
    this is found using the dispatcher's projection.__str__ method.""")

    def _getexpert(self):
        return self.__expert
    def _setexpert(self,value):
        self.__expert = value
    expert = property(_getexpert,_setexpert,
                      doc="used expert mode (True) or simple mode (False)?")

    def _getanaljac(self):
        return self.__analjac
    def _setanaljac(self,value):
        self.__analjac = value
    analjac = property(_getanaljac,_setanaljac,
                       doc="used analytic Jacobian or approximate? from options")

    def _gethascovx(self):
        return self.__havecovx
    def _sethascovx(self,value):
        self.__havecovx = value
    hascovx = property(_gethascovx,_sethascovx,
                        doc="has input image point covariances? from points")

    def _gethavedist(self):
        return self.__havedist
    def _sethavedist(self,value):
        self.__havedist = value
    hasDistortion = property(_gethavedist,_sethavedist,
                             doc="had distortion in camera model? from camera")
    
    def _getncdist(self):
        return self.__ncdist
    def _setncdist(self,value):
        self.__ncdist = value
    ncdist = property(_getncdist,_setncdist,
                      doc="if distortion was used, how many dist coeffs were fixed? from options")

    def _getfixedcal(self):
        return self.__fixedcal
    def _setfixedcal(self,value):
        self.__fixedcal = value
    hasFixedCal = property(_getfixedcal, _setfixedcal, 
                           doc="using fixed cal for all cameras? from options")

    def _getnccalib(self):
        return self.__nccalib
    def _setnccalib(self,value):
        self.__nccalib = value
    nccalib = property(_getnccalib, _setnccalib,
                       doc="if camera calibrations varied, how many intrinsics were fixed? from options")
    
    def _gett0(self):
        return self.__t0
    def _sett0(self,value):
        self.__t0 = value
    t0 = property(_gett0,_sett0,
                  doc="""start time of last sparse bundle adjustment.
    from time.clock() at start of SparseBundleAdjust""")

    def _gettf(self):
        return self.__tf
    def _settf(self,value):
        self.__tf = value
    tf = property(_gettf,_settf,
                  doc="""end time of last sparse bundle adjustment. 
    from time.clock() at end of SparseBundleAdjust""")

    def _fromInput(cls,cameras,points,options):
        """
        Alternate constructor (a class method) to create default
        Information object with values drawn from input objects
        cameras, points, and options.
        """
        newinfo = cls()

        newinfo.n = int(points.n)
        newinfo.ncameras = int(cameras.ncameras)
        newinfo.nprojs = int(np.sum(np.sum(points.vmask)))

        if options.motstruct == sba.options.OPTS_MOTSTRUCT:
            newinfo.nvars = int(cameras.ncameras*options.cnp+points.n*options.pnp)
        elif options.motstruct == sba.options.OPTS_MOT:
            newinfo.nvars = int(cameras.ncameras*options.cnp)
        elif options.motstruct == sba.options.OPTS_STRUCT:
            newinfo.nvars = int(points.n*options.pnp)

        newinfo.method = sba.routines.whichRoutine(options).__str__
        # do this during SparseBundleAdjust?

        newinfo.expert = options.expert
        newinfo.analjac = options.analjac
        newinfo.hascovx = points.hascovx
        newinfo.hasDistortion = options.hasDistortion
        newinfo.ncdist = options.ncdist
        newinfo.hasFixedCal = options.hasFixedCal
        newinfo.nccalib = options.nccalib
        return newinfo
    fromInput = classmethod(_fromInput)


    # pretty print methods
    def printResults(self):
        """Print results to stdout"""
        print ""
        print "SBA using {0} 3D pts, {1} frames and {2} image projections, {3} variables".format(self.n, self.ncameras, self.nprojs, self.nvars)
        print ""
        print "Method {0}, {expert} driver, {analjac} Jacobian, {hascovx} covariances, {hasDistortion} distortion{ncdist}, {fixedcal} intrinsics{nccalib}".format(self.method,expert="expert" if self.expert else "simple",analjac="analytic" if self.analjac else "approximate",hascovx="with" if self.hascovx else "without",hasDistortion="variable" if self.hasDistortion else "without",ncdist=" ({0} fixed)".format(self.ncdist) if self.hasDistortion else "",fixedcal="fixed" if self.hasFixedCal else "variable",nccalib=" ({0} fixed)".format(self.nccalib) if not self.hasFixedCal else "")
        print ""
        print "SBA returned {0} in {1} iter, reason {2}, error {3:f} [initial {4:f}, {5}/{6} func/fjac evals, {7} lin. systems".format(self.n,long(self.iterations),int(self.reason),self.epsF/float(self.nprojs),self.eps0/float(self.nprojs),long(self.fevals),long(self.jevals),long(self.nsystems))
        print "(Reason {0} is {1})".format(int(self.reason),self.why())
        print "Elapsed time: {0} seconds, {1} msecs".format(self.tf-self.t0,(self.tf-self.t0)*1000.)


        
