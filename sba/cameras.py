#!/usr/bin/env/python
"""
cameras.py
For dealing with cameras in sba only
"""


import logging
import numpy as np
import ctypes
import sba.quaternions as quaternions # need this for composition of rotations

CAMS_CNP_ALL = 5+5+3+3 # cnp for cameras with distortion - use imag q only
CAMS_CNP_NODIST = 5+3+3 # cnp for cameras without distortion - use imag q only
CAMS_CNP_INTR_FIXED = 3+3 # cnp for cameras with fixed intrinsics

class Cameras(object):
    """
    Cameras for use in sba routines are an array of ncamera rows,
    with cnp columns (number of camera parameters). The camera 
    parameters are (one row for each camera):

    fx cx cy AR s r2 r4 t1 t2 r6 qi qj qk tx ty tz
   
    note q0 is NOT included since q0 = sqrt(1-qi**2-qj**2-qk**2)!!!
    
    where: 
    fx - focal length for x in pixels
    cx, cy - principal point in pixels, typically in middle of image
    AR - aspect ratio, typically 1.
    s - skew, typically 0.
    r2,r4,r6 - radial distortion coeffs according to Bourguet (dimensionless)
    t1,t2 - tangential distortion coeffs (dimensionless)
    q0,qi,qj,qk - unit quaternion rotation of camera (dimensionless)
    tx,ty,tz - translation of camera in real world units

    The class can be asked for the array, the raveled array, 
    the number of cameras (ncameras) and the parameters for one camera (cnp).
    """
    def __init__(self,camarray=None):
        if camarray is None:
            self.__camarray = np.zeros((0,0),dtype=np.double)
        else:
            self.__camarray = np.array(camarray,dtype=np.double)
        self.__ncameras = self.__camarray.shape[0]
        self.__cnp = self.__camarray.shape[1]
        self.__rot0params = np.zeros((self.__ncameras,4),dtype=np.double)

    def _getcamarray(self):
        return self.__camarray
    def _setcamarray(self,newcamarray):
        if newcamarray is None:
            self.__camarray = None
            self.__initrot = None
            self.__ncameras = 0
            self.__cnp = 0
        else:
            self.__camarray = np.array(newcamarray,dtype=np.double)
            self.__ncameras = self.__camarray.shape[0]
            self.__cnp = self.__camarray.shape[1]
    camarray = property(_getcamarray,_setcamarray,doc="""
Camera array for sba
Note that it is of the format for full problem (cnp=5+5+3+3)
fx1 cx1 cy1 ar1 s1 r21 r41 t11 t21 r61 qi1 qj1 qk1 tx1 ty1 tz1

q0 is not carried around because it is not a free, independent variable
q is a unit quaternion so q0 must equal 1-(length of imag part)
""")

    def _getrot0params(self):
        return self.__rot0params
    def _setrot0params(self,newrot0params):
        self.__rot0params = np.array(newrot0params,dtype=np.double)
        assert self.__rot0params.shape == (self.__ncameras,4)
    rot0params = property(_getrot0params,_setrot0params,doc="""
Camera initial rotation for sba as quaternions
q0 qi qj qk
q0 is carried around here. ?!?
""")

    def _getcnp(self):
        return self.__cnp
    cnp = property(_getcnp,doc="cnp = number of parameters for 1 camera")
    
    @property
    def hasFixedIntrinsics(self):
        return self.cnp == CAM_CNP_INTR_FIXED
    
    @property
    def hasDistortion(self):
        return self.cnp == CAM_CNP_ALL

    def _getncameras(self):
        return self.__ncameras
    ncameras = property(_getncameras,doc="ncameras = m = number of cameras")

    def zeroLocalRotationEstimates(self):
        """This is because sba is crazy."""
        A = self.camarray
        A[:,-6:-3]=0.
        self.camarray = A
        
    def combineLocalAndInitialRotation(self,initial):
        """Not sure what this does but Lourakis has it. It is used to
        combine local (optimization result?) rotations and initial rotations
        in cases that allow camera motion."""
        for cam in xrange(self.ncameras):
            v = self.camarray[cam,-6:-3] # retrieve vector part from self
            v0 = 1-np.sqrt(v[0]**2.+v[1]**2.+v[2]**2) # compute v0

            # turn it into quaternions
            qlocal = quaternions.Quaternion(v0,v[0],v[1],v[2])
            qinitial = quaternions.Quaternion(initial[cam,0],initial[cam,1],
                                              initial[cam,2],initial[cam,3])
            # combine(compose) rotations by quaternion rotation
            qresult = qlocal.mul(qinitial)
            result = qresult.asVector()
            
            # Lourakis likes positive vectors i guess?
            self.camarray[cam,-6:-3] = np.copysign(result[1:],result[0])

    def Aravel(self):
        """
        Use this to ravel the camera array as A when providing to sba
        It will have to be combined with B when providing to sba C routines

        It will provide a one dimensional Numpy array like this:
        fx0 cx0 cy0 AR0 s0 r20 r40 t10 t20 r60 qi0 qj0 qk0 tx0 ty0 tz0...
        fx1 cx1 cy1 AR1 s1 r21 r41 t11 t21 r61 qi1 qj1 qk1 tx1 ty1 tz1...
        etc
        fxN cxN cyN ARN sN r2N r4N t1N t2N r6N q0N qiN qjN qkN txN tyN tzN...

        The length of this array should be ncameras*cnp
        """
        return self.__camarray.ravel()

    def cnpToC(self):
        """Provides cnp as a C int"""
        return ctypes.c_int(self.cnp)

    def ncamerasToC(self):
        """Provides ncameras as a C int"""
        return ctypes.c_int(self.ncameras)

    def __repr__(self):
        """Provide readable representation of Cameras"""
        return "<sba.Cameras() with {0} cameras, each with {1} parameters (q0 not independent)>".format(self.ncameras,self.cnp)

    def _fromTxt(cls,camfile):
        """
        Alternative constructor to read cameras from a text file of the format:

        #comment
        fx cx cy AR s r2 r4 t1 t2 r6 q0 qi qj qk tx ty tz
        
        fx - focal length for x in pixels
        cx, cy - principal point in pixels, typically in middle of image
        AR - aspect ratio, typically 1.
        s - skew, typically 0.
        r2,r4,r6 - radial distortion coeffs according to Bourguet (dimeless)
        t1,t2 - tangential distortion coeffs (dimensionless)
        q0,qi,qj,qk - unit quaternion rotation of camera (dimensionless)
        tx,ty,tz - translation of camera in real world units
        """
        #logging.debug("camerasFromTxt() reading cameras from {0}".format(camfile))
        raw = np.genfromtxt(camfile)
        if raw.shape[1]==17:
            # got a file with distortion coeffs in middle
            # cnp = 5+5+3+3
            camarray = np.hstack((raw[:,0:10],raw[:,-6:])) 
        else:
            # got a file with no distortion coeffs in middle
            # cnp = 5+3+3
            camarray = np.hstack((raw[:,0:5],raw[:,-6:]))
        #print raw,raw.shape
        #print camarray,camarray.shape
        result = Cameras(camarray)
        result.rot0params = raw[:,-7:-3]
        return result
    fromTxt = classmethod(_fromTxt)

    def toTxt(self,camoname):
        """
        Save cameras to a text file of the format:
        fx cx cy AR s r2 r4 t1 t2 r6 q0 qi qj qk tx ty tz
        
        fx - focal length for x in pixels
        cx, cy - principal point in pixels, typically in middle of image
        AR - aspect ratio, typically 1.
        s - skew, typically 0.
        r2,r4,r6 - radial distortion coeffs according to Bourguet (dimeless)
        t1,t2 - tangential distortion coeffs (dimensionless)
        q0,qi,qj,qk - unit quaternion rotation of camera (dimensionless)
        tx,ty,tz - translation of camera in real world units
        """
        #print self.camarray,self.camarray.shape
        outarray = np.zeros((self.ncameras,self.cnp+1),dtype=np.double)

        if self.cnp == 5+5+3+3:
            outarray[:,0:10] = self.camarray[:,0:10]
            outarray[:,-6:] = self.camarray[:,-6:]
        elif self.cnp == 5+3+3:
            outarray[:,0:5] = self.camarray[:,0:5]
            outarray[:,-6:] = self.camarray[:,-6:]
        elif self.cnp == 3+3:
            logging.warning("FIXED INTRINSICS NOT FULLY IMPLEMENTED")
        else:
            logging.warning("HUH????")

        # construct missing q0 real part of unit quaternions
        for cam in xrange(self.ncameras):
            outarray[cam,-7] = np.sqrt(1-self.camarray[cam,-6]**2.-self.camarray[cam,-5]**2.-self.camarray[cam,-4]**2.)

        # and save it
        np.savetxt(camoname,outarray,"%4.5f")

    def _fromTxtWithIntr(cls,camfile,intrfile):
        """
        Alternative constructor to read cameras from text files of the format:

        #comment for intrfile
        fx cx cy AR s r2 r4 t1 t2 r6 

        #comment for camfile
        q0 qi qj qk tx ty tz
        
        fx - focal length for x in pixels
        cx, cy - principal point in pixels, typically in middle of image
        AR - aspect ratio, typically 1.
        s - skew, typically 0.
        r2,r4,r6 - radial distortion coeffs according to Bourguet (dimeless)
        t1,t2 - tangential distortion coeffs (dimensionless)
        q0,qi,qj,qk - unit quaternion rotation of camera (dimensionless)
        tx,ty,tz - translation of camera in real world units
        """
        #logging.debug("camerasFromTxt() reading cameras from {0}".format(camfile))
        raw = np.genfromtxt(camfile)
        rawi = np.genfromtxt(intrfile)
        intr = np.array([rawi,]*raw.shape()[0]) # copy intrinsics to array

        camarray = np.hstack((intr,raw)) # and stack as before
        result = Cameras(camarray)
        result.rot0params = raw[:,-7:-3]
        return result
    fromTxtWithIntr = classmethod(_fromTxtWithIntr)



    # added for Dylan Ray to use in Python gui version
    def _fromDylan(cls,dylan_array):
        """
        Alternate constructor to read cameras from a numpy array of shape
        (n,17):

        each row is a different camera, the columns are:
        fx cx cy AR s r2 r4 t1 t2 r6 q0 qi qj qk tx ty tz
        
        fx - focal length for x in pixels
        cx, cy - principal point in pixels, typically in middle of image
        AR - aspect ratio, typically 1.
        s - skew, typically 0.
        r2,r4,r6 - radial distortion coeffs according to Bourguet (dimeless)
        t1,t2 - tangential distortion coeffs (dimensionless)
        q0,qi,qj,qk - unit quaternion rotation of camera (dimensionless)
        tx,ty,tz - translation of camera in real world units

        dtype is np.double
        """
        #logging.debug("sba.Cameras.fromDylan() reading cameras from array")
        #logging.warning("this feature not yet tested!")

        try: 
            assert(dylan_array.shape[1]==17)
        except:
            logging.error("bad array passed in sba.Cameras.fromDylan()")

        dylan_array = dylan_array.astype(np.double)
        camarray = np.hstack((dylan_array[:,0:10],dylan_array[:,-6:])) 
        result = Cameras(camarray)
        result.rot0params = dylan_array[:,-7:-3]
        return result
    fromDylan = classmethod(_fromDylan)



    def toDylan(self):
        """
        Returns a numpy array of shape (n,17) containing the cameras:
        each row is a different camera, the columns are:
        fx cx cy AR s r2 r4 t1 t2 r6 q0 qi qj qk tx ty tz
        dtype is np.double
        """
        result = np.array((self.ncameras, 17),dtype=np.double)
        result[:,0:10] = self.camarray[:,0:10]
        result[:,-6:] = self.camarray[:,-6:]

        # construct missing q0 real part of unit quaternions
        for cam in xrange(self.ncameras):
            result[cam,-7] = np.sqrt(1-self.camarray[cam,-6]**2.
                                     -self.camarray[cam,-5]**2.
                                     -self.camarray[cam,-4]**2.)
        return(result)

