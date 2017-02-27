#!/usr/bin/env python
"""
points.py
Dennis Evangelista, 2013
"""

import logging
import numpy as np
import ctypes

class Points(object):
    def __init__(self,B,X,vmask):
        self.__B=np.array(B,dtype=np.double)
        self.__n = self.__B.shape[0]
        self.__pnp = self.__B.shape[1]
        assert(self.__pnp==3)
        
        self.__X=np.array(X,dtype=np.double)
        self.__ncameras = self.__X.shape[1]
        self.__mnp = self.__X.shape[2]
        assert(self.__mnp==2)
        assert(self.__X.shape[0]==self.__n)

        self.__vmask=np.array(vmask,dtype=np.byte)
        assert(self.__vmask.shape==(self.__n,self.__ncameras))

        self.__covx = None

    def _getB(self):
        return self.__B
    def _setB(self,value):
        self.__B=np.array(value,dtype=np.double)
        self.__n = self.__B.shape[0]
        self.__pnp = self.__B.shape[1]
        assert(self.__pnp==3)
    B = property(_getB,_setB,doc="B is the 3D worldPoints array")

    def _getX(self):
        return self.__X
    def _setX(self,value):
        self.__X=np.array(value,dtype=np.double)
        self.__ncameras = self.__X.shape[1]
        self.__mnp = self.__X.shape[2]
        assert(self.__mnp==2)
    X = property(_getX,_setX,doc="X is the 2D imagePoints nonmasked array")
    
    def _getcovx(self):
        return self.__covx
    def _setcovx(self,value):
        self.__covx = value
    covx = property(_getcovx,_setcovx,doc="covx is not used and should be None")
    def covxToC(self):
        if self.__covx == None:
            return ctypes.POINTER(ctypes.c_double)()
        else:
            return self.__covx.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        
    def _getvmask(self):
        return self.__vmask
    def _setvmask(self,value):
        self.__vmask=np.array(value,dtype=np.byte)
    vmask = property(_getvmask,_setvmask,doc="vmask is the 2D imagePoints mask")

    @property
    def n(self):
        """n - number of points"""
        return self.__n

    @property
    def pnp(self):
        """pnp - number of parameters per world point, should be 3"""
        return self.__pnp

    @property
    def mnp(self):
        """mnp - number of parameters per image point, should be 2"""
        return self.__mnp

    @property
    def hascovx(self):
        """True if it has a covx that is not None"""
        if self.covx==None:
            return False
        else:
            return True

    @property
    def ncameras(self):
        """ncameras = m = number of parameters for 1 cameras"""
        return self.__ncameras

    def Bravel(self):
        """Use this to ravel B for combination with A in preparing sba inputs"""
        return self.B.ravel()

    def vmaskToC(self):
        """Converts vmask to C char ``*``vmask"""
        ct = self.__vmask.ravel()
        return ct.ctypes.data_as(ctypes.c_char_p)
    
    def XtoC(self):
        """Converts X to C double ``*``x
        Rewritted Jan 2015 by Nick Gravish to address strange bugs seen
        with large initial error; In eucsbademo.c from Lourakis, the
        points are read in this manner in readparams.c so this is correct.
        Points are raveled all packed together with the vmask holding 
        where they actually are (odd representation of sparse array done
        for reasons of Lourakis' C-ness.)
        -DE
        """
        temp = np.zeros(self.n*self.ncameras*self.mnp,dtype=np.double)
        idx = 0
        for b in xrange(self.n):
            for a in xrange(self.ncameras):
                if(self.vmask[b,a] == 1):
                    for c in xrange(2):
                        temp[idx] = self.__X[b,a,c]
                        idx += 1
        #logging.debug("Created temp {0},{1}".format(temp,temp.shape))
        return temp#.ctypes.data_as(ctypes.POINTER(ctypes.c_double))

    def __str__(self):
        """Provides readable representation of Points"""
        return("<sba.Points() with {0} world points (R{2}) in {1} cameras (R{3})>".format(self.n,self.ncameras,self.pnp,self.mnp))

    def _fromTxt(cls,ptsfile,ncameras):
        """Alternate constructor
        Construct array of points from sba formatted text file"""
        
        #logging.debug("pointsFromTxt() reading points from {0}".format(ptsfile))
        ifile = file(ptsfile,"r")
        blist = []
        for line in ifile:
            if line[0] != "#":
                splitted = line.split()
                blist.append([float(splitted[0]),
                              float(splitted[1]),
                              float(splitted[2])])
        ifile.close()
        b = np.array(blist,dtype=np.double)
        npoints = len(blist)

        ifile = file(ptsfile,"r")
        x = np.zeros((npoints,ncameras,2),dtype=np.double)
        vmask = np.zeros((npoints,ncameras),dtype=np.byte)
        pt = 0

        for line in ifile:
            if line[0] != "#":
                splitted = line.split()
                nvisible = int(float(splitted[3]))
                for a in xrange(nvisible):
                    cam = int(float(splitted[3+1+a*3]))
                    x[pt,cam,0] = float(splitted[3+2+a*3])
                    x[pt,cam,1] = float(splitted[3+3+a*3])
                    vmask[pt,cam] = 1
                pt += 1
        #print b,x,vmask
        result = cls(b,x,vmask)
        return result
    fromTxt = classmethod(_fromTxt)

    def toTxt(self,ptsoname):
        np.savetxt(ptsoname,self.B,"%4.6f")





    # added for Dylan Ray to use in Python gui version
    def _fromDylan(cls,dylan_array):
        """
        Alternate constructor to read points from a numpy array of shape
        (npoints,3+2*ncameras).

        Each row is a point + camera observations of the point.  
        The columns are as follows:
        X Y Z u1 v1 u2 v2...
        
        If a point is not observed in camera n, the value of un,vn is np.nan
        The shape of the array should be (Npoints,2*Ncameras+3)

        dtype is np.double
        """
        
        #logging.debug("sba.Points.fromDylan() reading points from numpy array")
        #logging.warning("this feature not yet tested!")

        # straight copy the points' X,Y,Z from dylan_array
        b = dylan_array[:,:3].astype(np.double)
        npoints = len(dylan_array)
        ncameras = (dylan_array.shape[1]-3)/2

        #ifile = file(ptsfile,"r")
        # to get x and vmask requires some reshaping
        x = np.zeros((npoints,ncameras,2),dtype=np.double)
        vmask = np.zeros((npoints,ncameras),dtype=np.byte)
        pt = 0

        for aa in xrange(npoints):
            for bb in xrange(ncameras):
                u = dylan_array[aa,3+2*bb]
                v = dylan_array[aa,3+2*bb+1]
                if (not(np.isnan(u))) and (not(np.isnan(v))):
                    x[aa,bb,0] = u
                    x[aa,bb,1] = v
                    vmask[aa,bb] = 1
                # otherwise vmask is already 0
 
        #print "b",b.shape,b
        #print "x",x.shape,x
        #print "vmask",vmask.shape,vmask

        result = cls(b,x,vmask)
        return result
    fromDylan = classmethod(_fromDylan)


    def toDylan(self):
        """
        Returns a single numpy array of shape (n,3+2*ncameras) 
        with np.nan if not visible in cameras. 
        """
        result = np.zeros((self.__n, 
                           3+self.__ncameras*2))
        result[:,:3] = self.X
        result[:,3:] = self.B
        r,c = self.vmask.shape
        for eachrow in xrange(self.__n):
            for eachcol in xrange(self.__ncameras):
                if not(self.vmask[eachrow,eachcol]):
                    result[eachrow,3+2*eachcol] = np.nan
                    result[eachrow,3+2*eachcol+1] = np.nan

        return(result)
