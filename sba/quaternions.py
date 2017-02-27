#!/usr/bin/env python
"""
quaternions.py
(c) 2013 Dennis Evangelista
Some helper functions that do quaternion rotations
"""

import logging
import numpy as np

def quatFromArray(a):
    """
    Construct a quaternion from four components
    a = something with four numerical components
    """
    assert a[0]**2+a[1]**2+a[2]**2+a[3]**2 == 1
    return Quaternion(re=a[0],i=a[1],j=a[2],k=a[3])

def quatFromAa(aa):
    """
    Construct a unit quaternion from an axis-angle representation
    aa = something with 3 numerical components that can have norm taken
    """
    theta = np.linalg.norm(aa)
    if (theta!=0.):
        uvec = aa/theta
    else:
        uvec = np.array([0.,0.,0.])
    re = np.cos(theta/2.)
    i = uvec[0]*np.sin(theta/2.)
    j = uvec[1]*np.sin(theta/2.)
    k = uvec[2]*np.sin(theta/2.)
    return Quaternion(re,i,j,k)

def quatFromRotationMatrix(R):
    """
    Construct a unit quaternion from a rotation matrix
    """
    a = 0.5*np.sqrt(1+R[0,0]+R[1,1]+R[2,2])
    b = np.copysign(0.5*np.sqrt(1+R[0,0]-R[1,1]-R[2,2]),R[2,1]-R[1,2])
    c = np.copysign(0.5*np.sqrt(1-R[0,0]+R[1,1]-R[2,2]),R[0,2]-R[2,0])
    d = np.copysign(0.5*np.sqrt(1-R[0,0]-R[1,1]+R[2,2]),R[1,0]-R[0,1])
    return Quaternion(a,b,c,d)





    
class Quaternion():
    """
    Class and methods for handling rotations using (unit) quaternions
    """
    def __init__(self,re=1.,i=0.,j=0.,k=0.):
        """
        Constructs a quaternion.  Doesn't currently check or adjust if
        the components are not for a unit quaternion, which shouldn't be
        a problem so long as they are always constructed using fromaa or
        from a rotation matrix. 
        """
        self.re=re
        self.i=i
        self.j=j
        self.k=k
    def abs(self):
        """
        Returns magnitude (norm) of a quaternion. 
        This should always be 1??? May not be if not unit quaternion. 
        """
        return np.linalg.norm([self.re,self.i,self.j,self.k])
    def imag(self):
        """
        Returns imaginary part of a quaternion as 3 element Numpy array 
        """
        return np.array([self.i,self.j,self.k])
    def asAa(self):
        """
        Returns the quaternion in axis-angle representation, as 
        a numpy array of 3 components. This part assumes it is a 
        unit quaternion!
        """
        theta = 2.*np.arctan2(np.linalg.norm(self.imag()),self.re)
        if (theta!=0.):
            return self.imag()/np.sin(theta/2.)
        else:
            return np.array([0.,0.,0.])
    def asVector(self):
        """
        Returns quaternion as a 4 component vector numpy array. 
        """
        return np.array([self.re,self.i,self.j,self.k])
    def mul(self,q2):
        """
        Multiplication by another Quaternion q2.
        """
        s = self.re
        v = self.imag()
        t = q2.re
        w = q2.imag()
        result = Quaternion()
        result.re = s*t-v.dot(w)
        vecpart = s*w+t*v+np.cross(v,w)
        result.i = vecpart[0]
        result.j = vecpart[1]
        result.k = vecpart[2]
        return result
    def rotate(self,vec):
        """
        Rotates vector according to the Hamiltonian product q v q*
        vec = numpy array, vector with 3 components
        Returns a numpy array vector with 3 components
        """
        result = vec+np.cross(2*self.imag(),(np.cross(self.imag(),vec)+self.re*vec))
        return result
    def asRotationMatrix(self):
        a,b,c,d = self.asVector().tolist()
        return np.matrix([[a**2.+b**2.-c**2.-d**2.,2.*b*c-2.*a*d,2.*b*d+2.*a*d],
                          [2.*b*c+2.*a*d,a**2.-b**2.+c**2.-d**2.,2.*c*d-2.*a*b],
                          [2.*b*d-2.*a*c,2.*c*d+2.*a*b,a**2.-b**2.-c**2.+d**2.]])

    def __repr__(self):
        """
        For printout to screen of quaternion as R4. 
        """
        return self.asVector().__repr__()

    
