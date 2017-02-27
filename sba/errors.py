#!/usr/bin/env python
"""
errors.py
Dennis Evangelista 2013

For handling errors in wrapped sba routines
"""


class SbaError(Exception):
    """Base class for exceptions in python-sba code"""
    def __init__(self,value):
        """Value should be the cause of the SbaError"""
        self.value = value
        #super(Exception,self).__init__()
    def __str__(self):
        """Returns string representation of SbaError"""
        return "{0}".format(self.value)
