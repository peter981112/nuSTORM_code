#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Class traceSpace:
=================

  Models tracespace co-ordinates

  Class attributes:

  --np       : numpy class
      
  Instance attributes:
  --------------------
  _TrcSpcCrd: Trace space (s, x, y, z, x', y') in numpy array - what are s and z?
    
  Methods:
  --------
  Built-in methods __new__, __repr__ and __str__.
      __init__ : Creates trace Space instance
      __repr__ : One liner with call.
      __str__  : Dump of values in the object
      __eq__   : two objects equal
      __ne__   : two objects not euqal

  Get/set methods:
        s  : get value of s (m)
        x  : get value of x (m)
        y  : get value of y (m)
        z  : get value of z (m)
        xp : get value of x'
        yp : get value of y'
  General methods:

Created on Tue 20Aug12: Version history:
----------------------------------------------
 1.0: 20Aug21: First implementation
 1.1: 22Aug21: Add __eq__ and __ne__ methods
 2.0: 18Jul22: Add a default constructor

@author: PaulKyberd
"""

from copy import deepcopy
import sys

class traceSpace:

    __Debug  = False
    __version = 2.0
    
#--------  "Built-in methods":
    def __init__(self, s=0.0, x=0.0, y=0.0, z=0.0, xp=0.0, yp=0.0):

        self._s  = s
        self._x  = x
        self._y  = y
        self._z  = z
        self._xp = xp
        self._yp = yp
        return

# Override the __eq__ method
    def __eq__(self, comp):
        if isinstance(comp, self.__class__):
            return self.__dict__ == comp.__dict__
        else:
            return False

# Override the __ne__ method
    def __ne__(self, comp):
        return not self.__eq__(comp)


    def __repr__(self):
        return "traceSpace co-ordinates"

    def __str__(self):
        return "s = %gm, x = %gm, y = %gm, z = %gm, x' = %g, y' = %g " % \
                  (self._s, self._x, self._y, self._z, self._xp, self._yp  )

    def __eq__(self, comp):
        if isinstance(comp, self.__class__):
            return self.__dict__ == comp.__dict__
        else:
            return False

    def __ne__(self,comp):
        return not self.__eq__(comp)
    
#--------  get/set methods:
    def s(self):
        return deepcopy(self._s)

    def x(self):
        return deepcopy(self._x)

    def y(self):
        return deepcopy(self._y)

    def z(self):
        return deepcopy(self._z)

    def xp(self):
        return deepcopy(self._xp)

    def yp(self):
        return deepcopy(self._yp)

    def Version(self):
        return deepcopy(self.__version)
    
