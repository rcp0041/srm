#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 26 13:04:02 2025

@author: ray
"""

import numpy as np
from srm import cpgrain, stargrain

NUM_SEGMENTS = 1000

class SegmentedGrain():
    """ Multiple grain segments composed together into a single grain """
    def __init__(self,segments: list):
        """ If you inherit this class, you must give it a list of segments """
        self.segments = segments
        pass
    
    def sum_segments(self,method,argument):
        """ Allows querying the segmented grain as a single grain object """
        results = []
        for segment in self.segments:
            results.append(getattr(segment,method)(argument))
        return sum(results)
    
    def burn_area(self,x):
        return self.sum_segments('burn_area',x)

class Tapered(SegmentedGrain):
    class CP(SegmentedGrain):
        # def __init__(self,Ro: float,Ri: tuple,length: float):
        def __init__(self,Ro:float,Ri:float,Ri_scale_factor:float,length:float):
            # assert len(Ri) == 2, "Ri must have dimension 2."
            self.Ro = Ro
            self.length = length
            # self.Ri_start = Ri[0]
            # self.Ri_end = Ri[1]
            self.Ri_start = Ri
            self.Ri_end = self.Ri_start*Ri_scale_factor
            self.delta_Ri = (self.Ri_end - self.Ri_start)/NUM_SEGMENTS
            self.delta_length = self.length/NUM_SEGMENTS
            
            segments = []
            for i in range(0,NUM_SEGMENTS):
                if i == 0:
                    Ro = self.Ro
                    Ri = self.Ri_start
                    L = self.delta_length
                else:
                    Ro = self.Ro
                    Ri = segments[i-1].Ri + self.delta_Ri
                    L = self.delta_length
                segments.append(cpgrain(Ro,Ri,L))
            
            self.segments = segments
        
        def __str__(self):
            a = "Tapered CP grain:\n"
            a = a + f"Ro: {self.Ro}\n"
            a = a + f"Ri_start = {self.Ri_start}\n"
            a = a + f"Ri_end = {self.Ri_end}\n"
            a = a + f"delta_Ri = {self.delta_Ri}\n"
            a = a + f"length = {self.length}\n"
            a = a + f"delta_length = {self.delta_length}"
            return a
            
Ro = 6
Ri = 5
L = 12

test_grain = cpgrain(Ro,Ri,L)

segments = []
for i in range(0,NUM_SEGMENTS):
    segments.append(cpgrain(Ro,Ri,L/NUM_SEGMENTS))

def segs_Ab(segments):
    Ab = 0
    for segment in segments:
        Ab = Ab + segment.burn_area(0)
    return Ab

def cone_SA(r,h):
    return np.pi*r*np.sqrt(r**2 + h**2)

def exact_taperCP(taperedcp):
    R1 = taperedcp.Ri_start
    R2 = taperedcp.Ri_end
    L = taperedcp.length
    C1 = 2*np.pi*R1
    C2 = 2*np.pi*R2
    return 0.5*(C1+C2)*L
    
# foo = Tapered.CP(Ro,(Ri,Ri),L)
foo = Tapered.CP(Ro,Ri,1.3,L)

# See if grain is neutral
for x in np.arange(0,3,0.5):
    print(foo.burn_area(x))

# """
# Setting Ri = Rp, f = 0, epsilon = 1e-6, and programmatically determining
# halfTheta from H1 and epsilonAngle very closely approximates a CP grain.
# """
# N=3
# epsilon=1e-6
# pinAngle = np.pi/N
# epsilonAngle = np.pi*epsilon/N
# # epsilonAngle = pinAngle
# tanEpsilon = np.tan(epsilonAngle)
# H1 = Ri*np.sin(np.pi*epsilon/N)
# halfTheta = np.arctan((H1*tanEpsilon)/(H1-Ri*tanEpsilon))
# bar = stargrain(N=N,Ro=Ro,Ri=Ri,Rp=Ri,f=0,epsilon=epsilon,halfTheta=np.deg2rad(60),length=L)