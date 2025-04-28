#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 26 13:04:02 2025

@author: ray
"""

import numpy as np
import srm

NUM_SEGMENTS = 1000

class SegmentedGrain():
    """ Multiple grain segments composed together into a single grain """
    def __init__(self,segments):
        """ You must give this class a numpy array of segments """
        self.segments = segments
        pass
    
    def sum_segments(self,method,argument):
        """ Allows querying the segmented grain as a single grain object """
        results = []
        for segment in self.segments:
            results.append(getattr(segment,method)(argument))
        return sum(results)
    
    def burn_area(self,x):
        # return self.sum_segments('burn_area',x)
        # def Ab(grain,x):
            # return grain.burn_area(x)
        Ab = np.vectorize(lambda grain,x: grain.burn_area(x))
        return Ab(self.segments, x).sum()    
    
    def pressure(self,x,nozzle,propellant):
        return (((self.burn_area(x)/nozzle.throat_area)*propellant.a*propellant.density*propellant.cstar)/srm.g0)**(1/(1-propellant.n))
    
    def propellant_volume(self,x):
        return self.sum_segments('propellant_volume',x)

class Tapered(SegmentedGrain):
    class CP(SegmentedGrain):
        # def __init__(self,Ro: float,Ri: tuple,length: float):
        def __init__(self,Ro:float,Ri:float,Ri_scale_factor:float,length:float):
            # assert len(Ri) == 2, "Ri must have dimension 2."
            assert Ri*Ri_scale_factor > 0, "Requested scale factor is physically impossible."
            assert Ri*Ri_scale_factor <= Ro, "Requested scale factor is physically impossible."
            self.Ro = Ro
            self.length = length
            # self.Ri_start = Ri[0]
            # self.Ri_end = Ri[1]
            self.Ri_start = Ri
            self.Ri_end = self.Ri_start*Ri_scale_factor
            self.delta_Ri = (self.Ri_end - self.Ri_start)/NUM_SEGMENTS
            self.delta_length = self.length/NUM_SEGMENTS
            self.graintype = 'tapered_CP'
            
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
                segments.append(srm.cpgrain(Ro,Ri,L))
            
            self.segments = segments
        
        def max_yf(self):
            yf = []
            for segment in self.segments:
                yf.append(segment.yf)
            return max(yf)
        
        def __str__(self):
            a = "Tapered CP grain:\n"
            a = a + f"Ro: {self.Ro}\n"
            a = a + f"Ri_start = {self.Ri_start}\n"
            a = a + f"Ri_end = {self.Ri_end}\n"
            a = a + f"delta_Ri = {self.delta_Ri}\n"
            a = a + f"length = {self.length}\n"
            a = a + f"delta_length = {self.delta_length}"
            return a
