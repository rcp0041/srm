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
    
    def pressure(self,x,nozzle,propellant):
        return (((self.burn_area(x)/nozzle.throat_area)*propellant.a*propellant.density*propellant.cstar)/srm.g0)**(1/(1-propellant.n))
    
    def propellant_volume(self,x):
        return self.sum_segments('propellant_volume',x)

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
            self.graintype = 'CP'
            
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
Ro = 6
Ri = 5
L = 12
foo = Tapered.CP(Ro,Ri,0.5,L)
foo.graintype = 'tapered_CP'

motor = srm.motor(
    propellant=srm.materials.propellants['PBAN_AP_AL'],
    grain=foo,
    nozzle=srm.nozzle(
        exit_pressure = 14.7,
        ambient_pressure = 14.7,
        exit_diameter = foo.Ro,
        expansion_ratio = 8
        )
    )

