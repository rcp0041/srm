#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 26 13:04:02 2025

@author: ray
"""

import numpy as np
from srm import cpgrain

class TaperedGrain():
    """ Grain with a smoothly tapered cross-sectional area """
    def __init__(self):
        pass
    
    def burn_area(self):
        pass

Ro = 6
Ri = 4
L = 12

test_grain = cpgrain(Ro,Ri,L)

N_segments = 10
segments = []
for i in range(0,N_segments):
    segments.append(cpgrain(Ro,Ri,L/N_segments))

def segs_Ab(segments):
    Ab = 0
    for segment in segments:
        Ab = Ab + segment.burn_area(0)
    return Ab

print(f"Burn area of test_grain: {test_grain.burn_area(0)}\nBurn area of test segmented grain: {segs_Ab(segments)}")