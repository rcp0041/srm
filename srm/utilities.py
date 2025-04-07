#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  7 14:58:24 2025

@author: ray
"""

import numpy as np
from matplotlib import pyplot as plt
from srm import g0

class data():
    def __init__(self,motor,timestep=0.1):
        self.motor = motor
        self.timestep = timestep
        self.Y = motor.burn_vector(timestep)
        self.time = self.Y[:,0]
        self.burn_distance = self.Y[:,1]
        self.pressure = self.Y[:,2]
        self.thrust = self.Y[:,3]
        self.mass_flow_rate = np.zeros_like(self.time)
        self.mass_flow_rate[0] = self.motor.mass_flow_rate(0)
        self.impulse = np.zeros_like(self.time)
        self.burn_area = np.zeros_like(self.time)
        self.specific_impulse = np.zeros_like(self.time)
        self.specific_impulse[0] = self.thrust[0]/(self.motor.mass_flow_rate(0)*g0)
        self.burn_area[0] = self.motor.grain.burn_area(0)
        for i in range(1,len(self.time)):
            self.impulse[i] = self.impulse[i-1] + self.thrust[i]*self.timestep
            self.burn_area[i] = self.motor.grain.burn_area(self.time[i])
            self.mass_flow_rate[i] = self.motor.mass_flow_rate(self.time[i])
            self.specific_impulse[i] = self.thrust[i]/(self.mass_flow_rate[i]*g0)
    
    def plot(self, items: list):
    # TODO: Use a dictionary or something to streamline this plot decoration
            for item in items:
                if item == 'burn_distance':
                    plot_title = "Burn Distance"
                    plot_label = "in"
                
                if item == 'pressure':
                    plot_title = "Chamber Pressure"
                    plot_label = "psi"
                    
                if item == 'thrust':
                    plot_title = "Thrust"
                    plot_label = "lbf"
                    
                if item == 'mass_flow_rate':
                    plot_title = "Mass flow rate"
                    plot_label = "slug/s"
                
                if item == 'impulse':
                    plot_title = "Impulse"
                    plot_label = "lbf-s"
                    
                if item == 'specific_impulse':
                    plot_title = "Specific Impulse"
                    plot_label = "Isp (s)"
                    
                if item == 'burn_area':
                    plot_title = "Burn Area"
                    plot_label = "in^2"
                
                fig, ax = plt.subplots(figsize=(7, 7), dpi=96)
                plt.title(plot_title)
                plt.xlabel("Burn Time (s)")
                plt.ylabel(plot_label)
                ax.plot(self.time,getattr(self,item))
                plt.show()