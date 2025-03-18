#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
Solid Rocket Motor library
Created on Wed Mar 12 22:43:15 2025

@author: Ray Patrick
"""

import numpy as np                    # Numerical operations
from scipy import optimize            # Root finder

# Global variables:
Ru = 1545
g0 = 32.174

class motor:
    def __init__(self,propellant,grain,nozzle):
        self.propellant = propellant
        self.grain = grain
        self.nozzle = nozzle
        
        # Isentropic relations
        if self.nozzle.P0 != None and self.nozzle.throat_area == None:
            self.nozzle.throat_area = (self.grain.Ab(0)*self.propellant.a*self.propellant.density*self.propellant.cstar)/(g0*self.nozzle.P0**(1-self.propellant.n))
        
        if self.nozzle.expansion_ratio != None and self.nozzle.throat_area != None:
            def f(Me):
                gamma = self.propellant.specific_heat_ratio
                return (1/Me)*((2/(gamma+1))*(1 + ((gamma-1)/2)*Me**2))**((gamma+1)/(2*(gamma-1))) - nozzle.expansion_ratio
            self.nozzle.Me = optimize.root_scalar(f, method='newton',x0=2).root
        
        # if self.nozzle.expansion ratio != None and self.nozzle.Pe != None and self.nozzle.Pa != None and self.nozzle.P0 != None:
            # self.nozzle.thrust_coefficient = self.nozzle.expansion_ratio*(self.nozzle.Pe-self.nozzle.Pa)/self.nozzle.P0
            
        if self.nozzle.Me == None and self.nozzle.P0 != None and self.nozzle.Pe != None:
            def f(Me):
                return (1 + 0.5*(self.propellant.specific_heat_ratio-1)*Me**2)**(self.propellant.specific_heat_ratio/(self.propellant.specific_heat_ratio-1))
            self.nozzle.Me = optimize.root_scalar(f, method='newton',x0=1).root

        if self.nozzle.Me != None:
            self.nozzle.pressure_ratio = (1 + ((self.propellant.specific_heat_ratio-1)/2)*self.nozzle.Me**2)**(self.propellant.specific_heat_ratio/(self.propellant.specific_heat_ratio-1))
            self.nozzle.temperature_ratio = self.nozzle.pressure_ratio**((self.propellant.specific_heat_ratio-1)/self.propellant.specific_heat_ratio)
            self.nozzle.area_ratio = (1/self.nozzle.Me)*((2/(self.propellant.specific_heat_ratio+1))*(1+((self.propellant.specific_heat_ratio-1)/2)*self.nozzle.Me**2))**((self.propellant.specific_heat_ratio+1)/(2*(self.propellant.specific_heat_ratio-1)))
        
            self.nozzle.Te = self.propellant.temperature/self.nozzle.temperature_ratio
            self.nozzle.Ae = self.nozzle.throat_area * self.nozzle.area_ratio
            self.nozzle.ue = self.nozzle.Me * np.sqrt(self.propellant.specific_heat_ratio*self.propellant.R*g0*self.nozzle.Te)
            
        # Grain-specific stuff
        if self.grain.graintype == "CP":
            # Burn rate
            self.grain.K1 = (2*np.pi*self.grain.length*self.propellant.a*((self.propellant.a*self.propellant.cstar*self.propellant.density)/(g0*self.nozzle.throat_area))**(self.propellant.n/(1-self.propellant.n)))
            self.grain.K2 = self.grain.K1*((2*self.propellant.n-1)/(self.propellant.n-1))
            self.grain.Abi = 2*np.pi*self.grain.Ri*self.grain.length
            self.grain.Abf = 2*np.pi*self.grain.Ro*self.grain.length
            self.grain.p = (2*self.propellant.n - 1)/(self.propellant.n-1)
            if self.propellant.n != 0.5:
                self.grain.q = 1/self.grain.p
                self.grain.tb = (self.grain.Abf**((2*self.propellant.n-1)/(self.propellant.n-1)) - self.grain.Abi**((2*self.propellant.n-1)/(self.propellant.n-1)))/self.grain.K2
            elif self.propellant.n == 0.5:
                self.grain.tb = np.log(self.grain.Abf)/(np.log(self.grain.Abi)*self.grain.K1)
        elif self.grain.graintype == "star":
            # Burn rate
            if self.nozzle.P0 != None:
                self.grain.r0 = self.propellant.a*self.nozzle.P0**self.propellant.n
                if self.grain.neutrality == 0:
                    self.grain.tb_phase1 = self.grain.web1/self.grain.r0
    
    def burn_vector(self,timestep,initial_time=0):
        """ Returns a state vector """
        Y = np.array([initial_time,0,self.pressure(0),self.thrust(0),self.specific_impulse(0)])
        t,y,r = initial_time,0,self.burn_rate(0)
        if self.grain.graintype == "CP":
            yf = self.grain.Ro - self.grain.Ri
            while t < self.grain.tb:
                t = t + timestep
                y = y + self.burn_rate(t)*timestep
                r = self.burn_rate(t)
                Y = np.vstack((Y,np.array([t,y,self.pressure(t),self.thrust(t),self.specific_impulse(t)])))    
        elif self.grain.graintype != "CP":
            yf = self.grain.Ro - self.grain.Rp - self.grain.f
            while y < yf:
                t = t + timestep
                y = y + r*timestep
                r = self.burn_rate(y)
                Y = np.vstack((Y,np.array([t,y,self.pressure(y),self.thrust(y),self.specific_impulse(y)])))
        return Y
    
    def burn_rate(self,x):
        """ For star grains, 'x' is burn distance """
        return self.propellant.a * self.pressure(x)**self.propellant.n
    
    def burn_area(self,x):
        """ For CP grains, 'x' is time. For star grains, 'x' is burn distance."""
        if self.grain.graintype == "CP":
            return self.grain.Ab(x,self.propellant.n)
        else:
            return self.grain.Ab(x)
    
    def pressure(self,x):
        """ For CP grains, 'x' is time. For star grains, 'x' is burn distance."""
        return self.grain.pressure(x,self.nozzle,self.propellant)
    
    def mdot(self,x):
        """ For CP grains, 'x' is time. For star grains, 'x' is burn distance."""
        return (self.pressure(x)*self.nozzle.throat_area)/self.propellant.cstar
    
    def prop_mass(self,x):
        """ 'x' is burn distance for all grain types """
        return ((np.pi*self.grain.Ro**2) - self.grain.Ap(x))*self.grain.length*self.propellant.density

    def thrust(self,x):
        """ For CP grains, 'x' is time. For star grains, 'x' is burn distance."""
        return self.mdot(x)*self.nozzle.ue + self.nozzle.Ae*(self.nozzle.Pe - self.nozzle.Pa)
    
    def specific_impulse(self,x):
        """ For CP grains, 'x' is time. For star grains, 'x' is burn distance."""
        return self.thrust(x)/(self.mdot(x)*g0)

    class propellant():
        """ Solid rocket propellant """
        def __init__(self,a,n,density,temperature,specific_heat_ratio,cstar,MW):
            self.a = a
            self.n = n
            self.density = density
            self.temperature = temperature
            self.specific_heat_ratio = specific_heat_ratio
            self.cstar = cstar
            self.MW = MW
            self.R = Ru/self.MW
    
    class grain():
        """ Solid rocket fuel grain """
        def pressure(self,x,nozzle,propellant):
            """ For CP grains, 'x' is time. For star grains, 'x' is burn distance."""
            if self.graintype == "CP":
                return (((self.Ab(x,propellant.n)/nozzle.throat_area)*propellant.a*propellant.density*propellant.cstar)/g0)**(1/(1-propellant.n))
            else:
                return (((self.Ab(x)/nozzle.throat_area)*propellant.a*propellant.density*propellant.cstar)/g0)**(1/(1-propellant.n))
        
    class cpgrain(grain):
        """ Circular-perforated grain """
        def __init__(self,Ro,Ri,length):
            self.Ro = Ro
            self.Ri = Ri
            self.length = length
            self.graintype = "CP"
        
        def Ab(self,t,n):
            if n != 0.5:
                return (self.K2*t + self.Abi**self.p)**self.q
            elif n == 0.5:
                return self.Abi * np.exp(self.K1*t)
        
        def Ap(self,y):
            return np.pi*(self.Ri + y)**2 
        
    class stargrain(grain):
        """ Internal burning star grain """
        def __init__(self,N,Ro,Ri,Rp,f,epsilon,halfTheta,length):
            self.N = N
            self.Ro = Ro
            self.Ri = Ri
            self.Rp = Rp
            self.f = f
            self.epsilon = epsilon
            self.halfTheta = halfTheta
            self.length = length
            self.graintype = "star"
            
            self.epsilonAngle = np.pi*self.epsilon/self.N
            self.H1 = self.Rp*np.sin(self.epsilonAngle)
            
            if self.Ri == None:
                self.Ri = self.find_neutral_Ri()
                
            if self.halfTheta == None:
                self.halfTheta = np.arctan((self.H1*np.tan(self.epsilonAngle))/(self.H1 - self.Ri*np.tan(self.epsilonAngle)))
            
            self.web1 = (self.Rp*np.sin(self.epsilonAngle)/np.cos(self.halfTheta))-self.f
            
            self.beta = np.pi/2 - self.halfTheta + self.epsilonAngle
        
            coeff = np.pi/2 - self.halfTheta + np.pi/self.N - 1/np.tan(self.halfTheta)
            if coeff == 0:
                self.neutrality = 0
            elif coeff > 1:
                self.neutrality = 1
            else:
                self.neutrality = -1
        
        def phase(self,y):
            y0 = self.H1/np.cos(self.halfTheta)
            if (y + self.f) < y0:
                return 1
            else:
                return 2
            
        def gamma(self,y):
            return np.arctan((np.sqrt((y+self.f)**2 - self.H1**2))/self.H1)-self.halfTheta
        
        def burn_perimeter(self,y,S=None):
            if self.phase(y) == 1:
                S1 = self.H1/np.sin(self.halfTheta) - (y+self.f)/np.tan(self.halfTheta)
                S2 = (y+self.f)*self.beta
                S3 = (self.Rp + y + self.f)*(np.pi/self.N - self.epsilonAngle)
            elif self.phase(y) == 2:
                S1 = 0
                S2 = (y+self.f)*(self.beta-self.gamma(y))
                S3 = (self.Rp+y+self.f)*(np.pi/self.N - self.epsilonAngle)
            if S == None:
                return 2*self.N*(S1+S2+S3)
            elif S == 1:
                return S1
            elif S == 2:
                return S2
            elif S == 3:
                return S3
                
        def Ab(self, y):
            burn_perimeter = self.burn_perimeter(y)
            return burn_perimeter*self.length
        
        def Ap(self, y):
            if self.phase(y) == 1:
                S1 = self.burn_perimeter(y,1)
                A1 = 0.5*self.H1*(self.Rp*np.cos(self.epsilonAngle)+self.H1*np.tan(self.halfTheta))-0.5*(S1**2)*np.tan(self.halfTheta)
                A2 = 0.5*self.beta*(y+self.f)**2
                A3 = 0.5*(np.pi/self.N - self.epsilonAngle)*(self.Rp+y+self.f)**2
            elif self.phase(y) == 2:
                A1 = 0.5*self.H1*(self.Rp*np.cos(self.epsilonAngle)+np.sqrt((y-self.f)**2 - self.H1**2))
                A2 = 0.5*(self.beta - self.gamma(y))*(y+self.f)**2
                A3 = 0.5*(np.pi/self.N - self.epsilonAngle)*(self.Rp + y + self.f)**2
            return 2*self.N*(A1+A2+A3)
        
        def find_neutral_Ri(self):
            def f(Ri):
                return np.pi/2 + np.pi/self.N - np.arctan((self.H1*np.tan(self.epsilonAngle))/(self.H1 - Ri*np.tan(self.epsilonAngle))) - (self.H1 - Ri*np.tan(self.epsilonAngle))/(self.H1*np.tan(self.epsilonAngle))
            return optimize.root_scalar(f, method='newton',x0=1).root
    
    class wagonwheelgrain(grain):
        """ Short- or long-spoke wagon wheel grain """
        def __init__(self,N,Ro,Ri,Rp,f,epsilon,halfTheta,length):
            self.N = N
            self.Ro = Ro
            self.Ri = Ri
            self.Rp = Rp
            self.f = f
            self.epsilon = epsilon
            self.halfTheta = halfTheta
            self.length = length
            self.epsilonAngle = np.pi*self.epsilon/self.N
            
            self.delta = np.arctan(
                (self.Rp*np.sin(self.epsilonAngle)-self.f)/
                (self.Ri + (self.f/np.sin(self.halfTheta)) + (self.Rp*np.sin(self.epsilonAngle)-self.f)/np.tan(self.halfTheta))
                )
            if self.delta >= np.pi/self.N:
                print("Maximum spoke length exceeded.")
                exit(1)
            
            self.h = (self.Rp*np.cos(self.epsilonAngle)-((self.Rp*np.sin(self.epsilonAngle)/(np.tan(self.halfTheta))))-self.Ri)*np.sin(self.halfTheta)
            if self.h > 0:
                # A wagon wheel grain
                if self.Rp*np.cos(self.epsilonAngle)-self.Ri-((self.Rp*np.sin(self.epsilonAngle))/(np.sin(self.halfTheta))) > 0:
                    self.graintype = "longspokewagonwheel"
                else:
                    self.graintype = "shortspokewagonwheel"
            else:
                print("The wagon wheel generator was called, but geometrically this is a star grain.")
                exit(1)
                
        def beta(self,y):
            return np.arccos((self.Rp*np.sin(self.epsilonAngle))/(y+self.f))
        
        def phi(self,y):
            return np.pi/2 + self.epsilonAngle - self.beta(y)
        
        def burn_perimeter(self,y,S=None):
            if self.phase(y) == 1:
                H = self.Rp*np.sin(self.epsilonAngle)-y-self.f
                x = ((y+self.f)/np.sin(self.halfTheta)) + (H/np.tan(self.halfTheta))
                S1 = (self.Rp + self.f + y)*(np.pi/self.N - self.epsilonAngle)
                S2 = (y+self.f)*(np.pi/2 + self.epsilonAngle)
                S3 = self.Rp*np.cos(self.epsilonAngle)-self.Ri-x
                S4 = H/np.sin(self.halfTheta)
            elif self.phase(y) == 2:
                if self.graintype == "longspokewagonwheel":
                    S1 = (self.Rp+self.f+y)*(np.pi/self.N - self.epsilonAngle)
                    S2 = (y+self.f)*self.phi(y)
                    S3 = 0
                    S4 = 0
                elif self.graintype == "shortspokewagonwheel":
                    alpha = np.arccos(1 - (self.h/(y+self.f)))
                    S1 = (self.Rp + self.f + y)*(np.pi/self.N - self.epsilonAngle)
                    S2 = (y+self.f)*(np.pi/2 + self.epsilonAngle - (self.halfTheta - alpha))
                    S3 = 0
                    S4 = (self.Rp*np.sin(self.epsilonAngle)-(y+self.f)*np.cos(self.halfTheta-alpha))/np.sin(self.halfTheta)
            if S == None:
                return 2*self.N*(S1+S2+S3+S4)
            elif S == 1:
                return S1
            elif S == 2:
                return S2
            elif S == 3:
                return S3
            elif S == 4:
                return S4
            elif S == 'H':
                return H
        
        def Ab(self, y):
            burn_perimeter = self.burn_perimeter(y)
            return burn_perimeter*self.length
        
        def Ap(self, y, S=None):
            if self.phase(y) == 1:
                H = self.burn_perimeter(y,'H')
                S3 = self.burn_perimeter(y,3)
                A1 = 0.5*(self.Rp + self.f + y)**2 * (np.pi/self.N - self.epsilonAngle)
                A2 = 0.5*(np.pi/2 + self.epsilonAngle)*(self.f+y)**2
                A3 = 0.5*self.Rp**2 * np.sin(self.epsilonAngle)*np.cos(self.epsilonAngle)
                A4 = -(H*S3 + 0.5*(H**2/np.tan(self.halfTheta)))
                A5 = 0
            elif self.phase(y) == 2:
                if self.graintype == "longspokewagonwheel":
                    A1 = 0.5*(np.pi/self.N - self.epsilonAngle)*(self.Rp+self.f+y)**2
                    A2 = 0.5*(np.pi/2 + self.epsilonAngle - self.beta(y))*(self.f+y)**2
                    A3 = 0.5*self.Rp*np.sin(self.epsilonAngle)*(self.Rp*np.cos(self.epsilonAngle)+(y+self.f)*np.sin(self.beta(y)))
                    A4 = 0
                    A5 = 0
                elif self.graintype == "shortspokewagonwheel":
                    alpha = np.arccos(1 - (self.h/(y+self.f)))
                    beta = self.halfTheta - alpha
                    phi = np.pi/2 + self.epsilonAngle - beta
                    A1 = 0.5*(np.pi/self.N - self.epsilonAngle)*(self.Rp+self.f+y)**2
                    A2 = 0.5*phi*(self.f+y)**2
                    A3 = 0.5*self.Rp**2 * np.sin(self.epsilonAngle) * np.cos(self.epsilonAngle)
                    A4 = 0.5*np.cos(beta)*np.sin(beta)*(y+self.f)**2
                    S4 = self.burn_perimeter(y,4)
                    A5 = (self.Rp * np.sin(self.epsilonAngle) - (y+self.f)*np.cos(beta))*((y+self.f)*np.sin(beta))-0.5*S4**2 * np.sin(self.halfTheta)*np.cos(self.halfTheta)
            if S == None:
                return 2*self.N*(A1+A2+A3+A4+A5)
            elif S == 1:
                return A1
            elif S == 2:
                return A2
            elif S == 3:
                return A3
            elif S == 4:
                return A4
            elif S == 5:
                return A5
            elif S == 'H':
                return H
        
        def phase(self,y):
            y0 = self.Rp*np.sin(self.epsilonAngle)
            if (y + self.f) < y0:
                return 1
            else:
                return 2
    
    class nozzle():
        """ Rocket nozzle """
        def __init__(self,throat_area=None,throat_diameter=None,exit_area=None,Me=None,P0=None,Pe=None,Pa=None,Te=None,expansion_ratio=None,specific_heat_ratio=None,ue=None):
            self.throat_area = throat_area
            if self.throat_area == None and throat_diameter != None:
                self.throat_area = np.pi*(throat_diameter/2)**2
                self.throat_diameter = throat_diameter
            self.exit_area = exit_area
            self.Me = Me
            self.P0 = P0
            self.Pe = Pe
            self.Pa = Pa
            self.Te = Te
            self.expansion_ratio = expansion_ratio
            self.specific_heat_ratio = specific_heat_ratio
            self.ue = ue