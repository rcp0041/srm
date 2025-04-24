#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
Solid Rocket Motor library
Created on Wed Mar 12 22:43:15 2025
@author: Ray Patrick
"""

# TODO: Add bool phase1_neutral to star and WW grains
# TODO: Implement grain type checking at grain rather than motor declaration
# TODO: Implement neutral WW theta table (slide 4:71)
# TODO: Implement slotted, spherical, tubular, and rod-in-tube grains
# TODO: Have a "design_altitude" attribute for the nozzle that queries the standard atmosphere! (And allows to set exit and ambient pressure dynamically)

import numpy as np                    # Numerical operations
from scipy import optimize            # Root finder
from srm.utilities import data

# Global variables:
Ru = 1545
g0 = 32.174
default_timestep = 0.01

def create_dict(obj):
    keys = [x for x in dir(obj) if not (x.startswith('__') or x.startswith('propellant'))]
    values = []
    for key in keys:
        values.append(getattr(obj,key))
    dictionary = {k: v for k, v in zip(keys, values)}
    return dictionary

def print_dict(dictionary):
    a = ""
    for i in dictionary:
        a = a + f"{i}: {dictionary[i]}\n"
    return a

class dummy_updater(object):
    def __init__(self, iterable=(), **kwargs):
        self.__dict__.update(iterable, **kwargs)

class motor:
    def __init__(self,propellant,grain,nozzle,case=None):
        self.propellant = propellant
        self.nozzle = nozzle
        # self.case = case

        # Do star/wagon wheel check before assigning self.grain        
        if grain.graintype == "star" or grain.graintype == "wagonwheel":
            # Perform spoke-adjacency check
            grain.delta = np.arctan(
                (grain.Rp*np.sin(grain.epsilonAngle)-grain.f)/
                (grain.Ri + (grain.f/np.sin(grain.halfTheta)) + (grain.Rp*np.sin(grain.epsilonAngle)-grain.f)/np.tan(grain.halfTheta))
                )
            if grain.delta >= np.pi/grain.N:
                # print("Maximum spoke length exceeded. Recheck your design parameters!")
                grain.spoke_collision = True
                # exit(1)
            # Detect star (h < 0) or wagon wheel (h > 0)
            grain.h = (grain.Rp*np.cos(grain.epsilonAngle)-((grain.Rp*np.sin(grain.epsilonAngle)/(np.tan(grain.halfTheta))))-grain.Ri)*np.sin(grain.halfTheta)
            if grain.h > 0:
                # Detect type of wagon wheel
                if grain.Rp*np.cos(grain.epsilonAngle)-grain.Ri-((grain.Rp*np.sin(grain.epsilonAngle))/(np.sin(grain.halfTheta))) > 0:
                    actualgraintype = "longspokewagonwheel"
                else:
                    actualgraintype = "shortspokewagonwheel"       
            else:
                actualgraintype = "star"
        # Now change the graintype (if required)
            if grain.graintype != actualgraintype:
                if grain.graintype == "wagonwheel" and (actualgraintype == "shortspokewagonwheel" or actualgraintype == "longspokewagonwheel"):
                    pass
                else:
                    # print("You defined it as a {}, but this grain is a {}.".format(grain.graintype,actualgraintype))
                    # print("Reinitializing the grain as a {}.".format(actualgraintype))
                    pass
                if actualgraintype == "star":
                    newgrain = stargrain(grain.N,grain.Ro,grain.Ri,grain.Rp,grain.f,grain.epsilon,grain.halfTheta,grain.length)
                elif actualgraintype == "longspokewagonwheel" or actualgraintype == "shortspokewagonwheel":
                    newgrain = wagonwheelgrain(grain.N,grain.Ro,grain.Ri,grain.Rp,grain.f,grain.epsilon,grain.halfTheta,grain.length)
                    newgrain.graintype = actualgraintype
                self.grain = newgrain
            elif grain.graintype == actualgraintype:
                self.grain = grain
        else:
            self.grain = grain
        
        if self.grain.graintype == "endburning":
            if hasattr(self.nozzle,'throat_area') and not hasattr(self.nozzle,'chamber_pressure'):
                self.nozzle.chamber_pressure = self.pressure(0)
                self.nozzle.__recalculate__()
            self.grain.burn_rate = self.propellant.a * self.nozzle.chamber_pressure ** self.propellant.n
            self.grain.burn_time = self.grain.length / self.grain.burn_rate
        if self.grain.graintype == "CP":
            self.grain.n = self.propellant.n
            # Burn rate
            self.grain.K1 = (2*np.pi*self.grain.length*self.propellant.a*((self.propellant.a*self.propellant.cstar*self.propellant.density)/(g0*self.nozzle.throat_area))**(self.propellant.n/(1-self.propellant.n)))
            self.grain.K2 = self.grain.K1*((2*self.propellant.n-1)/(self.propellant.n-1))
            self.grain.Abi = 2*np.pi*self.grain.Ri*self.grain.length
            self.grain.Abf = 2*np.pi*self.grain.Ro*self.grain.length
            self.grain.p = (2*self.propellant.n - 1)/(self.propellant.n-1)
            if self.propellant.n != 0.5:
                self.grain.q = 1/self.grain.p
                self.grain.burn_time = (self.grain.Abf**((2*self.propellant.n-1)/(self.propellant.n-1)) - self.grain.Abi**((2*self.propellant.n-1)/(self.propellant.n-1)))/self.grain.K2
            elif self.propellant.n == 0.5:
                self.grain.burn_time = np.log(self.grain.Abf)/(np.log(self.grain.Abi)*self.grain.K1)
        
        """ Case stuff """
        case.wall_radius = self.grain.Ro
        self.case = case
        
        """ Nozzle stuff """
        if hasattr(self.propellant,'specific_heat_ratio'):
            self.nozzle.k = self.propellant.specific_heat_ratio
            self.nozzle.__recalculate__()
        
        if hasattr(self.propellant,'temperature'):
            self.nozzle.__recalculate__(chamber_temperature=self.propellant.temperature)
            self.nozzle.__recalculate__()
            
        if hasattr(self.nozzle,'throat_area') and not hasattr(self.nozzle,'chamber_pressure'):
            self.nozzle.chamber_pressure = self.pressure(0)
            self.nozzle.__recalculate__()
            
        if hasattr(self.nozzle,'exit_temperature') and hasattr(self.nozzle,'chamber_temperature') and hasattr(self.propellant,'R'):
            # print(f"gamma = {self.propellant.specific_heat_ratio}, sqrt(2/gamma-1) = sqrt(2/{self.propellant.specific_heat_ratio-1})")
            self.nozzle.exit_velocity = np.sqrt(2 * (self.propellant.specific_heat_ratio)/(self.propellant.specific_heat_ratio-1) * self.propellant.R * g0 * (self.nozzle.chamber_temperature - self.nozzle.exit_temperature))
            self.nozzle.__recalculate__()
        
    def pressure(self,x):
        """ For CP grains, 'x' is time. For star grains, 'x' is burn distance."""
        return self.grain.pressure(x,self.nozzle,self.propellant)
    
    def mass_flow_rate(self,x):
        """Slug/s. For CP grains, 'x' is time. For star grains, 'x' is burn distance."""
        return (self.pressure(x)*self.nozzle.throat_area)/self.propellant.cstar
    
    def propellant_mass(self,x):
        """ 'x' is burn distance for all grain types """
        return ((np.pi*self.grain.Ro**2) - self.grain.port_area(x))*self.grain.length*self.propellant.density
    
    def burn_rate(self,x):
        """ For star grains, 'x' is burn distance """
        return self.propellant.a * self.pressure(x)**self.propellant.n
    
    def thrust(self,x):
        """ For CP grains, 'x' is time. For star grains, 'x' is burn distance."""
        if hasattr(self.nozzle,'exit_area'):
            return self.mass_flow_rate(x)*self.nozzle.exit_velocity + self.nozzle.exit_area*(self.nozzle.exit_pressure - self.nozzle.ambient_pressure)
        elif hasattr(self.nozzle,'exit_pressure') and hasattr(self.nozzle,'ambient_pressure') and hasattr(self.nozzle,'exit_velocity'):
            if self.nozzle.exit_pressure == self.nozzle.ambient_pressure:
                return self.mass_flow_rate(x)*self.nozzle.exit_velocity
    
    def specific_impulse(self,x):
        """ For CP grains, 'x' is time. For star grains, 'x' is burn distance."""
        return self.thrust(x)/(self.mass_flow_rate(x)*g0)
    
    def burn_vector(self,timestep,initial_time=0):
        """
        Returns a state vector of the form:
        [t, y, pressure, thrust, specific_impulse]
        """
        Y = np.array([initial_time,0,self.pressure(0),self.thrust(0),self.specific_impulse(0)])
        t,y,r = initial_time,0,self.burn_rate(0)
        if self.grain.graintype == "CP":
            yf = self.grain.Ro - self.grain.Ri
            while t < self.grain.burn_time:
                t = t + timestep
                y = y + self.burn_rate(t)*timestep
                r = self.burn_rate(t)
                Y = np.vstack((Y,np.array([t,y,self.pressure(t),self.thrust(t),self.specific_impulse(t)])))    
        elif self.grain.graintype != "CP":
            if self.grain.graintype == "endburning":
                yf = self.grain.Ro
            else:
                yf = self.grain.Ro - self.grain.Rp - self.grain.f
            while y < yf:
                t = t + timestep
                y = y + r*timestep
                r = self.burn_rate(y)
                Y = np.vstack((Y,np.array([t,y,self.pressure(y),self.thrust(y),self.specific_impulse(y)])))
        return Y
    
    def max_pressure(self,timestep=default_timestep):
        """ Computes peak chamber pressure """
        Y = self.burn_vector(timestep)
        """ Test for IndexError and print diagnostic """
        if Y.ndim == 1:
            # print("Burn vector is 1D! I wonder why?")
            # print(self)
            return 0
        else:
            return max(Y[:,2])
    
    def avg_pressure(self,timestep=default_timestep):
        """ Computes average chamber pressure """
        Y = self.burn_vector(timestep)
        if Y.ndim == 1:
            # print("Burn vector is 1D! I wonder why?")
            # print(self)
            return 0
        else:
            return np.average(Y[:,2])

    def max_thrust(self,timestep=default_timestep):
        """ Computes peak thrust """
        Y = self.burn_vector(timestep)
        if Y.ndim == 1:
            # print("Burn vector is 1D! I wonder why?")
            # print(self)
            return 0
        else:
            return max(Y[:,3])
    
    def min_thrust(self,timestep=default_timestep):
        """ Computes minimum thrust """
        Y = self.burn_vector(timestep)
        if Y.ndim == 1:
            # print("Burn vector is 1D! I wonder why?")
            # print(self)
            return 0
        else:
            return min(Y[:,3])

    def avg_thrust(self,timestep=default_timestep):
        """ Computes average thrust """
        Y = self.burn_vector(timestep)
        if Y.ndim == 1:
            # print("Burn vector is 1D! I wonder why?")
            # print(self)
            return 0
        else:
            return np.average(Y[:,3])
    
    def payload_mass(self):
        """ Computes payload mass (M_inert - M_case) for a given rocket motor """
        # Get case properties
        SF = self.case.safety_factor
        rho = self.case.density
        sigma = self.case.yield_strength
        MF = self.case.mass_fraction
        
        # Get motor properties
        Rc = self.grain.Ro
        L = self.grain.length
        Pc = self.max_pressure()
        Dt = self.nozzle.throat_diameter
        
        t = (SF*Rc*Pc)/sigma
        self.case.wall_thickness = t
        M_cyl = 2*np.pi*Rc*t*L*rho
        t_end = SF*(Rc*Pc)/(2*sigma)
        self.case.end_thickness = t_end
        M_nose = 2*np.pi*Rc**2 * t_end * rho
        phi = np.arctan(Dt/(2*Rc))
        M_aft = M_nose*np.cos(phi)
        t_nozzle = SF*((Dt*Pc)/(2*sigma))
        Rt = Dt/2
        Re = Rt*np.sqrt(self.nozzle.expansion_ratio)
        L_nozzle = (Re-Rt)/np.tan(phi)
        M_nozzle = rho*((np.pi*L_nozzle)/3) * (((Re+t_nozzle)**2 + (Re+t_nozzle)*(Rt+t_nozzle) + (Rt+t_nozzle)**2)-(Re**2+Re*Rt+Rt**2))
        M_inert = ((1/MF) - 1)*self.propellant_mass(0)
        M_case = M_cyl + M_nose + M_aft + M_nozzle
        # M_payload = M_inert - M_cyl - M_nose - M_aft - M_nozzle
        M_payload = M_inert - M_case
        self.dry_mass = M_inert + M_case
        self.case.mass = M_case
        return M_payload
    
    def plot(self, items: list):
        motordata = data(self)
        motordata.plot(items)
    
    def __str__(self):
        a = self.propellant.__str__() + "\n"
        a = a + self.grain.__str__() + "\n\n"
        a = a + self.nozzle.__str__() + "\n"
        # self.propellant.__str__()
        # self.grain.__str__()
        # self.nozzle.__str__()
        return a

class nozzle(dummy_updater):
    def __init__(self, iterable=(), **kwargs):
        super().__init__(iterable, **kwargs)
        
        # Compute expansion ratio if possible
        if hasattr(self, 'exit_diameter') and hasattr(self, 'throat_diameter') and not hasattr(self,'expansion_ratio'):
            self.exit_area = (np.pi/4)*self.exit_diameter**2
            self.throat_area = (np.pi/4)*self.throat_diameter**2
            self.expansion_ratio = self.exit_area/self.throat_area
        
        # Compute pressure ratio if possible
        if hasattr(self,'exit_pressure') and hasattr(self,'chamber_pressure'):
            self.pressure_ratio = self.exit_pressure/self.chamber_pressure
            
        # Compute throat area if possible (assumes circular cross-section)
        if not hasattr(self,'throat_area'):
            if hasattr(self,'throat_diameter'):
                self.throat_area = np.pi * (self.throat_diameter/2)**2
            elif hasattr(self,'propellant') and hasattr(self,'thrust') and hasattr(self,'chamber_pressure') and hasattr(self,'exit_pressure'):
                self.thrust_coefficient = self.propellant.cap_gamma *np.sqrt(((2*self.propellant.specific_heat_ratio)/(self.propellant.specific_heat_ratio-1))*(1-(self.exit_pressure/self.chamber_pressure)**((self.propellant.specific_heat_ratio-1)/self.propellant.specific_heat_ratio)))
                self.throat_area = self.thrust/(self.chamber_pressure * self.thrust_coefficient)
                self.throat_diameter = np.sqrt((4*self.throat_area)/np.pi)
          
        # Compute exit temperature if possible
        if hasattr(self,'chamber_temperature') and hasattr(self,'pressure_ratio') and hasattr(self,'k'):
            self.exit_temperature = self.chamber_temperature * (self.pressure_ratio)**((self.k-1)/self.k)
        
        # Compute exit Mach number if possible
        if hasattr(self,'k') and hasattr(self,'pressure_ratio'):
            # print(f"k = {self.k}, sqrt(2/k-1) = sqrt(2/{self.k-1})")
            self.exit_mach = np.sqrt((2/(self.k-1)) * (((1/self.pressure_ratio) ** ((self.k-1)/self.k))-1))
        
        # Compute exit area if possible
        if hasattr(self,'throat_area') and hasattr(self,'exit_mach') and hasattr(self,'k'):
            self.exit_area = (self.throat_area/self.exit_mach)*((2/(self.k+1))*(1+(0.5*(self.k-1))*self.exit_mach**2))**((self.k+1)/(2*(self.k-1)))
            
    def __str__(self):
        a = "Nozzle properties:\n"
        a = a + print_dict(create_dict(self))
        return a
    
    def __recalculate__(self,**kwargs):
        """ Add nozzle parameters, then recompute all nozzle properties """
        for k, val in kwargs.items():
            setattr(self,k,val)
        self.__init__(create_dict(self))

class case:
    def __init__(self,material,safety_factor,mass_fraction):
        self.material = material
        self.density = self.material.density
        self.yield_strength = self.material.yield_strength
        self.safety_factor = safety_factor
        self.mass_fraction = mass_fraction

    def yield_pressure(self):
        return self.yield_strength*(self.wall_thickness/self.wall_radius)

    def __str__(self):
        a = "Case properties:\n"
        if hasattr(self.material,'name'):
            a = a + f"Material: {self.material.name}\n"
        a = a + f"Density: {self.density} lbm/in^3\n"
        a = a + f"Yield strength: {self.yield_strength} psi\n"
        a = a + f"Safety factor: {self.safety_factor}\n"
        a = a + f"Mass fraction: {self.mass_fraction}\n"
        if hasattr(self,'wall_thickness'):
            a = a + f"Wall thickness: {self.wall_thickness} in\n"
        if hasattr(self,'end_thickness'):
            a = a + f"End thickness: {self.end_thickness} in\n"
        return a

class propellant:
    """ Solid rocket propellant """
    def __init__(self,a,n,density,temperature,specific_heat_ratio,cstar,MW,name=None):
        self.name = name
        self.a = a
        self.n = n
        self.density = density
        self.temperature = temperature
        self.specific_heat_ratio = specific_heat_ratio
        self.cap_gamma = np.sqrt(self.specific_heat_ratio*(2/(self.specific_heat_ratio+1))**((self.specific_heat_ratio+1)/(self.specific_heat_ratio-1)))
        self.cstar = cstar
        self.MW = MW
        # Compute specific gas constant from molecular weight or vice versa
        if hasattr(self,'MW') and not hasattr(self,'R'):
            self.R = Ru/self.MW
        if hasattr(self,'R') and not hasattr(self,'MW'):
            self.MW = Ru/self.R
        # Now compute specific heats
        self.cv = self.R/(self.specific_heat_ratio - 1)
        self.cp = self.cv * self.specific_heat_ratio
    
    def __str__(self):
        a = "Propellant properties:\n"
        if self.name != None:
            a = a + f"Name: {self.name}\n"
        a = a + f"a: {self.a}\n"
        a = a + f"n: {self.n}\n"
        a = a + f"density: {self.density}\n"
        a = a + f"Combustion temperature: {self.temperature}\n"
        a = a + f"Specific heat ratio, gamma: {self.specific_heat_ratio}\n"
        a = a + f"'Cap Gamma': {self.cap_gamma}\n"
        a = a + f"Characteristic velocity, c*: {self.cstar}\n"
        a = a + f"Molecular weight: {self.MW}\n"
        a = a + f"Specific gas constant, R: {self.R}\n"
        a = a + f"cv: {self.cv}\n"
        a = a + f"cp: {self.cp}\n"
        return a

class grain():
    """ Solid rocket fuel grain """
    def __init__(self):
        pass
    
    def pressure(self,x,nozzle,propellant):
        """ For CP grains, 'x' is time. For star grains, 'x' is burn distance."""
        if self.graintype == "CP":
            return (((self.burn_area(x)/nozzle.throat_area)*propellant.a*propellant.density*propellant.cstar)/g0)**(1/(1-propellant.n))
        else:
            return (((self.burn_area(x)/nozzle.throat_area)*propellant.a*propellant.density*propellant.cstar)/g0)**(1/(1-propellant.n))
    
    def __str__(self):
        a = "Grain properties:\n"
        if self.graintype == "shortspokewagonwheel" or self.graintype == "longspokewagonwheel" or self.graintype == "star":
            a = a + f"Type: {self.graintype}\n"
            a = a + f"Number of points, N: {self.N}\n"
            a = a + f"Outer radius, Ro: {self.Ro}\n"
            a = a + f"Inner radius, Ri: {self.Ri}\n"
            a = a + f"Radius to curvature center, Rp: {self.Rp}\n"
            a = a + f"Fillet radius, f: {self.f}\n"
            a = a + f"Epsilon: {self.epsilon}\n"
            a = a + f"Theta/2: {self.halfTheta} rad\n"
            a = a + f"Length: {self.length}\n"
            a = a + f"Phase I burning: {self.neutrality} (neutrality coefficient: {self.neutrality_coefficient})\n"
        if self.graintype == "endburning":
            a = a + f"Type: {self.graintype}\n"
            a = a + f"Outer radius, Ro: {self.Ro}\n"
            a = a + f"Length: {self.length}"
        if self.graintype == "CP":
            a = a + f"Type = {self.graintype}\n"
            a = a + f"Outer radius, Ro = {self.Ro}\n"
            a = a + f"Inner radius, Ri = {self.Ri}\n"
            a = a + f"Length = {self.length}\n"
        return a
    
class endburninggrain(grain):
    """ End-burning grain """
    def __init__(self,Ro,length):
        self.Ro = Ro
        self.length = length
        self.graintype = "endburning"
    
    def burn_area(self, x):
        return np.pi*self.Ro**2
    
    def port_area(self, x):
        return 0

class cpgrain(grain):
    """ Circular-perforated grain """
    def __init__(self,Ro,Ri,length):
        self.Ro = Ro
        self.Ri = Ri
        self.length = length
        self.graintype = "CP"
    
    def burn_area(self,t):
        if self.n != 0.5:
            return (self.K2*t + self.Abi**self.p)**self.q
        elif self.n == 0.5:
            return self.Abi * np.exp(self.K1*t)
    
    def port_area(self,y):
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
        self.spoke_collision = False
        
        self.epsilonAngle = np.pi*self.epsilon/self.N
        self.H1 = self.Rp*np.sin(self.epsilonAngle)
        
        if self.Ri == None:
            self.Ri = self.find_neutral_Ri()
            
        if self.halfTheta == None:
            self.halfTheta = np.arctan((self.H1*np.tan(self.epsilonAngle))/(self.H1 - self.Ri*np.tan(self.epsilonAngle)))
        
        self.web1 = (self.Rp*np.sin(self.epsilonAngle)/np.cos(self.halfTheta))-self.f
        
        self.beta = np.pi/2 - self.halfTheta + self.epsilonAngle
    
        self.neutrality_coefficient = np.pi/2 - self.halfTheta + np.pi/self.N - 1/np.tan(self.halfTheta)
        if self.neutrality_coefficient < 0:
            self.neutrality = 'Regressive'
        # elif self.neutrality_coefficient >= 1e-17:
        elif self.neutrality_coefficient == 0:
            self.neutrality = 'Neutral'
        else:
            self.neutrality = 'Progressive'
    
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
            
    def burn_area(self, y):
        burn_perimeter = self.burn_perimeter(y)
        return burn_perimeter*self.length
    
    def port_area(self, y):
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
        self.h = (self.Rp*np.cos(self.epsilonAngle)-((self.Rp*np.sin(self.epsilonAngle))/(np.tan(self.halfTheta)))-self.Ri)*np.sin(self.halfTheta)
        self.graintype = "wagonwheel"
        self.spoke_collision = False
        
        self.neutrality_coefficient = np.pi/self.N + np.pi/2 - (2/np.sin(self.halfTheta)) + 1/np.tan(self.halfTheta)
        if self.neutrality_coefficient < 0:
            self.neutrality = 'Regressive'
        # elif self.neutrality_coefficient >= 1e-17:
        elif self.neutrality_coefficient == 0:
            self.neutrality = 'Neutral'
        else:
            self.neutrality = 'Progressive'
            
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
    
    def burn_area(self, y):
        burn_perimeter = self.burn_perimeter(y)
        return burn_perimeter*self.length
    
    def port_area(self, y, S=None):
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

class tubulargrain(grain):
    def __init__(self,Rc,Ro,Ri,length):
        self.Rc = Rc
        self.Ro = Ro # Ra in slides
        self.Ri = Ri # Rb in slides
        self.length = length
        self.graintype = "tubular"
        
        def burn_area(self,y):
            #return 2*np.pi*(self.Ro-self.Ri)*self.length
            return 2*np.pi*(self.Ro - y)*self.length + 2*np.pi*(self.Ri + y)*self.length
        
        def port_area(self,y):
            inner_disk = np.pi * (self.Ri + y)**2
            outer_annulus = np.pi*self.Rc**2 - np.pi*(self.Ro - y)**2
            return inner_disk + outer_annulus