#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
GENOME:
    [N, Ro, Ri, Rp, f, epsilon, theta/2, grain_length, prop_i, material_i]
"""

import numpy as np                     # Math operations
import srm                             # Solid rocket motor performance library
from leap_ec.simple import ea_solve    # Evolutionary algorithm
import toml                            # For reading config file

# Create lists that allow numeric indexing of propellant/material dictionaries
propellants_list = list(srm.materials.propellants)
steel_list = list(srm.materials.steel)
num_propellants = len(propellants_list) - 1
num_steel = len(steel_list) - 1

def max_pressure(motor,timestep=0.1):
    """ Computes maximum chamber pressure for a given rocket motor """
    Y = motor.burn_vector(timestep)
    return max(Y[:,2])

def max_thrust(motor,timestep=0.1):
    """ Computes peak thrust for a given rocket motor """
    Y = motor.burn_vector(timestep)
    return max(Y[:,3])

def avg_thrust(motor,timestep=0.1):
    """ Computes average thrust for a given rocket motor """
    Y = motor.burn_vector(timestep)
    return np.average(Y[:,3])

def payload_mass(case,motor):
    """ Computes payload mass (M_inert - M_case) for a given rocket motor """
    # Get case properties
    SF = case.safety_factor
    rho = case.density
    sigma = case.yield_strength
    MF = case.mass_fraction
    
    # Get motor properties
    Rc = motor.grain.Ro
    L = motor.grain.length
    Pc = max_pressure(motor)
    Dt = motor.nozzle.throat_diameter
    
    t = (SF*Rc*Pc)/sigma
    M_cyl = 2*np.pi*Rc*t*L*rho
    t_end = SF*(Rc*Pc)/(2*sigma)
    M_nose = 2*np.pi*Rc**2 * t_end * rho
    phi = np.arctan(Dt/(2*Rc))
    M_aft = M_nose*np.cos(phi)
    t_nozzle = SF*((Dt*Pc)/(2*sigma))
    Rt = Dt/2
    Re = Rt*np.sqrt(motor.nozzle.expansion_ratio)
    L_nozzle = (Re-Rt)/np.tan(phi)
    M_nozzle = rho*((np.pi*L_nozzle)/3) * (((Re+t_nozzle)**2 + (Re+t_nozzle)*(Rt+t_nozzle) + (Rt+t_nozzle)**2)-(Re**2+Re*Rt+Rt**2))
    M_inert = ((1/MF) - 1)*motor.propellant_mass(0)
    M_payload = M_inert - M_cyl - M_nose - M_aft - M_nozzle
    return M_payload

def construct_motor(genome):
    """ Constructs a motor object from a given genome """
    N = int(genome[0])
    Ro = genome[1]
    Ri = genome[2]
    Rp = genome[3]
    f = genome[4]
    epsilon = genome[5]
    halfTheta = genome[6]
    length = genome[7]
    propellant_index = int(genome[8])
    propellant = srm.materials.propellants[propellants_list[propellant_index]]
    grain = srm.stargrain(N,Ro,Ri,Rp,f,epsilon,halfTheta,length)
    nozzle = srm.nozzle(expansion_ratio=8,
                        ambient_pressure=14.7,
                        exit_pressure=14.7,
                        throat_diameter=2.0)
    return srm.motor(propellant,grain,nozzle)

def construct_case(genome):
    """ Constructs a case object from a given genome """
    material_index = int(genome[9])
    return srm.case(material=srm.materials.steel[steel_list[material_index]],safety_factor=1.3,mass_fraction=0.85)

class fitness:
    def __init__(self):
        pass
    
    def payload_mass(genome):
        return payload_mass(construct_case(genome),construct_motor(genome))

    def max_pressure(genome):
        motor = construct_motor(genome)
        return max_pressure(motor)
    
    def avg_thrust(genome):
        motor = construct_motor(genome)
        return avg_thrust(motor)
    
    def max_thrust(genome):
        motor = construct_motor(genome)
        return max_thrust(motor)

def print_results(results):
    motor = construct_motor(results) # motor.propellant.name
    case = construct_case(results) # case.material.name
    M_payload = payload_mass(case,motor)
    print("Payload mass: {} lbm\nCase material: {}\nPropellant: {}\n{}".format(M_payload,case.material.name,motor.propellant.name,motor.grain))
    
def evolve(fitness_function,generations=50,pop_size=20,maximize=True):
    data = toml.load("./config.toml")
    def get_bounds(key_string):
        key_string = key_string + "_bounds"
        return (data[key_string][0],data[key_string][1])

    N_bounds = get_bounds('N')
    Ro_bounds = get_bounds('Ro')
    Ri_bounds = get_bounds('Ri')
    Rp_bounds = get_bounds('Rp')
    f_bounds = get_bounds('f')
    epsilon_bounds = get_bounds('epsilon')
    halfTheta_bounds = get_bounds('halfTheta')
    length_bounds = get_bounds('length')
    propellant_index_bounds = (0,num_propellants)
    material_index_bounds = (0,num_steel)
    
    bounds = [N_bounds,Ro_bounds,Ri_bounds,Rp_bounds,f_bounds,epsilon_bounds,
              halfTheta_bounds,length_bounds,
              propellant_index_bounds,material_index_bounds]
    
    results = ea_solve(fitness_function,
                       generations=generations,
                       pop_size=pop_size,
                       viz=True,
                       bounds=bounds,
                       hard_bounds=True,
                       maximize=maximize)
    print_results(results)
    return [construct_motor(results),construct_case(results)]

""" Examples """
evolve(fitness.payload_mass,50,20) # Evolve for maximum payload mass
# evolve(fitness.max_pressure,50,20,maximize=False) # Evolve for lowest max pressure
# evolve(fitness.avg_thrust,50,20) # Evolve for maximum average thrust
# evolve(fitness.max_thrust,50,20) # Evolve for maximum peak thrust
