#!/usr/bin/env python3

import numpy as np
from matplotlib import pyplot as plt
import srm.srm as srm

""" Tubular grain example from slide 4:106 """

# First define a propellant object

the_propellant = srm.propellant(
    a = 0.0563,
    n = 0.33,
    cstar = 5343,
    density = 0.064,
    temperature = 6260,
    MW = 25.3,
    specific_heat_ratio = 1.24,
)

# The 'nozzle' class can accept a propellant as an argument. The nozzle 
# __init__ method will then have access to the gas properties, with which it
# can follow the math on slide 4:107 to deduce the thrust coefficient and
# therefore the throat area. (You need to provide design thrust as well as
# chamber pressure and exit pressure)

the_nozzle = srm.nozzle(
    thrust = 3000,
    chamber_pressure = 600,
    exit_pressure = 14.7,
    ambient_pressure=14.7,
    propellant=the_propellant,
)

# Now for the tubular grain as on slide 4:106

the_grain = srm.tubulargrain(
    Ro = 2.499,
    Rc = 3,
    Ri = 1.104,
    length = 17.75
)

# Now assemble these pieces into a "motor" object

motor = srm.motor(
    propellant = the_propellant,
    grain = the_grain,
    nozzle = the_nozzle
)

# Motor, propellant, grain, and nozzle objects list their attributes to the
# console if you print them:
    
print(motor)