#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from srm import motor

motor = motor(
  motor.propellant(0.0563,0.33,0.064,6260,1.24,5343,25.3),
  motor.cpgrain(8,4,60),
  motor.nozzle(throat_diameter=5,Me=2.5,Pe=14.7,Pa=14.7),
)
timestep = 0.01
Y = motor.burn_vector(timestep)
time = Y[:,0]
thrust = Y[:,3]

from matplotlib import pyplot as plt
fig, ax = plt.subplots(figsize=(7, 7), dpi=96)
ax.plot(time,thrust)
plt.title("Thrust")
plt.xlabel("Burn Time (s)")
plt.ylabel("Thrust (lbf)")
# plt.savefig('thrust.png',dpi=100)
plt.show()

prop_mass = motor.prop_mass(Y[:,1])
from matplotlib import pyplot as plt
fig, ax = plt.subplots(figsize=(7, 7), dpi=96)
ax.plot(time,prop_mass)
plt.title("Propellant Mass")
plt.xlabel("Burn Time (s)")
plt.ylabel("Mass (lbm)")
# plt.savefig('prop_mass.png',dpi=100)
plt.show()
