#!/usr/bin/env python3

import numpy as np
from matplotlib import pyplot as plt
import srm.srm as srm

""" End burning grain example from slide 4:8 """
print("End-burning grain example from slide 4:8\n")

eb = srm.motor(
    propellant=srm.propellant(
        a = 0.040,
        n = 0.4,
        density = 0.064,
        temperature = 6059,
        cstar = 5364,
        MW = 24.3,
        specific_heat_ratio = 1.24,
    ),
    grain=srm.endburninggrain(Ro=3,length=36),
    nozzle=srm.nozzle(
        ambient_pressure = 14.7,
        exit_pressure = 14.7,
        throat_diameter = 1,
    )
)

print("Chamber pressure: {} psi".format(eb.nozzle.chamber_pressure))
print("Mass flow rate: {} lbm/s".format(srm.g0*eb.mass_flow_rate(0)))
print("Exit temperature: {} R".format(eb.nozzle.exit_temperature))
print("Exit velocity: {} ft/s".format(eb.nozzle.exit_velocity))
print("Thrust: {} lbf".format(eb.thrust(0)))
print("Isp: {} s".format(eb.thrust(0)/(srm.g0*eb.mass_flow_rate(0))))
print("Burn time: {} s".format(eb.grain.burn_time))
print("Total impulse: {} lbf-s".format(eb.thrust(0)*eb.grain.burn_time))
print("Exit Mach number: {}".format(eb.nozzle.exit_mach))
print("Exit area: {} in^2\n".format(eb.nozzle.exit_area))

# Printing a 'propellant' object will echo its properties to the console
print(eb.propellant)

timestep = 0.01
Y = eb.burn_vector(timestep)
time = Y[:,0]
burn_distance = Y[:,1]
pressure = Y[:,2]
thrust = Y[:,3]
Isp = Y[:,4]
impulse = np.zeros_like(time)
for i in range(1,len(time)):
    impulse[i] = impulse[i-1] + thrust[i]*timestep
    
fig, ax = plt.subplots(figsize=(7, 7), dpi=96)
ax.plot(time,burn_distance,label="Burn distance")
plt.hlines(eb.grain.Ro,color='red',linestyles=':',xmin=time[0],xmax=time[-1],label="Burnout distance")
plt.title("Burn distance (end-burning grain)")
plt.xlabel("Burn Time (s)")
plt.ylabel("Distance (in)")
plt.legend(loc='lower right')
plt.show()

fig, ax = plt.subplots(figsize=(7, 7), dpi=96)
ax.plot(time,pressure)
plt.title("Chamber pressure (end-burning grain)")
plt.xlabel("Burn Time (s)")
plt.ylabel("Pressure (psi)")
plt.show()

fig, ax = plt.subplots(figsize=(7, 7), dpi=96)
ax.plot(time,thrust)
plt.title("Thrust (end-burning grain)")
plt.xlabel("Burn Time (s)")
plt.ylabel("Thrust (lbf)")
plt.show()

fig, ax = plt.subplots(figsize=(7, 7), dpi=96)
ax.plot(time,Isp)
plt.title("Specific Impulse (end-burning grain)")
plt.xlabel("Burn Time (s)")
plt.ylabel("Isp (s)")
plt.show()

fig, ax = plt.subplots(figsize=(7, 7), dpi=96)
ax.plot(time,impulse)
plt.title("Total Impulse (end-burning grain)")
plt.xlabel("Burn Time (s)")
plt.ylabel("I (lbf-s)")
plt.show()

