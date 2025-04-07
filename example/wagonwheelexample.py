#!/usr/bin/env python3

import numpy as np
from matplotlib import pyplot as plt
import srm.srm as srm

# Declaring a star grain with wagon-wheel type parameters will result in an
# error message and the correct grain type will be automatically reloaded:
rocket = srm.motor(
    srm.propellant( a=0.0563,n=0.33,density=0.064,temperature=6260,specific_heat_ratio=1.24,cstar=5343,MW=25.3),
    srm.stargrain(
        N = 7,
        Ro = 6,
        Ri = 1.4,
        Rp = 3,
        f = 0.4,
        epsilon = 0.8,
        halfTheta = np.deg2rad(50),
        length = 36
    ),
    srm.nozzle(
        ambient_pressure = 14.7,
        exit_pressure = 14.7,
        throat_diameter = 5,
        exit_mach = 2.5
    )
)
print(rocket.grain)

timestep = 0.01
Y = rocket.burn_vector(timestep)
time = Y[:,0]
burn_distance = Y[:,1]
pressure = Y[:,2]
thrust = Y[:,3]
# Isp = Y[:,4]
mdot = np.zeros_like(time)
mdot[0] = rocket.mass_flow_rate(0)
impulse = np.zeros_like(time)
burn_area = np.zeros_like(time)
Isp = np.zeros_like(time)
Isp[0] = thrust[0]/(rocket.mass_flow_rate(0)*srm.g0)
burn_area[0] = rocket.grain.burn_area(0)
for i in range(1,len(time)):
    impulse[i] = impulse[i-1] + thrust[i]*timestep
    burn_area[i] = rocket.grain.burn_area(time[i])
    mdot[i] = rocket.mass_flow_rate(time[i])
    Isp[i] = thrust[i]/(mdot[i]*srm.g0)

fig, ax = plt.subplots(figsize=(7, 7), dpi=96)
ax.plot(time,burn_distance,label="Burn distance")
plt.hlines(rocket.grain.Ro-rocket.grain.Rp-rocket.grain.f,color='red',linestyles=':',xmin=time[0],xmax=time[-1],label="Burnout distance")
plt.title("Burn distance (CP grain)")
plt.xlabel("Burn Time (s)")
plt.ylabel("Distance (in)")
plt.legend(loc='lower right')
plt.show()

fig, ax = plt.subplots(figsize=(7, 7), dpi=96)
ax.plot(time,burn_area,label="Burn area")
# plt.hlines(cp.grain.Ro-cp.grain.Ri,color='red',linestyles=':',xmin=time[0],xmax=time[-1],label="Burnout distance")
plt.title("Burn area (CP grain)")
plt.xlabel("Burn Time (s)")
plt.ylabel("Area (in^2)")
# plt.legend(loc='lower right')
plt.show()

fig, ax = plt.subplots(figsize=(7, 7), dpi=96)
ax.plot(time,pressure)
plt.title("Chamber pressure (CP grain)")
plt.xlabel("Burn Time (s)")
plt.ylabel("Pressure (psi)")
plt.show()

fig, ax = plt.subplots(figsize=(7, 7), dpi=96)
ax.plot(time,thrust)
plt.title("Thrust (CP grain)")
plt.xlabel("Burn Time (s)")
plt.ylabel("Thrust (lbf)")
plt.show()

fig, ax = plt.subplots(figsize=(7, 7), dpi=96)
ax.plot(time,Isp)
plt.title("Specific Impulse (CP grain)")
plt.xlabel("Burn Time (s)")
plt.ylabel("Isp (s)")
plt.show()

fig, ax = plt.subplots(figsize=(7, 7), dpi=96)
ax.plot(time,impulse)
plt.title("Total Impulse (CP grain)")
plt.xlabel("Burn Time (s)")
plt.ylabel("I (lbf-s)")
plt.show()

fig, ax = plt.subplots(figsize=(7, 7), dpi=96)
ax.plot(time,mdot)
plt.title("Mass flow rate (CP grain)")
plt.xlabel("Burn Time (s)")
plt.ylabel("Mdot (slug/s)")
plt.show()