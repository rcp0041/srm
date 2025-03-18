#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from srm.srm import motor

motor = motor(
  motor.propellant(a=0.0563,                  # (in/s) / psi^n
                   n=0.33,                    # Burn index
                   density=0.064,             # lbm/in^3
                   temperature=6260,          # Combustion temperature (R)
                   specific_heat_ratio=1.24,  # "Gamma" or "k"
                   cstar=5343,                # Characteristic velocity
                   MW=25.3                    # Molecular weight
  ),
  motor.cpgrain(Ro=8,                         # Outer radius (in)
                Ri=4,                         # Inner (perforation) radius (in)
                length=60                     # Length of grain (in)
  ),
  motor.nozzle(throat_diameter=5,
               Me=2.5,                        # Exit Mach number
               Pe=14.7,                       # Exit pressure (psi)
               Pa=14.7                        # Ambient pressure (psi)
  ),
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
# plt.savefig('thrust.png',dpi=50)
plt.show()

# prop_mass = motor.prop_mass(Y[:,1])
# from matplotlib import pyplot as plt
# fig, ax = plt.subplots(figsize=(7, 7), dpi=96)
# ax.plot(time,prop_mass)
# plt.title("Propellant Mass")
# plt.xlabel("Burn Time (s)")
# plt.ylabel("Mass (lbm)")
# # plt.savefig('prop_mass.png',dpi=50)
# plt.show()
