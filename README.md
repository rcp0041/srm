# srm
Tool for modeling solid rocket motor performance.

## Getting Started

1. Clone the repository and navigate to the root directory.
   1. `git clone git@github.com:rcp0041/srm.git`
2. (Recommended) Create a virtual environment for the project and activate it.
   1. `python3 -m venv .venv`
   2. `source .venv/bin/activate` (Linux) or
   `.venv\Scripts\activate` (Windows)
3. Install the project requirements.
   1. `pip install -e .`
4. Run the example script (and then read below to see how it works):
   1. `python example/cpgrain.py`

***

### Example Problem
Consider a rocket with the propellant AP/PBAN having the following characteristics:
- a: 0.0563 (in/s)/psi^n
- n: 0.33
- c* = 5343 ft/s
- density: 0.064 lbm/in^3
- T: 6260 R
- MW: 25.3
- k: 1.24

The rocket has a cylindrically-perforated grain with the following characteristics:
- Grain performation diameter: 8 in
- Outer grain diameter: 16 in
- Grain length: 60 in

The nozzle has the following characteristics:
- Exit Mach number: 2.5
- Throat diameter: 5 in
- Design altitude: Sea level

Create a thrust-time plot for this motor.

The `motor` class is a composite object consisting of a propellant, a grain, and a nozzle. You can define them all at once:
```
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
```
The `motor` and `grain` objects have built-in functions (methods) you can call directly. For instance, you can query the `grain` object for the burnout time:
```
> motor.grain.tb
8.95849341578629
```
Or you can query the `motor` object for its initial thrust:
```
> motor.thrust(0)
8211.251083700156
```
For plotting motor performance, you can request a *state vector* with the `motor.burn_vector()` method. This method accepts a timestep and returns a state vector of the form `[ time, burn distance, pressure, thrust, specific impulse]`. You can then use `matplotlib` or your preferred plotting library to display the results:

```
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
```
![A thrust-time plot generated with matplotlib](https://github.com/rcp0041/srm/blob/main/example/thrust.png?raw=true)

***

## Project Goals

A modular object-based approach to solid rocket performance modeling. Propellants, grain designs, and nozzles can be defined as discrete objects, combined, and analyzed.

## To do

- [ ] Integrate [AstroPy](https://docs.astropy.org/en/stable/units/) for handling dimensional quantities. (The current codebase only works with US units.)
- [ ] Provide a standard library of propellants.
- [ ] Create full documentation on functions.

### Expected Outputs

- Time
- Burn distance
- Burn area
- Chamber pressure
- Thrust
- Total impulse
- Specific impulse
- Propellant mass

## Prerequisites

- Python 3
- NumPy >= 1.0
- SciPy >= 1.0

## Author

- Ray Patrick

## Assumptions

- Unidirectional and isentropic flow
- Homogeneous combustion products
- No heat exchange with chamber walls (adiabatic)
- No shockwaves or discontinuities in nozzle
- Erosive burning is neglected
- Some others
