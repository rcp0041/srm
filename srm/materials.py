#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from srm import propellant

rho_water = 0.03613 # lbm/in^3

class metal:
    def __init__(self,specific_gravity=None,density=None,yield_strength=None,name=None):
        self.name = name
        if specific_gravity != None:
            self.specific_gravity = specific_gravity
            self.density = specific_gravity*rho_water
        elif density != None:
            self.density = density
            self.specific_gravity = density/rho_water
        self.yield_strength = yield_strength
    
    def __str__(self):
        a = "Metal properties:\n"
        a = a + f"Name: {self.name}\n"
        a = a + f"Density: {self.density} lbm/in^3\n"
        a = a + f"Yield strength: {self.yield_strength} psi"
        return a

steel = {
    # "common" : metal(name="Common steel",specific_gravity=7.84,yield_strength=40000),
    # "low alloy" : metal(specific_gravity=7.84,yield_strength=190000,name="Low alloy steel"),
    # "maraging" : metal(specific_gravity=7.84,yield_strength=270000,name="Maraging steel"),
    # "Al alloy" : metal(specific_gravity=2.77,yield_strength=85000,name="Aluminum alloy steel"),
    # "Ti alloy" : metal(specific_gravity=4.43,yield_strength=95000,name="Titanium alloy steel"),
    # "Fiberglass-plastic" : metal(specific_gravity=1.75,yield_strength=65000,name="Fiberglass-plastic"),
    # These are from NASA SP-80025 Table 1
    # TODO: Implement yield strength bounds/uncertainties
    'HY-80'  : metal(name='HY-80',density=0.285,yield_strength=80000),
    'HY-130' : metal(name='HY-130',density=0.285,yield_strength=130000),
    '4130' : metal(name='4130',density=0.283,yield_strength=165000),
    '4335V' : metal(name='4335V',density=0.283,yield_strength=190000),
    'D6aC' : metal(name='D6aC',density=0.283,yield_strength=190000),
    'Grade 200' : metal(name='Grade 200',density=0.289,yield_strength=200000),
    'Grade 250' : metal(name='Grade 250',density=0.289,yield_strength=240000),
    'Grade 300' : metal(name='Grade 300',density=0.289,yield_strength=280000),
    '9 Ni-4 Co-0.250' : metal(name='9 Ni-4 Co-0.250',density=0.28,yield_strength=200000),
    '9 Ni-4 Co-0.450' : metal(name='9 Ni-4 Co-0.450',density=0.28,yield_strength=280000),
    'Ti-6A1-4V' : metal(name='Ti-6A1-4V',density=0.167,yield_strength=150000),
    '2000 Series' : metal(name='2000 Series',density=0.10,yield_strength=50000),
    '5000 Series' : metal(name='5000 Series',density=0.10,yield_strength=35000),
    '6000 Series' : metal(name='6000 Series',density=0.10,yield_strength=42000),
    '7000 Series' : metal(name='7000 Series',density=0.10,yield_strength=64000),
    }

propellants = {
    "DB" : propellant(
        name="DB",
        a=0.0567,
        n=0.30,
        density=0.058,
        temperature=4100,
        specific_heat_ratio=1.22,
        cstar=4658,
        MW=25.5,
    ),
    "DB_AP_AL" :  propellant(
        name="DB_AP_AL",
        a=0.0492,
        n=0.40,
        density=0.065,
        temperature=6500,
        specific_heat_ratio=1.24,
        cstar=5364,
        MW=28.0
    ),
    "PVC_AP_AL" : propellant(
        name="PVC_AP_AL",
        a=0.0401,
        n=0.35,
        density=0.064,
        temperature=5600,
        specific_heat_ratio=1.24,
        cstar=5364,
        MW=24.3
    ),
    "PS_AP_AL" : propellant(
        name="PS_AP_AL",
        a=0.0317,
        n=0.33,
        density=0.062,
        temperature=5000,
        specific_heat_ratio=1.24,
        cstar=5006,
        MW=25.2
    ),
    "PBAN_AP_AL" : propellant(
        name="PBAN_AP_AL",
        a=0.0563,
        n=0.33,
        density=0.064,
        temperature=5800,
        specific_heat_ratio=1.24,
        cstar=5343,
        MW=25.3
    ),
    "HTPB_AP_AL" : propellant(
        name="HTPB_AP_AL",
        a=0.025,
        n=0.40,
        density=0.064,
        temperature=5700,
        specific_heat_ratio=1.24,
        cstar=5364,
        MW=24.7
    ),
    "PBAA_AP_AN" : propellant(
        name="PBAA_AP_AN",
        a=0.0285,
        n=0.35,
        density=0.064,
        temperature=5700,
        specific_heat_ratio=1.24,
        cstar=5364,
        MW=24.7
    ),
    "AP_Polymer" : propellant(
        name="AP_Polymer",
        a=0.00475,
        n=0.60,
        density=0.063,
        temperature=2300,
        specific_heat_ratio=1.26,
        cstar=3803,
        MW=21.8
    ),
    }