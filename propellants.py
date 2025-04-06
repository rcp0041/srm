#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from srm.srm import propellant

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