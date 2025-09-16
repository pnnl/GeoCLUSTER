#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
UI to SI Conversion Helper for GeoCLUSTER
Converts current UI slider values (in the active unit system) back to SI for computation/HDF5 lookup
"""

from unit_conversions import unit_converter

def ui_to_SI(values):
    """
    Convert current UI slider values (in the *active* unit system)
    back to SI for computation / HDF5 lookup.
    `values` is a dict with raw UI numbers.
    Returns a dict of SI numbers.
    """
    prefs = unit_converter.user_preferences  # e.g. {'length':'ft', 'temperature':'F', ...}

    # Temperature (UI is C or F; model expects Kelvin)
    if prefs.get('temperature','C') == 'F':
        Tinj_K = unit_converter.convert_temperature(values['Tinj'], 'F', 'C') + 273.15
    else:
        Tinj_K = values['Tinj'] + 273.15  # UI shows °C

    # Length (m vs ft)
    L1_m = unit_converter.convert_length(values['L1'], prefs.get('length','m'), 'm')
    L2_m = unit_converter.convert_length(values['L2'], prefs.get('length','m'), 'm')
    D_m  = unit_converter.convert_length(values['D'],  prefs.get('length','m'), 'm')

    # Mass flow (kg/s vs lb/s)
    mdot_kg_s = unit_converter.convert_mass_flow(values['mdot'], prefs.get('mass_flow','kg/s'), 'kg/s')

    # Rock k (W/m-K vs Btu/ft-h-F)
    k_w_m_k = unit_converter.convert_thermal_conductivity(values['k'], prefs.get('thermal_conductivity','W/m-K'), 'W/m-K')

    # Gradient (K/m vs F/ft)
    grad_K_m = unit_converter.convert_geothermal_gradient(values['grad'], prefs.get('geothermal_gradient','K/m'), 'K/m')

    return dict(
        Tinj=Tinj_K, L1=L1_m, L2=L2_m, D=D_m, mdot=mdot_kg_s, k=k_w_m_k, grad=grad_K_m
    )
