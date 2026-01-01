#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test script to compare original SBT code results with tool results for coaxial SBT V2.0 H2O.
This runs the original SBT code directly to see if values match what's showing in the tool.
"""

import sys
import numpy as np
from sbt_v27 import run_sbt as run_sbt_final

print("=" * 80)
print("COAXIAL SBT V2.0 H2O - ORIGINAL CODE TEST")
print("=" * 80)
print("\nTesting with parameters that match the tool's defaults for coaxial H2O...")
print("(This should match what the tool shows for subsurface results)\n")

# Parameters matching tool defaults for coaxial SBT V2.0 H2O
# Based on app.py defaults and clgs.py processing
sbt_version = 2
mesh = 0
accuracy = 1
fluid = 1  # H2O
case = 1  # coaxial

# Default parameters from tool (coaxial defaults - matching tool UI exactly as shown in screenshot)
mdot = 21.0  # kg/s (tool default - slider shows 21)
Tinj = 50.0  # °C (tool default)
L1 = 3.5  # km (tool default - 3500 m)
L2 = 10.0  # km (tool default - 10000 m)
grad = 0.065  # °C/m (tool default - slider shows 0.065)
Tsurf = 25.0  # °C (tool default)
k_m = 3.0  # W/m-K (tool default)
c_m = 790.0  # J/kg-K (tool default)
rho_m = 2800.0  # kg/m³ (tool default - slider shows 2800)

# Geometry parameters (coaxial defaults from tool UI - matching screenshot exactly)
# Note: clgs.py converts Diameter2 from radius to diameter (line 328: Diameter2 = Diameter2 * 2)
Diameter1 = 0.445  # Wellbore diameter (m) - tool default (slider shows 0.445)
Diameter2_radius = 0.101  # Center pipe radius from slider - tool default (slider shows 0.101)
Diameter2 = Diameter2_radius * 2  # Convert radius to diameter (as clgs.py does) = 0.202 m
PipeParam3 = 0.0127  # thicknesscenterpipe (m) - tool default (matches start_vals_sbt["thicknesscenterpipe"])
PipeParam4 = 0.006  # k_center_pipe (insulation thermal conductivity)
PipeParam5 = 1  # Inject in Annulus

# SBT V2.0 hyperparameters
HYPERPARAM1 = 100  # Inlet pressure (bar)
HYPERPARAM2 = 1e-6  # Pipe roughness
HYPERPARAM3 = 1  # Variable fluid properties
HYPERPARAM4 = 1e-5  # Relative tolerance
HYPERPARAM5 = 15  # Max iterations

print("Parameters:")
print(f"  Model: SBT V2.0")
print(f"  Configuration: Coaxial")
print(f"  Fluid: H2O")
print(f"  Mass flow rate: {mdot} kg/s")
print(f"  Injection temperature: {Tinj} °C")
print(f"  L1 (drilling depth): {L1} km")
print(f"  L2 (horizontal extent): {L2} km")
print(f"  Geothermal gradient: {grad} °C/m")
print(f"  Surface temperature: {Tsurf} °C")
print(f"  Thermal conductivity: {k_m} W/m-K")
print(f"  Diameter1 (wellbore): {Diameter1} m")
print(f"  Diameter2_radius (from slider): {Diameter2_radius} m")
print(f"  Diameter2 (converted to diameter): {Diameter2} m")
print(f"  PipeParam3 (thickness): {PipeParam3} m")
print(f"  HYPERPARAM1 (inlet pressure): {HYPERPARAM1} bar")
print()

try:
    print("Running SBT calculation...")
    times, Tout, Pout = run_sbt_final(
        ## Model Specifications 
        sbt_version=sbt_version, 
        mesh_fineness=mesh, 
        HYPERPARAM1=HYPERPARAM1, 
        HYPERPARAM2=HYPERPARAM2, 
        HYPERPARAM3=HYPERPARAM3, 
        HYPERPARAM4=HYPERPARAM4, 
        HYPERPARAM5=HYPERPARAM5, 
        accuracy=accuracy,

        ## Operations
        clg_configuration=case, 
        mdot=mdot, 
        Tinj=Tinj,  # Note: SBT expects Kelvin, but clgs.py converts it
        fluid=fluid,
        DrillingDepth_L1=L1, 
        HorizontalExtent_L2=L2,
        Diameter1=Diameter1, 
        Diameter2=Diameter2, 
        PipeParam3=PipeParam3, 
        PipeParam4=PipeParam4, 
        PipeParam5=PipeParam5,

        ## Geologic Properties
        Tsurf=Tsurf, 
        GeoGradient=grad, 
        k_m=k_m, 
        c_m=c_m, 
        rho_m=rho_m,
    )
    
    print("\n" + "=" * 80)
    print("✓ SUCCESS! Simulation completed.")
    print("=" * 80)
    print(f"\nResults:")
    print(f"  Time points: {len(times)}")
    print(f"  Time range: {times.min():.2f} - {times.max():.2f} years")
    
    if len(Tout) > 0:
        print(f"  Tout range: {Tout.min():.2f} - {Tout.max():.2f} K")
        print(f"  Tout range: {Tout.min()-273.15:.2f} - {Tout.max()-273.15:.2f} °C")
        print(f"  Tout at t=0: {Tout[0]:.2f} K ({Tout[0]-273.15:.2f} °C)")
        if len(Tout) > 14:
            print(f"  Tout at t=0 (after slicing [14:]): {Tout[14]:.2f} K ({Tout[14]-273.15:.2f} °C)")
    
    if len(Pout) > 0:
        print(f"  Pout range: {Pout.min()/1e6:.2f} - {Pout.max()/1e6:.2f} MPa")
        print(f"  Pout range: {Pout.min()/1e5:.2f} - {Pout.max()/1e5:.2f} bar")
        print(f"  Pout at t=0: {Pout[0]/1e6:.2f} MPa ({Pout[0]/1e5:.2f} bar)")
        if len(Pout) > 14:
            print(f"  Pout at t=0 (after slicing [14:]): {Pout[14]/1e6:.2f} MPa ({Pout[14]/1e5:.2f} bar)")
    
    # Show first 20 values for comparison
    print(f"\nFirst 20 time points (after slicing [14:] if applicable):")
    print("  Index | Time (yr) | Tout (K) | Tout (°C) | Pout (MPa)")
    print("  " + "-" * 60)
    start_idx = 14 if len(times) > 14 else 0
    end_idx = min(start_idx + 20, len(times))
    for i in range(start_idx, end_idx):
        print(f"  {i:5d} | {times[i]:8.4f} | {Tout[i]:8.2f} | {Tout[i]-273.15:8.2f} | {Pout[i]/1e6:8.4f}")
    
    print("\n" + "=" * 80)
    print("KEY VALUES FOR COMPARISON WITH TOOL:")
    print("=" * 80)
    if len(Tout) > 14 and len(Pout) > 14:
        print(f"  Tout at t=0 (index 14): {Tout[14]-273.15:.2f} °C ({Tout[14]:.2f} K)")
        print(f"  Pout at t=0 (index 14): {Pout[14]/1e6:.2f} MPa ({Pout[14]/1e5:.2f} bar)")
    print(f"  Tout at t=0 (index 0): {Tout[0]-273.15:.2f} °C ({Tout[0]:.2f} K)")
    print(f"  Pout at t=0 (index 0): {Pout[0]/1e6:.2f} MPa ({Pout[0]/1e5:.2f} bar)")
    print(f"\nCompare these values with what shows in the tool's Subsurface Results tab")
    print("=" * 80)
    
except Exception as e:
    print(f"\n✗ FAILED: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)

