#!/usr/bin/env python3
"""
Quick test script to check CO2 LCOE calculation without full HDF5 initialization
"""

import sys
import numpy as np
from sbt_v27 import run_sbt

# Test parameters
mdot = 20.0
L1 = 5000
L2 = 10000
grad = 0.095
Tinj = 60.0
D = 0.2286
k = 3.0
Tsurf = 25
c_m = 790.0
rho_m = 2750

print(f"Quick CO2 LCOE test:")
print(f"  mdot={mdot}, L1={L1}, grad={grad}, Tinj={Tinj}")
print("Running SBT2 simulation...")

try:
    # Run SBT simulation
    times, Toutput, Poutput = run_sbt(
        sbt_version=2,
        clg_configuration=1,  # coaxial
        mdot=mdot,
        Tinj=Tinj + 273.15,  # Convert to Kelvin
        fluid=2,  # CO2
        DrillingDepth_L1=L1/1000,  # Convert to km
        HorizontalExtent_L2=L2/1000,  # Convert to km
        Diameter1=D,
        Diameter2=0.127,
        PipeParam3=0.0127,
        PipeParam4=0.006,
        PipeParam5=1,
        GeoGradient=grad,
        Tsurf=Tsurf,
        k_m=k,
        c_m=c_m,
        rho_m=rho_m,
        accuracy=1,
        mesh_fineness=0,
        HYPERPARAM1=100,  # 10 MPa = 100 bar
        HYPERPARAM2=1e-6,
        HYPERPARAM3=1,
        HYPERPARAM4=1e-5,
        HYPERPARAM5=15,
    )
    
    print(f"✅ Simulation complete!")
    print(f"   Times shape: {times.shape}")
    print(f"   Toutput shape: {Toutput.shape}")
    print(f"   Poutput shape: {Poutput.shape}")
    
    # Skip first 14 points
    times_skip = times[14:]
    Toutput_skip = Toutput[14:]
    Poutput_skip = Poutput[14:]
    
    print(f"\nAfter skipping first 14 points:")
    print(f"   Toutput range: {np.min(Toutput_skip)-273.15:.1f} to {np.max(Toutput_skip)-273.15:.1f} °C")
    print(f"   Poutput range: {np.min(Poutput_skip)/1e5:.1f} to {np.max(Poutput_skip)/1e5:.1f} bar")
    print(f"   Time range: {np.min(times_skip):.2f} to {np.max(times_skip):.2f} years")
    
    # Check if we have valid data
    if len(times_skip) == 0:
        print("❌ ERROR: No time points after skipping first 14!")
        sys.exit(1)
    
    if np.any(np.isnan(Toutput_skip)) or np.any(np.isnan(Poutput_skip)):
        print("❌ ERROR: NaN values in output!")
        sys.exit(1)
    
    print(f"\n✅ Basic SBT2 simulation check passed!")
    print(f"   Ready for TEA object creation (would need HDF5 initialization)")
    
except Exception as e:
    print(f"❌ Error: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)

