#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Quick test to validate coaxial calculations - matches existing test format
"""

import sys
from sbt_v27 import run_sbt as run_sbt_final

print("Testing coaxial with CORRECT unit conversions (as in clgs.py)...")
print("Diameter1=0.4572 m (diameter), Diameter2=0.254 m (diameter, converted from 0.127 radius)")

try:
    times, Tout, Pout = run_sbt_final(
        sbt_version=2, 
        mesh_fineness=0, 
        HYPERPARAM1=100,
        HYPERPARAM2=1e-6,
        HYPERPARAM3=1,
        HYPERPARAM4=1e-5,
        HYPERPARAM5=15,
        accuracy=1,
        clg_configuration=1,  # coaxial
        mdot=30.0, 
        Tinj=55.0 + 273.15,  # Convert to Kelvin
        fluid=1,  # H2O
        DrillingDepth_L1=4000.0/1000,  # Convert to km
        HorizontalExtent_L2=10.0,  # km
        Diameter1=0.4572,  # Wellbore diameter (correct - matches slider output)
        Diameter2=0.254,  # Center pipe diameter (converted from 0.127 radius, as clgs.py does)
        PipeParam3=0.0127, 
        PipeParam4=0.006,
        PipeParam5=1,
        Tsurf=25, 
        GeoGradient=0.080, 
        k_m=3.0, 
        c_m=790.0, 
        rho_m=2750,
    )
    print(f"\n✓ SUCCESS!")
    print(f"  Time points: {len(times)}")
    print(f"  Tout range: {Tout.min():.2f} - {Tout.max():.2f} K")
    print(f"  Pout range: {Pout.min()/1e6:.2f} - {Pout.max()/1e6:.2f} MPa")
except Exception as e:
    print(f"\n✗ FAILED: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)

