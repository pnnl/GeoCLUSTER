#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Simple test to validate coaxial SBT calculations with correct unit conversions.
Tests a single case to verify the calculation works.
"""

import sys
import numpy as np
from sbt_v27 import run_sbt as run_sbt_final

def test_single_case():
    """Test a single coaxial case with known parameters"""
    
    print("Testing coaxial SBT calculation...")
    print("Parameters:")
    print("  Fluid: H2O")
    print("  mdot: 30 kg/s")
    print("  Diameter1 (wellbore): 0.4572 m (diameter from slider)")
    print("  Diameter2 (center pipe): 0.127 m (radius from slider) → 0.254 m (diameter)")
    
    try:
        times, Tout, Pout = run_sbt_final(
            ## Model Specifications 
            sbt_version=2, 
            mesh_fineness=0, 
            HYPERPARAM1=100,  # 100 bar inlet pressure
            HYPERPARAM2=1e-6,  # pipe roughness
            HYPERPARAM3=1,  # variable fluid properties
            HYPERPARAM4=1e-5,  # rel tolerance
            HYPERPARAM5=15,  # max iterations
            accuracy=1,

            ## Operations
            clg_configuration=1,  # coaxial
            mdot=30.0, 
            Tinj=55.0, 
            fluid=1,  # H2O
            DrillingDepth_L1=4000.0, 
            HorizontalExtent_L2=10000.0,
            Diameter1=0.4572,  # Wellbore diameter (already a diameter)
            Diameter2=0.254,  # Center pipe diameter (converted from 0.127 radius)
            PipeParam3=0.0127,  # thicknesscenterpipe
            PipeParam4=0.006,  # k_center_pipe
            PipeParam5=1,  # Inject in Annulus

            ## Geologic Properties
            Tsurf=25.0, 
            GeoGradient=0.080, 
            k_m=3.0, 
            c_m=790.0, 
            rho_m=2800.0,
        )
        
        print("\n✓ SUCCESS! Simulation completed.")
        print(f"  Tout range: {Tout.min():.2f} - {Tout.max():.2f} K ({Tout.min()-273.15:.2f} - {Tout.max()-273.15:.2f} °C)")
        print(f"  Pout range: {Pout.min()/1e6:.2f} - {Pout.max()/1e6:.2f} MPa")
        print(f"  Time range: {times.min():.2f} - {times.max():.2f} years")
        print(f"  Number of time points: {len(times)}")
        return True
        
    except Exception as e:
        print(f"\n✗ FAILED: {e}")
        import traceback
        traceback.print_exc()
        return False

if __name__ == "__main__":
    success = test_single_case()
    sys.exit(0 if success else 1)

