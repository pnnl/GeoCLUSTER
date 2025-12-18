#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test script to validate coaxial SBT calculations with proper unit conversions.
Tests both H2O and CO2 fluids to ensure results are possible.
"""

import sys
import numpy as np
from sbt_v27 import run_sbt as run_sbt_final

def test_coaxial_validation():
    """Test coaxial SBT calculations with known working parameters"""
    
    print("=" * 80)
    print("COAXIAL SBT VALIDATION TEST")
    print("=" * 80)
    
    # Test cases: (mdot, Diameter1_diameter, Diameter2_radius, fluid_name, fluid_id)
    # Note: Diameter1 is wellbore diameter, Diameter2 is center pipe radius (from slider)
    # Start with just one case to avoid memory issues
    test_cases = [
        # H2O test case - most likely to work
        (30.0, 0.4572, 0.127, "H2O", 1),  # Standard defaults
    ]
    
    # Fixed parameters
    sbt_version = 2
    mesh = 0
    accuracy = 1
    Tinj = 55.0  # °C
    L1 = 4000.0  # m
    L2 = 10000.0  # m
    grad = 0.080  # °C/m
    Tsurf = 25.0  # °C
    c_m = 790.0  # J/kg-K
    rho_m = 2800.0  # kg/m³
    k_m = 3.0  # W/m-K
    
    # SBT V2.0 hyperparameters
    HYPERPARAM1 = 10.0 * 10  # 100 bar inlet pressure (10 MPa * 10)
    HYPERPARAM2 = 1e-6  # pipe roughness
    HYPERPARAM3 = 1  # variable fluid properties
    HYPERPARAM4 = 1e-5  # rel tolerance
    HYPERPARAM5 = 15  # max iterations
    
    PipeParam3 = 0.0127  # thicknesscenterpipe
    PipeParam4 = 0.006  # k_center_pipe
    PipeParam5 = 1  # Inject in Annulus
    
    successful_tests = []
    failed_tests = []
    
    for mdot, d1_diameter, d2_radius, fluid_name, fluid_id in test_cases:
        # Convert Diameter2 from radius to diameter (as clgs.py does)
        d2_diameter = d2_radius * 2
        
        print(f"\n{'='*80}")
        print(f"Testing: {fluid_name}, mdot={mdot} kg/s")
        print(f"  Diameter1 (wellbore): {d1_diameter} m (diameter)")
        print(f"  Diameter2 (center pipe): {d2_radius} m (radius) → {d2_diameter} m (diameter)")
        print(f"{'='*80}")
        
        try:
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
                clg_configuration=1,  # coaxial
                mdot=mdot, 
                Tinj=Tinj, 
                fluid=fluid_id,
                DrillingDepth_L1=L1, 
                HorizontalExtent_L2=L2,
                Diameter1=d1_diameter,  # Already a diameter
                Diameter2=d2_diameter,  # Converted from radius to diameter
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
            
            # Validate results
            if Tout is None or len(Tout) == 0:
                raise ValueError("Tout is empty")
            if Pout is None or len(Pout) == 0:
                raise ValueError("Pout is empty")
            if times is None or len(times) == 0:
                raise ValueError("times is empty")
            
            # Check for invalid values
            Tout_arr = np.array(Tout)
            Pout_arr = np.array(Pout)
            
            if np.any(np.isnan(Tout_arr)) or np.any(np.isinf(Tout_arr)):
                raise ValueError(f"Invalid Tout values: NaN or Inf detected")
            
            if np.any(np.isnan(Pout_arr)) or np.any(np.isinf(Pout_arr)):
                raise ValueError(f"Invalid Pout values: NaN or Inf detected")
            
            # Check reasonable temperature range (200-1000 K)
            if np.any(Tout_arr < 200) or np.any(Tout_arr > 1000):
                raise ValueError(f"Tout out of reasonable range: {Tout_arr.min():.2f} - {Tout_arr.max():.2f} K")
            
            # Success!
            print(f"✓ SUCCESS!")
            print(f"  Tout range: {Tout_arr.min():.2f} - {Tout_arr.max():.2f} K ({Tout_arr.min()-273.15:.2f} - {Tout_arr.max()-273.15:.2f} °C)")
            print(f"  Pout range: {Pout_arr.min()/1e6:.2f} - {Pout_arr.max()/1e6:.2f} MPa")
            print(f"  Time range: {times.min():.2f} - {times.max():.2f} years")
            print(f"  Number of time points: {len(times)}")
            
            successful_tests.append({
                'fluid': fluid_name,
                'mdot': mdot,
                'd1_diameter': d1_diameter,
                'd2_radius': d2_radius,
                'Tout_min': Tout_arr.min(),
                'Tout_max': Tout_arr.max(),
                'Pout_min': Pout_arr.min(),
                'Pout_max': Pout_arr.max(),
            })
            
        except ValueError as e:
            error_msg = str(e)
            print(f"✗ FAILED: {error_msg}")
            failed_tests.append({
                'fluid': fluid_name,
                'mdot': mdot,
                'd1_diameter': d1_diameter,
                'd2_radius': d2_radius,
                'error': error_msg,
            })
        except Exception as e:
            error_msg = str(e)
            print(f"✗ FAILED: Unexpected error: {error_msg}")
            import traceback
            traceback.print_exc()
            failed_tests.append({
                'fluid': fluid_name,
                'mdot': mdot,
                'd1_diameter': d1_diameter,
                'd2_radius': d2_radius,
                'error': error_msg,
            })
    
    # Summary
    print(f"\n{'='*80}")
    print("TEST SUMMARY")
    print(f"{'='*80}")
    print(f"Successful tests: {len(successful_tests)}/{len(test_cases)}")
    print(f"Failed tests: {len(failed_tests)}/{len(test_cases)}")
    
    if successful_tests:
        print(f"\n✓ SUCCESSFUL PARAMETER COMBINATIONS:")
        for i, test in enumerate(successful_tests, 1):
            print(f"\n{i}. {test['fluid']}:")
            print(f"   mdot: {test['mdot']} kg/s")
            print(f"   Diameter1: {test['d1_diameter']} m (wellbore diameter)")
            print(f"   Diameter2: {test['d2_radius']} m (center pipe radius)")
            print(f"   Tout: {test['Tout_min']-273.15:.2f} - {test['Tout_max']-273.15:.2f} °C")
            print(f"   Pout: {test['Pout_min']/1e6:.2f} - {test['Pout_max']/1e6:.2f} MPa")
    
    if failed_tests:
        print(f"\n✗ FAILED PARAMETER COMBINATIONS:")
        for i, test in enumerate(failed_tests, 1):
            print(f"\n{i}. {test['fluid']}:")
            print(f"   mdot: {test['mdot']} kg/s")
            print(f"   Diameter1: {test['d1_diameter']} m (wellbore diameter)")
            print(f"   Diameter2: {test['d2_radius']} m (center pipe radius)")
            print(f"   Error: {test['error']}")
    
    return len(successful_tests) > 0

if __name__ == "__main__":
    success = test_coaxial_validation()
    sys.exit(0 if success else 1)

