#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test script to find parameter combinations for coaxial SBT models that avoid laminar flow.
"""

import sys
import numpy as np
from sbt_v27 import run_sbt as run_sbt_final

def test_coaxial_params():
    """Test different parameter combinations for coaxial SBT models"""
    
    # Test a single high-value combination first
    test_cases = [
        {'mdot': 200, 'radius': 0.4445, 'radiuscenterpipe': 0.127},
        {'mdot': 250, 'radius': 0.4445, 'radiuscenterpipe': 0.127},
        {'mdot': 300, 'radius': 0.4445, 'radiuscenterpipe': 0.127},
    ]
    
    # Fixed parameters
    sbt_version = 2
    mesh = 0
    accuracy = 1
    fluid = 2  # CO2
    Tinj = 50.0
    L1 = 4000.0
    L2 = 10000.0
    grad = 0.080
    Tsurf = 25.0
    c_m = 790.0
    rho_m = 2800.0
    k_m = 3.0
    
    # SBT V2.0 hyperparameters
    HYPERPARAM1 = 10.0 * 10  # 100 bar inlet pressure
    HYPERPARAM2 = 1e-6  # pipe roughness
    HYPERPARAM3 = 1  # variable fluid properties
    HYPERPARAM4 = 1e-5  # rel tolerance
    HYPERPARAM5 = 15  # max iterations
    
    PipeParam3 = 0.0127  # thicknesscenterpipe
    PipeParam4 = 0.006  # k_center_pipe
    PipeParam5 = 1  # Inject in Annulus
    
    successful_params = []
    
    print("Testing coaxial SBT parameters to find combinations that avoid laminar flow...")
    print("=" * 80)
    
    for test_case in test_cases:
        mdot = test_case['mdot']
        radius = test_case['radius']
        radiuscenterpipe = test_case['radiuscenterpipe']
        
        try:
            print(f"\nTesting: mdot={mdot} kg/s, radius={radius} m, radiuscenterpipe={radiuscenterpipe} m")
            
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
                fluid=fluid,
                DrillingDepth_L1=L1, 
                HorizontalExtent_L2=L2,
                Diameter1=radius * 2,  # wellbore diameter (convert radius to diameter)
                Diameter2=radiuscenterpipe * 2,  # center pipe diameter (convert radius to diameter)
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
            
            # If we get here, simulation succeeded!
            print(f"✓ SUCCESS! Parameters work: mdot={mdot} kg/s, radius={radius} m")
            successful_params.append({
                'mdot': mdot,
                'radius': radius,
                'radiuscenterpipe': radiuscenterpipe,
                'Tout_range': (Tout.min(), Tout.max()),
                'Pout_range': (Pout.min(), Pout.max()),
            })
            
        except ValueError as e:
            if "laminar flow" in str(e):
                print(f"  ✗ Laminar flow error")
            else:
                print(f"  ✗ Other error: {e}")
        except Exception as e:
            print(f"  ✗ Unexpected error: {e}")
            import traceback
            traceback.print_exc()
    
    print("\n" + "=" * 80)
    print("SUMMARY OF SUCCESSFUL PARAMETER COMBINATIONS:")
    print("=" * 80)
    
    if successful_params:
        for i, params in enumerate(successful_params, 1):
            print(f"\n{i}. Mass flow rate: {params['mdot']} kg/s")
            print(f"   Wellbore radius: {params['radius']} m")
            print(f"   Center pipe radius: {params['radiuscenterpipe']} m")
            print(f"   Tout range: {params['Tout_range'][0]:.2f} - {params['Tout_range'][1]:.2f} K")
            print(f"   Pout range: {params['Pout_range'][0]/1e6:.2f} - {params['Pout_range'][1]/1e6:.2f} MPa")
        
        # Find minimum working mdot for each radius
        print("\n" + "=" * 80)
        print("MINIMUM WORKING PARAMETERS BY RADIUS:")
        print("=" * 80)
        radius_groups = {}
        for params in successful_params:
            radius = params['radius']
            if radius not in radius_groups or params['mdot'] < radius_groups[radius]['mdot']:
                radius_groups[radius] = params
        
        for radius in sorted(radius_groups.keys()):
            params = radius_groups[radius]
            print(f"\nRadius {radius} m: Minimum mdot = {params['mdot']} kg/s")
    else:
        print("\nNo successful parameter combinations found!")
        print("May need to test even higher mass flow rates or different parameters.")

if __name__ == "__main__":
    test_coaxial_params()

