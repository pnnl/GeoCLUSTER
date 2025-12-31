#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import numpy as np
from sbt_v27 import run_sbt as run_sbt_final
from econ_support import create_teaobject, fixed_economic_inputs
from reader import initialize_data
from paths import inpath_dict
import os

# Set up paths
properties_H2O_pathname = inpath_dict["properties_H2O_pathname"]
properties_CO2v2_pathname = inpath_dict["properties_CO2v2_pathname"]
additional_properties_CO2v2_pathname = inpath_dict["additional_properties_CO2v2_pathname"]

# Initialize HDF5 data objects (needed for TEA object)
print("Initializing HDF5 data objects (this may take a moment)...")
u_sCO2, u_H2O, c_sCO2, c_H2O = initialize_data()
print("HDF5 data initialized.")

def test_parameters(mdot, Tinj, L1, grad, k=3.0):
    """Test if a parameter combination produces valid CO2 LCOE"""
    print(f"\n{'='*60}")
    print(f"Testing: mdot={mdot} kg/s, Tinj={Tinj}°C, L1={L1}m, grad={grad}°C/m")
    print(f"{'='*60}")
    
    try:
        # Run SBT2 simulation
        times, Tout, Pout = run_sbt_final(
            sbt_version=2, 
            mesh_fineness=0, 
            HYPERPARAM1=100,  # Pin = 10 MPa = 100 bar
            HYPERPARAM2=1e-6,  # pipe roughness
            HYPERPARAM3=1,  # variable fluid properties
            HYPERPARAM4=1e-5,  # rel tolerance
            HYPERPARAM5=15,  # max iterations
            accuracy=1,
            clg_configuration=1,  # coaxial
            mdot=mdot, 
            Tinj=Tinj + 273.15,  # Convert to Kelvin
            fluid=2,  # CO2
            DrillingDepth_L1=L1/1000,  # Convert to km
            HorizontalExtent_L2=10.0,  # km
            Diameter1=0.2286, 
            Diameter2=0.127, 
            PipeParam3=0.0127, 
            PipeParam4=0.006,
            PipeParam5=1,
            Tsurf=25, 
            GeoGradient=grad, 
            k_m=k, 
            c_m=790.0, 
            rho_m=2750,
        )
        
        print(f"Simulation completed: {len(times)} time points")
        print(f"Tout range: {np.min(Tout)-273.15:.1f} to {np.max(Tout)-273.15:.1f} °C")
        print(f"Pout range: {np.min(Pout)/1e5:.1f} to {np.max(Pout)/1e5:.1f} bar")
        
        # Convert times from seconds to years (skip first 14 points like the app does)
        times_seconds = times[14:] * 365.25 * 24 * 3600
        Tout_skip = Tout[14:]
        Pout_skip = Pout[14:]
        
        # Create TandP_dict (using seconds for time, matching app behavior)
        TandP_dict = {
            "time": times_seconds.tolist(),
            "sCO2_Tout": (Tout_skip).tolist(),  # Already in Kelvin
            "sCO2_Pout": (Pout_skip).tolist(),  # Already in Pascal
        }
        
        # Economic parameters (use initialized data objects)
        case = "coaxial"
        end_use = "Electricity"
        fluid = "sCO2"
        model = "SBT V2.0"
        Flow_user = mdot
        Hor_length_user = 10000  # m
        Depth_user = L1  # m
        Gradient_user = grad * 1000  # Convert to K/km
        Diameter_user = 0.2286  # m
        Tin_user = Tinj  # °C
        krock_user = k  # W/m-K
        Drilling_cost_per_m = 500
        Discount_rate = 7.0
        Lifetime = 30
        Direct_use_heat_cost_per_kWth = 0.05
        Power_plant_cost_per_kWe = 3000
        Pre_Cooling_Delta_T = 5.0
        Turbine_outlet_pressure = 80  # bar (was 7.4 MPa = 74 bar, which is below minimum of 75 bar)
        
        # Create TEA object
        print("\nCreating TEA object...")
        try:
            teaobj = create_teaobject(
                TandP_dict,
                u_sCO2, u_H2O, c_sCO2, c_H2O,
                case, end_use, fluid, model,
                Flow_user, Hor_length_user, Depth_user, Gradient_user, Diameter_user, Tin_user, krock_user,
                Drilling_cost_per_m, Discount_rate, Lifetime,
                Direct_use_heat_cost_per_kWth, Power_plant_cost_per_kWe, Pre_Cooling_Delta_T, Turbine_outlet_pressure,
                properties_H2O_pathname,
                properties_CO2v2_pathname,
                additional_properties_CO2v2_pathname,
            )
        except Exception as e:
            print(f"❌ Error creating TEA object: {type(e).__name__}: {e}")
            import traceback
            traceback.print_exc()
            return False
        
        if teaobj is None:
            print("❌ TEA object creation failed")
            return False
        
        # Check conditions
        print(f"\nT_in: {teaobj.T_in:.2f}°C (need >32°C)")
        if hasattr(teaobj, 'T_turbine_out_actual') and len(teaobj.T_turbine_out_actual) > 0:
            min_turbine_out = np.min(teaobj.T_turbine_out_actual)
            print(f"Min turbine outlet temp: {min_turbine_out:.2f}°C (need >37°C)")
        else:
            print("❌ T_turbine_out_actual not available")
            return False
        
        # Check LCOE
        if teaobj.LCOE is None or teaobj.LCOE >= 9999:
            print("❌ LCOE is None or >= 9999")
            return False
        
        print(f"✅ LCOE: ${teaobj.LCOE:.2f}/MWh")
        
        # Check electricity production
        if hasattr(teaobj, 'Inst_Net_Electricity_production'):
            elec_prod = teaobj.Inst_Net_Electricity_production
            if elec_prod is not None and len(elec_prod) > 0:
                max_elec = np.max(elec_prod)
                avg_elec = np.mean(elec_prod[elec_prod > 0]) if np.any(elec_prod > 0) else 0
                print(f"Max electricity: {max_elec:.2f} kW")
                print(f"Avg electricity (positive): {avg_elec:.2f} kW")
        
        return True
        
    except Exception as e:
        print(f"❌ Error: {type(e).__name__}: {e}")
        import traceback
        traceback.print_exc()
        return False

if __name__ == "__main__":
    print("Testing parameter combinations for CO2 LCOE in SBT2")
    print("="*60)
    
    # Test different parameter combinations
    test_cases = [
        # (mdot, Tinj, L1, grad)
        (20.0, 60.0, 5000.0, 0.095),  # Current defaults
        (25.0, 60.0, 5000.0, 0.095),  # Higher flow
        (20.0, 70.0, 5000.0, 0.095),  # Higher injection temp
        (20.0, 60.0, 5500.0, 0.095),  # Deeper well
        (20.0, 60.0, 5000.0, 0.100),  # Higher gradient
        (15.0, 50.0, 4500.0, 0.090),  # Lower values
        (30.0, 50.0, 4000.0, 0.080),  # High flow, lower temp
        (12.0, 60.0, 5000.0, 0.095),  # Lower flow
    ]
    
    working_params = []
    
    for mdot, Tinj, L1, grad in test_cases:
        if test_parameters(mdot, Tinj, L1, grad):
            working_params.append((mdot, Tinj, L1, grad))
    
    print(f"\n{'='*60}")
    print("SUMMARY")
    print(f"{'='*60}")
    print(f"Working parameter combinations: {len(working_params)}/{len(test_cases)}")
    for mdot, Tinj, L1, grad in working_params:
        print(f"  ✅ mdot={mdot} kg/s, Tinj={Tinj}°C, L1={L1}m, grad={grad}°C/m")
    
    if len(working_params) == 0:
        print("❌ No working parameter combinations found!")
        sys.exit(1)
    else:
        print(f"\n✅ Found {len(working_params)} working combination(s)")
