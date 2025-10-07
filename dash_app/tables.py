#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# data manipulation libraries
import pandas as pd
import plotly.graph_objects as go

# source script
from paths import inpath_dict

# -----------------------
# Global colors.
# -----------------------

greyyellow = "#f2f2f0"
greygrey = "#ededed"
greyblue = "#e6eff5"
white = '#ffffff'
greybluel = "#f5fbff" # "#f2f5f7"

# -----------------------
# Read in data.
# -----------------------

summary_tbl_pathname = inpath_dict["summary_tbl_in"] # "data/summary_tbl_in.csv"
# summary_tbl_pathname = "/var/www/html/dash_app/data/summary_tbl_in.csv"
# summary_tbl_pathname = "/www/GeoCLUSTER/dash_app/data/summary_tbl_in.csv"
summary_tbl = pd.read_csv(summary_tbl_pathname) 

def generate_summary_table(mdot, L2, L1, grad, D, Tinj, k, Drilling_cost_per_m, Discount_rate, Lifetime, 
                             Direct_use_heat_cost_per_kWth, Power_plant_cost_per_kWe, Pre_Cooling_Delta_T, Turbine_outlet_pressure, 
                             interp_time, case, fluid, model,
                             thermal_dict, econ_dict, units="metric", **kwargs):

    # -----------------------------------------------------------------------------------------------------------------
    # Creates Plotly a plotly table. Then, returns the figure and data as a dictionary of lits for data downloading.
    # -----------------------------------------------------------------------------------------------------------------

    # Create variable input values list - only include the 20 parameters that match the CSV template rows
    # Convert values to imperial units for display if imperial is selected
    if units == "imperial":
        # Convert metric values to imperial for display
        mdot_display = round(mdot * 2.20462, 4)  # kg/s to lb/s
        L2_display = round(L2 * 3.28084, 4)  # m to ft
        L1_display = round(L1 * 3.28084, 4)  # m to ft
        grad_display = round(grad * 0.3048 * (9.0/5.0), 4)  # K/m to °F/ft
        D_display = round(D * 3.28084, 4)  # m to ft
        Tinj_display = round((Tinj * 9.0/5.0) + 32.0, 4)  # °C to °F
        k_display = round(k * 0.577789, 4)  # W/m·K to BTU/(hr·ft·°F)
        Drilling_cost_display = round(Drilling_cost_per_m * 0.3048, 4)  # $/m to $/ft
        Pre_Cooling_Delta_T_display = round(Pre_Cooling_Delta_T * (9.0/5.0), 4)  # °C to °F
        Turbine_outlet_pressure_display = round(Turbine_outlet_pressure * 145.038, 4)  # MPa to psi
    else:
        # Metric values - pass through as-is with 4 decimal places
        mdot_display = round(mdot, 4)
        L2_display = round(L2, 4)
        L1_display = round(L1, 4)
        grad_display = round(grad, 4)
        D_display = round(D, 4)
        Tinj_display = round(Tinj, 4)
        k_display = round(k, 4)
        Drilling_cost_display = round(Drilling_cost_per_m, 4)
        Pre_Cooling_Delta_T_display = round(Pre_Cooling_Delta_T, 4)
        Turbine_outlet_pressure_display = round(Turbine_outlet_pressure, 4)
    
    variable_input_values = [mdot_display, L2_display, L1_display, grad_display, D_display, Tinj_display, k_display, 
                             Drilling_cost_display, Discount_rate, Lifetime, 
                             Direct_use_heat_cost_per_kWth, Power_plant_cost_per_kWe, Pre_Cooling_Delta_T_display, Turbine_outlet_pressure_display, 
                             interp_time, case, fluid]
    
    # Note: SBT-specific parameters are handled separately and don't get added to the main summary table
    # to avoid creating extra rows that don't match the CSV template structure

    # Handle different data structures for HDF5 vs SBT models
    if model == "HDF5" or model == "CovHDF5":
        # HDF5 and CovHDF5 data structure
        print(f"DEBUG: Model = {model}")
        print(f"DEBUG: thermal_dict keys = {thermal_dict.keys()}")
        print(f"DEBUG: econ_dict keys = {econ_dict.keys()}")
        print(f"DEBUG: thermal_dict = {thermal_dict}")
        print(f"DEBUG: econ_dict = {econ_dict}")
        
        results = [econ_dict.get('LCOH H2O', '-'), econ_dict.get('LCOE H2O', '-'), 
                   thermal_dict.get('Mean H2O Tout', '-'), thermal_dict.get('Mean H2O Pout', '-'), 
                   econ_dict.get('Mean H2O Net HProd', '-'), econ_dict.get('Mean H2O Net EProd', '-'),
                   econ_dict.get('LCOH sCO2', '-'), econ_dict.get('LCOE sCO2', '-'), 
                   thermal_dict.get('Mean sCO2 Tout', '-'), thermal_dict.get('Mean sCO2 Pout', '-'),
                   econ_dict.get('Mean sCO2 Net HProd', '-'), econ_dict.get('Mean sCO2 Net EProd', '-')
                   ]
    else:
        # SBT model data structure - calculate mean values from arrays
        import numpy as np
        
        # Process TandP data for SBT models
        tandp_data = thermal_dict.get('TandP-data', None)
        if model != "HDF5" and tandp_data:
            # Calculate mean values from arrays for SBT models
            if 'H2O_Tout' in tandp_data and tandp_data['H2O_Tout'] is not None:
                h2o_tout_array = np.array(tandp_data['H2O_Tout'])
                # Check if array is not empty and has valid data
                if len(h2o_tout_array) > 0 and not (h2o_tout_array == None).all():
                    try:
                        if units == "imperial":
                            mean_h2o_tout = round((np.mean(h2o_tout_array - 273.15) * 9.0/5.0) + 32.0, 4)  # Convert from K to °F
                        else:
                            mean_h2o_tout = round(np.mean(h2o_tout_array - 273.15), 4)  # Convert from K to °C
                    except Exception as e:
                        mean_h2o_tout = '-'
                else:
                    mean_h2o_tout = '-'
            else:
                mean_h2o_tout = '-'

            if 'H2O_Pout' in tandp_data and tandp_data['H2O_Pout'] is not None:
                h2o_pout_array = np.array(tandp_data['H2O_Pout'])
                if len(h2o_pout_array) > 0 and not (h2o_pout_array == None).all():
                    try:
                        if units == "imperial":
                            mean_h2o_pout = round(np.mean(h2o_pout_array) * 0.145038, 4)  # Convert from Pa to psi
                        else:
                            mean_h2o_pout = round(np.mean(h2o_pout_array / 1e6), 4)  # Convert from Pa to MPa
                    except Exception as e:
                        mean_h2o_pout = '-'
                else:
                    mean_h2o_pout = '-'
            else:
                mean_h2o_pout = '-'

            # For SBT V2.0, also calculate sCO2 values
            if model == "SBT V2.0":
                if 'sCO2_Tout' in tandp_data and tandp_data['sCO2_Tout'] is not None:
                    sco2_tout_array = np.array(tandp_data['sCO2_Tout'])
                    if len(sco2_tout_array) > 0 and not (sco2_tout_array == None).all():
                        try:
                            if units == "imperial":
                                mean_sco2_tout = round((np.mean(sco2_tout_array - 273.15) * 9.0/5.0) + 32.0, 4)  # Convert from K to °F
                            else:
                                mean_sco2_tout = round(np.mean(sco2_tout_array - 273.15), 4)  # Convert from K to °C
                        except Exception as e:
                            mean_sco2_tout = '-'
                    else:
                        mean_sco2_tout = '-'
                else:
                    mean_sco2_tout = '-'

                if 'sCO2_Pout' in tandp_data and tandp_data['sCO2_Pout'] is not None:
                    sco2_pout_array = np.array(tandp_data['sCO2_Pout'])
                    if len(sco2_pout_array) > 0 and not (sco2_pout_array == None).all():
                        try:
                            if units == "imperial":
                                mean_sco2_pout = round(np.mean(sco2_pout_array) * 0.145038, 4)  # Convert from Pa to psi
                            else:
                                mean_sco2_pout = round(np.mean(sco2_pout_array / 1e6), 4)  # Convert from Pa to MPa
                        except Exception as e:
                            mean_sco2_pout = '-'
                    else:
                        mean_sco2_pout = '-'
                else:
                    mean_sco2_pout = '-'
            else:
                # SBT V1.0 doesn't support sCO2
                mean_sco2_tout = '-'
                mean_sco2_pout = '-'
        else:
            # For HDF5 models, use existing thermal_dict values and convert if needed
            if units == "imperial":
                # Convert from metric to imperial
                h2o_tout_c = thermal_dict.get('Mean H2O Tout', '-')
                if h2o_tout_c != '-':
                    mean_h2o_tout = round((h2o_tout_c * 9.0/5.0) + 32.0, 4)  # °C to °F
                else:
                    mean_h2o_tout = '-'
                
                h2o_pout_mpa = thermal_dict.get('Mean H2O Pout', '-')
                if h2o_pout_mpa != '-':
                    mean_h2o_pout = round(h2o_pout_mpa * 145.038, 4)  # MPa to psi
                else:
                    mean_h2o_pout = '-'
                
                sco2_tout_c = thermal_dict.get('Mean sCO2 Tout', '-')
                if sco2_tout_c != '-':
                    mean_sco2_tout = round((sco2_tout_c * 9.0/5.0) + 32.0, 4)  # °C to °F
                else:
                    mean_sco2_tout = '-'
                
                sco2_pout_mpa = thermal_dict.get('Mean sCO2 Pout', '-')
                if sco2_pout_mpa != '-':
                    mean_sco2_pout = round(sco2_pout_mpa * 145.038, 4)  # MPa to psi
                else:
                    mean_sco2_pout = '-'
            else:
                # Metric values - pass through as-is
                mean_h2o_tout = thermal_dict.get('Mean H2O Tout', '-')
                mean_h2o_pout = thermal_dict.get('Mean H2O Pout', '-')
                mean_sco2_tout = thermal_dict.get('Mean sCO2 Tout', '-')
                mean_sco2_pout = thermal_dict.get('Mean sCO2 Pout', '-')

        # Get economic values
        lcoh_h2o = econ_dict.get('LCOH H2O', '-')
        lcoe_h2o = econ_dict.get('LCOE H2O', '-')
        lcoh_sco2 = econ_dict.get('LCOH sCO2', '-')
        lcoe_sco2 = econ_dict.get('LCOE sCO2', '-')
        mean_h2o_net_hprod = econ_dict.get('Mean H2O Net HProd', '-')
        mean_h2o_net_eprod = econ_dict.get('Mean H2O Net EProd', '-')
        mean_sco2_net_hprod = econ_dict.get('Mean sCO2 Net HProd', '-')
        mean_sco2_net_eprod = econ_dict.get('Mean sCO2 Net EProd', '-')
        

        # Check if economic calculations failed but thermal results are available
        if 'error_codes' in econ_dict and econ_dict['error_codes']:
            # Economic calculations failed, but we can still show thermal results
            pass

        
        results = [lcoh_h2o, lcoe_h2o, mean_h2o_tout, mean_h2o_pout, 
                   mean_h2o_net_hprod, mean_h2o_net_eprod,
                   lcoh_sco2, lcoe_sco2, mean_sco2_tout, mean_sco2_pout,
                   mean_sco2_net_hprod, mean_sco2_net_eprod]
        
    result_names = summary_tbl['Result'].to_list()

    # Convert fixed values to imperial units if needed
    if units == "imperial":
        # CovHDF5 model has different fixed parameters
        if model == "CovHDF5":
            # Hard-coded CovHDF5 fixed values converted to imperial
            fixed_values = [
                20 * 145.038,  # Inlet Pressure: 20 MPa to psi
                round((300.0 - 273.15) * (9.0/5.0) + 32.0, 2),  # Ambient Temperature: 300K to °F
                round((298.15 - 273.15) * (9.0/5.0) + 32.0, 2),  # Surface Temperature: 298.15K to °F
                0.025 * 0.0393701,  # Pipe Roughness: 0.025 mm to inches
                2750 * 0.062428,  # Rock Density: 2750 kg/m³ to lb/ft³
                790 * 0.238846,  # Rock Specific Heat: 790 J/kg·K to BTU/lb·°F
                0.4445 * 3.28084,  # Borehole Diameter: 0.4445 m to ft
                0.1,  # Porosity (dimensionless)
                3.05 * 0.577789,  # Rock Thermal Conductivity: 3.05 W/m·K to BTU/(hr·ft·°F)
                0.015,  # Operation and Maintenance Cost (dimensionless)
                0.8,  # Pump Efficiency (dimensionless)
                0.9,  # Turbine Isentropic Efficiency (dimensionless)
                0.98,  # Generator Conversion Efficiency (dimensionless)
                (293.15 - 273.15) * (9.0/5.0) + 32.0,  # Dead-State Temperature: 293.15K to °F
                100000 * (9.0/5.0) + 32.0,  # Dead-State Temperature: 100000°C to °F (note: this seems like an error in original)
                0.1,  # Electricity Rate in Direct-Use (dimensionless)
                '-',  # Empty
                '-'   # Empty
            ]
            
            # Update parameter names for CovHDF5 (no thermal conductivity slider, has permeability slider)
            param_names = [
                "Mass Flow Rate (lb/s)",
                "Horizontal Extent (ft)", 
                "Drilling Depth (ft)",
                "Geothermal Gradient (°F/ft)",
                "Permeability (HWR)",
                "Injection Temperature (°F)",
                "-",  # No Rock Thermal Conductivity slider for CovHDF5
                "Drilling Cost ($/ft)",
                "Discount Rate (%)",
                "Lifetime (years)",
                "Plant CAPEX ($/kWt)",
                "Plant CAPEX ($/kWe)",
                "Pre-Cooling (°F)",
                "Turbine Outlet Pressure (psi)",
                "Select Interpolation",
                "Select Heat-Exchanger Design", 
                "Select Fluid",
                "-"
            ]
            
            fixed_param_names = [
                "Inlet Pressure (psi)",
                "Ambient Temperature (°F)",
                "Surface Temperature (°F)", 
                "Pipe Roughness (in)",
                "Rock Density (lb/ft³)",
                "Rock Specific Heat (BTU/lb·°F)",
                "Borehole Diameter (ft)",
                "Porosity",
                "Rock Thermal Conductivity (BTU/(hr·ft·°F))",
                "Operation and Maintence Cost of Plant",
                "Pump Efficiency For Circulation Pump",
                "Turbine Isentropic Efficency",
                "Generator Conversion Efficiency",
                "Dead-State Temperature (°F)",
                "Dead-State Temperature (°F)",
                "Electricity Rate in Direct-Use",
                "-",
                "-"
            ]
        else:
            # Standard HDF5 and SBT fixed values converted to imperial
            fixed_values = [
                20 * 145.038,  # Inlet Pressure: 20 MPa to psi
                26.85 * (9.0/5.0) + 32.0,  # Ambient Temperature: 26.85°C to °F
                25 * (9.0/5.0) + 32.0,  # Surface Temperature: 25°C to °F
                0.025 * 0.0393701,  # Pipe Roughness: 0.025 mm to inches
                2750 * 0.062428,  # Rock Density: 2750 kg/m³ to lb/ft³
                790 * 0.238846,  # Rock Specific Heat: 790 J/kg·K to BTU/lb·°F
                1,  # Inner Pipe / Annulus Area Ratio (dimensionless)
                0.0192 * 3.28084,  # Inner Pipe Thickness: 0.0192 m to ft
                0.06 * 0.577789,  # Inner Pipe Thermal Conductivity: 0.06 W/m·K to BTU/(hr·ft·°F)
                0.015,  # Operation and Maintenance Cost (dimensionless)
                0.8,  # Pump Efficiency (dimensionless)
                0.9,  # Turbine Isentropic Efficiency (dimensionless)
                0.98,  # Generator Conversion Efficiency (dimensionless)
                293.15 * (9.0/5.0) + 32.0,  # Dead-State Temperature: 293.15°C to °F
                100000 * (9.0/5.0) + 32.0,  # Dead-State Temperature: 100000°C to °F
                0.1,  # Electricity Rate in Direct-Use (dimensionless)
                '-',  # Empty
                '-'   # Empty
            ]
            
            # Update parameter names to show imperial units
            param_names = [
                "Mass Flow Rate (lb/s)",
                "Horizontal Extent (ft)", 
                "Drilling Depth (ft)",
                "Geothermal Gradient (°F/ft)",
                "Borehole Diameter (ft)",
                "Injection Temperature (°F)",
                "Rock Thermal Conductivity (BTU/(hr·ft·°F))",
                "Drilling Cost ($/ft)",
                "Discount Rate (%)",
                "Lifetime (years)",
                "Plant CAPEX ($/kWt)",
                "Plant CAPEX ($/kWe)",
                "Pre-Cooling (°F)",
                "Turbine Outlet Pressure (psi)",
                "Select Interpolation",
                "Select Heat-Exchanger Design", 
                "Select Fluid",
                "-"
            ]
            
            fixed_param_names = [
                "Inlet Pressure (psi)",
                "Ambient Temperature (°F)",
                "Surface Temperature (°F)", 
                "Pipe Roughness (in)",
                "Rock Density (lb/ft³)",
                "Rock Specific Heat (BTU/lb·°F)",
                "Inner Pipe / Annulus Area Ratio",
                "Inner Pipe Thickness (ft)",
                "Inner Pipe Thermal Conductivity (BTU/(hr·ft·°F))",
                "Operation and Maintence Cost of Plant",
                "Pump Efficiency For Circulation Pump",
                "Turbine Isentropic Efficency",
                "Generator Conversion Efficiency",
                "Dead-State Temperature (°F)",
                "Dead-State Temperature (°F)",
                "Electricity Rate in Direct-Use",
                "-",
                "-"
            ]
        
        # Update result names to show imperial units
        result_names = [
            "H2O LCOH ($/MWh th)",
            "H2O LCOE ($/MWh e)",
            "H2O 40-Year Average Outlet Temperature (°F)",
            "H2O 40-Year Average Outlet Pressure (psi)",
            f"H2O {Lifetime}-Year Average Heat (MWt)",
            f"H2O {Lifetime}-Year Average Electricity Production (MWe)",
            "sCO2 LCOH ($/MWh th)",
            "sCO2 LCOE ($/ MWh e)",
            "sCO2 40-Year Average Outlet Temperature (°F)",
            "sCO2 40-Year Average Outlet Pressure (psi)",
            f"sCO2 {Lifetime}-Year Average Heat (MWt)",
            f"sCO2 {Lifetime}-Year Average Electricity Production (MWe)",
            "-",
            "-",
            "-",
            "-",
            "-",
            "-"
        ]
        
        # Round fixed values to 4 decimal places
        fixed_values = [round(val, 4) if isinstance(val, (int, float)) else val for val in fixed_values]
        
    else:
        # Metric values
        if model == "CovHDF5":
            # CovHDF5 metric fixed values
            fixed_values = [
                20,  # Inlet Pressure: 20 MPa
                round(300.0 - 273.15, 2),  # Ambient Temperature: 300K to °C
                round(298.15 - 273.15, 2),  # Surface Temperature: 298.15K to °C
                0.025,  # Pipe Roughness: 0.025 mm
                2750,  # Rock Density: 2750 kg/m³
                790,  # Rock Specific Heat: 790 J/kg·K
                0.4445,  # Borehole Diameter: 0.4445 m
                0.1,  # Porosity (dimensionless)
                3.05,  # Rock Thermal Conductivity: 3.05 W/m·K
                0.015,  # Operation and Maintenance Cost (dimensionless)
                0.8,  # Pump Efficiency (dimensionless)
                0.9,  # Turbine Isentropic Efficiency (dimensionless)
                0.98,  # Generator Conversion Efficiency (dimensionless)
                293.15,  # Dead-State Temperature: 293.15°C
                100000,  # Dead-State Temperature: 100000°C
                0.1,  # Electricity Rate in Direct-Use (dimensionless)
                '-',  # Empty
                '-'   # Empty
            ]
            
            param_names = [
                "Mass Flow Rate (kg/s)",
                "Horizontal Extent (m)", 
                "Drilling Depth (m)",
                "Geothermal Gradient (K/m)",
                "Permeability (HWR)",
                "Injection Temperature (˚C)",
                "-",  # No Rock Thermal Conductivity slider for CovHDF5
                "Drilling Cost ($/m)",
                "Discount Rate (%)",
                "Lifetime (years)",
                "Plant CAPEX ($/kWt)",
                "Plant CAPEX ($/kWe)",
                "Pre-Cooling (˚C)",
                "Turbine Outlet Pressure (bar)",
                "Select Interpolation",
                "Select Heat-Exchanger Design", 
                "Select Fluid",
                "-"
            ]
            
            fixed_param_names = [
                "Inlet Pressure (MPa)",
                "Ambient Temperature (˚C)",
                "Surface Temperature (˚C)", 
                "Pipe Roughness (mm)",
                "Rock Density (kg/m3)",
                "Rock Specific Heat (J/kg-K)",
                "Borehole Diameter (m)",
                "Porosity",
                "Rock Thermal Conductivity (W/m-K)",
                "Operation and Maintence Cost of Plant",
                "Pump Efficiency For Circulation Pump",
                "Turbine Isentropic Efficency",
                "Generator Conversion Efficiency",
                "Dead-State Temperature (˚C)",
                "Dead-State Temperature (˚C)",
                "Electricity Rate in Direct-Use",
                "-",
                "-"
            ]
        else:
            # Standard HDF5 and SBT - use original CSV values
            param_names = summary_tbl['Variable Parameter'].to_list()
            fixed_param_names = summary_tbl['Fixed Parameter'].to_list()
            fixed_values = summary_tbl['Fixed Value'].to_list()
    
    result_names = [name.replace('X', str(Lifetime)) for name in result_names]
    
    # Fix column lengths to match across all columns
    row_count = max(
        len(param_names),
        len(variable_input_values),
        len(fixed_param_names),
        len(fixed_values),
        len(result_names),
        len(results),
    )

    def _pad(lst, target, fill='-'):
        return list(lst) + [fill] * (target - len(lst))

    param_names = _pad(param_names, row_count)
    variable_input_values = _pad(variable_input_values, row_count)
    fixed_param_names = _pad(fixed_param_names, row_count)
    fixed_values = _pad(fixed_values, row_count)
    result_names = _pad(result_names, row_count)
    results = _pad(results, row_count)

    # Create the table
    list_of_lists = [param_names, variable_input_values, 
                     fixed_param_names, fixed_values, 
                     result_names, results]

    fig_table = go.Figure(data=[go.Table(
                                columnwidth = [200,75,200,75,200,75],
                                header=dict(values=summary_tbl.columns,
                                            line_color=[greyyellow, greyyellow, greyyellow, greyyellow, greyblue, greyblue],
                                            fill_color=[greyyellow, greyyellow, greyyellow, greyyellow, greyblue, greyblue],
                                            font=dict(color='black', size=14)), 
                                cells=dict(values=list_of_lists,
                                           line_color=[greygrey, greygrey, greygrey, greygrey, greyblue, greyblue],
                                           fill_color=[white, white, white, white, greybluel, greybluel],
                                           font=dict(color='black', size=12),
                                           height=30))])
    
    fig_table.update_layout(
        margin=dict(l=0, r=0, t=0, b=0),
        height=1100
    )

    data_dict_of_lists = {summary_tbl.columns[0]: summary_tbl['Variable Parameter'].to_list(),
	                      "Variable Value": variable_input_values,
	                      summary_tbl.columns[2]: summary_tbl['Fixed Parameter'].to_list(),
	                      "Fixed Value": summary_tbl['Fixed Value'].to_list(),
	                      summary_tbl.columns[4]: result_names,
	                      "Result Value": results
                    }

    return fig_table, data_dict_of_lists



