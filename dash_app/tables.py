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
                             thermal_dict, econ_dict):

    # -----------------------------------------------------------------------------------------------------------------
    # Creates Plotly a plotly table. Then, returns the figure and data as a dictionary of lits for data downloading.
    # -----------------------------------------------------------------------------------------------------------------

    variable_input_values = [mdot, L2, L1, grad, D, Tinj, k, Drilling_cost_per_m, Discount_rate, Lifetime, 
                             Direct_use_heat_cost_per_kWth, Power_plant_cost_per_kWe, Pre_Cooling_Delta_T, Turbine_outlet_pressure, 
                             interp_time, case, fluid]

    # Handle different data structures for HDF5 vs SBT models
    if model == "HDF5":
        # Original HDF5 data structure
        results = [econ_dict['LCOH H2O'], econ_dict['LCOE H2O'], thermal_dict['Mean H2O Tout'], thermal_dict['Mean H2O Pout'], 
                   econ_dict['Mean H2O Net HProd'], econ_dict['Mean H2O Net EProd'],
                   econ_dict['LCOH sCO2'], econ_dict['LCOE sCO2'], thermal_dict['Mean sCO2 Tout'], thermal_dict['Mean sCO2 Pout'],
                   econ_dict['Mean sCO2 Net HProd'], econ_dict['Mean sCO2 Net EProd']
                   ]
    else:
        # SBT model data structure - calculate mean values from arrays
        import numpy as np
        
        # Get economic values - use '-' if economic calculations failed
        lcoh_h2o = econ_dict.get('LCOH H2O', '-')
        lcoe_h2o = econ_dict.get('LCOE H2O', '-')
        mean_h2o_net_hprod = econ_dict.get('Mean H2O Net HProd', '-')
        mean_h2o_net_eprod = econ_dict.get('Mean H2O Net EProd', '-')
        
        lcoh_sco2 = econ_dict.get('LCOH sCO2', '-')
        lcoe_sco2 = econ_dict.get('LCOE sCO2', '-')
        mean_sco2_net_hprod = econ_dict.get('Mean sCO2 Net HProd', '-')
        mean_sco2_net_eprod = econ_dict.get('Mean sCO2 Net EProd', '-')
        
        # If economic calculations failed, show thermal results anyway
        if 'error_codes' in econ_dict and econ_dict['error_codes']:
            print(f"DEBUG TABLES: Economic calculations failed, but showing thermal results")
        
        # Calculate mean thermal values from arrays
        mean_h2o_tout = '-'
        mean_h2o_pout = '-'
        mean_sco2_tout = '-'
        mean_sco2_pout = '-'
        
        # Check if we have TandP data for SBT models
        if 'TandP-data' in thermal_dict:
            tandp_data = thermal_dict['TandP-data']
            print(f"DEBUG TABLES: Processing TandP data for {model}")
            print(f"DEBUG TABLES: H2O_Tout type: {type(tandp_data.get('H2O_Tout'))}, value: {tandp_data.get('H2O_Tout')}")
            
            # Calculate mean temperature and pressure for H2O
            if 'H2O_Tout' in tandp_data and tandp_data['H2O_Tout'] is not None:
                h2o_tout_array = np.array(tandp_data['H2O_Tout'])
                print(f"DEBUG TABLES: H2O_Tout array length: {len(h2o_tout_array)}, first few values: {h2o_tout_array[:5] if len(h2o_tout_array) > 0 else 'empty'}")
                # Check if array is not empty and has valid data
                if len(h2o_tout_array) > 0 and not (h2o_tout_array == None).all():
                    try:
                        mean_h2o_tout = round(np.mean(h2o_tout_array - 273.15), 2)  # Convert from K to °C
                        print(f"DEBUG TABLES: Calculated mean H2O Tout: {mean_h2o_tout}")
                    except Exception as e:
                        print(f"DEBUG TABLES: Error calculating H2O Tout mean: {e}")
                        mean_h2o_tout = '-'
                else:
                    print(f"DEBUG TABLES: H2O_Tout array is empty or all None")
                    mean_h2o_tout = '-'
            else:
                print(f"DEBUG TABLES: H2O_Tout not found or is None")
                mean_h2o_tout = '-'
            
            if 'H2O_Pout' in tandp_data and tandp_data['H2O_Pout'] is not None:
                h2o_pout_array = np.array(tandp_data['H2O_Pout'])
                # Check if array is not empty and has valid data
                if len(h2o_pout_array) > 0 and not (h2o_pout_array == None).all():
                    try:
                        mean_h2o_pout = round(np.mean(h2o_pout_array / 1000000), 2)  # Convert from Pa to MPa
                    except:
                        mean_h2o_pout = '-'
            
            # Calculate mean temperature and pressure for sCO2
            if 'sCO2_Tout' in tandp_data and tandp_data['sCO2_Tout'] is not None:
                sco2_tout_array = np.array(tandp_data['sCO2_Tout'])
                # Check if array is not empty and has valid data
                if len(sco2_tout_array) > 0 and not (sco2_tout_array == None).all():
                    try:
                        mean_sco2_tout = round(np.mean(sco2_tout_array - 273.15), 2)  # Convert from K to °C
                    except:
                        mean_sco2_tout = '-'
            
            if 'sCO2_Pout' in tandp_data and tandp_data['sCO2_Pout'] is not None:
                sco2_pout_array = np.array(tandp_data['sCO2_Pout'])
                # Check if array is not empty and has valid data
                if len(sco2_pout_array) > 0 and not (sco2_pout_array == None).all():
                    try:
                        mean_sco2_pout = round(np.mean(sco2_pout_array / 1000000), 2)  # Convert from Pa to MPa
                    except:
                        mean_sco2_pout = '-'
        else:
            # Fallback to thermal_dict if TandP-data is not available
            mean_h2o_tout = thermal_dict.get('Mean H2O Tout', '-')
            mean_h2o_pout = thermal_dict.get('Mean H2O Pout', '-')
            mean_sco2_tout = thermal_dict.get('Mean sCO2 Tout', '-')
            mean_sco2_pout = thermal_dict.get('Mean sCO2 Pout', '-')
        
        results = [lcoh_h2o, lcoe_h2o, mean_h2o_tout, mean_h2o_pout, 
                   mean_h2o_net_hprod, mean_h2o_net_eprod,
                   lcoh_sco2, lcoe_sco2, mean_sco2_tout, mean_sco2_pout,
                   mean_sco2_net_hprod, mean_sco2_net_eprod]
        
        print(f"DEBUG TABLES: Final results list: {results}")
        print(f"DEBUG TABLES: mean_h2o_tout = {mean_h2o_tout}, mean_h2o_pout = {mean_h2o_pout}")
        print(f"DEBUG TABLES: mean_sco2_tout = {mean_sco2_tout}, mean_sco2_pout = {mean_sco2_pout}")

    result_names = summary_tbl['Result'].to_list()
    result_names = [name.replace('X', str(Lifetime)) for name in result_names]

    list_of_lists = [summary_tbl['Variable Parameter'].to_list(), 
                       variable_input_values,
                       summary_tbl['Fixed Parameter'].to_list(),
                       summary_tbl['Fixed Value'].to_list(),
                       result_names,
                       results,
                       ]

    fig_table = go.Figure(data=[go.Table(
                                columnwidth = [200,75,200,75,200,75],
                                header=dict(values=summary_tbl.columns,
                                            line_color=[greyyellow, greyyellow, greyyellow, greyyellow, greyblue, greyblue],
                                            fill_color=[greyyellow, greyyellow, greyyellow, greyyellow, greyblue, greyblue],
                                            font=dict(color='black', size=14)), 
                                cells=dict(values=list_of_lists,
                                           line_color=[greygrey, greygrey, greygrey, greygrey, greyblue, greyblue],
                                           fill_color=[white, white, white, white, greybluel, greybluel],
                                           align='center',
                                           height=22,
                                           font=dict(color='black', size=13) )
                                )
                     ])
    fig_table.update_layout(height=1100, margin=dict(r=10, l=10, t=0, b=0))

    print(f"DEBUG TABLES: Table figure created successfully")
    print(f"DEBUG TABLES: list_of_lists structure: {len(list_of_lists)} lists")
    print(f"DEBUG TABLES: Results list (6th element): {list_of_lists[5] if len(list_of_lists) > 5 else 'Not found'}")
    print(f"DEBUG TABLES: First few results: {list_of_lists[5][:6] if len(list_of_lists) > 5 and len(list_of_lists[5]) > 6 else 'Not available'}")

    data_dict_of_lists = {summary_tbl.columns[0]: summary_tbl['Variable Parameter'].to_list(),
	                      "Variable Value": variable_input_values,
	                      summary_tbl.columns[2]: summary_tbl['Fixed Parameter'].to_list(),
	                      "Fixed Value": summary_tbl['Fixed Value'].to_list(),
	                      summary_tbl.columns[4]: result_names,
	                      "Result Value": results
                    }

    return fig_table, data_dict_of_lists



