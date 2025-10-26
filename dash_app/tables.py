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

    variable_input_values = [mdot, L2, L1, grad, D, Tinj, k, Drilling_cost_per_m, Discount_rate, Lifetime, 
                             Direct_use_heat_cost_per_kWth, Power_plant_cost_per_kWe, Pre_Cooling_Delta_T, Turbine_outlet_pressure, 
                             interp_time, case, fluid]

    # Handle different data structures for HDF5 vs SBT models
    if model == "HDF5" or model == "CovHDF5":
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
                # SBT V1.0 doesn't provide pressure data, so set to '-'
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
            # For SBT models without TandP-data or HDF5 models, use existing thermal_dict values and convert if needed
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

        # Get economic values and handle None/9999 values
        def format_value(value):
            # Check for various 9999 formats
            if (value is None or 
                value == 9999 or 
                value == 9999.0 or 
                value == "9999.00" or 
                value == "9999" or 
                str(value) == "9999.00" or
                str(value) == "9999" or
                (isinstance(value, (int, float)) and abs(value - 9999) < 0.01)):
                return '-'
            return value
            
        lcoh_h2o = format_value(econ_dict.get('LCOH H2O', '-'))
        lcoe_h2o = format_value(econ_dict.get('LCOE H2O', '-'))
        lcoh_sco2 = format_value(econ_dict.get('LCOH sCO2', '-'))
        lcoe_sco2 = format_value(econ_dict.get('LCOE sCO2', '-'))
        mean_h2o_net_hprod = format_value(econ_dict.get('Mean H2O Net HProd', '-'))
        mean_h2o_net_eprod = format_value(econ_dict.get('Mean H2O Net EProd', '-'))
        mean_sco2_net_hprod = format_value(econ_dict.get('Mean sCO2 Net HProd', '-'))
        mean_sco2_net_eprod = format_value(econ_dict.get('Mean sCO2 Net EProd', '-'))
        
        # Also format thermal values
        mean_h2o_tout = format_value(mean_h2o_tout)
        mean_h2o_pout = format_value(mean_h2o_pout)
        mean_sco2_tout = format_value(mean_sco2_tout)
        mean_sco2_pout = format_value(mean_sco2_pout)
        
        results = [lcoh_h2o, lcoe_h2o, mean_h2o_tout, mean_h2o_pout, 
                   mean_h2o_net_hprod, mean_h2o_net_eprod,
                   lcoh_sco2, lcoe_sco2, mean_sco2_tout, mean_sco2_pout,
                   mean_sco2_net_hprod, mean_sco2_net_eprod]

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

    data_dict_of_lists = {summary_tbl.columns[0]: summary_tbl['Variable Parameter'].to_list(),
	                      "Variable Value": variable_input_values,
	                      summary_tbl.columns[2]: summary_tbl['Fixed Parameter'].to_list(),
	                      "Fixed Value": summary_tbl['Fixed Value'].to_list(),
	                      summary_tbl.columns[4]: result_names,
	                      "Result Value": results
                    }

    return fig_table, data_dict_of_lists



