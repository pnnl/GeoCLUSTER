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
                             interp_time, case, fluid,
                             thermal_dict, econ_dict):

    # -----------------------------------------------------------------------------------------------------------------
    # Creates Plotly a plotly table. Then, returns the figure and data as a dictionary of lits for data downloading.
    # -----------------------------------------------------------------------------------------------------------------

    variable_input_values = [mdot, L2, L1, grad, D, Tinj, k, Drilling_cost_per_m, Discount_rate, Lifetime, 
                             Direct_use_heat_cost_per_kWth, Power_plant_cost_per_kWe, Pre_Cooling_Delta_T, Turbine_outlet_pressure, 
                             interp_time, case, fluid]

    results = [econ_dict['LCOH H2O'], econ_dict['LCOE H2O'], thermal_dict['Mean H2O Tout'], thermal_dict['Mean H2O Pout'], 
               econ_dict['Mean H2O Net HProd'], econ_dict['Mean H2O Net EProd'],
               econ_dict['LCOH sCO2'], econ_dict['LCOE sCO2'], thermal_dict['Mean sCO2 Tout'], thermal_dict['Mean sCO2 Pout'],
               econ_dict['Mean sCO2 Net HProd'], econ_dict['Mean sCO2 Net EProd']
               ]

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



