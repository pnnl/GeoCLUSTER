#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# standard system libraries
import sys

# data manipulation libraries
import numpy as np
from numpy import genfromtxt
from scipy.interpolate import interpn

# plotting libraries
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots

# sourced scripts
from reader import initialize_data, data_dict
from econ_support import create_teaobject
from plots_support import * 
import traceback

# -----------------------
# Read in data.
# -----------------------
u_sCO2, u_H2O, c_sCO2, c_H2O = initialize_data() # 3 GB of memory
param_dict = data_dict(u_sCO2, u_H2O, c_sCO2, c_H2O)

# -----------------------
# Global properties.
# -----------------------
to_kelvin_factor = 273.15

labels_cat = ['Working Fluid2', 'Working Fluid1'] # legend visibility and syncing 
labels = ['H2O','sCO2'] # legend visibility and syncing 
lw = 1.4 # line width

# --------------------------------------------------------------------------------------------------------------------
# THERMAL PERFORMANCE PLOTS
# --------------------------------------------------------------------------------------------------------------------

def param_nearest_init(arg_mdot, arg_L2, arg_L1, arg_grad, arg_D, arg_Tinj, arg_k):

    arg_mdot_v, arg_mdot_i = find_nearest(u_sCO2.mdot, arg_mdot) 
    arg_L2_v, arg_L2_i = find_nearest(u_sCO2.L2, arg_L2)
    arg_L1_v, arg_L1_i = find_nearest(u_sCO2.L1, arg_L1)
    arg_grad_v, arg_grad_i = find_nearest(u_sCO2.grad, arg_grad)
    arg_D_v, arg_D_i = find_nearest(u_sCO2.D, arg_D)
    arg_Tinj_v, arg_Tinj_i = find_nearest(u_sCO2.Tinj, arg_Tinj + to_kelvin_factor) # to kelvin
    arg_k_v, arg_k_i = find_nearest(u_sCO2.k, arg_k)

    return arg_mdot_v, arg_mdot_i, arg_L2_v, arg_L2_i, arg_L1_v, arg_L1_i, arg_grad_v, arg_grad_i, arg_D_v, arg_D_i, \
                arg_Tinj_v, arg_Tinj_i, arg_k_v, arg_k_i


def get_kWe_kWt_over_mass_or_time(case, fluid, point, arg_L2_i, arg_L1_i, arg_grad_i, arg_D_i, arg_Tinj_i, arg_k_i):

    # --------------------------------------------------------------------------------------------------------------
    # For both H2O and sCO2, calculate:
    #       - kWe and kWt averages over mass flow rate & 
    #       - instantaneous kWe and kWt overtime 
    #
    # Outputs can be used for when interp_time == True, the following exists:
    #           sCO2_kWe_avg, sCO2_kWt_avg, H2O_kWe_avg, H2O_kWt_avg
    # Then, the following also get overriden when interp_time == True:
    #           sCO2_kWe, sCO2_kWt, H2O_kWe, H2O_kWt
    # 
    # When interp_time == False then all the output will exist and be plotted.
    # --------------------------------------------------------------------------------------------------------------

    # need initialization here
    sCO2_kWe_avg = None
    sCO2_kWt_avg = None 
    H2O_kWe_avg = None
    H2O_kWt_avg = None
    error_messages_d = {}

    if case == "utube":


        if fluid == "H2O" or fluid == "All":
            H2O_kWe_avg = u_H2O.kWe_avg[:, arg_L2_i, arg_L1_i, arg_grad_i, arg_D_i, arg_Tinj_i, arg_k_i]
            H2O_kWt_avg = u_H2O.kWt_avg[:, arg_L2_i, arg_L1_i, arg_grad_i, arg_D_i, arg_Tinj_i, arg_k_i]

        if fluid == "sCO2" or fluid == "All":
            sCO2_kWe_avg = u_sCO2.kWe_avg[:, arg_L2_i, arg_L1_i, arg_grad_i, arg_D_i, arg_Tinj_i, arg_k_i]
            sCO2_kWt_avg = u_sCO2.kWt_avg[:, arg_L2_i, arg_L1_i, arg_grad_i, arg_D_i, arg_Tinj_i, arg_k_i]

    if case == "coaxial":


        if fluid == "H2O" or fluid == "All":
            H2O_kWe_avg = c_H2O.kWe_avg[:, arg_L2_i, arg_L1_i, arg_grad_i, arg_D_i, arg_Tinj_i, arg_k_i]
            H2O_kWt_avg = c_H2O.kWt_avg[:, arg_L2_i, arg_L1_i, arg_grad_i, arg_D_i, arg_Tinj_i, arg_k_i]

        if fluid == "sCO2" or fluid == "All":
            sCO2_kWe_avg = c_sCO2.kWe_avg[:, arg_L2_i, arg_L1_i, arg_grad_i, arg_D_i, arg_Tinj_i, arg_k_i]
            sCO2_kWt_avg = c_sCO2.kWt_avg[:, arg_L2_i, arg_L1_i, arg_grad_i, arg_D_i, arg_Tinj_i, arg_k_i]

    return sCO2_kWe_avg, sCO2_kWt_avg, H2O_kWe_avg, H2O_kWt_avg, error_messages_d


def generate_subsurface_lineplots(interp_time, fluid, case, arg_mdot, arg_L2, arg_L1, arg_grad, arg_D, arg_Tinj, arg_k, scale, model,
            Tsurf, c_m, rho_m, 
            # radius_vertical, radius_lateral, n_laterals, lateral_flow, lateral_multiplier,
            Diameter1, Diameter2, PipeParam3, PipeParam4, PipeParam5,
            mesh, accuracy, mass_mode, temp_mode):

    # -----------------------------------------------------------------------------------------------------------------
    # Creates Plotly with 5 subplots:
    #
    #       kWe vs. mass flow rate
    #       kWt vs. mass flow rate
    #       Temperature vs. Time (years)
    #       Pressure vs. Time (years)
    #       kWe vs. time
    #
    # -----------------------------------------------------------------------------------------------------------------

    sCO2_Tout = sCO2_Pout = H2O_Tout = H2O_Pout = sCO2_kWe = sCO2_kWt = H2O_kWe = H2O_kWt = None

    if model == "HDF5":
        sbt_version = 0
    elif model == "SBT V1.0":
        sbt_version = 1
    else:
        sbt_version = 0

    mean_H2O_Tout = "-"
    mean_H2O_Pout = "-"
    mean_sCO2_Tout = "-"
    mean_sCO2_Pout = "-"

    error_messages_dict = {}

    time = u_sCO2.time
    m_dot = u_sCO2.mdot

    mass_flow_rates_dict = {"Mass Flow Rate (kg/s)": m_dot}
    time_dict = {"Time (year)": time}

    arg_mdot_v, arg_mdot_i, arg_L2_v, arg_L2_i, arg_L1_v, arg_L1_i, arg_grad_v, arg_grad_i, arg_D_v, arg_D_i, \
                arg_Tinj_v, arg_Tinj_i, arg_k_v, arg_k_i = param_nearest_init(arg_mdot, arg_L2, arg_L1, arg_grad, arg_D, arg_Tinj, arg_k)

    point = (arg_mdot, arg_L2, arg_L1, arg_grad, arg_D, arg_Tinj + to_kelvin_factor, arg_k) # to kelvin

    # ** Average calculations
    sCO2_kWe_avg, sCO2_kWt_avg, H2O_kWe_avg, H2O_kWt_avg, error_messages_d = \
                            get_kWe_kWt_over_mass_or_time(case, fluid, point, arg_L2_i, arg_L1_i, arg_grad_i, arg_D_i, arg_Tinj_i, arg_k_i)

    error_messages_dict.update(error_messages_d)

    # print(sCO2_kWe_avg, sCO2_kWt_avg, H2O_kWe_avg, H2O_kWt_avg)

    is_blank_data = False

    if interp_time == "False":

        # this doesn't impact the kWe over time

        if case == "utube":

            if fluid == "H2O" or fluid == "All":
                H2O_Tout = u_H2O.Tout[arg_mdot_i, arg_L2_i, arg_L1_i, arg_grad_i, arg_D_i, arg_Tinj_i, arg_k_i, :]
                H2O_Pout = u_H2O.Pout[arg_mdot_i, arg_L2_i, arg_L1_i, arg_grad_i, arg_D_i, arg_Tinj_i, arg_k_i, :]

            if fluid == "sCO2" or fluid == "All":
                sCO2_Tout = u_sCO2.Tout[arg_mdot_i, arg_L2_i, arg_L1_i, arg_grad_i, arg_D_i, arg_Tinj_i, arg_k_i, :]
                sCO2_Pout = u_sCO2.Pout[arg_mdot_i, arg_L2_i, arg_L1_i, arg_grad_i, arg_D_i, arg_Tinj_i, arg_k_i, :]

        if case == "coaxial":

            if fluid == "H2O" or fluid == "All":
                H2O_Tout = c_H2O.Tout[arg_mdot_i, arg_L2_i, arg_L1_i, arg_grad_i, arg_D_i, arg_Tinj_i, arg_k_i, :]
                H2O_Pout = c_H2O.Pout[arg_mdot_i, arg_L2_i, arg_L1_i, arg_grad_i, arg_D_i, arg_Tinj_i, arg_k_i, :]

            if fluid == "sCO2" or fluid == "All":
                sCO2_Tout = c_sCO2.Tout[arg_mdot_i, arg_L2_i, arg_L1_i, arg_grad_i, arg_D_i, arg_Tinj_i, arg_k_i, :]
                sCO2_Pout = c_sCO2.Pout[arg_mdot_i, arg_L2_i, arg_L1_i, arg_grad_i, arg_D_i, arg_Tinj_i, arg_k_i, :]

    # print("interp time", interp_time)
    if interp_time == "True":

        # this doesn't impact the kWe and kWt averages

        if case == "utube":

            try:
                if fluid == "sCO2" or fluid == "All":
                    sCO2_Tout, sCO2_Pout, time = u_sCO2.interp_outlet_states(point, sbt_version,
                                                        Tsurf, c_m, rho_m, 
                                                        # radius_vertical, radius_lateral, n_laterals, lateral_flow, lateral_multiplier,
                                                        Diameter1, Diameter2, PipeParam3, PipeParam4, PipeParam5,
                                                        mesh, accuracy, mass_mode, temp_mode)
                    sCO2_kWe, sCO2_kWt = u_sCO2.interp_kW(point, sCO2_Tout, sCO2_Pout)
                if fluid == "H2O" or fluid == "All":
                    H2O_Tout, H2O_Pout, time = u_H2O.interp_outlet_states(point, sbt_version,
                                                        Tsurf, c_m, rho_m, 
                                                        # radius_vertical, radius_lateral, n_laterals, lateral_flow, lateral_multiplier,
                                                        Diameter1, Diameter2, PipeParam3, PipeParam4, PipeParam5,
                                                        mesh, accuracy, mass_mode, temp_mode)
                    H2O_kWe, H2O_kWt = u_H2O.interp_kW(point, H2O_Tout, H2O_Pout)

            except ValueError as e:
                sCO2_Tout, sCO2_Pout, H2O_Tout, H2O_Pout, sCO2_kWe, sCO2_kWt, H2O_kWe, H2O_kWt = blank_data()
                error_message = parse_error_message(e=e, e_name='Err SubRes3')
                error_messages_dict['Err SubRes3'] = error_message
                is_blank_data = True
                
        if case == "coaxial":

            try:
                if fluid == "sCO2" or fluid == "All":
                    sCO2_Tout, sCO2_Pout, time = c_sCO2.interp_outlet_states(point, sbt_version,
                                                        Tsurf, c_m, rho_m, 
                                                        # radius_vertical, radius_lateral, n_laterals, lateral_flow, lateral_multiplier,
                                                        Diameter1, Diameter2, PipeParam3, PipeParam4, PipeParam5,
                                                        mesh, accuracy, mass_mode, temp_mode)
                    sCO2_kWe, sCO2_kWt = c_sCO2.interp_kW(point, sCO2_Tout, sCO2_Pout)
                if fluid == "H2O" or fluid == "All":
                    H2O_Tout, H2O_Pout, time = c_H2O.interp_outlet_states(point, sbt_version,
                                                        Tsurf, c_m, rho_m, 
                                                        # radius_vertical, radius_lateral, n_laterals, lateral_flow, lateral_multiplier,
                                                        Diameter1, Diameter2, PipeParam3, PipeParam4, PipeParam5,
                                                        mesh, accuracy, mass_mode, temp_mode)
                    H2O_kWe, H2O_kWt = c_H2O.interp_kW(point,H2O_Tout, H2O_Pout )

            except ValueError as e:
                sCO2_Tout, sCO2_Pout, H2O_Tout, H2O_Pout, sCO2_kWe, sCO2_kWt, H2O_kWe,H2O_kWt = blank_data()
                error_message = parse_error_message(e=e, e_name='Err SubRes4')
                error_messages_dict['Err SubRes4'] = error_message
                is_blank_data = True
                


    # Figure generation 

    fig = make_subplots(rows=2, cols=3,
                        specs=[[{}, {}, None],
                        [{}, {}, {}]],
                        horizontal_spacing = 0.11,
                        vertical_spacing = 0.21
                        )

    if fluid == "sCO2" or fluid == "All":

        # kWe_avg and kWt_avg
        fig.add_trace(go.Scatter(x=m_dot, y=sCO2_kWe_avg,
                      hovertemplate='<b>Mass Flow Rate (kg/s)</b>: %{x:.1f}<br><b>Average kWt</b>: %{y:.3f} ',
                      line = dict(color='#c26219', width=lw),
                      legendgroup=labels_cat[1], name=labels[1], showlegend=False),
                      # margin=dict(pad=20),
                      row=1, col=1)
        fig.add_trace(go.Scatter(x=m_dot, y=sCO2_kWt_avg,
                      hovertemplate='<b>Mass Flow Rate (kg/s)</b>: %{x:.1f}<br><b>Average kWt</b>: %{y:.3f} ',
                      line = dict(color='orange', width=lw),
                      legendgroup=labels_cat[1], name=labels[1], showlegend=False),
                      row=1, col=2)

        # kWe vs. time
        fig.add_trace(go.Scatter(x=time, y=sCO2_kWe,
                      hovertemplate='<b>Time (year)</b>: %{x:.1f}<br><b>kWe</b>: %{y:.3f} ',
                      line = dict(color='black', width=lw), # #379dbf
                      legendgroup=labels_cat[1], name=labels[1], showlegend=True),
                      row=2, col=1)
        # temperature
        fig.add_trace(go.Scatter(x=time, y=sCO2_Tout - to_kelvin_factor, # to celsius
                      hovertemplate='<b>Time (year)</b>: %{x:.1f}<br><b>Outlet Temperature (˚C)</b>: %{y:.3f} ',
                      line = dict(color='royalblue', width=lw),
                      legendgroup=labels_cat[1], name=labels[1], showlegend=False),
                      row=2, col=2)
        # pressure
        fig.add_trace(go.Scatter(x=time, y=sCO2_Pout / 1000000,
                      hovertemplate='<b>Time (year)</b>: %{x:.1f}<br><b>Outlet Pressure (MPa)</b>: %{y:.3f} ',
                      line = dict(color='#16b8a2', width=lw),
                      legendgroup=labels_cat[1], name=labels[1], showlegend=False),
                      row=2, col=3)

        if is_blank_data:
            blank_canvas(fig=fig, row_n=2, col_n=1)
            blank_canvas(fig=fig, row_n=2, col_n=2)
            blank_canvas(fig=fig, row_n=2, col_n=3)

        mean_sCO2_Tout = round(np.mean(sCO2_Tout - to_kelvin_factor),2) # to celsius
        mean_sCO2_Pout = round(np.mean(sCO2_Pout / 1000000),2)

        if str(mean_sCO2_Tout) == "nan":
            mean_sCO2_Tout = "-"
        if str(mean_sCO2_Pout) == "nan":
            mean_sCO2_Pout = "-"

        mass_flow_rates_dict["sCO2 40-Year Average of Exergy (kWe)"] = sCO2_kWe_avg
        mass_flow_rates_dict["sCO2 40-Year Average of Available Thermal Output (kWt)"] = sCO2_kWt_avg

        time_dict["sCO2 Exergy (kWe)"] = sCO2_kWe
        time_dict["sCO2 Outlet Temperature (˚C)"] = sCO2_Tout - to_kelvin_factor # to celsius
        time_dict["sCO2 Outlet Pressure (MPa)"] = sCO2_Pout / 1000000

    if fluid == "H2O" or fluid == "All":

        # kWe_avg and kWt_avg
        fig.add_trace(go.Scatter(x=m_dot, y=H2O_kWe_avg,
                      hovertemplate='<b>Mass Flow Rate (kg/s)</b>: %{x:.1f}<br><b>Average kWe</b>: %{y:.3f} ',
                      line = dict(color='#c26219', width=lw, dash='dash'),
                      legendgroup=labels_cat[0], name=labels[0], showlegend=False),
                      # margin=dict(pad=20),
                      row=1, col=1)

        fig.add_trace(go.Scatter(x=m_dot, y=H2O_kWt_avg,
                      hovertemplate='<b>Mass Flow Rate (kg/s)</b>: %{x:.1f}<br><b>Average kWt</b>: %{y:.3f} ',
                      line = dict(color='orange', width=lw, dash='dash'),
                      legendgroup=labels_cat[0], name=labels[0], showlegend=False),
                      row=1, col=2)

        # kWe vs. time
        fig.add_trace(go.Scatter(x=time, y=H2O_kWe, 
                      hovertemplate='<b>Time (year)</b>: %{x:.1f}<br><b>kWe</b>: %{y:.3f} ',
                      line = dict(color='black', width=lw, dash='dash'), # #379dbf
                      legendgroup=labels_cat[0], name=labels[0], showlegend=True),
                      row=2, col=1)
        # temperature
        fig.add_trace(go.Scatter(x=time, y=H2O_Tout - to_kelvin_factor, # to celsius
                      hovertemplate='<b>Time (year)</b>: %{x:.1f}<br><b>Outlet Temperature (˚C)</b>: %{y:.3f} ',
                      line = dict(color='royalblue', width=lw, dash='dash'),
                      legendgroup=labels_cat[0], name=labels[0], showlegend=False),
                      row=2, col=2)
        # pressure
        fig.add_trace(go.Scatter(x=time, y=H2O_Pout / 1000000,
                      hovertemplate='<b>Time (year)</b>: %{x:.1f}<br><b>Outlet Pressure (MPa)</b>: %{y:.3f} ',
                      line = dict(color='#16b8a2', width=lw, dash='dash'),
                      legendgroup=labels_cat[0], name=labels[0], showlegend=False),
                      row=2, col=3)

        if is_blank_data:
            blank_canvas(fig=fig, row_n=2, col_n=1)
            blank_canvas(fig=fig, row_n=2, col_n=2)
            blank_canvas(fig=fig, row_n=2, col_n=3)

        mean_H2O_Tout = round(np.mean(H2O_Tout - to_kelvin_factor),2) # to celsius
        mean_H2O_Pout = round(np.mean(H2O_Pout / 1000000),2)

        if str(mean_H2O_Tout) == "nan":
            mean_H2O_Tout = "-"
        if str(mean_H2O_Pout) == "nan":
            mean_H2O_Pout = "-"

        mass_flow_rates_dict["H2O 40-Year Average of Exergy (kWe)"] = H2O_kWe_avg
        mass_flow_rates_dict["H2O 40-Year Average of Available Thermal Output (kWt)"] = H2O_kWt_avg

        time_dict["H2O Exergy (kWe)"] = H2O_kWe
        time_dict["H2O Outlet Temperature (˚C)"] = H2O_Tout - to_kelvin_factor # to celsius
        time_dict["H2O Outlet Pressure (MPa)"] = H2O_Pout / 1000000

    
    fig = update_layout_properties_subsurface_results(fig=fig, m_dot=m_dot, time=time, plot_scale=scale)

    forty_yr_means_dict = {'Mean H2O Tout': mean_H2O_Tout, 
                            'Mean H2O Pout': mean_H2O_Pout,
                            'Mean sCO2 Tout': mean_sCO2_Tout,
                            'Mean sCO2 Pout': mean_sCO2_Pout}
   
    TandP_dict = {"sCO2_Tout": sCO2_Tout,
                            "sCO2_Pout": sCO2_Pout,
                            "H2O_Tout": H2O_Tout,
                            "H2O_Pout": H2O_Pout,
                            "sCO2_kWe": sCO2_kWe,
                            "sCO2_kWt": sCO2_kWt,
                            "H2O_kWe": H2O_kWe,
                            "H2O_kWt": H2O_kWt,
                            "time": time,
                            "mdot": m_dot,
                            }

    
    return fig, forty_yr_means_dict, mass_flow_rates_dict, time_dict, error_messages_dict, TandP_dict



def generate_subsurface_contours(interp_time, fluid, case, param, arg_mdot, arg_L2, arg_L1, arg_grad, arg_D, arg_Tinj, arg_k):

    # -----------------------------------------------------------------------------------------------------------------
    # Creates Plotly with 4 subplot contours:
    #
    #       <parameter selected> vs. mass flow rate vs. kWe
    #       <parameter selected> vs. mass flow rate vs. kWt
    #       <parameter selected> vs. mass flow rate vs. Tout
    #       <parameter selected> vs. mass flow rate vs. Pout
    #
    # -----------------------------------------------------------------------------------------------------------------

    error_messages_dict = {}

    # initializer to prevent errors
    if fluid == "All":
        fluid = "H2O"

    # contour elements
    param_y = param_dict[(case, fluid, param)]
    mdot_x = param_dict[(case, fluid, "mdot")]
    mdot_ij, param_ij = np.meshgrid(mdot_x, param_y, indexing='ij') # returns coordinate matrices from coordinate vectors

    arg_mdot_v, arg_mdot_i, arg_L2_v, arg_L2_i, arg_L1_v, arg_L1_i, arg_grad_v, arg_grad_i, arg_D_v, arg_D_i, \
                arg_Tinj_v, arg_Tinj_i, arg_k_v, arg_k_i = param_nearest_init(arg_mdot, arg_L2, arg_L1, arg_grad, arg_D, arg_Tinj, arg_k)

    # initial time conditions
    # arg_i = 160 # 4
    # arg_v = 40 # 20
    arg_i = 160 - (1*101) 
    arg_v = 40 - (0.25*101)
    if interp_time == "False":

        # print('FALSE')
        # print(param)
        # print(param_y.shape) # (20,)
        # print(mdot_x.shape) # (26,)
        # print(mdot_ij[:,0].shape) # (26,)
        # print(param_ij[0,:].shape) # (20,)
        
        if param == "Horizontal Extent (m)":
            arg_L2_i = slice(None)
        if param == "Vertical Extent (m)":
            arg_L1_i = slice(None)
        if param == "Geothermal Gradient (K/m)":
            arg_grad_i = slice(None)
        if param == "Borehole Diameter (m)":
            arg_D_i = slice(None)
        if param == "Injection Temperature (˚C)":
            arg_Tinj_i = slice(None)
        if param == "Rock Thermal Conductivity (W/m-K)":
            arg_k_i = slice(None)

        if fluid == "sCO2" and case == "utube":

            kWe_avg_flipped = np.transpose(u_sCO2.kWe_avg[:, arg_L2_i, arg_L1_i, arg_grad_i, arg_D_i, arg_Tinj_i, arg_k_i])
            kWt_avg_flipped = np.transpose(u_sCO2.kWt_avg[:, arg_L2_i, arg_L1_i, arg_grad_i, arg_D_i, arg_Tinj_i, arg_k_i])
            Tout_flipped = np.transpose(u_sCO2.Tout[:, arg_L2_i, arg_L1_i, arg_grad_i, arg_D_i, arg_Tinj_i, arg_k_i, arg_i])
            Pout_flipped = np.transpose(u_sCO2.Pout[:, arg_L2_i, arg_L1_i, arg_grad_i, arg_D_i, arg_Tinj_i, arg_k_i, arg_i])

        if fluid == "sCO2" and case == "coaxial":

            kWe_avg_flipped = np.transpose(c_sCO2.kWe_avg[:, arg_L2_i, arg_L1_i, arg_grad_i, arg_D_i, arg_Tinj_i, arg_k_i])
            kWt_avg_flipped = np.transpose(c_sCO2.kWt_avg[:, arg_L2_i, arg_L1_i, arg_grad_i, arg_D_i, arg_Tinj_i, arg_k_i])
            Tout_flipped = np.transpose(c_sCO2.Tout[:, arg_L2_i, arg_L1_i, arg_grad_i, arg_D_i, arg_Tinj_i, arg_k_i, arg_i])
            Pout_flipped = np.transpose(c_sCO2.Pout[:, arg_L2_i, arg_L1_i, arg_grad_i, arg_D_i, arg_Tinj_i, arg_k_i, arg_i])

        if fluid == "H2O" and case == "utube":

            kWe_avg_flipped = np.transpose(u_H2O.kWe_avg[:, arg_L2_i, arg_L1_i, arg_grad_i, arg_D_i, arg_Tinj_i, arg_k_i])
            kWt_avg_flipped = np.transpose(u_H2O.kWt_avg[:, arg_L2_i, arg_L1_i, arg_grad_i, arg_D_i, arg_Tinj_i, arg_k_i])
            Tout_flipped = np.transpose(u_H2O.Tout[:, arg_L2_i, arg_L1_i, arg_grad_i, arg_D_i, arg_Tinj_i, arg_k_i, arg_i])
            Pout_flipped = np.transpose(u_H2O.Pout[:, arg_L2_i, arg_L1_i, arg_grad_i, arg_D_i, arg_Tinj_i, arg_k_i, arg_i])

            # print(kWe_avg_flipped.shape) # (20, 26)
            # print('\n')

        if fluid == "H2O" and case == "coaxial":

            kWe_avg_flipped = np.transpose(c_H2O.kWe_avg[:, arg_L2_i, arg_L1_i, arg_grad_i, arg_D_i, arg_Tinj_i, arg_k_i])
            kWt_avg_flipped = np.transpose(c_H2O.kWt_avg[:, arg_L2_i, arg_L1_i, arg_grad_i, arg_D_i, arg_Tinj_i, arg_k_i])
            Tout_flipped = np.transpose(c_H2O.Tout[:, arg_L2_i, arg_L1_i, arg_grad_i, arg_D_i, arg_Tinj_i, arg_k_i, arg_i])
            Pout_flipped = np.transpose(c_H2O.Pout[:, arg_L2_i, arg_L1_i, arg_grad_i, arg_D_i, arg_Tinj_i, arg_k_i, arg_i])

    if interp_time == "True": # would this even make sense? It doesn't because can manipulate the other factors, when they should be fixed?

        # print('TRUE')
        # print(param)
        # print(param_y.shape) # (20,)
        # print(mdot_x.shape) # (26,)
        # print(mdot_ij[:,0].shape) # (26,)
        # print(param_ij[0,:].shape) # (20,)

        point = (arg_L2, arg_L1, arg_grad, arg_D, arg_Tinj + to_kelvin_factor, arg_k, arg_v) # to kelvin
        # print("point for contour", point)
        # point2 = (arg_mdot, arg_L2, arg_L1, arg_grad, arg_D, arg_Tinj + to_kelvin_factor, arg_k) # to kelvin

        if fluid == "sCO2" and case == "utube":
            try:
                sCO2_Tout, sCO2_Pout = u_sCO2.interp_outlet_states_contour(param, point) 
                sCO2_kWe, sCO2_kWt = u_sCO2.interp_kW_contour(param, point, sCO2_Tout, sCO2_Pout)
                kWe_avg_flipped, kWt_avg_flipped, Tout_flipped, Pout_flipped = rename_for_contour(sCO2_kWe, sCO2_kWt, sCO2_Tout, sCO2_Pout)

            except ValueError as e:
                kWe_avg_flipped, kWt_avg_flipped, Tout_flipped, Pout_flipped = blank_data_kW() 
                error_message = parse_error_message(e=e, e_name='Err SubContour1')
                error_messages_dict['Err SubContour1'] = error_message
        
        if fluid == "sCO2" and case == "coaxial":
            try:
                sCO2_Tout, sCO2_Pout = c_sCO2.interp_outlet_states_contour(param, point)
                sCO2_kWe, sCO2_kWt = c_sCO2.interp_kW_contour(param, point, sCO2_Tout, sCO2_Pout)
                kWe_avg_flipped, kWt_avg_flipped, Tout_flipped, Pout_flipped = rename_for_contour(sCO2_kWe, sCO2_kWt, sCO2_Tout, sCO2_Pout)

            except ValueError as e:
                kWe_avg_flipped, kWt_avg_flipped, Tout_flipped, Pout_flipped = blank_data_kW() 
                error_message = parse_error_message(e=e, e_name='Err SubContour2')
                error_messages_dict['Err SubContour2'] = error_message
        
        if fluid == "H2O" and case == "utube":
            try:
                H2O_Tout, H2O_Pout = u_H2O.interp_outlet_states_contour(param, point)
                H2O_kWe, H2O_kWt = u_H2O.interp_kW_contour(param, point, H2O_Tout, H2O_Pout)
                kWe_avg_flipped, kWt_avg_flipped, Tout_flipped, Pout_flipped = rename_for_contour(H2O_kWe, H2O_kWt, H2O_Tout, H2O_Pout)

            except ValueError as e:
                kWe_avg_flipped, kWt_avg_flipped, Tout_flipped, Pout_flipped = blank_data_kW() 
                error_message = parse_error_message(e=e, e_name='Err SubContour3')
                error_messages_dict['Err SubContour3'] = error_message
        
            # print(kWe_avg_flipped.shape) # should have an x and y right for the z component?
            # print('\n')

        if fluid == "H2O" and case == "coaxial":
            try:
                H2O_Tout, H2O_Pout = c_H2O.interp_outlet_states_contour(param, point)
                H2O_kWe, H2O_kWt = c_H2O.interp_kW_contour(param, point, H2O_Tout, H2O_Pout)
                kWe_avg_flipped, kWt_avg_flipped, Tout_flipped, Pout_flipped = rename_for_contour(H2O_kWe, H2O_kWt, H2O_Tout, H2O_Pout)

            except ValueError as e:
                kWe_avg_flipped, kWt_avg_flipped, Tout_flipped, Pout_flipped = blank_data_kW() 
                error_message = parse_error_message(e=e, e_name='Err SubContour4')
                error_messages_dict['Err SubContour4'] = error_message
        

    fig = make_subplots(rows=2, cols=2, start_cell="top-left",
                        horizontal_spacing = 0.12,
                        vertical_spacing = 0.12,
                        subplot_titles=('40-Year Average Exergy (kWe)', 
                                       '40-Year Average Thermal Output (kWt)',
                                       'Outlet Temperature (°C)',
                                       'Outlet Pressure (MPa)')
           )


    contour_formatting = dict(showlabels = True, #coloring ='heatmap'
                              labelfont = dict( size = 14,
                                                color = 'white'
                                                )
                              )
    fig.add_trace( 
        go.Contour(showscale=False,
            colorbar=dict(title='Average Exergy (kWe)'), #titleside='right'
            colorscale=colorscaleR,
            name="",
            hovertemplate='<b>Mass Flow Rate (kg/s)</b>: %{x:.1f}<br><b>Y</b>: %{y:.4f}<br><b>Average Exergy (kWe)</b>: %{z:.5f} ',
            contours=contour_formatting,
            z=kWe_avg_flipped, y=param_ij[0,:], x=mdot_ij[:,0]
        ), row=1, col=1 )
                  

    fig.add_trace( 
        go.Contour(showscale=False,
            colorbar=dict(title='Average Thermal Output (kWe)'),
            colorscale=colorscaleY,
            name="",
            hovertemplate='<b>Mass Flow Rate (kg/s)</b>: %{x:.1f}<br><b>Y</b>: %{y:.4f}<br><b>Average Thermal Output (kWt)</b>: %{z:.5f} ',
            contours=contour_formatting,
            z=kWt_avg_flipped, y=param_ij[0,:], x=mdot_ij[:,0]
        ), row=1, col=2 )


    fig.add_trace( 
        go.Contour(showscale=False,
            colorbar=dict(title='Outlet Temperature (°C)'),
            colorscale=colorscaleB,
            name="",
            hovertemplate='<b>Mass Flow Rate (kg/s)</b>: %{x:.1f}<br><b>Y</b>: %{y:.4f}<br><b>Outlet Temperature (°C)</b>: %{z:.5f} ',
            contours=contour_formatting,
            z=Tout_flipped - to_kelvin_factor, y=param_ij[0,:], x=mdot_ij[:,0] # to celsius
        ), row=2, col=1 )

    fig.add_trace( 
        go.Contour(showscale=False,
            colorbar=dict(title='Outlet Pressure (MPa)'),
            colorscale=colorscaleG,
            name="",
            hovertemplate='<b>Mass Flow Rate (kg/s)</b>: %{x:.1f}<br><b>Y</b>: %{y:.4f}<br><b>Outlet Pressure (MPa)</b>: %{z:.0f} ',
            contours=contour_formatting,
            z=Pout_flipped / 1000000, y=param_ij[0,:], x=mdot_ij[:,0]
        ), row=2, col=2 )

    fig = update_layout_properties_subsurface_contours(fig, param)
    
    return fig, error_messages_dict



# --------------------------------------------------------------------------------------------------------------------
# ECONOMIC PERFORMANCE PLOTS
# --------------------------------------------------------------------------------------------------------------------


def generate_econ_lineplots(TandP_dict,
                            interp_time, case, end_use, fluid,
                            mdot, L2, L1, grad, D, Tinj, k,
                            Drilling_cost_per_m, Discount_rate, Lifetime, 
                            Direct_use_heat_cost_per_kWth, Power_plant_cost_per_kWe, Pre_Cooling_Delta_T, Turbine_outlet_pressure,
                            scale,  
                            properties_H2O_pathname, 
                            properties_CO2v2_pathname, 
                            additional_properties_CO2v2_pathname,
                            tmatrix_pathname,
                            model,
                            is_plot_ts_check
                            ):

    # -----------------------------------------------------------------------------------------------------------------
    # Creates Plotly with 4 subplot lineplots:
    #
    #       heat prod vs. time
    #       annual heat prod vs. time
    #       electricity prod vs. time
    #       annual electricity prod vs. time
    #
    # -----------------------------------------------------------------------------------------------------------------

    if model == "HDF5":
        sbt_version = 0
    elif model == "SBT V1.0":
        sbt_version = 1
    else:
        sbt_version = 0

    lcoh_sCO2 = '-'
    lcoh_H2O = '-'
    lcoe_sCO2 = '-'
    lcoe_H2O = '-'
    
    mean_H2O_Net_HProd = '-'
    mean_H2O_Net_EProd = '-'
    mean_sCO2_Net_HProd = '-'
    mean_sCO2_Net_EProd = '-'

    econ_values_dict = {}
    error_messages_dict = {}

    fig = make_subplots(rows=2, cols=5,
                        specs=[[{'colspan': 2}, None, {'colspan': 2}, None, {"type": "table"}],
                                [{'colspan': 2}, None, {'colspan': 2}, None, {"type": "table"}]],
                        horizontal_spacing = 0.11
                        )

    fig = make_subplots(rows=3, cols=5,
                        specs=[[{'colspan': 2}, None, {'colspan': 2}, None, {"type": "table"}],
                                [{'colspan': 2}, None, {'colspan': 2}, None, {"type": "table"}],
                                [{'colspan': 2}, None, None, None, None]
                                ],
                        horizontal_spacing = 0.11
                        )

    # ts_fig = make_subplots(rows=1, cols=1,
    #                         horizontal_spacing = 0.11
    #                         )

    if end_use == "Heating" or end_use == "All":

        if end_use == "Heating":
            is_display_legend = True
        if end_use == "All":
            is_display_legend = True

        if fluid == "sCO2" or fluid == "All":

            try:
                teaobj_sCO2 = create_teaobject(TandP_dict,
                                                u_sCO2, u_H2O, c_sCO2, c_H2O,
                                                case, end_use, fluid, sbt_version,
                                                mdot, L2, L1, grad, D, Tinj, k,
                                                Drilling_cost_per_m, Discount_rate, Lifetime, 
                                                Direct_use_heat_cost_per_kWth, Power_plant_cost_per_kWe, Pre_Cooling_Delta_T, Turbine_outlet_pressure,
                                                properties_H2O_pathname, 
                                                properties_CO2v2_pathname, 
                                                additional_properties_CO2v2_pathname,
                                                is_heating=True)
                lcoh_sCO2 = format( teaobj_sCO2.LCOH, '.2f')

                # Heat Production 
                fig.add_trace(go.Scatter(x=teaobj_sCO2.Linear_time_distribution, y=teaobj_sCO2.Instantaneous_heat_production/1e3,
                              hovertemplate='<b>Time (year)</b>: %{x:.1f}<br><b>Heat Production (MWt)</b>: %{y:.3f} ',
                              line = dict(color='black', width=lw),
                              legendgroup=labels_cat[1], name=labels[1], showlegend=is_display_legend
                              ),
                              row=1, col=1)


                # Anuual Heat Production
                fig.add_trace(go.Bar(x=np.arange(1,teaobj_sCO2.Lifetime+1), y=teaobj_sCO2.Annual_heat_production/1e6,
                                    # hovertemplate='<b>Time (year)</b>: %{x:.1f}<br><b>Annual Heat Production (GWh)</b>: %{y:.3f} ',
                                    name=labels[1],
                                    showlegend=False, hovertemplate='<b>Annual Heat Production (GWh)</b>: %{y:.3f}'), 
                                row=1, col=3)

                fig.add_trace(go.Scatter(x=np.arange(1,teaobj_sCO2.Lifetime+1), y=teaobj_sCO2.Annual_heat_production/1e6,
                              hovertemplate='<b>Time (year)</b>: %{x:.1f}<br>',
                              line = dict(color='black', width=lw),
                              legendgroup=labels_cat[1], name=labels[1], showlegend=False
                              ),
                              row=1, col=3)

                # Table Data
                mean_sCO2_Net_HProd = round(np.mean(teaobj_sCO2.Instantaneous_heat_production/1e3),2)

                time_dict = {"Time (year)": teaobj_sCO2.Linear_time_distribution,
                             "sCO2 Heat Production (MWt)": teaobj_sCO2.Instantaneous_heat_production/1e3,
                             "sCO2 Annual Heat Production (GWh)": teaobj_sCO2.Annual_heat_production/1e6}
                econ_values_dict.update(time_dict)

            except ValueError as e:
                fig, lcoh_sCO2 = update_blank_econ2(fig=fig, nrow1=1, ncol1=1, nrow2=1, ncol2=3)
                error_message = parse_error_message(e=e, e_name='Err Econ1a')
                error_messages_dict['Err Econ1a'] = error_message
            
            except AttributeError as e:
                fig, lcoh_sCO2 = update_blank_econ2(fig=fig, nrow1=1, ncol1=1, nrow2=1, ncol2=3)
                error_message = parse_error_message(e=e, e_name='Err Econ1b')
                error_messages_dict['Err Econ1b'] = error_message

        if fluid == "H2O" or fluid == "All":

            # print(" ...... H20 HEATING LCOH .... ")

            try:
                # TODO: update D ... based on radial
                teaobj_H2O = create_teaobject(TandP_dict,
                                                u_sCO2, u_H2O, c_sCO2, c_H2O,
                                                case, end_use, fluid, sbt_version,
                                                mdot, L2, L1, grad, D, Tinj, k,
                                                Drilling_cost_per_m, Discount_rate, Lifetime, 
                                                Direct_use_heat_cost_per_kWth, Power_plant_cost_per_kWe, Pre_Cooling_Delta_T, Turbine_outlet_pressure,
                                                properties_H2O_pathname, 
                                                properties_CO2v2_pathname, 
                                                additional_properties_CO2v2_pathname,
                                                is_H20=True, is_heating=True
                                                )
                lcoh_H2O = format(teaobj_H2O.LCOH, '.2f')
                # print(lcoh_H2O)
                # print("Error on LCOH ... ")
                # print(teaobj_H2O)
                # print(lcoh_H2O)

                # HERE !!!!! "'TEA' object has no attribute 'LCOH'"
                # print('here')
                # print(teaobj_H2O.Linear_time_distribution)
 
                # Heat Production 
                fig.add_trace(go.Scatter(x=teaobj_H2O.Linear_time_distribution, y=teaobj_H2O.Instantaneous_heat_production/1e3,
                              hovertemplate='<b>Time (year)</b>: %{x:.1f}<br><b>Heat Production (MWt)</b>: %{y:.3f} ',
                              line = dict(color='black', width=lw, dash='dash'),
                              legendgroup=labels_cat[0], name=labels[0], showlegend=is_display_legend
                              ),
                              row=1, col=1)

                # Anuual Heat Production
                fig.add_trace(go.Bar(x=np.arange(1,teaobj_H2O.Lifetime+1), y=teaobj_H2O.Annual_heat_production/1e6,
                                    # hovertemplate='<b>Time (year)</b>: %{x:.1f}<br><b>Annual Heat Production (GWh)</b>: %{y:.3f} ',
                                    name=labels[0],
                                    showlegend=False, hovertemplate='<b>Annual Heat Production (GWh)</b>: %{y:.3f}'), 
                                row=1, col=3)


                fig.add_trace(go.Scatter(x=np.arange(1,teaobj_H2O.Lifetime+1), y=teaobj_H2O.Annual_heat_production/1e6,
                              hovertemplate='<b>Time (year)</b>: %{x:.1f}<br>',
                              line = dict(color='black', width=lw, dash='dash'),
                              legendgroup=labels_cat[0], name=labels[0], showlegend=False
                              ),
                              row=1, col=3)

                # Table Data
                mean_H2O_Net_HProd = round(np.mean(teaobj_H2O.Instantaneous_heat_production/1e3),2)

                time_dict = {"Time (year)": teaobj_H2O.Linear_time_distribution,
                             "H2O Heat Production (MWt)": teaobj_H2O.Instantaneous_heat_production/1e3,
                             "H2O Annual Heat Production (GWh)": teaobj_H2O.Annual_heat_production/1e6}
                econ_values_dict.update(time_dict)
            
            except ValueError as e:
                fig, lcoh_H2O = update_blank_econ2(fig=fig, nrow1=1, ncol1=1, nrow2=1, ncol2=3)
                error_message = parse_error_message(e=e, e_name='Err Econ2a')
                error_messages_dict['Err Econ2a'] = error_message
            
            except AttributeError as e:
                fig, lcoh_H2O = update_blank_econ2(fig=fig, nrow1=1, ncol1=1, nrow2=1, ncol2=3)
                error_message = parse_error_message(e=e, e_name='Err Econ2b')
                error_messages_dict['Err Econ2b'] = error_message

    if end_use == "Electricity" or end_use == "All":

        if end_use == "Electricity":
            is_display_legend = True
            row_num = 1
        else:
            is_display_legend = False
            row_num = 2

        if fluid == "sCO2" or fluid == "All":

            try:
                teaobj_sCO2 = create_teaobject(TandP_dict, 
                                                u_sCO2, u_H2O, c_sCO2, c_H2O,
                                                case, end_use, fluid, sbt_version,
                                                mdot, L2, L1, grad, D, Tinj, k,
                                                Drilling_cost_per_m, Discount_rate, Lifetime, 
                                                Direct_use_heat_cost_per_kWth, Power_plant_cost_per_kWe, Pre_Cooling_Delta_T, Turbine_outlet_pressure,
                                                properties_H2O_pathname, 
                                                properties_CO2v2_pathname, 
                                                additional_properties_CO2v2_pathname)
                
                # convert any negative value to 0
                teaobj_sCO2.Inst_Net_Electricity_production[teaobj_sCO2.Inst_Net_Electricity_production<0] = 0

                lcoe_sCO2 = format( teaobj_sCO2.LCOE, '.2f')

                # Electricity 
                fig.add_trace(go.Scatter(x=teaobj_sCO2.Linear_time_distribution, y=teaobj_sCO2.Inst_Net_Electricity_production/1e3,
                              hovertemplate='<b>Time (year)</b>: %{x:.1f}<br><b>Electricity Production (MWe)</b>: %{y:.3f} ',
                              line = dict(color='black', width=lw),
                              legendgroup=labels_cat[1], name=labels[1], showlegend=is_display_legend
                              ),
                              row=row_num, col=1)

                
                # Annual Electricity
                fig.add_trace(go.Bar(x=np.arange(1,teaobj_sCO2.Lifetime+1), y=teaobj_sCO2.Annual_electricity_production/1e6,
                                # hovertemplate='<b>Time (year)</b>: %{x:.1f}<br><b>Annual Electricity Production (GWe)</b>: %{y:.3f} ',
                                name=labels[1],
                                showlegend=False, hovertemplate='<b>Annual Electricity Production (GWe)</b>: %{y:.3f}'), 
                                row=row_num, col=3)

                fig.add_trace(go.Scatter(x=np.arange(1,teaobj_sCO2.Lifetime+1), y=teaobj_sCO2.Annual_electricity_production/1e6,
                              hovertemplate='<b>Time (year)</b>: %{x:.1f}<br>',
                              line = dict(color='black', width=lw),
                              legendgroup=labels_cat[1], name=labels[1], showlegend=False
                              ),
                              row=row_num, col=3)

                # Table Data
                time_dict = {"Time (year)": teaobj_sCO2.Linear_time_distribution,
                             "sCO2 Electricity Production (MWe)": teaobj_sCO2.Inst_Net_Electricity_production/1e3,
                             "sCO2 Annual Electricity Production (GWe)": teaobj_sCO2.Annual_electricity_production/1e6}
                econ_values_dict.update(time_dict)

                if is_plot_ts_check:
                    # if fluid == 2 and np.in1d(3000,teaobj_sCO2.error_codes) == False and np.in1d(4000,teaobj_sCO2.error_codes) == False:
                    if end_use == "Electricity":
                        # ts_fig = get_Ts_diagram(fig=ts_fig, teaobj=teaobj_sCO2, nrow=1, ncol=1)
                        get_Ts_diagram(fig=fig, teaobj=teaobj_sCO2, nrow=2, ncol=1, tmatrix_pathname=tmatrix_pathname)
                    if end_use == "All":
                        get_Ts_diagram(fig=fig, teaobj=teaobj_sCO2, nrow=3, ncol=1, tmatrix_pathname=tmatrix_pathname)

                    # ths throws an error sometimes ... to see why and where
                    # else:
                    #     ts_fig = blank_canvas(fig=ts_fig, row_n=1, col_n=1)

                # mean_sCO2_Net_HProd = round(np.mean(teaobj_sCO2.Instantaneous_heat_production/1e3),2)
                mean_sCO2_Net_EProd = round(np.mean(teaobj_sCO2.Inst_Net_Electricity_production/1e3),2)

            except ValueError as e:
                fig, lcoe_sCO2 = update_blank_econ2(fig=fig, nrow1=row_num, ncol1=1, nrow2=row_num, ncol2=3)
                error_message = parse_error_message(e=e, e_name='Err Econ3a')
                error_messages_dict['Err Econ3a'] = error_message

            except AttributeError as e:
                fig, lcoe_sCO2 = update_blank_econ2(fig=fig, nrow1=row_num, ncol1=1, nrow2=row_num, ncol2=3)
                error_message = parse_error_message(e=e, e_name='Err Econ3b')
                error_messages_dict['Err Econ3b'] = error_message

        if fluid == "H2O" or fluid == "All":

            try:
                teaobj_H2O = create_teaobject(TandP_dict,
                                                u_sCO2, u_H2O, c_sCO2, c_H2O,
                                                case, end_use, fluid, sbt_version,
                                                mdot, L2, L1, grad, D, Tinj, k,
                                                Drilling_cost_per_m, Discount_rate, Lifetime, 
                                                Direct_use_heat_cost_per_kWth, Power_plant_cost_per_kWe, Pre_Cooling_Delta_T, Turbine_outlet_pressure,
                                                properties_H2O_pathname, 
                                                properties_CO2v2_pathname, 
                                                additional_properties_CO2v2_pathname,
                                                is_H20=True
                                                )
                
                # convert any negative value to 0
                teaobj_H2O.Inst_Net_Electricity_production[teaobj_H2O.Inst_Net_Electricity_production<0] = 0

                lcoe_H2O = format(teaobj_H2O.LCOE, '.2f')

                # print(" ********* ")
                # print(teaobj_H2O.Inst_Net_Electricity_production/1e3)
                # print(teaobj_H2O.Linear_time_distribution)
                # print("\n")

                # Electricity 
                fig.add_trace(go.Scatter(x=teaobj_H2O.Linear_time_distribution, y=teaobj_H2O.Inst_Net_Electricity_production/1e3,
                              hovertemplate='<b>Time (year)</b>: %{x:.1f}<br><b>Electricity Production (MWe)</b>: %{y:.3f} ',
                              line = dict(color='black', width=lw, dash='dash'),
                              legendgroup=labels_cat[0], name=labels[0], showlegend=is_display_legend
                              ),
                              row=row_num, col=1)
                
                # Annual Electricity
                fig.add_trace(go.Bar(x=np.arange(1,teaobj_H2O.Lifetime+1), y=teaobj_H2O.Annual_electricity_production/1e6,
                                # hovertemplate='<b>Time (year)</b>: %{x:.1f}<br><b>Annual Electricity Production (GWe)</b>: %{y:.3f} ',
                                name=labels[0],
                                showlegend=False, hovertemplate='<b>Annual Electricity Production (GWe)</b>: %{y:.3f}'), 
                                row=row_num, col=3)

                fig.add_trace(go.Scatter(x=np.arange(1,teaobj_H2O.Lifetime+1), y=teaobj_H2O.Annual_electricity_production/1e6,
                              hovertemplate='<b>Time (year)</b>: %{x:.1f}<br>',
                              line = dict(color='black', width=lw, dash='dash'),
                              legendgroup=labels_cat[0], name=labels[0], showlegend=False
                              ),
                              row=row_num, col=3)

                # Table Data
                time_dict = {"Time (year)": teaobj_H2O.Linear_time_distribution,
                             "H2O Electricity Production (MWe)": teaobj_H2O.Inst_Net_Electricity_production/1e3,
                             "H2O Annual Electricity Production (GWe)": teaobj_H2O.Annual_electricity_production/1e6}
                econ_values_dict.update(time_dict)


                # mean_H2O_Net_HProd = round(np.mean(teaobj_H2O.Instantaneous_heat_production/1e3),2)
                mean_H2O_Net_EProd = round(np.mean(teaobj_H2O.Inst_Net_Electricity_production/1e3),2)

            except ValueError as e:
                fig, lcoe_H2O = update_blank_econ2(fig=fig, nrow1=row_num, ncol1=1, nrow2=row_num, ncol2=3)
                error_message = parse_error_message(e=e, e_name='Err Econ4a')
                error_messages_dict['Err Econ4a'] = error_message

            except AttributeError as e:
                fig, lcoe_H2O = update_blank_econ2(fig=fig, nrow1=row_num, ncol1=1, nrow2=row_num, ncol2=3)
                error_message = parse_error_message(e=e, e_name='Err Econ4b')
                error_messages_dict['Err Econ4b'] = error_message

    fig = update_layout_properties_econ_results(fig, end_use, scale)
    fig = update_lcoh_lcoe_table(fig, fluid, end_use, lcoh_sCO2, lcoh_H2O, lcoe_sCO2, lcoe_H2O) # table
    fig.update_traces(cells_font=dict(size = 13), row=1, col=5)
    fig.update_traces(cells_font=dict(size = 13), row=2, col=5)

    econ_data_dict = {'LCOH sCO2': lcoh_sCO2, 
                        'LCOH H2O': lcoh_H2O, 
                        'LCOE sCO2': lcoe_sCO2,
                        'LCOE H2O': lcoe_H2O,
                        'Mean H2O Net HProd': mean_H2O_Net_HProd,
                        'Mean H2O Net EProd': mean_H2O_Net_EProd,
                        'Mean sCO2 Net HProd': mean_sCO2_Net_HProd,
                        'Mean sCO2 Net EProd': mean_sCO2_Net_EProd}


    fig.update_layout(paper_bgcolor='rgba(255,255,255,0.10)', # or 0.40
                      plot_bgcolor='rgba(255,255,255,0)')

    # print(error_messages_dict)

    return fig, econ_data_dict, econ_values_dict, error_messages_dict



def CO2_isobaric_lines(tmatrix_pathname):

    # CO2 isobaric lines data for Ts diagram
    pvectorforTS = np.array([40,50,60,70,80,90,100,120,150,190,240,300]) # bar
    svectorforTS = np.arange(1000,2000,20) # J/kg/K                              
    tmatrix = genfromtxt(tmatrix_pathname, delimiter=',')

    # dome data for Ts diagram
    sdome = np.array([1001.37898699, 1046.06882698, 1088.21416035, 1128.86481841, 1169.09596369, 1210.24526707, 
                      1254.66210399, 1309.32653659, 1323.58218729, 1340.74189756, 1364.61655928, 1400.86587985, 
                      1468.75018408, 1516.27446567, 1546.66849166, 1567.63482054, 1584.41961907, 1643.56442582, 
                      1686.22125974, 1722.06785103, 1754.37017763, 1784.81415615, 1814.51645442, 1844.39848637])
    Tdome = np.array([273.31081538, 278.44972407, 283.13044128, 287.43392381, 291.41872472, 295.12790099, 
                      298.59249754, 301.8325153,  302.45422394, 303.06687319, 303.6699029,  304.08533585, 
                      304.08533585, 303.6699029,  303.06687319, 302.45422394, 301.8325153,  298.59249754, 
                      295.12790099, 291.41872472, 287.43392381, 283.13044128, 278.44972407, 273.31081538])

    return pvectorforTS, svectorforTS, tmatrix, sdome, Tdome

def get_Ts_diagram(fig, teaobj, nrow, ncol, tmatrix_pathname):

    error_dict = {}

    # in case of CO2 power production
    mint = "#bad9c2"
    orange = "#ff9500"
    darkred = "#940704"
    darkblue = "#04007a"
    pvectorforTS, svectorforTS, tmatrix, sdome, Tdome = CO2_isobaric_lines(tmatrix_pathname)

    subplot_len = len(fig.data[0].__dict__['_parent'].__dict__['_data'])

    # plot static parts
    for i in reversed(range(len(pvectorforTS))):
        fig.add_trace(go.Scatter(x=svectorforTS, y=tmatrix[i,:]-to_kelvin_factor,
                          hovertemplate='<b>Constant Entropy Line</b>: %{y:.1f}', # Isentropic
                          line = dict(color=orange, width=1), 
                          showlegend=False,
                          name="",
                          ),
                          row=nrow, col=ncol)

    template_i = 12 # last line

    for i in range(len(pvectorforTS)):

        template_i -= 1

        fig.add_trace(go.Scatter(
                            x=[svectorforTS[-1]],
                            y=[tmatrix[i,:][-1]-to_kelvin_factor],
                            hovertemplate=None,
                            hoverinfo='skip',
                            showlegend=False,
                            line = dict(color=orange, width=1),
                            mode="lines+markers+text",
                            name="point",
                            text=[f'P{i}'],
                            textposition="middle right"
                            ), row=nrow, col=ncol)
        
        # adaptable to new situations?
        hover_replace = fig.data[0].__dict__['_parent'].__dict__['_data'][subplot_len + template_i]['hovertemplate'].replace("Constant Entropy Line", f"Constant Entropy Line, P{i}")
        fig.data[0].__dict__['_parent'].__dict__['_data'][subplot_len + template_i]['hovertemplate'] = hover_replace

    fig.add_trace(go.Scatter(
                        x=[1100],
                        y=[200],
                        hovertemplate=None,
                        hoverinfo='skip',
                        showlegend=False,
                        line = dict(color="white", width=1),
                        mode="lines+markers+text",
                        name="point",
                        text=['<i>sCO2 only</i>'],
                        textposition="middle right"
                        ), row=nrow, col=ncol)

    fig.add_trace(go.Scatter(x=sdome, y=Tdome-to_kelvin_factor,
                              hovertemplate='<b>Saturation Curve</b>: %{y:.1f}',
                              line = dict(color="black", width=1.5), showlegend=False,
                              name="",
                              ),
                              row=nrow, col=ncol)

    
    # Temperature-entropy (Ts) at turbine inlet
    turbin_inletT = teaobj.Linear_production_temperature[-1]
    turbin_inletS = teaobj.s_prod[-1]

    # Ts at turbine outlet
    turbine_outletT = teaobj.T_turbine_out_actual[-1]
    turbine_outletS = interpn((teaobj.Pvector, teaobj.Tvector), teaobj.entropy, 
                                np.array([teaobj.Turbine_outlet_pressure*1e5, turbine_outletT+to_kelvin_factor]))[0]

    
    try:
        # Ts at compressor inlet
        compressor_inletT = teaobj.Pre_cooling_temperature
        compressor_inletS = interpn((teaobj.Pvector, teaobj.Tvector), teaobj.entropy, 
                                    np.array([teaobj.Turbine_outlet_pressure*1e5, compressor_inletT+to_kelvin_factor]))[0]
    except Exception as e:
        error_dict['TS-error1'] = e

   
    try:
         # Ts at compressor outlet
        compressor_outletT = teaobj.Post_compressor_T_actual[0]
        compressor_outletS = interpn((teaobj.Pvector, teaobj.Tvector), teaobj.entropy, 
                                    np.array([teaobj.P_in, compressor_outletT+to_kelvin_factor]))[0]
    except Exception as e:
        error_dict['TS-error2'] = e

    try:
        # Ts at re-injection
        reinjectionT = teaobj.T_in
        reinjectionS = interpn((teaobj.Pvector, teaobj.Tvector), teaobj.entropy, 
                                    np.array([teaobj.P_in, reinjectionT+to_kelvin_factor]))[0]

        cycleT = np.array([turbin_inletT, turbine_outletT, compressor_inletT, compressor_outletT, reinjectionT, turbin_inletT])
        cycleS = np.array([turbin_inletS, turbine_outletS, compressor_inletS, compressor_outletS, reinjectionS, turbin_inletS])

        fig.add_trace(go.Scatter(x=cycleS, y=cycleT,
                                  hovertemplate='<b>Direct Turbine Expansion Cycle</b>: %{y:.1f}<br>',
                                  line = dict(color="grey", width=1.5), showlegend=False,
                                  name="",
                                  ),
                                  row=nrow, col=ncol)
        # turbine injection (state 1) to turbine expansion ('trubine outlet state', state 2)
        fig.add_trace(go.Scatter(x=cycleS[:2], y=cycleT[:2],
                                  hovertemplate=None,
                                  hoverinfo='skip',
                                  line = dict(color=darkred, width=1.5), showlegend=False,
                                  name="",
                                  ),
                                  row=nrow, col=ncol)
        # pre-cooling
        fig.add_trace(go.Scatter(x=cycleS[1:3], y=cycleT[1:3],
                                  hovertemplate=None,
                                  hoverinfo='skip',
                                  line = dict(color=darkblue, width=1.5), showlegend=False,
                                  name="",
                                  ),
                                  row=nrow, col=ncol)
        # compression, depends on turbine outlet pressure
        fig.add_trace(go.Scatter(x=cycleS[2:4], y=cycleT[2:4],
                                  hovertemplate=None,
                                  hoverinfo='skip',
                                  line = dict(color="#5192f5", width=1.5), showlegend=False,
                                  name="",
                                  ),
                                  row=nrow, col=ncol)
        # secondary cooling
        fig.add_trace(go.Scatter(x=cycleS[3:], y=cycleT[3:],
                                  hovertemplate=None,
                                  hoverinfo='skip',
                                  line = dict(color=darkred, width=1.5), showlegend=False,
                                  name="",
                                  ),
                                  row=nrow, col=ncol)
    except Exception as e:
        error_dict['TS-error3'] = e


    fig.update_layout(title_text=f'<b>CO2 Temperature-entropy (T-s) Diagram</b>', 
                        title_x=0.35, title_y=0.99,
                        font=dict(size=10)
                        )
    fig.update_yaxes(title_text="Temperature (°C)", 
                            row=nrow, col=ncol,
                            tickfont = dict(size=12), title_font=dict(size=14))
    fig.update_xaxes(title_text="Entropy (J/kg/°C)", range=[1000, 2300],
                            row=nrow, col=ncol,
                            tickfont = dict(size=12), title_font=dict(size=14))
    return fig


