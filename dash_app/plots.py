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
# HELPER FUNCTIONS FOR FLUID EXPANSION
# --------------------------------------------------------------------------------------------------------------------

def fluids_to_run(model, case, fluid):
    """
    Expand "All" fluid selection to concrete fluids for simulator.
    For simulator+fluid="All", returns ["H2O", "sCO2"].
    Otherwise returns [fluid].
    """
    is_sim = model in ("SBT V1.0", "SBT V2.0")
    if is_sim and fluid == "All":
        return ["H2O", "sCO2"]
    return [fluid]


def model_for_fluid(model, case, fluid):
    """
    Get the correct model version for a specific fluid.
    For simulator+utube: H2O uses SBT V1.0, sCO2 uses SBT V2.0.
    Otherwise returns the original model.
    """
    if model in ("SBT V1.0", "SBT V2.0") and case == "utube":
        if fluid == "H2O":
            return "SBT V1.0"
        if fluid == "sCO2":
            return "SBT V2.0"
    return model

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
            mesh, accuracy, 
            # mass_mode, temp_mode
            HyperParam1, HyperParam3, HyperParam5
            ):

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
    elif model == "SBT V2.0":
        sbt_version = 2
    else:
        sbt_version = 0
    
    sbt_version_sco2 = sbt_version
    sbt_version_h2o = sbt_version
    
    # For simulator (SBT models) only, not database (HDF5)
    # Override sbt_version for each fluid based on fluid selection
    if model in ("SBT V1.0", "SBT V2.0"):
        if fluid == "H2O":
            sbt_version_h2o = 1
            # print(f"[DEBUG generate_subsurface_lineplots] model={model}, case={case}, fluid={fluid} -> Running H2O only with SBT V{sbt_version_h2o}.0", flush=True)
        elif fluid == "sCO2":
            sbt_version_sco2 = 2
            # print(f"[DEBUG generate_subsurface_lineplots] model={model}, case={case}, fluid={fluid} -> Running sCO2 only with SBT V{sbt_version_sco2}.0", flush=True)
        elif fluid == "All":
            sbt_version_h2o = 1
            sbt_version_sco2 = 2
            # print(f"[DEBUG generate_subsurface_lineplots] model={model}, case={case}, fluid={fluid} -> Running H2O: SBT V{sbt_version_h2o}.0, sCO2: SBT V{sbt_version_sco2}.0", flush=True)
    
    # sCO2 is supported if using SBT V2.0 or HDF5 database
    sCO2_supported = (sbt_version_sco2 == 2) or (model == "HDF5")

    mean_H2O_Tout = "-"
    mean_H2O_Pout = "-"
    mean_sCO2_Tout = "-"
    mean_sCO2_Pout = "-"

    error_messages_dict = {}

    # Initialize time and m_dot based on case - use appropriate data source
    if case == "coaxial":
        time = c_sCO2.time if hasattr(c_sCO2, 'time') else u_sCO2.time
        m_dot = c_sCO2.mdot if hasattr(c_sCO2, 'mdot') else u_sCO2.mdot
    else:  # utube
        time = u_sCO2.time
        m_dot = u_sCO2.mdot

    mass_flow_rates_dict = {"Mass Flow Rate (kg/s)": m_dot}
    time_dict = {"Time (year)": time}

    arg_mdot_v, arg_mdot_i, arg_L2_v, arg_L2_i, arg_L1_v, arg_L1_i, arg_grad_v, arg_grad_i, arg_D_v, arg_D_i, \
                arg_Tinj_v, arg_Tinj_i, arg_k_v, arg_k_i = param_nearest_init(arg_mdot, arg_L2, arg_L1, arg_grad, arg_D, arg_Tinj, arg_k)

    point = (arg_mdot, arg_L2, arg_L1, arg_grad, arg_D, arg_Tinj + to_kelvin_factor, arg_k) # to kelvin

    # ** Average calculations
    # For SBT models, HDF5 array access may fail, so wrap in try-except
    try:
        sCO2_kWe_avg, sCO2_kWt_avg, H2O_kWe_avg, H2O_kWt_avg, error_messages_d = \
                                get_kWe_kWt_over_mass_or_time(case, fluid, point, arg_L2_i, arg_L1_i, arg_grad_i, arg_D_i, arg_Tinj_i, arg_k_i)
        error_messages_dict.update(error_messages_d)
    except (IndexError, ValueError, KeyError) as e:
        # For SBT models, HDF5 arrays may not have data for these parameter ranges
        # Set averages to None and continue - SBT calculations will still work
        sCO2_kWe_avg = None
        sCO2_kWt_avg = None
        H2O_kWe_avg = None
        H2O_kWt_avg = None
        if model in ["SBT V1.0", "SBT V2.0"]:
            # Don't add error for SBT models - this is expected
            pass
        else:
            # For HDF5 models, this is a real error
            error_message = parse_error_message(e=e, e_name='Err SubResAvg', model=model)
            error_messages_dict['Err SubResAvg'] = error_message

    # print(sCO2_kWe_avg, sCO2_kWt_avg, H2O_kWe_avg, H2O_kWt_avg)

    is_blank_data = False
    
    # Initialize success flags for all cases (HDF5 always succeeds, SBT may fail)
    sCO2_success = True
    H2O_success = True

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

    if interp_time == "True":

        # this doesn't impact the kWe and kWt averages

        if case == "utube":

            # Track success/failure for each fluid separately
            sCO2_success = False
            H2O_success = False
            
            # Handle each fluid separately so one failure doesn't prevent the other from running
            if (fluid == "sCO2" or fluid == "All") and sCO2_supported:
                try:
                    sCO2_Tout, sCO2_Pout, time = u_sCO2.interp_outlet_states(point, sbt_version_sco2,
                                                        Tsurf, c_m, rho_m, 
                                                        # radius_vertical, radius_lateral, n_laterals, lateral_flow, lateral_multiplier,
                                                        Diameter1, Diameter2, PipeParam3, PipeParam4, PipeParam5,
                                                        mesh, accuracy, HyperParam1, HyperParam3, HyperParam5
                                                        )
                    # For SBT V2.0 sCO2, use HyperParam1 (in MPa) converted to Pa as inlet pressure
                    # For SBT V1.0, HyperParam1 is Mass Flow Rate Mode (string like 'Constant'), not inlet pressure
                    if sbt_version_sco2 == 2 and HyperParam1 is not None:
                        try:
                            Pinj_sco2 = float(HyperParam1) * 1e6
                        except (ValueError, TypeError):
                            # HyperParam1 is not a numeric value (e.g., 'Constant' for SBT V1.0)
                            Pinj_sco2 = 1e7
                    else:
                        # For SBT V1.0 or if HyperParam1 is None, use default 100 bar (1e7 Pa)
                        Pinj_sco2 = 1e7 if sbt_version_sco2 > 0 else None
                    sCO2_kWe, sCO2_kWt = u_sCO2.interp_kW(point, sCO2_Tout, sCO2_Pout, Pinj=Pinj_sco2)
                    sCO2_success = True
                except Exception as e:
                    sCO2_Tout, sCO2_Pout, _, _, sCO2_kWe, sCO2_kWt, _, _ = blank_data()
                    error_message = parse_error_message(e=e, e_name='Err SubRes3', model=model)
                    error_messages_dict['Err SubRes3'] = error_message
                    print(f"[ERROR] Exception in generate_subsurface_lineplots (utube, SBT{sbt_version}, sCO2): {e}")
                    print(f"  Parameters: case={case}, fluid={fluid}, PipeParam5={PipeParam5}, Diameter1={Diameter1}, Diameter2={Diameter2}")
                    import traceback
                    traceback.print_exc()
            elif (fluid == "sCO2" or fluid == "All") and not sCO2_supported:
                # SBT V1.0 doesn't support sCO2 - skip calculation and set to blank
                sCO2_Tout, sCO2_Pout, _, _, sCO2_kWe, sCO2_kWt, _, _ = blank_data()
                sCO2_success = False
                print(f"[INFO] Skipping sCO2 calculation: SBT V1.0 doesn't support sCO2. Only H2O calculations will be performed.", flush=True)
            
            if fluid == "H2O" or fluid == "All":
                try:
                    H2O_Tout, H2O_Pout, time = u_H2O.interp_outlet_states(point, sbt_version_h2o,
                                                        Tsurf, c_m, rho_m, 
                                                        # radius_vertical, radius_lateral, n_laterals, lateral_flow, lateral_multiplier,
                                                        Diameter1, Diameter2, PipeParam3, PipeParam4, PipeParam5,
                                                        mesh, accuracy, HyperParam1, HyperParam3, HyperParam5)
                    # For SBT V2.0 H2O, use HyperParam1 (in MPa) converted to Pa as inlet pressure
                    # For SBT V1.0, HyperParam1 is Mass Flow Rate Mode (string like 'Constant'), not inlet pressure
                    if sbt_version_h2o == 2 and HyperParam1 is not None:
                        try:
                            Pinj_h2o = float(HyperParam1) * 1e6
                        except (ValueError, TypeError):
                            # HyperParam1 is not a numeric value (e.g., 'Constant' for SBT V1.0)
                            Pinj_h2o = None
                    else:
                        Pinj_h2o = None
                    H2O_kWe, H2O_kWt = u_H2O.interp_kW(point, H2O_Tout, H2O_Pout, Pinj=Pinj_h2o)
                    H2O_success = True
                except Exception as e:
                    _, _, H2O_Tout, H2O_Pout, _, _, H2O_kWe, H2O_kWt = blank_data()
                    error_message = parse_error_message(e=e, e_name='Err SubRes3', model=model)
                    error_messages_dict['Err SubRes3'] = error_message
                    print(f"[ERROR] Exception in generate_subsurface_lineplots (utube, SBT{sbt_version}, H2O): {e}")
                    print(f"  Parameters: case={case}, fluid={fluid}, PipeParam5={PipeParam5}, Diameter1={Diameter1}, Diameter2={Diameter2}")
                    import traceback
                    traceback.print_exc()
            
            # Set is_blank_data only if all selected fluids failed
            if (fluid == "sCO2" and not sCO2_success) or \
               (fluid == "H2O" and not H2O_success) or \
               (fluid == "All" and not sCO2_success and not H2O_success):
                is_blank_data = True
                # Only set time to empty if no fluid succeeded
                if not sCO2_success and not H2O_success:
                    time = pd.Series(dtype=object)
                
        if case == "coaxial":

            # Track success/failure for each fluid separately
            sCO2_success = False
            H2O_success = False
            
            # Initialize variables to ensure they exist even if both fluids fail
            sCO2_Tout, sCO2_Pout, H2O_Tout, H2O_Pout, sCO2_kWe, sCO2_kWt, H2O_kWe, H2O_kWt = blank_data()
            
            # Initialize time variables for each fluid separately
            sCO2_time = None
            H2O_time = None
            
            # Handle each fluid separately so one failure doesn't prevent the other from running
            # Store geometry values for comparison
            geometry_before_sCO2 = {'Diameter1': Diameter1, 'Diameter2': Diameter2, 'PipeParam3': PipeParam3, 'PipeParam5': PipeParam5}
            
            if (fluid == "sCO2" or fluid == "All") and sCO2_supported:
                try:
                    sCO2_Tout, sCO2_Pout, sCO2_time = c_sCO2.interp_outlet_states(point, sbt_version_sco2,
                                                            Tsurf, c_m, rho_m, 
                                                            # radius_vertical, radius_lateral, n_laterals, lateral_flow, lateral_multiplier,
                                                            Diameter1, Diameter2, PipeParam3, PipeParam4, PipeParam5,
                                                            mesh, accuracy, HyperParam1, HyperParam3, HyperParam5)
                    
                    # For SBT V2.0 sCO2, use HyperParam1 (in MPa) converted to Pa as inlet pressure
                    # For SBT V1.0, HyperParam1 is Mass Flow Rate Mode (string like 'Constant'), not inlet pressure
                    if sbt_version_sco2 == 2 and HyperParam1 is not None:
                        try:
                            Pinj_sco2 = float(HyperParam1) * 1e6
                        except (ValueError, TypeError):
                            # HyperParam1 is not a numeric value (e.g., 'Constant' for SBT V1.0)
                            Pinj_sco2 = 1e7
                    else:
                        # For SBT V1.0 or if HyperParam1 is None, use default 100 bar (1e7 Pa)
                        Pinj_sco2 = 1e7 if sbt_version_sco2 > 0 else None
                    sCO2_kWe, sCO2_kWt = c_sCO2.interp_kW(point, sCO2_Tout, sCO2_Pout, Pinj=Pinj_sco2)
                    sCO2_success = True
                    # Store sCO2 time for later use
                    if sCO2_time is not None:
                        time = sCO2_time
                except Exception as e:
                    sCO2_Tout, sCO2_Pout, _, _, sCO2_kWe, sCO2_kWt, _, _ = blank_data()
                    sCO2_success = False
                    error_message = parse_error_message(e=e, e_name='Err SubRes4', model=model)
                    error_messages_dict['Err SubRes4'] = error_message
                    print(f"[ERROR] Exception in generate_subsurface_lineplots (coaxial, SBT{sbt_version}, sCO2): {type(e).__name__}: {e}", flush=True)
                    print(f"  Parameters: case={case}, fluid={fluid}, PipeParam5={PipeParam5}, Diameter1={Diameter1}, Diameter2={Diameter2}, mdot={arg_mdot}", flush=True)
                    import traceback
                    traceback.print_exc()
            elif (fluid == "sCO2" or fluid == "All") and not sCO2_supported:
                # SBT V1.0 doesn't support sCO2 - skip calculation and set to blank
                sCO2_Tout, sCO2_Pout, _, _, sCO2_kWe, sCO2_kWt, _, _ = blank_data()
                sCO2_success = False
                print(f"[INFO] Skipping sCO2 calculation: SBT V1.0 doesn't support sCO2. Only H2O calculations will be performed.", flush=True)
            
            # Store geometry values for comparison
            geometry_before_H2O = {'Diameter1': Diameter1, 'Diameter2': Diameter2, 'PipeParam3': PipeParam3, 'PipeParam5': PipeParam5}
            
            if fluid == "H2O" or fluid == "All":
                # Compare geometry values if both fluids are being run
                if fluid == "All" and 'geometry_before_sCO2' in locals():
                    if geometry_before_sCO2['Diameter1'] != geometry_before_H2O['Diameter1'] or geometry_before_sCO2['Diameter2'] != geometry_before_H2O['Diameter2']:
                        print(f"[WARNING] Geometry values differ between sCO2 and H2O!", flush=True)
                try:
                    print(f"[DEBUG] Attempting H2O calculation for coaxial SBT (sbt_version={sbt_version_h2o}, Diameter1={Diameter1}, Diameter2={Diameter2})", flush=True)
                    H2O_Tout, H2O_Pout, H2O_time = c_H2O.interp_outlet_states(point, sbt_version_h2o,
                                                            Tsurf, c_m, rho_m, 
                                                            # radius_vertical, radius_lateral, n_laterals, lateral_flow, lateral_multiplier,
                                                            Diameter1, Diameter2, PipeParam3, PipeParam4, PipeParam5,
                                                            mesh, accuracy, HyperParam1, HyperParam3, HyperParam5)
                    # For SBT V2.0 H2O, use HyperParam1 (in MPa) converted to Pa as inlet pressure
                    # For SBT V1.0, HyperParam1 is Mass Flow Rate Mode (string like 'Constant'), not inlet pressure
                    if sbt_version_h2o == 2 and HyperParam1 is not None:
                        try:
                            Pinj_h2o = float(HyperParam1) * 1e6
                        except (ValueError, TypeError):
                            # HyperParam1 is not a numeric value (e.g., 'Constant' for SBT V1.0)
                            Pinj_h2o = None
                    else:
                        Pinj_h2o = None
                    H2O_kWe, H2O_kWt = c_H2O.interp_kW(point, H2O_Tout, H2O_Pout, Pinj=Pinj_h2o)
                    H2O_success = True
                    print(f"[DEBUG] H2O calculation succeeded for coaxial SBT", flush=True)
                    # Use H2O time if sCO2 didn't succeed, otherwise keep sCO2 time
                    if not sCO2_success and H2O_time is not None:
                        time = H2O_time
                except Exception as e:
                    _, _, H2O_Tout, H2O_Pout, _, _, H2O_kWe, H2O_kWt = blank_data()
                    error_message = parse_error_message(e=e, e_name='Err SubRes4', model=model)
                    error_messages_dict['Err SubRes4'] = error_message
                    print(f"[ERROR] Exception in generate_subsurface_lineplots (coaxial, SBT{sbt_version}, H2O): {e}", flush=True)
                    print(f"  Parameters: case={case}, fluid={fluid}, PipeParam5={PipeParam5}, Diameter1={Diameter1}, Diameter2={Diameter2}, mdot={arg_mdot}", flush=True)
                    import traceback
                    traceback.print_exc()
            
            # Set is_blank_data only if all selected fluids failed
            if (fluid == "sCO2" and not sCO2_success) or \
               (fluid == "H2O" and not H2O_success) or \
               (fluid == "All" and not sCO2_success and not H2O_success):
                is_blank_data = True
                # Only set time to empty if no fluid succeeded
                if not sCO2_success and not H2O_success:
                    time = pd.Series(dtype=object)
            else:
                # At least one fluid succeeded, ensure time is set
                if time is None or (hasattr(time, '__len__') and len(time) == 0):
                    # Fallback to HDF5 time if SBT simulation didn't return time
                    time = c_sCO2.time if hasattr(c_sCO2, 'time') else u_sCO2.time
            
                


    # Figure generation 
    is_simulator = model in ("SBT V1.0", "SBT V2.0")
    
    if is_simulator:
        fig = make_subplots(rows=2, cols=3,
                            specs=[[{}, {}, None],
                            [{}, {}, {}]],
                            horizontal_spacing = 0.11,
                            vertical_spacing = 0.21,
                            row_heights=[0, 0.5]
                            )
    else:
        fig = make_subplots(rows=2, cols=3,
                            specs=[[{}, {}, None],
                            [{}, {}, {}]],
                            horizontal_spacing = 0.11,
                            vertical_spacing = 0.21
                            )

    if (fluid == "sCO2" or fluid == "All") and sCO2_success:

        if not is_simulator:
            if sCO2_kWe_avg is not None and sCO2_kWt_avg is not None:
                try:
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
                except (ValueError, TypeError) as e:
                    print(f"[WARNING] Could not plot sCO2 averages: {e}", flush=True)

        # Plot time-series data only if calculation succeeded
        if not is_blank_data:
            # kWe vs. time
            try:
                fig.add_trace(go.Scatter(x=time, y=sCO2_kWe,
                              hovertemplate='<b>Time (year)</b>: %{x:.1f}<br><b>kWe</b>: %{y:.3f} ',
                              line = dict(color='black', width=lw), # #379dbf
                              legendgroup=labels_cat[1], name=labels[1], showlegend=True),
                              row=2, col=1)
            except (ValueError, TypeError):
                pass  # Skip if data is empty
            # temperature
            try:
                fig.add_trace(go.Scatter(x=time, y=sCO2_Tout - to_kelvin_factor, # to celsius
                              hovertemplate='<b>Time (year)</b>: %{x:.1f}<br><b>Outlet Temperature (˚C)</b>: %{y:.3f} ',
                              line = dict(color='royalblue', width=lw),
                              legendgroup=labels_cat[1], name=labels[1], showlegend=False),
                              row=2, col=2)
            except (ValueError, TypeError):
                pass  # Skip if data is empty
            # pressure
            try:
                fig.add_trace(go.Scatter(x=time, y=sCO2_Pout / 1000000,
                              hovertemplate='<b>Time (year)</b>: %{x:.1f}<br><b>Outlet Pressure (MPa)</b>: %{y:.3f} ',
                              line = dict(color='#16b8a2', width=lw),
                              legendgroup=labels_cat[1], name=labels[1], showlegend=False),
                              row=2, col=3)
            except (ValueError, TypeError):
                pass  # Skip if data is empty
        else:
            blank_canvas(fig=fig, row_n=2, col_n=1)
            blank_canvas(fig=fig, row_n=2, col_n=2)
            blank_canvas(fig=fig, row_n=2, col_n=3)
            if not is_simulator:
                if sCO2_kWe_avg is None or sCO2_kWt_avg is None:
                    blank_canvas(fig=fig, row_n=1, col_n=1)
                    blank_canvas(fig=fig, row_n=1, col_n=2)
        
        if is_blank_data or not sCO2_success:
            mean_sCO2_Tout = "-"
            mean_sCO2_Pout = "-"
        else:
            try:
                mean_sCO2_Tout = round(np.mean(sCO2_Tout - to_kelvin_factor),2) # to celsius
                mean_sCO2_Pout = round(np.mean(sCO2_Pout / 1000000),2)
            except (ValueError, TypeError):
                mean_sCO2_Tout = "-"
                mean_sCO2_Pout = "-"

        if str(mean_sCO2_Tout) == "nan":
            mean_sCO2_Tout = "-"
        if str(mean_sCO2_Pout) == "nan":
            mean_sCO2_Pout = "-"

        # Only add to dicts if data is not blank and sCO2 succeeded
        if not is_blank_data and sCO2_success:
            if sCO2_kWe_avg is not None:
                mass_flow_rates_dict["sCO2 40-Year Average of Exergy (kWe)"] = sCO2_kWe_avg
            if sCO2_kWt_avg is not None:
                mass_flow_rates_dict["sCO2 40-Year Average of Available Thermal Output (kWt)"] = sCO2_kWt_avg

            try:
                time_dict["sCO2 Exergy (kWe)"] = sCO2_kWe
                time_dict["sCO2 Outlet Temperature (˚C)"] = sCO2_Tout - to_kelvin_factor # to celsius
                time_dict["sCO2 Outlet Pressure (MPa)"] = sCO2_Pout / 1000000
            except (ValueError, TypeError):
                pass  # Skip if data is empty

    if (fluid == "H2O" or fluid == "All") and H2O_success:

        if not is_blank_data:
            if not is_simulator:
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
            try:
                fig.add_trace(go.Scatter(x=time, y=H2O_Tout - to_kelvin_factor, # to celsius
                          hovertemplate='<b>Time (year)</b>: %{x:.1f}<br><b>Outlet Temperature (˚C)</b>: %{y:.3f} ',
                          line = dict(color='royalblue', width=lw, dash='dash'),
                          legendgroup=labels_cat[0], name=labels[0], showlegend=False),
                          row=2, col=2)
            except (ValueError, TypeError):
                pass  # Skip if data is empty
            # pressure
            try:
                fig.add_trace(go.Scatter(x=time, y=H2O_Pout / 1000000,
                          hovertemplate='<b>Time (year)</b>: %{x:.1f}<br><b>Outlet Pressure (MPa)</b>: %{y:.3f} ',
                          line = dict(color='#16b8a2', width=lw, dash='dash'),
                          legendgroup=labels_cat[0], name=labels[0], showlegend=False),
                          row=2, col=3)
            except (ValueError, TypeError):
                pass  # Skip if data is empty

        if is_blank_data:
            if not is_simulator:
                blank_canvas(fig=fig, row_n=1, col_n=1)
                blank_canvas(fig=fig, row_n=1, col_n=2)
            blank_canvas(fig=fig, row_n=2, col_n=1)
            blank_canvas(fig=fig, row_n=2, col_n=2)
            blank_canvas(fig=fig, row_n=2, col_n=3)
            mean_H2O_Tout = "-"
            mean_H2O_Pout = "-"
        else:
            try:
                mean_H2O_Tout = round(np.mean(H2O_Tout - to_kelvin_factor),2) # to celsius
                mean_H2O_Pout = round(np.mean(H2O_Pout / 1000000),2)
            except (ValueError, TypeError):
                mean_H2O_Tout = "-"
                mean_H2O_Pout = "-"

        if str(mean_H2O_Tout) == "nan":
            mean_H2O_Tout = "-"
        if str(mean_H2O_Pout) == "nan":
            mean_H2O_Pout = "-"

        # Only add to dicts if data is not blank
        if not is_blank_data:
            if H2O_kWe_avg is not None:
                mass_flow_rates_dict["H2O 40-Year Average of Exergy (kWe)"] = H2O_kWe_avg
            if H2O_kWt_avg is not None:
                mass_flow_rates_dict["H2O 40-Year Average of Available Thermal Output (kWt)"] = H2O_kWt_avg

            try:
                time_dict["H2O Exergy (kWe)"] = H2O_kWe
                time_dict["H2O Outlet Temperature (˚C)"] = H2O_Tout - to_kelvin_factor # to celsius
                time_dict["H2O Outlet Pressure (MPa)"] = H2O_Pout / 1000000
            except (ValueError, TypeError):
                pass  # Skip if data is empty

    
    fig = update_layout_properties_subsurface_results(fig=fig, m_dot=m_dot, time=time, plot_scale=scale, is_simulator=is_simulator)

    forty_yr_means_dict = {'Mean H2O Tout': mean_H2O_Tout, 
                            'Mean H2O Pout': mean_H2O_Pout,
                            'Mean sCO2 Tout': mean_sCO2_Tout,
                            'Mean sCO2 Pout': mean_sCO2_Pout}
   
    # Convert numpy arrays to lists for JSON serialization (needed for dcc.Store)
    def convert_to_list(value):
        if value is not None and hasattr(value, 'tolist'):
            return value.tolist()
        return value
    
    # Convert time to list - handle pandas Series and numpy arrays
    def convert_time_to_list(time_value):
        if time_value is None:
            return []
        if hasattr(time_value, 'tolist'):
            return time_value.tolist()
        if hasattr(time_value, 'values'):  # pandas Series
            return time_value.values.tolist()
        if isinstance(time_value, (list, tuple)):
            return list(time_value)
        return []
    
    TandP_dict = {"sCO2_Tout": convert_to_list(sCO2_Tout),
                            "sCO2_Pout": convert_to_list(sCO2_Pout),
                            "H2O_Tout": convert_to_list(H2O_Tout),
                            "H2O_Pout": convert_to_list(H2O_Pout),
                            "sCO2_kWe": convert_to_list(sCO2_kWe),
                            "sCO2_kWt": convert_to_list(sCO2_kWt),
                            "H2O_kWe": convert_to_list(H2O_kWe),
                            "H2O_kWt": convert_to_list(H2O_kWt),
                            "time": convert_time_to_list(time),
                            "mdot": convert_to_list(m_dot),
                            }
    
    if is_blank_data:
        pass

    
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
        if param == "Geothermal Gradient (°C/m)":
            arg_grad_i = slice(None)
        if param == "Borehole Diameter (m)":
            arg_D_i = slice(None)
        if param == "Injection Temperature (˚C)":
            arg_Tinj_i = slice(None)
        if param == "Rock Thermal Conductivity (W/m-°C)":
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
        

    # Fill NaN values in contour data to prevent cut-off corners
    # Use a combination of forward/backward fill and nearest neighbor interpolation
    def fill_contour_nans(data):
        """Fill NaN values in 2D contour data using interpolation and fallback methods"""
        if not isinstance(data, np.ndarray):
            data = np.array(data)
        
        # Convert to regular numpy array (not masked array or zarr)
        data = np.asarray(data, dtype=np.float64)
        
        # Replace Inf values with NaN so they get filled too
        data[np.isinf(data)] = np.nan
        
        # Check if there are any NaN values
        if not np.any(np.isnan(data)):
            return data
        
        # Create a copy to avoid modifying original
        filled_data = data.copy()
        
        # Get valid (non-NaN) indices
        valid_mask = ~np.isnan(filled_data)
        
        if not np.any(valid_mask):
            # All NaN - fill with zeros
            filled_data[:] = 0
            return filled_data
        
        # Strategy 1: Use forward/backward fill along both axes for edge cases
        # Forward fill along rows (axis 0)
        for i in range(1, filled_data.shape[0]):
            nan_mask = np.isnan(filled_data[i, :])
            if np.any(nan_mask):
                filled_data[i, nan_mask] = filled_data[i-1, nan_mask]
        
        # Backward fill along rows
        for i in range(filled_data.shape[0]-2, -1, -1):
            nan_mask = np.isnan(filled_data[i, :])
            if np.any(nan_mask):
                filled_data[i, nan_mask] = filled_data[i+1, nan_mask]
        
        # Forward fill along columns (axis 1)
        for j in range(1, filled_data.shape[1]):
            nan_mask = np.isnan(filled_data[:, j])
            if np.any(nan_mask):
                filled_data[nan_mask, j] = filled_data[nan_mask, j-1]
        
        # Backward fill along columns
        for j in range(filled_data.shape[1]-2, -1, -1):
            nan_mask = np.isnan(filled_data[:, j])
            if np.any(nan_mask):
                filled_data[nan_mask, j] = filled_data[nan_mask, j+1]
        
        # Strategy 2: Use scipy interpolation for better handling of large NaN regions
        if np.any(np.isnan(filled_data)):
            try:
                from scipy.interpolate import griddata
                
                # Get coordinates of all points
                rows, cols = np.meshgrid(np.arange(filled_data.shape[0]), 
                                         np.arange(filled_data.shape[1]), 
                                         indexing='ij')
                
                # Get valid (non-NaN) points
                valid_mask = ~np.isnan(filled_data)
                valid_points = np.column_stack((rows[valid_mask], cols[valid_mask]))
                valid_values = filled_data[valid_mask]
                
                # Get invalid (NaN) points
                invalid_mask = np.isnan(filled_data)
                invalid_points = np.column_stack((rows[invalid_mask], cols[invalid_mask]))
                
                if len(invalid_points) > 0 and len(valid_points) > 0:
                    # Use griddata to interpolate NaN values
                    # 'nearest' method is fast and works well for filling gaps
                    interpolated = griddata(valid_points, valid_values, invalid_points, 
                                           method='nearest', fill_value=np.nanmin(valid_values))
                    filled_data[invalid_mask] = interpolated
            except ImportError:
                # Fallback to manual nearest neighbor if scipy not available
                rows, cols = np.meshgrid(np.arange(filled_data.shape[0]), 
                                         np.arange(filled_data.shape[1]), 
                                         indexing='ij')
                
                valid_mask = ~np.isnan(filled_data)
                valid_rows = rows[valid_mask]
                valid_cols = cols[valid_mask]
                valid_values = filled_data[valid_mask]
                
                invalid_rows = rows[~valid_mask]
                invalid_cols = cols[~valid_mask]
                
                if len(invalid_rows) > 0 and len(valid_rows) > 0:
                    # Use vectorized nearest neighbor for remaining NaN values
                    for inv_row, inv_col in zip(invalid_rows, invalid_cols):
                        distances_sq = (valid_rows - inv_row)**2 + (valid_cols - inv_col)**2
                        nearest_idx = np.argmin(distances_sq)
                        filled_data[inv_row, inv_col] = valid_values[nearest_idx]
        
        # Strategy 3: If still any NaN, fill with minimum valid value or extrapolate from edges
        if np.any(np.isnan(filled_data)):
            # Try to extrapolate from valid edge values
            # Get edge values (first and last row/column)
            edge_values = []
            if filled_data.shape[0] > 0 and filled_data.shape[1] > 0:
                # Top row
                top_row = filled_data[0, :]
                edge_values.extend(top_row[~np.isnan(top_row)])
                # Bottom row
                bottom_row = filled_data[-1, :]
                edge_values.extend(bottom_row[~np.isnan(bottom_row)])
                # Left column
                left_col = filled_data[:, 0]
                edge_values.extend(left_col[~np.isnan(left_col)])
                # Right column
                right_col = filled_data[:, -1]
                edge_values.extend(right_col[~np.isnan(right_col)])
            
            if len(edge_values) > 0:
                fill_value = np.mean(edge_values)
            else:
                valid_min = np.nanmin(filled_data)
                fill_value = valid_min if not np.isnan(valid_min) else 0
            
            filled_data[np.isnan(filled_data)] = fill_value
        
        return filled_data
    
    # Fill NaN values in all contour arrays
    kWe_avg_flipped = fill_contour_nans(kWe_avg_flipped)
    kWt_avg_flipped = fill_contour_nans(kWt_avg_flipped)
    Tout_flipped = fill_contour_nans(Tout_flipped)
    Pout_flipped = fill_contour_nans(Pout_flipped)
    
    # Final verification: ensure no NaN or Inf values remain
    for name, data in [("kWe_avg", kWe_avg_flipped), ("kWt_avg", kWt_avg_flipped), 
                       ("Tout", Tout_flipped), ("Pout", Pout_flipped)]:
        if np.any(np.isnan(data)) or np.any(np.isinf(data)):
            # If still NaN/Inf, fill with nearest valid value or edge value
            valid_mask = ~(np.isnan(data) | np.isinf(data))
            if np.any(valid_mask):
                fill_value = np.nanmean(data[valid_mask])
                if np.isnan(fill_value) or np.isinf(fill_value):
                    fill_value = 0
                data[np.isnan(data) | np.isinf(data)] = fill_value
            else:
                data[:] = 0

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
            z=kWt_avg_flipped, y=param_ij[0,:], x=mdot_ij[:,0],
            connectgaps=True  # Connect across NaN gaps
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
                            is_plot_ts_check,
                            HyperParam1=None
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
    elif model == "SBT V2.0":
        sbt_version = 2
    else:
        sbt_version = 0
    
    sbt_version_sco2 = sbt_version
    sbt_version_h2o = sbt_version
    
    # For simulator (SBT models) only, not database (HDF5)
    # Override sbt_version for each fluid based on fluid selection
    if model in ("SBT V1.0", "SBT V2.0"):
        if fluid == "H2O":
            sbt_version_h2o = 1
            # print(f"[DEBUG generate_econ_lineplots] model={model}, case={case}, fluid={fluid} -> Running H2O only with SBT V{sbt_version_h2o}.0", flush=True)
        elif fluid == "sCO2":
            sbt_version_sco2 = 2
            # print(f"[DEBUG generate_econ_lineplots] model={model}, case={case}, fluid={fluid} -> Running sCO2 only with SBT V{sbt_version_sco2}.0", flush=True)
        elif fluid == "All":
            sbt_version_h2o = 1
            sbt_version_sco2 = 2
            # print(f"[DEBUG generate_econ_lineplots] model={model}, case={case}, fluid={fluid} -> Running H2O: SBT V{sbt_version_h2o}.0, sCO2: SBT V{sbt_version_sco2}.0", flush=True)
    
    # Convert numeric versions to model strings for create_teaobject
    model_sco2 = "HDF5" if sbt_version_sco2 == 0 else (f"SBT V{sbt_version_sco2}.0" if sbt_version_sco2 > 0 else model)
    model_h2o = "HDF5" if sbt_version_h2o == 0 else (f"SBT V{sbt_version_h2o}.0" if sbt_version_h2o > 0 else model)

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

    fig = make_subplots(rows=3, cols=5,
                        specs=[[{'colspan': 2}, None, {'colspan': 2}, None, {"type": "table"}],
                                [{'colspan': 2}, None, {'colspan': 2}, None, {"type": "table"}],
                                [{'colspan': 2}, None, None, None, None]
                                ],
                        horizontal_spacing = 0.11,
                        vertical_spacing = 0.18
                        )
    
    fig.data = []

    teaobj_sCO2 = None
    teaobj_H2O = None
    teaobj_H2O_electricity = None  # Separate variable for electricity section
    teaobj_sCO2_electricity = None  # Separate variable for electricity section

    # ts_fig = make_subplots(rows=1, cols=1,
    #                         horizontal_spacing = 0.11
    #                         )

    if end_use == "Heating" or end_use == "All":

        if end_use == "Heating":
            is_display_legend = True
        if end_use == "All":
            is_display_legend = True
            is_display_legend_electricity = False

        if fluid == "sCO2" or fluid == "All":

            try:
                required_keys = ["time", "sCO2_Tout", "sCO2_Pout"]
                if not TandP_dict or not all(key in TandP_dict for key in required_keys):
                    teaobj_sCO2 = None
                else:
                    time_arr = TandP_dict.get("time", [])
                    tout_arr = TandP_dict.get("sCO2_Tout", [])
                    pout_arr = TandP_dict.get("sCO2_Pout", [])
                    
                    if (not time_arr or not tout_arr or not pout_arr or 
                        len(time_arr) == 0 or len(tout_arr) == 0 or len(pout_arr) == 0 or
                        len(time_arr) != len(tout_arr) or len(time_arr) != len(pout_arr)):
                        teaobj_sCO2 = None
                    else:
                        print(f"[DEBUG generate_econ_lineplots] Creating sCO2 TEA object with model={model_sco2}", flush=True)
                        teaobj_sCO2 = create_teaobject(TandP_dict,
                                                        u_sCO2, u_H2O, c_sCO2, c_H2O,
                                                        case, end_use, "sCO2", model_sco2,
                                                        mdot, L2, L1, grad, D, Tinj, k,
                                                        Drilling_cost_per_m, Discount_rate, Lifetime, 
                                                        Direct_use_heat_cost_per_kWth, Power_plant_cost_per_kWe, Pre_Cooling_Delta_T, Turbine_outlet_pressure,
                                                        properties_H2O_pathname, 
                                                        properties_CO2v2_pathname, 
                                                        additional_properties_CO2v2_pathname,
                                                        is_heating=True,
                                                        HyperParam1=HyperParam1)
                
                if teaobj_sCO2 is None:
                    lcoh_sCO2 = "Insufficient Inputs"
                    mean_sCO2_Net_HProd = '-'
                    fig, lcoh_sCO2 = update_blank_econ2(fig=fig, nrow1=1, ncol1=1, nrow2=1, ncol2=3)
                else:
                    try:
                        if not hasattr(teaobj_sCO2, 'LCOH') or teaobj_sCO2.LCOH is None or teaobj_sCO2.LCOH >= 9999:
                            lcoh_sCO2 = "Insufficient Inputs"
                        else:
                            lcoh_sCO2 = format(teaobj_sCO2.LCOH, '.2f')
                        has_linear_time = hasattr(teaobj_sCO2, 'Linear_time_distribution') and teaobj_sCO2.Linear_time_distribution is not None
                        has_inst_heat = hasattr(teaobj_sCO2, 'Instantaneous_heat_production') and teaobj_sCO2.Instantaneous_heat_production is not None
                        has_annual_heat = hasattr(teaobj_sCO2, 'Annual_heat_production') and teaobj_sCO2.Annual_heat_production is not None
                        has_lifetime = hasattr(teaobj_sCO2, 'Lifetime') and teaobj_sCO2.Lifetime is not None
                        
                        if not (has_linear_time and has_inst_heat and has_annual_heat and has_lifetime):
                            raise AttributeError("TEA object missing required attributes for plotting")
                        
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
                    except Exception as e:
                        print(f"[ERROR] Economic plots (sCO2): Exception when creating plots: {type(e).__name__}: {e}", flush=True)
                        import traceback
                        traceback.print_exc()
                        lcoh_sCO2 = "Insufficient Inputs"
                        mean_sCO2_Net_HProd = '-'
                        fig, lcoh_sCO2 = update_blank_econ2(fig=fig, nrow1=1, ncol1=1, nrow2=1, ncol2=3)

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
                required_keys = ["time", "H2O_Tout", "H2O_Pout"]
                if not TandP_dict or not all(key in TandP_dict for key in required_keys):
                    teaobj_H2O = None
                else:
                    time_arr = TandP_dict.get("time") or []
                    tout_arr = TandP_dict.get("H2O_Tout") or []
                    pout_arr = TandP_dict.get("H2O_Pout") or []
                    
                    if (not time_arr or not tout_arr or not pout_arr or 
                        len(time_arr) == 0 or len(tout_arr) == 0 or len(pout_arr) == 0 or
                        len(time_arr) != len(tout_arr) or len(time_arr) != len(pout_arr)):
                        teaobj_H2O = None
                    else:
                        # TODO: update D ... based on radial
                        print(f"[DEBUG generate_econ_lineplots] Creating H2O TEA object with model={model_h2o}", flush=True)
                        teaobj_H2O = create_teaobject(TandP_dict,
                                                    u_sCO2, u_H2O, c_sCO2, c_H2O,
                                                    case, end_use, "H2O", model_h2o,
                                                    mdot, L2, L1, grad, D, Tinj, k,
                                                    Drilling_cost_per_m, Discount_rate, Lifetime, 
                                                    Direct_use_heat_cost_per_kWth, Power_plant_cost_per_kWe, Pre_Cooling_Delta_T, Turbine_outlet_pressure,
                                                    properties_H2O_pathname, 
                                                    properties_CO2v2_pathname, 
                                                    additional_properties_CO2v2_pathname,
                                                    is_H20=True, is_heating=True,
                                                    HyperParam1=HyperParam1
                                                    )
                
                if teaobj_H2O is None:
                    lcoh_H2O = "Insufficient Inputs"
                    # Add blank plots when TEA object is None
                    fig, lcoh_H2O = update_blank_econ2(fig=fig, nrow1=1, ncol1=1, nrow2=1, ncol2=3)
                else:
                    lcoh_H2O = "Insufficient Inputs" if (teaobj_H2O.LCOH is None or teaobj_H2O.LCOH >= 9999) else format(teaobj_H2O.LCOH, '.2f')

                    # Heat Production 
                    fig.add_trace(go.Scatter(x=teaobj_H2O.Linear_time_distribution, y=teaobj_H2O.Instantaneous_heat_production/1e3,
                              hovertemplate='<b>Time (year)</b>: %{x:.1f}<br><b>Heat Production (MWt)</b>: %{y:.3f} ',
                              line = dict(color='black', width=lw, dash='dash'),
                              legendgroup=labels_cat[0], name=labels[0], showlegend=is_display_legend
                              ),
                              row=1, col=1)

                    # Anuual Heat Production
                    fig.add_trace(go.Bar(x=np.arange(1,teaobj_H2O.Lifetime+1), y=teaobj_H2O.Annual_heat_production/1e6,
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
                required_keys = ["time", "sCO2_Tout", "sCO2_Pout"]
                if not TandP_dict or not all(key in TandP_dict for key in required_keys):
                    teaobj_sCO2_electricity = None
                else:
                    time_arr = TandP_dict.get("time", [])
                    tout_arr = TandP_dict.get("sCO2_Tout", [])
                    pout_arr = TandP_dict.get("sCO2_Pout", [])
                    
                    if (not time_arr or not tout_arr or not pout_arr or 
                        len(time_arr) == 0 or len(tout_arr) == 0 or len(pout_arr) == 0 or
                        len(time_arr) != len(tout_arr) or len(time_arr) != len(pout_arr)):
                        teaobj_sCO2_electricity = None
                    else:
                        print(f"[DEBUG generate_econ_lineplots] Creating sCO2 electricity TEA object with model={model_sco2}", flush=True)
                        teaobj_sCO2_electricity = create_teaobject(TandP_dict, 
                                                    u_sCO2, u_H2O, c_sCO2, c_H2O,
                                                    case, end_use, "sCO2", model_sco2,
                                                    mdot, L2, L1, grad, D, Tinj, k,
                                                    Drilling_cost_per_m, Discount_rate, Lifetime, 
                                                    Direct_use_heat_cost_per_kWth, Power_plant_cost_per_kWe, Pre_Cooling_Delta_T, Turbine_outlet_pressure,
                                                    properties_H2O_pathname, 
                                                    properties_CO2v2_pathname, 
                                                    additional_properties_CO2v2_pathname,
                                                    HyperParam1=HyperParam1)
                
                if teaobj_sCO2_electricity is None:
                    lcoe_sCO2 = "Insufficient Inputs"
                    # Add blank plots when TEA object is None
                    fig, lcoe_sCO2 = update_blank_econ2(fig=fig, nrow1=row_num, ncol1=1, nrow2=row_num, ncol2=3)
                else:
                    # convert any negative value to 0
                    if teaobj_sCO2_electricity.Inst_Net_Electricity_production is not None:
                        teaobj_sCO2_electricity.Inst_Net_Electricity_production[teaobj_sCO2_electricity.Inst_Net_Electricity_production<0] = 0

                    lcoe_sCO2 = "Insufficient Inputs" if (teaobj_sCO2_electricity.LCOE is None or teaobj_sCO2_electricity.LCOE >= 9999) else format(teaobj_sCO2_electricity.LCOE, '.2f')

                    # Electricity 
                    if teaobj_sCO2_electricity.Inst_Net_Electricity_production is not None:
                        fig.add_trace(go.Scatter(x=teaobj_sCO2_electricity.Linear_time_distribution, y=teaobj_sCO2_electricity.Inst_Net_Electricity_production/1e3,
                                      hovertemplate='<b>Time (year)</b>: %{x:.1f}<br><b>Electricity Production (MWe)</b>: %{y:.3f} ',
                                      line = dict(color='black', width=lw),
                                      legendgroup=labels_cat[1], name=labels[1], showlegend=is_display_legend
                                      ),
                                  row=row_num, col=1)

                    
                    # Annual Electricity
                    fig.add_trace(go.Bar(x=np.arange(1,teaobj_sCO2_electricity.Lifetime+1), y=teaobj_sCO2_electricity.Annual_electricity_production/1e6,
                                    # hovertemplate='<b>Time (year)</b>: %{x:.1f}<br><b>Annual Electricity Production (GWe)</b>: %{y:.3f} ',
                                    name=labels[1],
                                    showlegend=False, hovertemplate='<b>Annual Electricity Production (GWe)</b>: %{y:.3f}'), 
                                    row=row_num, col=3)

                    fig.add_trace(go.Scatter(x=np.arange(1,teaobj_sCO2_electricity.Lifetime+1), y=teaobj_sCO2_electricity.Annual_electricity_production/1e6,
                                  hovertemplate='<b>Time (year)</b>: %{x:.1f}<br>',
                                  line = dict(color='black', width=lw),
                                  legendgroup=labels_cat[1], name=labels[1], showlegend=False
                                  ),
                                  row=row_num, col=3)

                    # Table Data
                    if teaobj_sCO2_electricity.Inst_Net_Electricity_production is not None:
                        time_dict = {"Time (year)": teaobj_sCO2_electricity.Linear_time_distribution,
                                     "sCO2 Electricity Production (MWe)": teaobj_sCO2_electricity.Inst_Net_Electricity_production/1e3,
                                     "sCO2 Annual Electricity Production (GWe)": teaobj_sCO2_electricity.Annual_electricity_production/1e6}
                        econ_values_dict.update(time_dict)

            except ValueError as e:
                fig, lcoe_sCO2 = update_blank_econ2(fig=fig, nrow1=row_num, ncol1=1, nrow2=row_num, ncol2=3)
                error_message = parse_error_message(e=e, e_name='Err Econ3a')
                error_messages_dict['Err Econ3a'] = error_message

            except AttributeError as e:
                fig, lcoe_sCO2 = update_blank_econ2(fig=fig, nrow1=row_num, ncol1=1, nrow2=row_num, ncol2=3)
                error_message = parse_error_message(e=e, e_name='Err Econ3b')
                error_messages_dict['Err Econ3b'] = error_message

        if is_plot_ts_check and (fluid == "sCO2" or fluid == "All"):
            ts_teaobj = None
            if (end_use == "Electricity" or end_use == "All") and teaobj_sCO2_electricity is not None:
                ts_teaobj = teaobj_sCO2_electricity
            elif (end_use == "Heating" or end_use == "All") and teaobj_sCO2 is not None:
                ts_teaobj = teaobj_sCO2
            
            if ts_teaobj is not None:
                if end_use == "Electricity":
                    try:
                        get_Ts_diagram(fig=fig, teaobj=ts_teaobj, nrow=2, ncol=1, tmatrix_pathname=tmatrix_pathname)
                    except Exception as e:
                        print(f"[ERROR] Error drawing TS diagram (Electricity): {e}")
                        traceback.print_exc()
                if end_use == "All":
                    try:
                        get_Ts_diagram(fig=fig, teaobj=ts_teaobj, nrow=3, ncol=1, tmatrix_pathname=tmatrix_pathname)
                    except Exception as e:
                        print(f"[ERROR] Error drawing TS diagram (All): {e}")
                        traceback.print_exc()

                # mean_sCO2_Net_HProd = round(np.mean(teaobj_sCO2.Instantaneous_heat_production/1e3),2)
                mean_sCO2_Net_EProd = round(np.mean(teaobj_sCO2_electricity.Inst_Net_Electricity_production/1e3),2) if teaobj_sCO2_electricity is not None and hasattr(teaobj_sCO2_electricity, 'Inst_Net_Electricity_production') and teaobj_sCO2_electricity.Inst_Net_Electricity_production is not None else 0

        if fluid == "H2O" or fluid == "All":

            try:
                required_keys = ["time", "H2O_Tout", "H2O_Pout"]
                if not TandP_dict or not all(key in TandP_dict for key in required_keys):
                    teaobj_H2O_electricity = None
                else:
                    time_arr = TandP_dict.get("time") or []
                    tout_arr = TandP_dict.get("H2O_Tout") or []
                    pout_arr = TandP_dict.get("H2O_Pout") or []
                    
                    if (not time_arr or not tout_arr or not pout_arr or 
                        len(time_arr) == 0 or len(tout_arr) == 0 or len(pout_arr) == 0 or
                        len(time_arr) != len(tout_arr) or len(time_arr) != len(pout_arr)):
                        teaobj_H2O_electricity = None
                    else:
                        teaobj_H2O_electricity = create_teaobject(TandP_dict,
                                                        u_sCO2, u_H2O, c_sCO2, c_H2O,
                                                        case, end_use, "H2O", model_h2o,
                                                        mdot, L2, L1, grad, D, Tinj, k,
                                                        Drilling_cost_per_m, Discount_rate, Lifetime, 
                                                        Direct_use_heat_cost_per_kWth, Power_plant_cost_per_kWe, Pre_Cooling_Delta_T, Turbine_outlet_pressure,
                                                        properties_H2O_pathname, 
                                                        properties_CO2v2_pathname, 
                                                        additional_properties_CO2v2_pathname,
                                                        is_H20=True,
                                                        HyperParam1=HyperParam1
                                                        )
                
                if teaobj_H2O_electricity is None:
                    lcoe_H2O = "Insufficient Inputs"
                    fig, lcoe_H2O = update_blank_econ2(fig=fig, nrow1=row_num, ncol1=1, nrow2=row_num, ncol2=3)
                else:
                    if teaobj_H2O_electricity.Inst_Net_Electricity_production is not None:
                        teaobj_H2O_electricity.Inst_Net_Electricity_production[teaobj_H2O_electricity.Inst_Net_Electricity_production<0] = 0

                    lcoe_H2O = "Insufficient Inputs" if (teaobj_H2O_electricity.LCOE is None or teaobj_H2O_electricity.LCOE >= 9999) else format(teaobj_H2O_electricity.LCOE, '.2f')

                    # Electricity 
                    if teaobj_H2O_electricity.Inst_Net_Electricity_production is not None:
                        fig.add_trace(go.Scatter(x=teaobj_H2O_electricity.Linear_time_distribution, y=teaobj_H2O_electricity.Inst_Net_Electricity_production/1e3,
                                      hovertemplate='<b>Time (year)</b>: %{x:.1f}<br><b>Electricity Production (MWe)</b>: %{y:.3f} ',
                                      line = dict(color='black', width=lw, dash='dash'),
                                      legendgroup=labels_cat[0], name=labels[0], showlegend=is_display_legend
                                      ),
                                  row=row_num, col=1)
                    
                    # Annual Electricity
                    fig.add_trace(go.Bar(x=np.arange(1,teaobj_H2O_electricity.Lifetime+1), y=teaobj_H2O_electricity.Annual_electricity_production/1e6,
                                    name=labels[0],
                                    showlegend=False, hovertemplate='<b>Annual Electricity Production (GWe)</b>: %{y:.3f}'), 
                                    row=row_num, col=3)

                    fig.add_trace(go.Scatter(x=np.arange(1,teaobj_H2O_electricity.Lifetime+1), y=teaobj_H2O_electricity.Annual_electricity_production/1e6,
                                  hovertemplate='<b>Time (year)</b>: %{x:.1f}<br>',
                                  line = dict(color='black', width=lw, dash='dash'),
                                  legendgroup=labels_cat[0], name=labels[0], showlegend=False
                                  ),
                                  row=row_num, col=3)

                    # Table Data
                    if teaobj_H2O_electricity.Inst_Net_Electricity_production is not None:
                        time_dict = {"Time (year)": teaobj_H2O_electricity.Linear_time_distribution,
                                     "H2O Electricity Production (MWe)": teaobj_H2O_electricity.Inst_Net_Electricity_production/1e3,
                                     "H2O Annual Electricity Production (GWe)": teaobj_H2O_electricity.Annual_electricity_production/1e6}
                        econ_values_dict.update(time_dict)

                    mean_H2O_Net_EProd = round(np.mean(teaobj_H2O_electricity.Inst_Net_Electricity_production/1e3),2) if teaobj_H2O_electricity.Inst_Net_Electricity_production is not None else 0

            except ValueError as e:
                fig, lcoe_H2O = update_blank_econ2(fig=fig, nrow1=row_num, ncol1=1, nrow2=row_num, ncol2=3)
                error_message = parse_error_message(e=e, e_name='Err Econ4a')
                error_messages_dict['Err Econ4a'] = error_message

            except AttributeError as e:
                fig, lcoe_H2O = update_blank_econ2(fig=fig, nrow1=row_num, ncol1=1, nrow2=row_num, ncol2=3)
                error_message = parse_error_message(e=e, e_name='Err Econ4b')
                error_messages_dict['Err Econ4b'] = error_message

    fig = update_layout_properties_econ_results(fig, end_use, scale, is_plot_ts_check, fluid)
    fig = update_lcoh_lcoe_table(fig, fluid, end_use, lcoh_sCO2, lcoh_H2O, lcoe_sCO2, lcoe_H2O) # table
    fig.update_traces(cells_font=dict(size = 13), row=1, col=5)
    fig.update_traces(cells_font=dict(size = 13), row=2, col=5)

    error_codes = []
    #getting errors codes
    if teaobj_H2O is not None:
        error_codes += teaobj_H2O.error_codes.tolist()
    
    if teaobj_H2O_electricity is not None:
        error_codes += teaobj_H2O_electricity.error_codes.tolist()

    if teaobj_sCO2 is not None:
        error_codes += teaobj_sCO2.error_codes.tolist()
    #removing duplicates
    error_codes = list(set(error_codes))

    econ_data_dict = {'LCOH sCO2': lcoh_sCO2, 
                        'LCOH H2O': lcoh_H2O, 
                        'LCOE sCO2': lcoe_sCO2,
                        'LCOE H2O': lcoe_H2O,
                        'Mean H2O Net HProd': mean_H2O_Net_HProd,
                        'Mean H2O Net EProd': mean_H2O_Net_EProd,
                        'Mean sCO2 Net HProd': mean_sCO2_Net_HProd,
                        'Mean sCO2 Net EProd': mean_sCO2_Net_EProd,
                        "error_codes": error_codes}


    fig.update_layout(paper_bgcolor='rgba(255,255,255,0.10)',
                      plot_bgcolor='rgba(255,255,255,0)')
    
    import time
    # Include checkbox state in uirevision and set datarevision to force recalculation
    fig.update_layout(
        uirevision=f"model-ts={int(is_plot_ts_check)}-{time.time()}",
        datarevision=time.time()
    )

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

    # Check if teaobj has required attributes
    if teaobj is None:
        return fig
    
    required_attrs = ['Linear_production_temperature', 's_prod', 'T_turbine_out_actual', 
                      'Pvector', 'Tvector', 'entropy', 'Turbine_outlet_pressure']
    if not all(hasattr(teaobj, attr) for attr in required_attrs):
        print(f"Warning: teaobj missing required attributes for TS diagram")
        return fig

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


    fig.update_yaxes(title_text="Temperature (°C)", 
                            row=nrow, col=ncol,
                            tickfont = dict(size=12), title_font=dict(size=14))
    fig.update_xaxes(title_text="Entropy (J/kg/°C)", range=[1000, 2300],
                            row=nrow, col=ncol,
                            tickfont = dict(size=12), title_font=dict(size=14))
    return fig


