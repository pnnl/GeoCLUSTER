#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# sourced scripts
import clg_tea_module as clg_tea_module_v3

# -----------------------
# Global properties.
# -----------------------
to_kelvin_factor = 273.15

def fixed_economic_inputs():

    # --------------------------------------------------------------------------------
    # Inputs that don't get changed by the app.
    # --------------------------------------------------------------------------------

    O_and_M_cost_plant = 0.015                 # Operation & maintance cost of surface plant expressed as fraction of total surface plant capital cost [-] (between 0 and 0.2)
    Pump_efficiency = 0.8                      # Pump efficiency for cirulcation pump (if required) [-] (between 0.5 and 1)
    Turbine_isentropic_efficiency = 0.9        # Isentropic efficiency for turbine when CO2 is working fluid [-] (between 0.8 and 1)
    Generator_efficiency = 0.98                # Conversion efficiency from mechanical turbine work to electricity [-] (between 0.8 and 1)
    Compressor_isentropic_efficiency = 0.9     # Isentropic efficiency for compressor when CO2 is working fluid [-] (between 0.8 and 1)
    T0 = 20 + to_kelvin_factor # to kelvin     # Dead-state temperature [K] (between 278.15 and 303.15)
    P0 = 1e5                                   # Dead-state pressure [Pa] (between 0.8e5 and 1.1e5)
    Electricity_rate = 0.1                     # Electricity rate in direct-use for pumping power (if pumping is required) [$/kWh] (between 0 and 0.5)

    return O_and_M_cost_plant, Pump_efficiency, Turbine_isentropic_efficiency, Generator_efficiency, Compressor_isentropic_efficiency, T0, P0, Electricity_rate


def create_teaobject(TandP_dict, 
                    u_sCO2, u_H2O, c_sCO2, c_H2O,
                     case, end_use, fluid, model,
                     Flow_user, Hor_length_user, Depth_user, Gradient_user, Diameter_user, Tin_user, krock_user, 
                     Drilling_cost_per_m, Discount_rate, Lifetime, 
                     Direct_use_heat_cost_per_kWth, Power_plant_cost_per_kWe, Pre_Cooling_Delta_T, Turbine_outlet_pressure, 
                     properties_H2O_pathname, 
                     properties_CO2v2_pathname, 
                     additional_properties_CO2v2_pathname,
                     is_heating=False, is_H20=False
                     ):
    
    # ------------------------------------------
    # Create TEA object for economic results.
    # ------------------------------------------
    # print(Flow_user, Gradient_user, Diameter_user, Tin_user, krock_user)

    # Handle None values
    if Flow_user is None:
        Flow_user = 20.0
    if Hor_length_user is None:
        Hor_length_user = 10000.0
    if Depth_user is None:
        Depth_user = 3500.0
    if Gradient_user is None:
        Gradient_user = 0.05
    if Diameter_user is None:
        Diameter_user = 0.35
    if Tin_user is None:
        Tin_user = 30.0
    if krock_user is None:
        krock_user = 3.0
    if Drilling_cost_per_m is None:
        Drilling_cost_per_m = 1000.0
    if Discount_rate is None:
        Discount_rate = 7.0
    if Lifetime is None:
        Lifetime = 40
    if Direct_use_heat_cost_per_kWth is None:
        Direct_use_heat_cost_per_kWth = 100.0
    if Power_plant_cost_per_kWe is None:
        Power_plant_cost_per_kWe = 3000.0
    if Pre_Cooling_Delta_T is None:
        Pre_Cooling_Delta_T = 13.0
    if Turbine_outlet_pressure is None:
        Turbine_outlet_pressure = 80.0

    # Gradient_user = Gradient_user / 1000
    Tin_user = Tin_user + to_kelvin_factor # to kelvin
    Discount_rate = Discount_rate / 100 
    
    # Inputs that don't get changed by the GUI
    O_and_M_cost_plant, Pump_efficiency, Turbine_isentropic_efficiency, \
        Generator_efficiency, Compressor_isentropic_efficiency, T0, P0, Electricity_rate = fixed_economic_inputs()

    if case == 'utube':
        Configuration = 1  # Closed-loop configuration: 1 is U-loop; 2 is co-axial with injection in annulus
    else:
        Configuration = 2
    
    if end_use == 'Heating':
        End_use = 1        # End-use application: 1 is heating, 2 is electricity (must be 1 or 2)    
    else:
        End_use = 2 # 'Electricity or All', then get End_use == 2, but if is_heating can switch
        if is_heating:
            End_use = 1

    if fluid == 'H2O':
        Fluid = 1          # Heat transfer fluid: 1 is water, 2 is sCO2 (must be 1 or 2)
    if fluid == 'sCO2':
        Fluid = 2
    if fluid == 'All':
        if is_H20:
            Fluid = 1
        else:
            Fluid = 2


    # create object
    # print(u_sCO2, u_H2O, c_sCO2, c_H2O)
    teaobject = clg_tea_module_v3.TEA(u_sCO2, u_H2O, c_sCO2, c_H2O,
                                          Fluid, End_use, Configuration, # user
                                          Flow_user, Hor_length_user, Depth_user, Gradient_user, Diameter_user, # user  
                                          Tin_user, krock_user, # user 
                                          Drilling_cost_per_m, # user
                                          O_and_M_cost_plant, 
                                          Discount_rate, # user
                                          Pump_efficiency, 
                                          Lifetime, # user 
                                          Direct_use_heat_cost_per_kWth, # user
                                          Electricity_rate, 
                                          Power_plant_cost_per_kWe, # user
                                          T0, P0, # initialized
                                          Turbine_isentropic_efficiency, Generator_efficiency, # initialized
                                          Compressor_isentropic_efficiency, 
                                          Pre_Cooling_Delta_T, # user
                                          Turbine_outlet_pressure, # user
                                          properties_H2O_pathname, # read in 
                                          properties_CO2v2_pathname, # read in 
                                          additional_properties_CO2v2_pathname # read in
                                          )

    # load case and property data
    teaobject.initialize(TandP_dict, u_sCO2, u_H2O, c_sCO2, c_H2O, 
                                          properties_H2O_pathname, 
                                          properties_CO2v2_pathname, 
                                          additional_properties_CO2v2_pathname)

    # get interpolated temperature and pressure array
    teaobject.getTandP(u_sCO2, u_H2O, c_sCO2, c_H2O, model, TandP_dict)
    teaobject.calculateLC() # ERROR STARTS HERE
    # teaobject.printresults() # uncomment to debug
    
    return teaobject
    
