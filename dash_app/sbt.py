#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import time
#import tensorflow as tf
import scipy.special as sp
from scipy.special import erf,erfc, jv, yv,exp1
import pdb
import scipy.io
import math
from scipy.interpolate import RegularGridInterpolator
from traceback import print_stack
#has SBT v1 for co-axial and U-loop, SBT v2 for co-axial,as well as FMM algorithm

# sourced scripts
is_plot = False
from plot_sbt import plot_borehole_geometry, plot_final_fluid_temp_profile_v1, plot_final_fluid_temp_profile_v2
from plot_sbt import plot_heat_production, plot_production_temperature_linear, plot_production_tempterature_log

#%% -------
# 1. Input
# Generally, the user should only make changes to this section
#---------

def set_sbt_hyperparameters(sbt_version, clg_configuration, accuracy, mesh_fineness, fluid, HYPERPARAM1, HYPERPARAM2, HYPERPARAM3, HYPERPARAM4, HYPERPARAM5):

    """ 
        ## SBT Model Specifications
        sbt_version:                  # Must be 1 or 2. 1 means SBT v1 (Temperature only); 2 means SBT v2 (Temperature and Pressure including allowing TP dependent fluid properties)
        
    """
    # TODO: Really this function should be able the speed and accuracy of the SBT model

    ## FINETUNE ACROSS V1 AND V2 (simulation and SBT algorithm settings)

    # accuracy = 1                          # Must be 1,2,3,4 or 5 with 1 lowest accuracy and 5 highest accuracy. Lowest accuracy runs fastest. Accuracy level impacts number of discretizations for numerical integration and decision tree thresholds in SBT algorithm.
    FMM = 1                               # if 1, use fast multi-pole methold-like approach (i.e., combine old heat pulses to speed up simulation)
    FMMtriggertime = 3600*24*10           # threshold time beyond which heat pulses can be combined with others [s]
    # mesh_fineness = 0                   # AVB: less mesh (0), finer mesh (1), finest mesh (2): SEE DETAILS BELOW FOR KOENRAAD'S BREAKDOWN

    # load accuracy parameters    
    if accuracy == 1:
        NoArgumentsFinitePipeCorrection = 25
        NoDiscrFinitePipeCorrection = 200
        NoArgumentsInfCylIntegration = 25
        NoDiscrInfCylIntegration = 200
        LimitPointSourceModel = 1.5
        LimitCylinderModelRequired = 25
        LimitInfiniteModel = 0.05
        LimitNPSpacingTime = 0.1
        LimitSoverL = 1.5
        M = 3
    elif accuracy == 2:
        NoArgumentsFinitePipeCorrection = 50
        NoDiscrFinitePipeCorrection = 400
        NoArgumentsInfCylIntegration = 50
        NoDiscrInfCylIntegration = 400
        LimitPointSourceModel = 2.5
        LimitCylinderModelRequired = 50
        LimitInfiniteModel = 0.01
        LimitNPSpacingTime = 0.04
        LimitSoverL = 2
        M = 4
    elif accuracy == 3:
        NoArgumentsFinitePipeCorrection = 100
        NoDiscrFinitePipeCorrection = 500
        NoArgumentsInfCylIntegration = 100
        NoDiscrInfCylIntegration = 500
        LimitPointSourceModel = 5
        LimitCylinderModelRequired = 100
        LimitInfiniteModel = 0.004
        LimitNPSpacingTime = 0.02
        LimitSoverL = 3
        M = 5
    elif accuracy == 4:
        NoArgumentsFinitePipeCorrection = 200
        NoDiscrFinitePipeCorrection = 1000
        NoArgumentsInfCylIntegration = 200
        NoDiscrInfCylIntegration = 1000
        LimitPointSourceModel = 10
        LimitCylinderModelRequired = 200
        LimitInfiniteModel = 0.002
        LimitNPSpacingTime = 0.01
        LimitSoverL = 5
        M = 10
    elif accuracy == 5:
        NoArgumentsFinitePipeCorrection = 400
        NoDiscrFinitePipeCorrection = 2000
        NoArgumentsInfCylIntegration = 400
        NoDiscrInfCylIntegration = 2000
        LimitPointSourceModel = 20
        LimitCylinderModelRequired = 400
        LimitInfiniteModel = 0.001
        LimitNPSpacingTime = 0.005
        LimitSoverL = 9
        M = 20


    if mesh_fineness == 0:
        times = np.concatenate((np.linspace(0,9900,100), np.logspace(np.log10(100*100), np.log10(40*365*24*3600), 75))) # simulation times [s] (must start with 0; to obtain smooth results, abrupt changes in time step size should be avoided. logarithmic spacing is recommended)
        #times = np.concatenate((np.linspace(0,9900,100),  np.linspace(10000, int(20*3.1*10**7), num=(int(20*3.1*10**7) - 10000) // 3600 + 1)))
    elif mesh_fineness == 1:
        #Note 1: When providing a variable injection temperature or flow rate, a finer time grid may be required. Below is an example with long term time steps of about 36 days.
        times = [0] + list(range(100, 10000, 100)) + list(np.logspace(np.log10(100*100), np.log10(0.1*365*24*3600), 40)) + list(np.arange(0.2*365*24*3600, 20*365*24*3600, 0.1*365*24*3600))
    elif mesh_fineness == 2:
        #Note 2: To capture the start-up effects, several small time steps are taken during the first 10,000 seconds in the time vector considered. To speed up the simulation, this can be avoided with limited impact on the long-term results. For example, an alternative time vector would be:
        times = [0] + list(range(100, 1000, 100)) + list(range(1000, 10000, 1000)) + list(np.logspace(np.log10(100*100), np.log10(20*365*24*3600), 75))
    
    fullyimplicit = None
    if clg_configuration == 2:
        fullyimplicit = 1            # Should be between 0 and 1. Only required when clg_configuration is 2. Most stable is setting it to 1 which results in a fully implicit Euler scheme when calculting the fluid temperature at each time step. With a value of 0, the convective term is modelled using explicit Euler. A value of 0.5 would model the convective term 50% explicit and 50% implicit, which may be slightly more accurate than fully implicit.


    variablefluidproperties = None
    ## SBT ASSUMPTIONS FOR V1 v.s. V2
    if sbt_version == 1: ## SBT v1 specific settings (ONLY WATER)
        
        variableflowrate = HYPERPARAM1 #0                                            # Must be 0 or 1. "0" means the user provides a constant mass flow rate m. "1" means the user provides an excel file with a mass flow rate profile. [only works in SBT v1]
        flowratefilename = HYPERPARAM2 # 'MassFlowRate.xlsx'                         # Name of excel file with mass flow rate profile. Must be provided if user sets variableflowrate to 1. First column stores time in seconds, second column stores mass flow rate in kg/s. Time must start with 0 and end with the final simulation time chosen (as specified in the array times"). [only works in SBT v1]
        variableinjectiontemperature = HYPERPARAM3 # 0                               # Must be 0 or 1. "0" means the user provides a constant injection temperature Tin. "1" means the user provides an excel file with an injection temperature profile. [only works in SBT v1]
        injectiontemperaturefilename = HYPERPARAM4 # 'InjectionTemperatures.xlsx'    # Name of excel file with injection temperature profile. Must be provided if user sets variableinjectiontemperature to 1. First column stores time in seconds, second column stores injection temperature in degrees C. Time must start with 0 and end with the final simulation time chosen (as specified in the array times"). [only works in SBT v1]
        HYPERPARAM5 = None

        if fluid == 2:
            print("ERROR!! CAN'T DO sCO2 for SBT V1.0") #will always only be able to run H20 and NOT sCO2 || # Heat transfer fluid selection: 1 = H2O; 2 = CO2

    if sbt_version == 2: ## SBT v2 specific settings
        
        Pin = HYPERPARAM1 # 100                             # Fluid input pressure [bar]
        piperoughness = HYPERPARAM2 # 1e-6                  # Pipe/borehole roughness to calculate friction pressure drop [m]
        variablefluidproperties = HYPERPARAM3 # 1           # Must be 0 or 1. "0" means the fluid properties remain constant and are specified by cp_f, rho_f, k_f and mu_f. "1" means that the fluid properties (e.g., density) are calculated internally each time step and are a function of temperature and pressure. 
        
        # converge parameters
        reltolerance = HYPERPARAM4 # 1e-5                   # Target maximum acceptable relative tolerance each time step [-]. Lower tolerance will result in more accurate results but requires longer computational time
        maxnumberofiterations = HYPERPARAM5 #15             # Maximum number of iterations each time step [-]. Each time step, solution converges until relative tolerance criteria is met or maximum number of time steps is reached.

    return locals()


def set_wellbore_geometry(clg_configuration, DrillingDepth_L1, HorizontalExtent_L2, numberoflaterals):

    """ 
        ## Wellbore Geometry
        "vertical_depth" and Depth (L1): # L1 works for pure vertical wells as well (HorizontalExtent (L2) = 0)
        HorizontalExtent (L2): 

        #( x,y,z) geometry of CLG heat exchanger
        # The vectors storing the x-, y- and z-coordinates should be column vectors
        # To obtain smooth results, abrupt changes in segment lengths should be avoided.
    """

    x = y = z = None
    xinj = yinj = zinj = zprod = yprod = xprod = xlat = ylat = zlat = None

    if clg_configuration == 1: # co-axial geometry: (x,y,z)-coordinates of centerline of co-axial heat exchanger [m]
        # Example 1: 2 km vertical well
        # TODO: we need link this to depth
        # NOTE: option for speeding up, change step size here

        # NOTE: new code by Koenraad
        verticaldepthsection = np.arange(0, -DrillingDepth_L1*1000-1, -100)
        horizontalextentsection = np.arange(100, HorizontalExtent_L2*1000+1, 100)
        z = np.concatenate((verticaldepthsection, 
                            verticaldepthsection[-1]*np.ones(len(horizontalextentsection)))
                            ).reshape(-1, 1)
        x = np.concatenate((np.zeros(len(verticaldepthsection)),horizontalextentsection)).reshape(-1, 1)
        y = np.zeros((len(z), 1))

        # NOTE: old code by Koenraad (another way to get the geometry just on Vertical Depth alone)
        # vertical_depth = 2000
        # max_dim_1 = -1 * (vertical_depth + 1)
        # z = np.arange(0, max_dim_1, -100).reshape(-1, 1)
        # x = np.zeros((len(z), 1))
        # y = np.zeros((len(z), 1))

    elif clg_configuration == 2: # U-loop geometry: (x,y,z)-coordinates of centerline of injection well, production well and laterals
        # Coordinates of injection well (coordinates are provided from top to bottom in the direction of flow)
        # DrillingDepth_L1 == -2000 (default)
        # HorizontalExtent_L2 == 1000 (default)
        HorizontalExtent_L2_half = HorizontalExtent_L2/2
        zinj = np.arange(0, -DrillingDepth_L1 - 100, -100).reshape(-1, 1) # 2k down, 100 m horizontal, 100 m up
        yinj = np.zeros((len(zinj), 1))
        xinj = -HorizontalExtent_L2_half * np.ones((len(zinj), 1))
        
        # Coordinates of production well (coordinates are provided from bottom to top in the direction of flow)
        zprod = np.arange(-DrillingDepth_L1, 0 + 100, 100).reshape(-1, 1)
        yprod = np.zeros((len(zprod), 1))
        xprod = HorizontalExtent_L2_half * np.ones((len(zprod), 1))
        
        # (x, y, z)-coordinates of laterals are stored in three matrices (one each for the x, y, and z coordinates). 
        # The number of columns in each matrix corresponds to the number of laterals. The number of discretizations 
        # should be the same for each lateral. Coordinates are provided in the direction of flow; the first coordinate should match 
        # the last coordinate of the injection well, and the last coordinate should match the first coordinate of the 
        # production well
        # NOTE: laterals (for validation purposes, have one lateral that doesn't bend if it's 1 lateral)

        xlat = np.concatenate((np.array([-1000, -918, -814])*HorizontalExtent_L2_half/1000, # -1000 x coord
                                np.linspace(-706, 706, 14)*HorizontalExtent_L2_half/1000, 
                                    np.array([814, 918, 1000])*HorizontalExtent_L2_half/1000)).reshape(-1,1)

        ylat = np.concatenate((100 * np.cos(np.linspace(-np.pi/2, 0, 3)), 100 * np.ones(14), 100 *\
        np.cos(np.linspace(0, np.pi/2, 3)))).reshape(-1,1) # COS = curving the lateral
        
        zlat = (-DrillingDepth_L1 * np.ones((len(xlat)))).reshape(-1,1) # -2000 y coord
        
        if numberoflaterals == 1:
            xlat = np.hstack((np.linspace(-HorizontalExtent_L2_half,HorizontalExtent_L2_half,20).reshape(-1,1)))
            ylat = np.hstack((np.zeros(len(ylat)).reshape(-1,1)))
            zlat = np.hstack((zlat)) # same as line 208
            xlat = np.expand_dims(xlat, axis=1) 
            ylat = np.expand_dims(ylat, axis=1) 
            zlat = np.expand_dims(zlat, axis=1) 

        else:
            xlat = np.hstack((xlat,xlat, np.linspace(-HorizontalExtent_L2_half,HorizontalExtent_L2_half,20).reshape(-1,1)))
            ylat = np.hstack((ylat,-ylat, np.zeros(len(ylat)).reshape(-1,1)))
            zlat = np.hstack((zlat,zlat,zlat))
            
        # Merge x-, y-, and z-coordinates
        x = np.concatenate((xinj, xprod))
        y = np.concatenate((yinj, yprod))
        z = np.concatenate((zinj, zprod))
        
        for i in range(numberoflaterals):
            x = np.concatenate((x, xlat[:, i].reshape(-1,1)))
            y = np.concatenate((y, ylat[:, i].reshape(-1,1)))
            z = np.concatenate((z, zlat[:, i].reshape(-1,1)))

    return locals()
    

def set_tube_geometry(clg_configuration, Diameter1, Diameter2, PipeParam3, PipeParam4, PipeParam5):
   
    """ 
        ## Tube Geometry
        Diameter1                        # If coaxial: radius, radiuscenterpipe and If U-loop: radiusvertical and radiuscenterpipe
        Diameter2                        #
    """

    radius = radiuscenterpipe = thicknesscenterpipe = coaxialflowtype = None # coaxial
    numberoflaterals = radiuslateral = lateralflowallocation =  None # uloop
    
    if clg_configuration == 1:  # co-axial geometry (1)
        radius = Diameter1/2  # 0.2286/2                # Wellbore radius [m] (everything is assumed open-hole)
        radiuscenterpipe = Diameter2/2  # 0.127/2       # Inner radius of inner pipe [m]
        thicknesscenterpipe = PipeParam3 # 0.0127       # Thickness of inner pipe [m]
        k_center_pipe = PipeParam4 # 0.006                                           # Thermal conductivity of insulation of center pipe wall [W/m/K]
        coaxialflowtype = PipeParam5 # 1                                             # 1 = CXA (fluid injection in annulus); 2 = CXC (fluid injection in center pipe)

    elif clg_configuration == 2: # U-loop geometry (2)

        radiusvertical = Diameter1/2 # 0.15                         # Radius of "vertical" injection and production well (open hole assumed for heat transfer) [m] (it is labeled here as vertical but it is allowed to be deviated)
        radiuslateral = Diameter2/2 # 0.30                          # Radius of laterals (open hole assumed for heat transfer) [m]
        
        numberoflaterals = PipeParam3 # 3                                            # Number of laterals (must be integer) [-]
        lateralflowallocation = PipeParam4 # [1/3, 1/3, 1/3]                         # Distribution of flow accross laterals, must add up to 1 (it will get normalized below if sum does not equal to 1). Length of array must match number of laterals.
        lateralflowmultiplier = PipeParam5 # 1                                       # Multiplier to allow decreasing the lateral flow rate to account for other laterals not simulated. 

    return locals()

def admin_fluid_properties():

    """
        # NOT editable by the user 
        cp_f:                            # Fluid specific heat capacity [J/kgK]
        rho_f:                           # Fluid density [kg/m3]
        k_f:                             # Fluid heat conductivity [W/m.K]
        mu_f:                            # Fluid dynamic viscosity [Pa.s]
    """

    cp_f=4200
    rho_f=1000
    k_f=0.68
    mu_f=600*10**-6

    g = 9.81                       # Gravitational acceleration [m/s^2]
    gamma = 0.577215665            # Euler's constant
    alpha_f = k_f / rho_f / cp_f   # Fluid thermal diffusivity [m2/s]
    Pr_f = mu_f / rho_f / alpha_f  # Fluid Prandtl number [-]

    return locals()


def compute_tube_geometry(sbt_version, clg_configuration, radiuscenterpipe, thicknesscenterpipe, 
                                xinj, xprod, xlat, numberoflaterals, radiuslateral, lateralflowallocation):

    interconnections = None
    Deltaz = np.sqrt((x[1:] - x[:-1]) ** 2 + (y[1:] - y[:-1]) ** 2 + (z[1:] - z[:-1]) ** 2)  # Length of each segment [m]
    Deltaz = Deltaz.reshape(-1)

    if clg_configuration == 1: # COAXIAL

        outerradiuscenterpipe = radiuscenterpipe+thicknesscenterpipe    # Outer radius of inner pipe [m]
        A_flow_annulus = math.pi*(radius**2-outerradiuscenterpipe**2)   # Flow area of annulus pipe [m^2]
        A_flow_centerpipe = math.pi*radiuscenterpipe**2                 # Flow area of center pipe [m^2]
        Dh_annulus = 2*(radius-outerradiuscenterpipe)                   # Hydraulic diameter of annulus [m]    
        if sbt_version == 2:
            eps_annulus = Dh_annulus*piperoughness                      # Relative roughness annulus [-]
            eps_centerpipe = 2*radiuscenterpipe*piperoughness           # Relative roughness inner pipe [-]
        LoverR = Deltaz / radius                                        # Ratio of pipe segment length to radius along the wellbore [-]
        RelativeLengthChanges = (Deltaz[1:] - Deltaz[:-1]) / Deltaz[:-1]
        
    elif clg_configuration == 2: # U-LOOP

        interconnections = np.concatenate((np.array([len(xinj)],dtype=int), np.array([len(xprod)],dtype=int), (np.ones(numberoflaterals - 1, dtype=int) * len(xlat))))
        interconnections = np.cumsum(interconnections)  # lists the indices of interconnections between inj, prod, and laterals (this will used to take care of the duplicate coordinates of the start and end points of the laterals)
        radiusvector = np.concatenate([np.ones(len(xinj) + len(xprod) - 2) * radiusvertical, np.ones(numberoflaterals * len(xlat) - numberoflaterals) * radiuslateral])  # Stores radius of each element in a vector [m]
        Dvector = radiusvector * 2  # Diameter of each element [m]
        lateralflowallocation_norm = lateralflowallocation / np.sum(lateralflowallocation)  # Ensure the sum equals 1   
        Deltaz = np.delete(Deltaz, interconnections - 1)  # Removes the phantom elements due to duplicate coordinates
        LoverR = Deltaz / radiusvector  # Ratio of pipe segment length to radius along the wellbore [-]
        if numberoflaterals > 1:
            DeltazOrdered = np.concatenate((Deltaz[0:(interconnections[0]-1)], Deltaz[(interconnections[1]-2):(interconnections[2]-3)], Deltaz[(interconnections[0]-1):(interconnections[1]-2)]))
        else:
            DeltazOrdered = np.concatenate((Deltaz[0:interconnections[0] - 1], Deltaz[interconnections[1] - 1:-1], Deltaz[interconnections[0]:interconnections[1] - 2]))
        RelativeLengthChanges = (DeltazOrdered[1:] - DeltazOrdered[:-1]) / DeltazOrdered[:-1]

        # QUALITY CONTROL U-LOOP
        for dd in range(1, numberoflaterals + 1):
            if abs(xinj[-1] - xlat[0][dd - 1]) > 1e-12 or abs(yinj[-1] - ylat[0][dd - 1]) > 1e-12 or abs(zinj[-1] - zlat[0][dd - 1]) > 1e-12:
                print(f'Error: Coordinate mismatch between bottom of injection well and start of lateral #{dd}')
            if abs(xprod[0] - xlat[-1][dd - 1]) > 1e-12 or abs(yprod[0] - ylat[-1][dd - 1]) > 1e-12 or abs(zprod[0] - zlat[-1][dd - 1]) > 1e-12:
                print(f'Error: Coordinate mismatch between bottom of production well and end of lateral #{dd}')
        if len(lateralflowallocation_norm) != numberoflaterals:
            print('Error: Length of array "lateralflowallocation" does not match the number of laterals')

    TotalLength = np.sum(Deltaz)  # Total length of all elements (for informational purposes only) [m]
    smallestLoverR = np.min(LoverR)  # Smallest ratio of pipe segment length to pipe radius. This ratio should be larger than 10. [-]
    
    # QUALITY CONTROL
    if smallestLoverR < 10:
        print('Warning: smallest ratio of segment length over radius is less than 10. Good practice is to keep this ratio larger than 10.')
    if max(abs(RelativeLengthChanges)) > 0.5:
        print('Warning: abrupt change(s) in segment length detected, which may cause numerical instabilities. Good practice is to avoid abrupt length changes to obtain smooth results.')

    return locals()

def calc_tube_min_time_steps(clg_configuration, radius, radiusvector, alpha_m, Deltaz, LimitPointSourceModel, LimitCylinderModelRequired, LimitInfiniteModel):

    timeforpointssource = max(Deltaz)**2 / alpha_m * LimitPointSourceModel  # Calculates minimum time step size when point source model becomes applicable [s]

    if clg_configuration == 1: # co-axial geometry (1)
        timeforlinesource = radius**2 / alpha_m * LimitCylinderModelRequired  # Calculates minimum time step size when line source model becomes applicable [s]

    elif clg_configuration == 2: # U-loop geometry (2)
        timeforlinesource = max(radiusvector)**2 / alpha_m * LimitCylinderModelRequired  # Calculates minimum time step size when line source model becomes applicable [s]

    timeforfinitelinesource = max(Deltaz)**2 / alpha_m * LimitInfiniteModel  # Calculates minimum time step size when finite line source model should be considered [s]

    return timeforpointssource, timeforlinesource, timeforfinitelinesource,


def prepare_interpolators(sbt_version, variablefluidproperties, fluid, rho_f, cp_f, k_f, mu_f):

    interpolator_density = interpolator_enthalpy = interpolator_entropy = None
    interpolator_heatcapacity = interpolator_heatcapacity = interpolator_phase = None 
    interpolator_thermalconductivity = interpolator_thermalexpansion = interpolator_viscosity = None

    if sbt_version == 2:

        if variablefluidproperties == 0:  # For computational purposes, use constant fluid property tables
            # Define vectors for pressure and temperature
            Pvector = [1, 1e9]
            Tvector = [1, 1e4]
        
            # Create 2x2 arrays with constant fluid properties
            density = [[rho_f] * 2] * 2
            heatcapacity = [[cp_f] * 2] * 2
            thermalconductivity = [[k_f] * 2] * 2
            viscosity = [[mu_f] * 2] * 2
            thermalexpansion = [[0] * 2] * 2  # Incompressible fluid has zero thermal expansion coefficient
        else:  # If variable fluid properties, import pre-generated tables
            print('Loading fluid properties...')
            if fluid == 1:  # H2O
                try:
                    mat = scipy.io.loadmat('properties_H2O.mat') 
                    Pvector = mat['Pvector']
                    print("pVector", Pvector)
                    Tvector = mat['Tvector']
                    density = mat['density']
                    enthalpy = mat['enthalpy']
                    entropy = mat['entropy']
                    heatcapacity = mat['heatcapacity']
                    phase = mat['phase']
                    thermalconductivity = mat['thermalconductivity']
                    thermalexpansion = mat['thermalexpansion']
                    viscosity = mat['viscosity']
                    print('Fluid properties for water loaded successfully')
                except Exception as e:
                    print(f"Error loading properties for water: {e}")
                    raise
            elif fluid == 2:  #CO2
                try:
                    mat = scipy.io.loadmat('properties_CO2.mat') 
                    Pvector = mat['Pvector']
                    Tvector = mat['Tvector']
                    density = mat['density']
                    enthalpy = mat['enthalpy']
                    entropy = mat['entropy']
                    heatcapacity = mat['heatcapacity']
                    phase = mat['phase']
                    thermalconductivity = mat['thermalconductivity']
                    thermalexpansion = mat['thermalexpansion']
                    viscosity = mat['viscosity']                
                    print('Fluid properties for CO2 loaded successfully')
                except Exception as e:
                    print(f"Error loading properties for CO2: {e}")
                    raise
            else:
                print('No valid fluid selected')
                exit()
        
        # Prepare interpolators
        Pvector_1d = Pvector.ravel()
        Tvector_1d = Tvector.ravel()
        interpolator_density = RegularGridInterpolator((Pvector_1d, Tvector_1d), density)
        interpolator_enthalpy = RegularGridInterpolator((Pvector_1d, Tvector_1d), enthalpy)
        interpolator_entropy = RegularGridInterpolator((Pvector_1d, Tvector_1d), entropy)
        interpolator_heatcapacity = RegularGridInterpolator((Pvector_1d, Tvector_1d), heatcapacity)
        interpolator_phase = RegularGridInterpolator((Pvector_1d, Tvector_1d), phase)
        interpolator_thermalconductivity = RegularGridInterpolator((Pvector_1d, Tvector_1d), thermalconductivity)
        interpolator_thermalexpansion = RegularGridInterpolator((Pvector_1d, Tvector_1d), thermalexpansion)
        interpolator_viscosity = RegularGridInterpolator((Pvector_1d, Tvector_1d), viscosity)
    
    return interpolator_density, interpolator_enthalpy, \
                interpolator_entropy, interpolator_heatcapacity, interpolator_heatcapacity, \
                    interpolator_phase, interpolator_thermalconductivity, interpolator_thermalexpansion, interpolator_viscosity


def get_profiles(sbt_version, variableinjectiontemperature, variableflowrate, flowratefilename, Tinj, mdot):

    # Read injection temperature profile if provided
    Tinstore = np.zeros(len(times))
    if variableinjectiontemperature == 1 and sbt_version == 1:
        # User has provided injection temperature in an Excel spreadsheet. (can currently only be used with sbt version 1)
        num = pd.read_excel(injectiontemperaturefilename)
        Tintimearray = np.array(num.iloc[:, 0])
        Tintemperaturearray = np.array(num.iloc[:, 1])
        # Quality control
        if len(Tintimearray) < 2:
            print('Error: Provided injection temperature profile should have at least 2 values')
        
        if Tintimearray[0] != 0:
            print('Error: First time value in the user-provided injection temperature profile does not equal 0 s')
        
        if abs(Tintimearray[-1] - times[-1]) > 1e-5:
            print('Error: Last time value in the user-provided injection temperature profile does not equal the final value in the "times" array')
    
        else:
            Tintimearray[-1] = times[-1]  # Ensure final time values "exactly" match to prevent interpolation issues at the final time step
        Tinstore[0] = Tintemperaturearray[0]
    else:
        Tinstore[0] = Tinj

    # Read mass flow rate profile if provided
    mstore = np.zeros(len(times))  # The value for m used at each time step is stored in this array (is either constant or interpolated from a user-provided mass flow rate profile)

    if variableflowrate == 1 and sbt_version == 1:  # User has provided mass flow rate in an Excel spreadsheet. (can currently only be used with sbt version 1)
        data = pd.read_excel(flowratefilename)
        mtimearray = data.iloc[:, 0].values  # This array has the times provided by the user
        mflowratearray = data.iloc[:, 1].values  # This array has the injection temperatures provided by the user

        # Quality control
        if len(mtimearray) < 2:
            print('Error: Provided flow rate profile should have at least 2 values')
        
        if mtimearray[0] != 0:
            print('Error: First time value in user-provided flow rate profile does not equal to 0 s')
            
        if abs(mtimearray[-1] - times[-1]) > 1e-5:
            print('Error: Last time value in user-provided flow rate profile does not equal to final value in "times" array')
        
        else:
            mtimearray[-1] = times[-1]  # Ensure final time values "exactly" match to prevent interpolation issues at the final time step

        mstore[0] = mflowratearray[0]
    else:
        mstore[0] = mdot
    
    return Tinstore, mstore

def precalculations(clg_configuration, Deltaz, alpha_m, k_m, times, NoArgumentsFinitePipeCorrection, NoDiscrFinitePipeCorrection,
                    timeforlinesource, radius, radiusvector, interconnections):

    fpcminarg = min(Deltaz)**2 / (4 * alpha_m * times[-1])
    fpcmaxarg = max(Deltaz)**2 / (4 * alpha_m * (min(times[1:] - times[:-1])))
    Amin1vector = np.logspace(np.log10(fpcminarg) - 0.1, np.log10(fpcmaxarg) + 0.1, NoArgumentsFinitePipeCorrection)
    finitecorrectiony = np.zeros(NoArgumentsFinitePipeCorrection)

    for i, Amin1 in enumerate(Amin1vector):
        Amax1 = (16)**2
        if Amin1 > Amax1:
            Amax1 = 10 * Amin1
        Adomain1 = np.logspace(np.log10(Amin1), np.log10(Amax1), NoDiscrFinitePipeCorrection)
        finitecorrectiony[i] = np.trapz(-1 / (Adomain1 * 4 * np.pi * k_m) * erfc(1/2 * np.power(Adomain1, 1/2)), Adomain1)

    #precalculate besselintegration for infinite cylinder
    if clg_configuration == 1: # co-axial geometry (1)
        besselminarg = alpha_m * (min(times[1:] - times[:-1])) / radius**2
        besselmaxarg = alpha_m * timeforlinesource / radius**2
        
    elif clg_configuration == 2: # U-loop geometry (2)
        besselminarg = alpha_m * (min(times[1:] - times[:-1])) / max(radiusvector)**2
        besselmaxarg = alpha_m * timeforlinesource / min(radiusvector)**2
    
    deltazbessel = np.logspace(-10, 8, NoDiscrInfCylIntegration)
    argumentbesselvec = np.logspace(np.log10(besselminarg) - 0.5, np.log10(besselmaxarg) + 0.5, NoArgumentsInfCylIntegration)
    besselcylinderresult = np.zeros(NoArgumentsInfCylIntegration)

    for i, argumentbessel in enumerate(argumentbesselvec):
        besselcylinderresult[i] = 2 / (k_m * np.pi**3) * np.trapz((1 - np.exp(-deltazbessel**2 * argumentbessel)) / (deltazbessel**3 * (jv(1, deltazbessel)**2 + yv(1, deltazbessel)**2)), deltazbessel)

    N = len(Deltaz)  # Number of elements
    elementcenters = 0.5 * np.column_stack((x[1:], y[1:], z[1:])) + 0.5 * np.column_stack((x[:-1], y[:-1], z[:-1]))  # Matrix that stores the mid point coordinates of each element
    if clg_configuration == 2: # U-loop geometry (2)
        interconnections_new = interconnections - 1
        elementcenters = np.delete(elementcenters, interconnections_new.reshape(-1,1), axis=0)  # Remove duplicate coordinates

    SMatrix = np.zeros((N, N))  # Initializes the spacing matrix, which holds the distance between center points of each element [m]
    SoverL = np.zeros((N, N))  # Initializes the ratio of spacing to element length matrix
    for i in range(N):
        SMatrix[i, :] = np.sqrt((elementcenters[i, 0] - elementcenters[:, 0])**2 + (elementcenters[i, 1] - elementcenters[:, 1])**2 + (elementcenters[i, 2] - elementcenters[:, 2])**2)
        SoverL[i, :] = SMatrix[i, :] / Deltaz[i] #Calculates the ratio of spacing between two elements and element length
    
    return Amin1vector, argumentbesselvec, besselcylinderresult, elementcenters, SMatrix, SoverL, N, interconnections_new, finitecorrectiony





def run_sbt(
            ## Model Specifications 
            sbt_version=1, mesh_fineness=0, HYPERPARAM1=0, HYPERPARAM2="MassFlowRate.xlsx", 
            HYPERPARAM3=0, HYPERPARAM4="InjectionTemperatures.xlsx", HYPERPARAM5=None, 
            accuracy=1,

             ## Operations
            clg_configuration=2, mdot=20, Tinj=20, fluid=1, ## Operations
            DrillingDepth_L1=2, HorizontalExtent_L2=1, #BoreholeDiameter=1, ## Wellbore Geometry
            Diameter1=0.2286, Diameter2=0.2286, PipeParam3=3, PipeParam4=[1/3, 1/3, 1/3], PipeParam5=1, ## Tube Geometry

            ## Geologic Properties
            Tsurf=20, GeoGradient=90/1000, k_m=2.83, c_m=825, rho_m=2875, 
            ):
    """
    Runs the SBT Model
    EDITABLE BY A USER:

        ## Operations
        clg_configuration:               # Must be 1 or 2. "1" mean co-axial, "2" U-loop
        mdot:                            # Total fluid mass flow rate [kg/s]. m must be provided if the user sets variableflowrate to 0.
        Tinj:                            # Constant injection temperature [deg.C]
        
        ## Geologic Properties
        Tsurf:                           # Surface temperature [deg C]
        GeoGradient:                     # Geothermal gradient [C/m]
        k_m:                             # Rock thermal conductivity [W/m.K]
        c_m:                             # Rock specific heat capacity [J/kgK]
        rho_m:                           # Rock density [kg/m3]
    """

    HorizontalExtent_L2 = HorizontalExtent_L2*1000
    DrillingDepth_L1 = DrillingDepth_L1*1000
    # print("\n")
    # print(" -------------------------------- SBT USER INPUTS -------------------------------- ")
    # all input possibilities can be placed into a dataframe at some point ...
    # Geologic properties are set at the start of run_sbt() 
    # rest of variables are "set" below

    tube_geometry_dict = set_tube_geometry(clg_configuration=clg_configuration, 
                                            Diameter1=Diameter1, Diameter2=Diameter2, 
                                            PipeParam3=PipeParam3, PipeParam4=PipeParam4, PipeParam5=PipeParam5
                                            )
    globals().update(tube_geometry_dict)
    # print(tube_geometry_dict.keys())

    wellbore_geometry_dict = set_wellbore_geometry(clg_configuration=clg_configuration, 
                                                    DrillingDepth_L1=DrillingDepth_L1, HorizontalExtent_L2=HorizontalExtent_L2,
                                                    numberoflaterals=numberoflaterals
                                                    )
    globals().update(wellbore_geometry_dict)
    # print(wellbore_geometry_dict.keys())

    sbt_hyperparams_dict = set_sbt_hyperparameters(sbt_version=sbt_version, clg_configuration=clg_configuration, 
                                                accuracy=accuracy,
                                                mesh_fineness=mesh_fineness, fluid=fluid, 
                                                HYPERPARAM1=HYPERPARAM1, HYPERPARAM2=HYPERPARAM2, 
                                                HYPERPARAM3=HYPERPARAM3, HYPERPARAM4=HYPERPARAM4, HYPERPARAM5=HYPERPARAM5)
    globals().update(sbt_hyperparams_dict)
    # print(sbt_hyperparams_dict.keys())

    if is_plot:
        plot_borehole_geometry(clg_configuration=clg_configuration, numberoflaterals=numberoflaterals, 
                                x=x, y=y, z=z, 
                                xinj=xinj, yinj=yinj, zinj=zinj, xprod=xprod, yprod=yprod, zprod=zprod, xlat=xlat, ylat=ylat, zlat=zlat)

    #%% ----------------
    # 2. Pre-Processing
    # Generally, nothing should be changed by the user in this section
    #------------------

    # print("\n")
    # print(" -------------------------------- SBT ADMIN INPUTS -------------------------------- ")
    fluid_properties = admin_fluid_properties()
    globals().update(fluid_properties)
    # print(fluid_properties.keys())

    # print(" -------------------------------- COMPUTATIONS -------------------------------- ")

    ### COMPUTE TUBE GEOMETRY 
    tube_geometry_dict2 = compute_tube_geometry(sbt_version=sbt_version, clg_configuration=clg_configuration, 
                                                radiuscenterpipe=radiuscenterpipe, thicknesscenterpipe=thicknesscenterpipe, 
                                                xinj=xinj, xprod=xprod, xlat=xlat, numberoflaterals=numberoflaterals, 
                                                radiuslateral=radiuslateral, 
                                                lateralflowallocation=lateralflowallocation)
    globals().update(tube_geometry_dict2)
    # print(tube_geometry_dict2.keys())

    ### PREPARE TEMPERATURE AND PRESSURE INTERPOLATORS 
    interpolator_density, interpolator_enthalpy, \
                interpolator_entropy, interpolator_heatcapacity, interpolator_heatcapacity, \
                    interpolator_phase, interpolator_thermalconductivity, interpolator_thermalexpansion, interpolator_viscosity = \
                    prepare_interpolators(sbt_version=sbt_version, variablefluidproperties=variablefluidproperties, 
                                fluid=fluid, rho_f=rho_f, cp_f=cp_f, k_f=k_f, mu_f=mu_f)

    ### GET INJECTION TEMPERATURE as an array AND MASS FLOW RATE PROFILES as an array
    # e.g. [30.  0.  0.  0.  0. ....] or similar
    Tinstore, mstore = get_profiles(sbt_version=sbt_version, variableinjectiontemperature=variableinjectiontemperature, 
                                variableflowrate=variableflowrate, flowratefilename=flowratefilename, 
                                Tinj=Tinj, mdot=mdot)
    
    ### FLUID TEMPERATURE AND PRESSURE COMPUTATIONS
    Tw_down_previous = Tw_up_previous = Tfluiddownnodes = TwMatrix = None
    Pfluiddownnodes = Pfluidupnodes = None

    ### MINIMUM TIME STEPS (model accuracy parameters integrated here)
    alpha_m = k_m / rho_m / c_m    # Thermal diffusivity medium [m2/s]
    timeforpointssource, timeforlinesource, timeforfinitelinesource = calc_tube_min_time_steps(
                                                                            clg_configuration=clg_configuration, radius=radius, radiusvector=radiusvector,
                                                                            alpha_m=alpha_m, Deltaz=Deltaz, 
                                                                            LimitPointSourceModel=LimitPointSourceModel, 
                                                                            LimitCylinderModelRequired=LimitCylinderModelRequired, 
                                                                            LimitInfiniteModel=LimitInfiniteModel
                                                                            )
    ### PRECALCULATIONS (model accuracy parameters integrated here)
    # precalculate the thermal response with a line and cylindrical heat source. Precalculating allows to speed up the SBT algorithm.
    # precalculate finite pipe correction
    Amin1vector, argumentbesselvec, besselcylinderresult, elementcenters, SMatrix, SoverL, N, interconnections_new, finitecorrectiony = precalculations(
                                                                clg_configuration=clg_configuration, 
                                                                Deltaz=Deltaz, alpha_m=alpha_m, k_m=k_m, 
                                                                times=times, 
                                                                NoArgumentsFinitePipeCorrection=NoArgumentsFinitePipeCorrection, 
                                                                NoDiscrFinitePipeCorrection=NoDiscrFinitePipeCorrection,
                                                                timeforlinesource=timeforlinesource, radius=radius, radiusvector=radiusvector,
                                                                interconnections=interconnections
                                                                )

    # Element ranking based on spacinng is required for SBT algorithm as elements in close proximity to each other use different analytical heat transfer models than elements far apart
    # print(Deltaz) # Length of each segment [m]
    SortedIndices = np.argsort(SMatrix, axis=1, kind = 'stable') # Getting the indices of the sorted elements
    SMatrixSorted = np.take_along_axis(SMatrix, SortedIndices, axis=1)  # Sorting the spacing matrix 
    SoverLSorted = SMatrixSorted / Deltaz
    #filename = 'smatrixpython.mat'
    #scipy.io.savemat(filename, dict(SortedIndicesPython=SortedIndices,SMatrixSortedPython=SMatrixSorted))

    mindexNPCP = np.where(np.min(SoverLSorted, axis=0) < LimitSoverL)[0][-1]  # Finding the index where the ratio is less than the limit

    midpointsx = elementcenters[:, 0]  # x-coordinate of center of each element [m]
    midpointsy = elementcenters[:, 1]  # y-coordinate of center of each element [m]
    midpointsz = elementcenters[:, 2]  # z-coordinate of center of each element [m]
    BBinitial = Tsurf - GeoGradient * midpointsz  # Initial temperature at center of each element [degC]
    verticalchange = z[1:]-z[:-1]        #Vertical change between nodes to calculate impact of gravity on pressure [m]
    verticalchange = verticalchange.ravel()

    if clg_configuration == 2: #U-loop geometry
        previouswaterelements = np.zeros(N)
        previouswaterelements[0:] = np.arange(-1,N-1)
        
        for i in range(numberoflaterals):
            previouswaterelements[interconnections_new[i + 1] - i-1] = len(xinj) - 2
        
        previouswaterelements[len(xinj) - 1] = 0
        
        lateralendpoints = []
        for i in range(1,numberoflaterals+1):
            lateralendpoints.append(len(xinj) - 2 + len(xprod) - 1 + i * ((xlat[:, 0]).size- 1))
        lateralendpoints = np.array(lateralendpoints)

    if clg_configuration == 1 and sbt_version == 2: #co-axial sbt v2 calculates nodal fluid temperatures
        Tfluidupnodes =  Tsurf-GeoGradient*(z)     #initial temperature of upflowing fluid at nodes [degC]
        Tfluiddownnodes =  Tsurf-GeoGradient*(z)   #initial temperature of downflowing fluid at nodes [degC]
        Tfluiddownmidpoints = 0.5*Tfluiddownnodes[1:]+0.5*Tfluiddownnodes[:-1] #initial midpoint temperature of downflowing fluid [deg.C]
        Tfluidupmidpoints = 0.5*Tfluidupnodes[1:]+0.5*Tfluidupnodes[:-1] #initial midpoint temperature of upflowing fluid [deg.C]



    MaxSMatrixSorted = np.max(SMatrixSorted, axis=0)

    indicesyoucanneglectupfront = alpha_m * (np.ones((N-1, 1)) * times) / (MaxSMatrixSorted[1:].reshape(-1, 1) * np.ones((1, len(times))))**2 / LimitNPSpacingTime
    indicesyoucanneglectupfront[indicesyoucanneglectupfront > 1] = 1

    lastneighbourtoconsider = np.zeros(len(times))
    for i in range(len(times)):
        lntc = np.where(indicesyoucanneglectupfront[:, i] == 1)[0]
        if len(lntc) == 0:
            lastneighbourtoconsider[i] = 0
        else:
            lastneighbourtoconsider[i] = max(1, lntc[-1])

    distributionx = np.zeros((len(x) - 1, M + 1))
    distributiony = np.zeros((len(x) - 1, M + 1))
    distributionz = np.zeros((len(x) - 1, M + 1))

    for i in range(len(x) - 1):
        distributionx[i, :] = np.linspace(x[i], x[i + 1], M + 1).reshape(-1)
        distributiony[i, :] = np.linspace(y[i], y[i + 1], M + 1).reshape(-1)
        distributionz[i, :] = np.linspace(z[i], z[i + 1], M + 1).reshape(-1)

    if clg_configuration == 2: #U-loop geometry
        # Remove duplicates
        distributionx = np.delete(distributionx, interconnections_new, axis=0)
        distributiony = np.delete(distributiony, interconnections_new, axis=0)
        distributionz = np.delete(distributionz, interconnections_new, axis=0)
        
        

    if sbt_version == 2: #calculate and converge intial fluid properties in sbt v2
        if clg_configuration == 1: #co-axial system
            
            # Calculate initial pressure distribution
            if fluid == 1:  # H2O
                Pfluidupnodes = Pin * 1e5 - 1000 * g * z  # Initial guess for pressure distribution at nodes [Pa]
                Pfluiddownnodes = np.copy(Pfluidupnodes)  # As initial guess, assume downflowing and upflowing water have the same pressure at each depth [Pa]
                Pfluidupmidpoints = Pin * 1e5 - 1000 * g * midpointsz  # Initial guess for pressure distribution at midpoints [Pa]
                Pfluiddownmidpoints = np.copy(Pfluidupmidpoints) #As initial guess, assume downflowing and upflowing water have the same pressure at each depth [Pa]
            elif fluid == 2:  # CO2
                Pfluidupnodes = Pin * 1e5 - 500 * g * z  # Initial guess for pressure distribution at nodes [Pa]
                Pfluiddownnodes = np.copy(Pfluidupnodes) #As initial guess, assume downflowing and upflowing water have the same pressure at each depth [Pa]
                Pfluidupmidpoints = Pin * 1e5 - 500 * g * midpointsz  # Initial guess for pressure distribution at midpoints [Pa]
                Pfluiddownmidpoints = np.copy(Pfluidupmidpoints) #As initial guess, assume downflowing and upflowing water have the same pressure at each depth [Pa]
            
            kk = 1
            maxrelativechange = 1
            print("Calculating initial pressure field ... | Iteration = 1")
            while kk < maxnumberofiterations and maxrelativechange > reltolerance: #Iterate to converge to initial pressure distribution
                # Store old values
                Pfluidupmidpoints_old = np.copy(Pfluidupmidpoints) #Store current guess for upflowing pressure distribution at midpoints to previous guess [Pa]
                Pfluiddownmidpoints_old = np.copy(Pfluiddownmidpoints) #Store current guess for downflowing pressure distribution at midpoints to previous guess [Pa]
            
                # Calculate fluid density
                densityfluidupmidpoints = interpolator_density(np.array([[x, y] for x, y in zip(Pfluidupmidpoints, BBinitial + 273.15)])) #Calculate density distribution of upflowing fluid at midpoints [kg/m3]
                densityfluiddownmidpoints = np.copy(densityfluidupmidpoints) #At time 0 there is no flow yet so upflowing and downflowing fluid have same pressure, temperature and density distribution
            
                # Update pressure distributions
                np.cumsum(g * verticalchange * densityfluiddownmidpoints)
                Pfluiddownnodes = Pin * 1e5 - np.cumsum(np.append([0], g * verticalchange * densityfluiddownmidpoints)) #Calculate pressure distribution of downflowing fluid at nodes [Pa]
                Pfluidupnodes = np.copy(Pfluiddownnodes) #At time 0 there is no flow yet so upflowing and downflowing fluid have same pressure, temperature and density distribution
                Pfluiddownmidpoints = 0.5 * (Pfluiddownnodes[1:] + Pfluiddownnodes[:-1]) #Pressure at midpoints is calculated by interpolating between nodes
                Pfluidupmidpoints = np.copy(Pfluiddownmidpoints) #Upflowing and downflowing fluid have same initial pressure at time 0
            
                # Calculate maximum relative change
                maxrelativechange = np.max(np.abs((Pfluiddownmidpoints_old - Pfluiddownmidpoints) / Pfluiddownmidpoints_old))
            
                # Print iteration status
                print(f"Calculating initial pressure field ... | Iteration = {kk} | Max. Rel. change = {maxrelativechange}")
                kk += 1
            
            # Calculate initial density distribution
            densityfluiddownnodes = interpolator_density(np.array([[x, y] for x, y in zip(Pfluiddownnodes, Tfluiddownnodes + 273.15)])) #After initial pressure distribution converged, calculate initial density distribution [kg/m3]
            densityfluidupnodes = np.copy(densityfluiddownnodes) #Upflowing and downflowing fluid have the same initial density distribution at time 0
            
            if maxrelativechange < reltolerance:
                print("Initial pressure field calculated successfully")
            else:
                print("Initial pressure field calculated but maximum relative tolerance not met")
            
            # Calculate velocity field
            if coaxialflowtype == 1:  # CXA
                velocityfluiddownmidpoints = m / A_flow_annulus / densityfluiddownmidpoints #Downgoing fluid velocity at midpoints in annulus [m/s]
                velocityfluidupmidpoints = m / A_flow_centerpipe / densityfluidupmidpoints #Upgoing fluid velocity at midpoint in center pipe [m/s]
                velocityfluiddownnodes = m / A_flow_annulus / densityfluiddownnodes #Downgoing fluid velocity at nodes in annulus [m/s]
                velocityfluidupnodes = m / A_flow_centerpipe / densityfluidupnodes #Upgoing fluid velocity at nodes in center pipe [m/s]
            elif coaxialflowtype == 2:  # CXC
                velocityfluiddownmidpoints = m / A_flow_centerpipe / densityfluiddownmidpoints #Downgoing fluid velocity at midpoints in center pipe [m/s]
                velocityfluidupmidpoints = m / A_flow_annulus / densityfluidupmidpoints #Upgoing fluid velocity at midpoint in annulus [m/s]
                velocityfluiddownnodes = m / A_flow_centerpipe / densityfluiddownnodes #Downgoing fluid velocity at nodes in center pipe [m/s]
                velocityfluidupnodes = m / A_flow_annulus / densityfluidupnodes #Upgoing fluid velocity at nodes in annulus [m/s]
            
            # Obtain initial viscosity distribution [Pa*s]
            viscosityfluiddownmidpoints = interpolator_viscosity(np.array([[x, y] for x, y in zip(Pfluiddownmidpoints, BBinitial + 273.15)]))
            viscosityfluidupmidpoints = interpolator_viscosity(np.array([[x, y] for x, y in zip(Pfluidupmidpoints, BBinitial + 273.15)]))
            
            # Obtain initial specific heat capacity distribution [J/kg/K]
            heatcapacityfluiddownmidpoints = interpolator_heatcapacity(np.array([[x, y] for x, y in zip(Pfluiddownmidpoints, BBinitial + 273.15)]))
            heatcapacityfluidupmidpoints = interpolator_heatcapacity(np.array([[x, y] for x, y in zip(Pfluidupmidpoints, BBinitial + 273.15)]))
            
            # Obtain initial thermal conductivity distribution [W/m/K]
            thermalconductivityfluiddownmidpoints = interpolator_thermalconductivity(np.array([[x, y] for x, y in zip(Pfluiddownmidpoints, BBinitial + 273.15)]))
            thermalconductivityfluidupmidpoints = interpolator_thermalconductivity(np.array([[x, y] for x, y in zip(Pfluidupmidpoints, BBinitial + 273.15)]))
            
            # Obtain initial thermal diffusivity distribution [m2/s]
            alphafluiddownmidpoints = thermalconductivityfluiddownmidpoints / densityfluiddownmidpoints / heatcapacityfluiddownmidpoints
            alphafluidupmidpoints = thermalconductivityfluidupmidpoints / densityfluidupmidpoints / heatcapacityfluidupmidpoints
            
            # Obtain initial thermal expansion coefficient distribution [1/K]
            thermalexpansionfluiddownmidpoints = interpolator_thermalexpansion(np.array([[x, y] for x, y in zip(Pfluiddownmidpoints, BBinitial + 273.15)]))
            thermalexpansionfluidupmidpoints = interpolator_thermalexpansion(np.array([[x, y] for x, y in zip(Pfluidupmidpoints, BBinitial + 273.15)]))
            
            # Obtain initial Prandtl number distribution
            Prandtlfluiddownmidpoints = viscosityfluiddownmidpoints / densityfluiddownmidpoints / alphafluiddownmidpoints
            Prandtlfluidupmidpoints = viscosityfluidupmidpoints / densityfluidupmidpoints / alphafluidupmidpoints
            
            # Obtain initial Reynolds number distribution
            if coaxialflowtype == 1:  # CXA
                Refluiddownmidpoints = densityfluiddownmidpoints * velocityfluiddownmidpoints * Dh_annulus / viscosityfluiddownmidpoints
                Refluidupmidpoints = densityfluidupmidpoints * velocityfluidupmidpoints * (2 * radiuscenterpipe) / viscosityfluidupmidpoints
            elif coaxialflowtype == 2:  # CXC
                Refluiddownmidpoints = densityfluiddownmidpoints * velocityfluiddownmidpoints * (2 * radiuscenterpipe) / viscosityfluiddownmidpoints
                Refluidupmidpoints = densityfluidupmidpoints * velocityfluidupmidpoints * Dh_annulus / viscosityfluidupmidpoints
        

    # Initialize SBT algorithm linear system of equation matrices
    if clg_configuration == 1: #co-axial geometry
        L = np.zeros((4 * N, 4 * N))                # Will store the "left-hand side" of the system of equations
        R = np.zeros((4 * N, 1))                    # Will store the "right-hand side" of the system of equations
    elif clg_configuration == 2: #U-loop geometry
        L = np.zeros((3 * N, 3 * N))                # Will store the "left-hand side" of the system of equations
        R = np.zeros((3 * N, 1))                    # Will store the "right-hand side" of the system of equations
    Q = np.zeros((N, len(times)))                   # Initializes the heat pulse matrix, i.e., the heat pulse emitted by each element at each time step
    Toutput = np.zeros(len(times))                  # Initializes the production temperatures array
    if sbt_version == 2:
        Poutput = np.zeros(len(times))             # Initializes the production pressure array
    if clg_configuration == 1: #co-axial geometry
        Tw_up_previous = BBinitial                  # At time zero, the initial upflowing fluid temperature corresponds to the initial local rock temperature
        Tw_down_previous = BBinitial                # At time zero, the initial downflowing fluid temperature corresponds to the initial local rock temperature
    elif clg_configuration == 2: #U-loop geometry
        Twprevious = BBinitial                      # At time zero, the initial fluid temperature corresponds to the initial local rock temperature                          
    Toutput[0] = Tsurf                              # At time zero, the outlet temperature is the initial local fluid temperature at the surface, which corresponds to the surface temperature
    if sbt_version == 2:
        Poutput[0] = Pin                            # At time 0, the fluid flow is assumed zero and we assume the outlet pressure equals the inlet pressure
        
    if sbt_version == 1:    
        if clg_configuration == 1: #co-axial geometry
            Tw_up_Matrix = np.zeros((len(times), N))       # Initializes the matrix that holds the upflowing fluid temperature over time
            Tw_up_Matrix[0, :] = Tw_up_previous
            Tw_down_Matrix = np.zeros((len(times), N))     # Initializes the matrix that holds the downflowing fluid temperature over time
            Tw_down_Matrix[0, :] = Tw_down_previous
        elif clg_configuration == 2: #U-loop geometry
            TwMatrix = np.zeros((len(times), N))            # Initializes the matrix that holds the fluid temperature over time
            TwMatrix[0, :] = Twprevious
    elif sbt_version == 2: 
        if clg_configuration == 1: #co-axial geometry
            Tfluidupnodesstore = np.zeros((N + 1, len(times)))      #Initializes the matrix that stores the upflowing fluid temperature at the nodes
            Tfluiddownnodesstore = np.zeros((N + 1, len(times)))    #Initializes the matrix that stores the downflowing fluid temperature at the nodes
            Tfluidupmidpointsstore = np.zeros((N, len(times)))    #Initializes the matrix that stores the upflowing fluid temperature at the midpoints
            Tfluiddownmidpointsstore  =np.zeros((N, len(times)))  #Initializes the matrix that stores the downflowing fluid temperature at the midpoints
            Pfluidupnodesstore = np.zeros((N + 1, len(times)))     #Initializes the matrix that stores the upflowing fluid pressure at the nodes
            Pfluiddownnodesstore = np.zeros((N + 1, len(times)))    #Initializes the matrix that stores the downflowing fluid pressure at the nodes
            Pfluidupmidpointsstore = np.zeros((N, len(times)))   #Initializes the matrix that stores the upflowing fluid pressure at the midpoints
            Pfluiddownmidpointsstore = np.zeros((N, len(times)))  #Initializes the matrix that stores the downflowing fluid pressure at the midpoints
            
            Tfluidupnodesstore[:,0] = Tfluidupnodes.ravel()
            Tfluiddownnodesstore[:,0]  = Tfluiddownnodes.ravel()
            Tfluidupmidpointsstore[:,0]  = BBinitial.ravel()
            Tfluiddownmidpointsstore[:,0]  = BBinitial.ravel()
            Pfluidupnodesstore[:,0]  = Pfluidupnodes.ravel()
            Pfluiddownnodesstore[:,0]  = Pfluiddownnodes.ravel()
            Pfluidupmidpointsstore[:,0]  = Pfluidupmidpoints.ravel()
            Pfluiddownmidpointsstore[:,0]  = Pfluiddownmidpoints.ravel()
        
        
        
    #Initialize FMM arrays
    combinedtimes = np.array([0])
    combinedQ = np.zeros((N, 1))
    combinedtimes2ndlevel = np.array([0])
    combinedtimes3rdlevel = np.array([0])
    timesFMM = 0
    QFMM = 0

    #%% -------------
    # 3. Calculating
    # Generally, nothing should be changed by the user in this section
    #---------------

    for i in range(1, len(times)):
        Deltat = times[i] - times[i - 1]  # Current time step size [s]

        # If the user has provided an injection temperature profile, current value of Tinj is calculated (only allowed in sbt version 1)
        if variableinjectiontemperature == 1 and sbt_version == 1:
            Tinj = np.interp(times[i], Tintimearray, Tintemperaturearray)
        Tinstore[i] = Tinj  # Value that is used for Tinj at each time step gets stored for postprocessing purposes

        # If the user has provided a flow rate profile, current value of m is calculated (only allowed in sbt version 1)
        if variableflowrate == 1 and sbt_version == 1:
            m = np.interp(times[i], mtimearray, mflowratearray)
        mstore[i] = mdot  # Value that is used for m at each time step gets stored for postprocessing purposes


        # ------------------------------------
        # CPCP (= Current pipe, current pulse)
        # ------------------------------------
        if clg_configuration == 1: #co-axial geometry     
        
            if alpha_m * Deltat / radius**2 > LimitCylinderModelRequired:
                CPCP = np.ones(N) * 1 / (4 * np.pi * k_m) * exp1(radius**2 / (4 * alpha_m * Deltat))  # Use line source model if possible
            else:
                CPCP = np.ones(N) * np.interp(alpha_m * Deltat / radius**2, argumentbesselvec, besselcylinderresult)  # Use cylindrical source model if required
            
            if Deltat > timeforfinitelinesource:  # For long time steps, the finite length correction should be applied
                CPCP = CPCP + np.interp(Deltaz**2 / (4 * alpha_m * Deltat), Amin1vector, finitecorrectiony)
        
        elif clg_configuration == 2: #U-loop geometry  
            if alpha_m * Deltat / max(radiusvector)**2 > LimitCylinderModelRequired:
                CPCP = np.ones(N) * 1 / (4 * np.pi * k_m) * exp1(radiusvector**2 / (4 * alpha_m * Deltat))  # Use line source model if possible
            else:
                CPCP = np.ones(N) * np.interp(alpha_m * Deltat / radiusvector**2, argumentbesselvec, besselcylinderresult)  # Use cylindrical source model if required
            
            if Deltat > timeforfinitelinesource:  # For long time steps, the finite length correction should be applied
                CPCP = CPCP + np.interp(Deltaz**2 / (4 * alpha_m * Deltat), Amin1vector, finitecorrectiony)


        # ---------------------------------
        # CPOP (= Current pipe, old pulses)
        # ---------------------------------
        if i > 1:  # After the second time step, we need to keep track of previous heat pulses

            CPOP = np.zeros((N, len(timesFMM)-1))
            
            indexpsstart = 0

            indexpsend = np.where(timeforpointssource < (times[i] - timesFMM[1:]))[-1]
            if indexpsend.size > 0:
                indexpsend = indexpsend[-1] + 1
            else:
                indexpsend = indexpsstart - 1
            if indexpsend >= indexpsstart:  # Use point source model if allowed
            
                CPOP[:, 0:indexpsend] = Deltaz * np.ones((N, indexpsend)) / (4 * np.pi * np.sqrt(alpha_m * np.pi) * k_m) * (
                        np.ones(N) * (1 / np.sqrt(times[i] - timesFMM[indexpsstart + 1:indexpsend + 2]) -
                        1 / np.sqrt(times[i] - timesFMM[indexpsstart:indexpsend+1])))
            indexlsstart = indexpsend + 1
            indexlsend = np.where(timeforlinesource < (times[i] - timesFMM[1:]))[0]
            if indexlsend.size == 0:
                indexlsend = indexlsstart - 1
            else:
                indexlsend = indexlsend[-1]
            
            if indexlsend >= indexlsstart:  # Use line source model for more recent heat pulse events
                if clg_configuration == 1: #co-axial geometry 
                    CPOP[:, indexlsstart:indexlsend+1] = np.ones((N,1)) * 1 / (4*np.pi*k_m) * (exp1(radius**2\
                            / (4*alpha_m*(times[i]-timesFMM[indexlsstart:indexlsend+1])).reshape(1,len(4 * alpha_m * (times[i] - timesFMM[indexlsstart:indexlsend+1]))))-\
                        exp1((radius**2)\
                            / (4 * alpha_m * (times[i]-timesFMM[indexlsstart+1:indexlsend+2])).reshape(1,len(4 * alpha_m * (times[i] - timesFMM[indexlsstart+1:indexlsend+2])))))

                elif clg_configuration == 2: #U-loop geometry 
                    CPOP[:, indexlsstart:indexlsend+1] = np.ones((N,1)) * 1 / (4*np.pi*k_m) * (exp1((radiusvector**2).reshape(len(radiusvector ** 2),1)\
                            / (4*alpha_m*(times[i]-timesFMM[indexlsstart:indexlsend+1])).reshape(1,len(4 * alpha_m * (times[i] - timesFMM[indexlsstart:indexlsend+1]))))-\
                        exp1((radiusvector**2).reshape(len(radiusvector ** 2),1)\
                            / (4 * alpha_m * (times[i]-timesFMM[indexlsstart+1:indexlsend+2])).reshape(1,len(4 * alpha_m * (times[i] - timesFMM[indexlsstart+1:indexlsend+2])))))

            indexcsstart = max(indexpsend, indexlsend) + 1
            #indexcsend = i - 2
            indexcsend = len(timesFMM)-2

            if indexcsstart <= indexcsend:  # Use cylindrical source model for the most recent heat pulses
                CPOPPH = np.zeros((CPOP[:, indexcsstart:indexcsend+1].shape))   
                CPOPdim = CPOP[:, indexcsstart:indexcsend+1].shape
                CPOPPH = CPOPPH.T.ravel()
                if clg_configuration == 1: #co-axial geometry 
                    CPOPPH = (np.ones(N) * ( \
                                np.interp(alpha_m * (times[i] - timesFMM[indexcsstart:indexcsend+1]).reshape(len(times[i] - timesFMM[indexcsstart:indexcsend+1]),1) / (radius**2), argumentbesselvec, besselcylinderresult) - \
                                np.interp(alpha_m * (times[i] - timesFMM[indexcsstart+1: indexcsend+2]).reshape(len(times[i] - timesFMM[indexcsstart+1:indexcsend+2]),1) / (radius**2), argumentbesselvec, besselcylinderresult))).reshape(-1,1)
                
                elif clg_configuration == 2: #U-loop geometry 
                    CPOPPH = (np.ones(N) * ( \
                                np.interp(alpha_m * (times[i] - timesFMM[indexcsstart:indexcsend+1]).reshape(len(times[i] - timesFMM[indexcsstart:indexcsend+1]),1) / (radiusvector ** 2).reshape(len(radiusvector ** 2),1).T, argumentbesselvec, besselcylinderresult) - \
                                np.interp(alpha_m * (times[i] - timesFMM[indexcsstart+1: indexcsend+2]).reshape(len(times[i] - timesFMM[indexcsstart+1:indexcsend+2]),1) / (radiusvector ** 2).reshape(len(radiusvector ** 2),1).T, argumentbesselvec, besselcylinderresult))).reshape(-1,1)
                CPOPPH=CPOPPH.reshape((CPOPdim),order='F')
                CPOP[:, indexcsstart:indexcsend+1] = CPOPPH
    
            indexflsstart = indexpsend + 1
            indexflsend = np.where(timeforfinitelinesource < (times[i] - timesFMM[1:]))[-1]
            if indexflsend.size == 0:
                indexflsend = indexflsstart - 1
            else:
                indexflsend = indexflsend[-1] - 1
        
            if indexflsend >= indexflsstart:  # Perform finite length correction if needed
                CPOP[:, indexflsstart:indexflsend+2] = CPOP[:, indexflsstart:indexflsend+2] + (np.interp(np.matmul((Deltaz.reshape(len(Deltaz),1) ** 2),np.ones((1,indexflsend-indexflsstart+2))) / np.matmul(np.ones((N,1)),(4 * alpha_m * (times[i] - timesFMM[indexflsstart:indexflsend+2]).reshape(len(times[i] - timesFMM[indexflsstart:indexflsend+2]),1)).T), Amin1vector, finitecorrectiony) - \
                np.interp(np.matmul((Deltaz.reshape(len(Deltaz),1) ** 2),np.ones((1,indexflsend-indexflsstart+2))) / np.matmul(np.ones((N,1)),(4 * alpha_m * (times[i] - timesFMM[indexflsstart+1:indexflsend+3]).reshape(len(times[i] - timesFMM[indexflsstart:indexflsend+2]),1)).T), Amin1vector, finitecorrectiony))

        # ------------------------------------------
        # NPCP (= Neighbouring pipes, current pulse)
        # ------------------------------------------
        NPCP = np.zeros((N, N))
        np.fill_diagonal(NPCP, CPCP)
        

        spacingtest = alpha_m * Deltat / SMatrixSorted[:, 1:]**2 / LimitNPSpacingTime
        maxspacingtest = np.max(spacingtest,axis=0)

        
        if maxspacingtest[0] < 1:
            maxindextoconsider = 0  #KB 11/21/2024: maybe this should be -1 at no index should be consdiered?????
        else:
            maxindextoconsider = np.where(maxspacingtest > 1)[0][-1]+1 #KB 11/21/2024: is plus 1 needed here?????

        #Calculate and store neighbouring pipes for current pulse as point sources        
        if mindexNPCP < maxindextoconsider:      #KB 11/21/2024: is plus 1 needed here?????
            indicestocalculate = SortedIndices[:, mindexNPCP + 1:maxindextoconsider + 1]
            indicestocalculatetranspose = indicestocalculate.T
            indicestocalculatelinear = indicestocalculate.ravel()
            indicestostorematrix = (indicestocalculate - 1) * N + np.arange(1, N) * np.ones((1, maxindextoconsider - mindexNPCP + 1))
            indicestostorematrixtranspose = indicestostorematrix.T
            indicestostorelinear = indicestostorematrix.ravel()
            NPCP[indicestostorelinear] = Deltaz[indicestocalculatelinear] / (4 * np.pi * k_m * SMatrix[indicestostorelinear]) * erfc(SMatrix[indicestostorelinear] / np.sqrt(4 * alpha_m * Deltat))
        #Calculate and store neighbouring pipes for current pulse as set of line sources
        if mindexNPCP > 1 and maxindextoconsider > 0:
            lastindexfls = min(mindexNPCP, maxindextoconsider + 1)
            indicestocalculate = SortedIndices[:, 1:lastindexfls]
            indicestocalculatetranspose = indicestocalculate.T
            indicestocalculatelinear = indicestocalculate.ravel()
            indicestostorematrix = (indicestocalculate) * N + np.arange(N).reshape(-1,1) * np.ones((1, lastindexfls - 1), dtype=int)
        # #pdb.set_trace()
            indicestostorematrixtranspose = indicestostorematrix.T
            indicestostorelinear = indicestostorematrix.ravel()
            midpointindices = np.matmul(np.ones((lastindexfls - 1, 1)), np.arange( N ).reshape(1,N)).T
            midpointsindices = midpointindices.ravel().astype(int)
            rultimate = np.sqrt(np.square((midpointsx[midpointsindices].reshape(len(midpointsindices),1)*( np.ones((1, M + 1))) - distributionx[indicestocalculatelinear,:])) +
                                np.square((midpointsy[midpointsindices].reshape(len(midpointsindices),1)*( np.ones((1, M + 1))) - distributiony[indicestocalculatelinear,:])) +
                                np.square((midpointsz[midpointsindices].reshape(len(midpointsindices),1)*( np.ones((1, M + 1))) - distributionz[indicestocalculatelinear,:])))
            
            NPCP[np.unravel_index(indicestostorelinear, NPCP.shape, 'F')] =  Deltaz[indicestocalculatelinear] / M * np.sum((1 - erf(rultimate / np.sqrt(4 * alpha_m * Deltat))) / (4 * np.pi * k_m * rultimate) * np.matmul(np.ones((N*(lastindexfls-1),1)),np.concatenate((np.array([1/2]), np.ones(M-1), np.array([1/2]))).reshape(-1,1).T), axis=1)
            
        # ------------------------------------------
        # NPOP (= Neighbouring pipes, old pulses)
        # ------------------------------------------ 
        if i>1:
            NPOP = np.zeros((N * (int(lastneighbourtoconsider[i])+1), len(timesFMM) - 1))
        BB = np.zeros((N, 1))
        if i > 1 and lastneighbourtoconsider[i] > 0:
            SMatrixRelevant = SMatrixSorted[:, 1 : int(lastneighbourtoconsider[i] + 2)]
            SoverLRelevant = SoverLSorted[:, 1 : int(lastneighbourtoconsider[i]) + 2]
            SortedIndicesRelevant = SortedIndices[:, 1 : int(lastneighbourtoconsider[i]) + 2] 
            #maxtimeindexmatrix = alpha_m * np.ones((N * int(lastneighbourtoconsider[i]), 1)) * (times[i] - times[1:i]) / (SMatrixRelevant.ravel().reshape(-1,1) * np.ones((1,i-1)))**2
            maxtimeindexmatrix = alpha_m * np.ones((N * (int(lastneighbourtoconsider[i])+1), 1)) * (times[i] - timesFMM[1:]) / (SMatrixRelevant.T.ravel().reshape(-1,1) * np.ones((1,len(timesFMM)-1)))**2
            
            #allindices = np.arange(N * int(lastneighbourtoconsider[i]) * (i - 1))
            allindices = np.arange(N * (int(lastneighbourtoconsider[i])+1) * (len(timesFMM) - 1))
            #if (i>=154):
            #   
            pipeheatcomesfrom = np.matmul(SortedIndicesRelevant.T.ravel().reshape(len(SortedIndicesRelevant.ravel()),1), np.ones((1,len(timesFMM) - 1)))
            pipeheatgoesto = np.arange(N).reshape(N,1) * np.ones((1, int(lastneighbourtoconsider[i]+1)))
            pipeheatgoesto = pipeheatgoesto.transpose().ravel().reshape(len(pipeheatgoesto.ravel()),1) * np.ones((1, len(timesFMM) - 1))
            # Delete everything smaller than LimitNPSpacingTime
            # 
            indicestoneglect = np.where((maxtimeindexmatrix.transpose()).ravel() < LimitNPSpacingTime)[0]
            
            maxtimeindexmatrix = np.delete(maxtimeindexmatrix.transpose().ravel(), indicestoneglect)
            allindices = np.delete(allindices.transpose().ravel(), indicestoneglect)
            indicesFoSlargerthan = np.where(maxtimeindexmatrix.ravel() > 10)[0]
        #  
            indicestotakeforpsFoS = allindices[indicesFoSlargerthan]
    
            allindices2 = allindices.copy()

            allindices2[indicesFoSlargerthan] = []
            SoverLinearized = SoverLRelevant.transpose().ravel().reshape(len(SoverLRelevant.ravel()),1) * np.ones((1, len(timesFMM) - 1))
            indicestotakeforpsSoverL = np.where(SoverLinearized.transpose().ravel()[allindices2] > LimitSoverL)[0]
            overallindicestotakeforpsSoverL = allindices2[indicestotakeforpsSoverL]
            remainingindices = allindices2.copy() 
        
            remainingindices=np.delete(remainingindices,indicestotakeforpsSoverL)
        
            NPOP = np.zeros((N * (int(lastneighbourtoconsider[i])+1), len(timesFMM) - 1))
            
            # Use point source model when FoS is very large
            if len(indicestotakeforpsFoS) > 0:
                #deltatlinear1 = np.ones(N * int(lastneighbourtoconsider[i]), 1) * (times[i] - times[1:i-1])
                deltatlinear1 = np.ones((N * int(lastneighbourtoconsider[i]), 1)) * (times[i] - timesFMM[1:])
                deltatlinear1 = deltatlinear1.ravel()[indicestotakeforpsFoS]
                deltatlinear2 = np.ones((N * int(lastneighbourtoconsider[i]), 1)) * (times[i] - timesFMM[0:-1])
                deltatlinear2 = deltatlinear2[indicestotakeforpsFoS]
                deltazlinear = pipeheatcomesfrom[indicestotakeforpsFoS]
                SMatrixlinear = SMatrixRelevant.flatten(order='F')
                NPOPFoS = Deltaz[deltazlinear] / (4 * np.pi * k_m * SMatrixlinear[indicestotakeforpsFoS]) * (erfc(SMatrixlinear[indicestotakeforpsFoS] / np.sqrt(4 * alpha_m * deltatlinear2)) -
                    erfc(SMatrixlinear[indicestotakeforpsFoS] / np.sqrt(4 * alpha_m * deltatlinear1)))
            
                NPOP[indicestotakeforpsFoS] = NPOPFoS
            
            # Use point source model when SoverL is very large
            if len(overallindicestotakeforpsSoverL) > 0:
                deltatlinear1 = np.ones((N * (int(lastneighbourtoconsider[i])+1), 1)) * (times[i] - timesFMM[1:]).ravel()
                deltatlinear1 = deltatlinear1[overallindicestotakeforpsSoverL]
                deltatlinear2 = np.ones((N * int(lastneighbourtoconsider[i]), 1)) * (times[i] - timesFMM[0:-1]).ravel()
                deltatlinear2 = deltatlinear2[overallindicestotakeforpsSoverL]
                deltazlinear = pipeheatcomesfrom[overallindicestotakeforpsSoverL]
                SMatrixlinear = SMatrixRelevant.flatten(order='F')
                NPOPSoverL = Deltaz[deltazlinear] / (4 * np.pi * k_m * SMatrixlinear[overallindicestotakeforpsSoverL]) * (erfc(SMatrixlinear[overallindicestotakeforpsSoverL] / np.srt(4 * alpha_m * deltatlinear2)) -
                    erfc(SMatrixlinear[overallindicestotakeforpsSoverL] / np.sqrt(4 * alpha_m * deltatlinear1)))
            
                NPOP[overallindicestotakeforpsSoverL] = NPOPSoverL
            
            # Use finite line source model for remaining pipe segments
            if len(remainingindices) > 0:
            
                deltatlinear1 = np.ones((N * (int(lastneighbourtoconsider[i])+1), 1)) * (times[i] - timesFMM[1:])
                deltatlinear1 = (deltatlinear1.transpose()).ravel()[remainingindices]
                deltatlinear2 = np.ones((N * (int(lastneighbourtoconsider[i])+1), 1)) * (times[i] - timesFMM[0:-1])
                deltatlinear2 = (deltatlinear2.transpose()).ravel()[remainingindices]
                deltazlinear = (pipeheatcomesfrom.T).ravel()[remainingindices]
                midpointstuff = (pipeheatgoesto.transpose()).ravel()[remainingindices]
                rultimate = np.sqrt(np.square((midpointsx[midpointstuff.astype(int)].reshape(len(midpointsx[midpointstuff.astype(int)]),1)*( np.ones((1, M + 1))) - distributionx[deltazlinear.astype(int),:])) +
                                np.square((midpointsy[midpointstuff.astype(int)].reshape(len(midpointsy[midpointstuff.astype(int)]),1)*( np.ones((1, M + 1))) - distributiony[deltazlinear.astype(int),:])) +
                                np.square((midpointsz[midpointstuff.astype(int)].reshape(len(midpointsz[midpointstuff.astype(int)]),1)*( np.ones((1, M + 1))) - distributionz[deltazlinear.astype(int),:])))
            # #pdb.set_trace()
                NPOPfls = Deltaz[deltazlinear.astype(int)].reshape(len(Deltaz[deltazlinear.astype(int)]),1).T / M * np.sum((-erf(rultimate / np.sqrt(4 * alpha_m * np.ravel(deltatlinear2).reshape(len(np.ravel(deltatlinear2)),1)*np.ones((1, M + 1)))) + erf(rultimate / np.sqrt(4 * alpha_m * np.ravel(deltatlinear1).reshape(len(np.ravel(deltatlinear1)),1)*np.ones((1, M + 1))))) / (4 * np.pi * k_m * rultimate) *  np.matmul((np.ones((len(remainingindices),1))),(np.concatenate((np.array([1/2]),np.ones(M - 1),np.array([1/2])))).reshape(-1,1).T), axis=1)
                NPOPfls = NPOPfls.T
                dimensions = NPOP.shape
            #  #pdb.set_trace()
                NPOP=NPOP.T.ravel()
                NPOP[remainingindices.reshape((len(remainingindices),1))] = NPOPfls
                NPOP = NPOP.reshape((dimensions[1],dimensions[0])).T
                
                
                
            ##pdb.set_trace()
        # Put everything together and calculate BB (= impact of all previous heat pulses from old neighbouring elements on current element at current time)
        #  
            Qindicestotake = SortedIndicesRelevant.T.ravel().reshape((N * (int(lastneighbourtoconsider[i])+1), 1))*np.ones((1,len(timesFMM)-1)) + \
                            np.ones((N * (int(lastneighbourtoconsider[i])+1), 1)) * N * np.arange(len(timesFMM) - 1)
            Qindicestotake = Qindicestotake.astype(int)
            Qlinear = QFMM.T.ravel()[Qindicestotake]
        # Qlinear = Qlinear[:,:,0]
            BBPS = NPOP * Qlinear
            BBPS = np.sum(BBPS, axis=1)
            BBPSindicestotake = np.arange(N).reshape((N, 1)) + N * np.arange(int(lastneighbourtoconsider[i])+1).reshape((1, int(lastneighbourtoconsider[i])+1))
            BBPSMatrix = BBPS[BBPSindicestotake]
            BB = np.sum(BBPSMatrix, axis=1)

        if i > 1:
            #BBCPOP = np.sum(CPOP * Q[:, 1:i], axis=1)
            BBCPOP = np.sum(CPOP * QFMM[:,1:], axis=1)
        else:
            BBCPOP = np.zeros(N)

        if sbt_version == 1: #constnat fluid properties and no convergence needed
            # Velocities and thermal resistances are calculated each time step as the flow rate is allowed to vary each time step
            if clg_configuration == 1: #co-axial geometry 
                if coaxialflowtype == 1: #CXA
                    u_down = mdot/rho_f/A_flow_annulus               # Downgoing fluid velocity in annulus [m/s]
                    u_up = mdot/rho_f/A_flow_centerpipe              # Upgoing fluid velocity in center pipe [m/s]
                elif coaxialflowtype == 2: #CXC
                    u_up = mdot/rho_f/A_flow_annulus                 # Upgoing fluid velocity in annulus [m/s]
                    u_down = mdot/rho_f/A_flow_centerpipe            # Downgoing fluid velocity in center pipe [m/s]
        
        
                if coaxialflowtype == 1: #CXA (injection in annulus; production from center pipe)
                    #Thermal resistance in annulus (downflowing)
                    Re_down = rho_f*u_down*(Dh_annulus)/mu_f      # Reynolds Number [-]
                    Nuturb_down = 0.023*Re_down**(4/5)*Pr_f**(0.4)  # Turbulent Flow Nusselt Number for Annulus [-] (Dittus-Boelter equation for heating)
                    if (Re_down>2300):                              #Based on Section 8.6 in Bergman (2011), for annulus turbulent flow, the Nusselt numbers for the inner and outer wall can be assumed the same
                        Nu_down_o = Nuturb_down
                        Nu_down_i = Nuturb_down
                    else:
                        Nu_down_o = 5                             # Laminar flow annulus Nusselt number for outer wall (approximate; see Table 8.2 and 8.3 in Bergman (2011)
                        Nu_down_i = 6                             # Laminar flow annulus Nusselt number for inner wall (approximate; see Table 8.2 and 8.3 in Bergman (2011)
                    
                    if mdot < 0.1: # We assume at very low flow rates, we are actually simulating well shut-in. The Nusselt numbers get set to 1 to represent thermal conduction.
                        Nu_down_o = 1
                        Nu_down_i = 1
                    
                    h_down_o = Nu_down_o*k_f/Dh_annulus
                    h_down_i = Nu_down_i*k_f/Dh_annulus
                    Rt = 1/(np.pi*h_down_o*radius*2)                 # Thermal resistance between annulus flow and surrounding rock (open-hole assumed)
                    
                    #Thermal resistance in center pipe (upflowing)
                    Re_up = rho_f*u_up*(2*radiuscenterpipe)/mu_f  # Reynolds Number [-]
                    Nuturb_up = 0.023*Re_up**(4/5)*Pr_f**(0.4)    # Turbulent Flow Nusselt Number [-] (Dittus-Boelter equation for heating)
                    Nulam_up = 4                                  # Laminar Flow Nusselt Number [-] (this is approximate)
                    if Re_up>2300:
                        Nu_up = Nuturb_up
                    else:
                        Nu_up = Nulam_up
                    
                    
                    if mdot < 0.1: # We assume at very low flow rates, we are actually simulating well shut-in. The Nusselt number is set to 1 to represent thermal conduction.
                        Nu_up = 1
                    
                    h_up = Nu_up*k_f/(2*radiuscenterpipe)
                    R_cp = 1/(np.pi*h_up*2*radiuscenterpipe)+np.log((2*outerradiuscenterpipe)/(2*radiuscenterpipe))/(2*np.pi*k_center_pipe) + 1/(np.pi*h_down_i*2*outerradiuscenterpipe) # thermal resistance between annulus flow and center pipe flow
                    
                elif coaxialflowtype == 2: #CXC (injection in center pipe; production from annulus)
                    #Thermal resistance in annulus (upflowing)
                    Re_up = rho_f*u_up*(Dh_annulus)/mu_f          # Reynolds Number [-]
                    Nuturb_up = 0.023*Re_up**(4/5)*Pr_f**(0.4)    # Turbulent Flow Nusselt Number [-] (Dittus-Boelter equation for heating)
                    if Re_up>2300:                                # Based on Section 8.6 in Bergman (2011), for annulus turbulent flow, the Nusselt numbers for the inner and outer wall can be assumed the same
                        Nu_up_o = Nuturb_up
                        Nu_up_i = Nuturb_up
                    else:
                        Nu_up_o = 5                               # Laminar flow annulus Nusselt number for outer wall (approximate; see Table 8.2 and 8.3 in Bergman (2011)
                        Nu_up_i = 6                               # Laminar flow annulus Nusselt number for inner wall (approximate; see Table 8.2 and 8.3 in Bergman (2011)
                    
                    if mdot < 0.1: # We assume at very low flow rates, we are actually simulating well shut-in. The Nusselt numbers get set to 1 to represent thermal conduction.
                        Nu_up_o = 1
                        Nu_up_i = 1
                    
                    h_up_o = Nu_up_o*k_f/Dh_annulus
                    h_up_i = Nu_up_i*k_f/Dh_annulus
                    Rt = 1/(np.pi*h_up_o*radius*2)                   # Thermal resistance between annulus flow and surrounding rock (open-hole assumed)
                    
                    # Thermal resistance in center pipe (downflowing)
                    Re_down = rho_f*u_down*(2*radiuscenterpipe)/mu_f # Reynolds Number [-]
                    Nuturb_down = 0.023*Re_down**(4/5)*Pr_f**(0.4)   # Turbulent Flow Nusselt Number [-] (Dittus-Boelter equation for heating)
                    Nulam_down = 4                                   # Laminar Flow Nusselt Number [-] (this is approximate)
                    if Re_down>2300:
                        Nu_down = Nuturb_down
                    else:
                        Nu_down = Nulam_down
                    
                    if m < 0.1: # We assume at very low flow rates, we are actually simulating well shut-in. The Nusselt number is to 1 to represent thermal conduction.
                        Nu_down = 1
                    
                    h_down = Nu_down*k_f/(2*radiuscenterpipe)
                    R_cp = 1/(np.pi*h_down*2*radiuscenterpipe)+np.log((2*outerradiuscenterpipe)/(2*radiuscenterpipe))/(2*np.pi*k_center_pipe) + 1/(np.pi*h_up_i*2*outerradiuscenterpipe) # Thermal resistance between annulus flow and center pipe flow
        
            
            elif clg_configuration == 2: #U-loop geometry  
        
                uvertical = mdot / rho_f / (np.pi * radiusvertical ** 2)  # Fluid velocity in vertical injector and producer [m/s]
                ulateral = mdot / rho_f / (np.pi * radiuslateral ** 2) * lateralflowallocation_norm * lateralflowmultiplier  # Fluid velocity in each lateral [m/s]
                uvector = np.hstack((uvertical * np.ones(len(xinj) + len(xprod) - 2)))
            
                for dd in range(numberoflaterals):
                    uvector = np.hstack((uvector, ulateral[dd] * np.ones(len(xlat[:, 0]) - 1)))
            
                if mdot > 0.1:
                    Revertical = rho_f * uvertical * (2 * radiusvertical) / mu_f  # Fluid Reynolds number in injector and producer [-]
                    Nuvertical = 0.023 * Revertical ** (4 / 5) * Pr_f ** 0.4  # Nusselt Number in injector and producer (we assume turbulent flow) [-]
                else:
                    Nuvertical = 1  # At low flow rates, we assume we are simulating the condition of well shut-in and set the Nusselt number to 1 (i.e., conduction only) [-]
                    print(f'Vertical flow shut-in assumed')
            
                hvertical = Nuvertical * k_f / (2 * radiusvertical)  # Heat transfer coefficient in injector and producer [W/m2/K]
                Rtvertical = 1 / (np.pi * hvertical * 2 * radiusvertical)  # Thermal resistance in injector and producer (open-hole assumed)
            
                if mdot > 0.1:
                    Relateral = rho_f * ulateral * (2 * radiuslateral) / mu_f  # Fluid Reynolds number in lateral [-]
                    Nulateral = 0.023 * Relateral ** (4 / 5) * Pr_f ** 0.4  # Nusselt Number in lateral (we assume turbulent flow) [-]
                else:
                    Nulateral = np.ones(numberoflaterals)  # At low flow rates, we assume we are simulating the condition of well shut-in and set the Nusselt number to 1 (i.e., conduction only) [-]
            
                hlateral = Nulateral * k_f / (2 * radiuslateral)  # Heat transfer coefficient in lateral [W/m2/K]
                Rtlateral = 1 / (np.pi * hlateral * 2 * radiuslateral)  # Thermal resistance in lateral (open-hole assumed)
            
            
                Rtvector = Rtvertical * np.ones(len(radiusvector))  # Store thermal resistance of each element in a vector
                
                for dd in range(1, numberoflaterals + 1):
                    if dd < numberoflaterals:
                        Rtvector[interconnections_new[dd] - dd : interconnections_new[dd + 1] - dd] = Rtlateral[dd - 1] * np.ones(len(xlat[:, 0]))
                    else:
                        Rtvector[interconnections_new[dd] - numberoflaterals:] = Rtlateral[dd - 1] * np.ones(len(xlat[:, 0]) - 1)
        
        
        
            
            if clg_configuration == 1: #co-axial geometry 
                if coaxialflowtype == 1: #CXA    
                    #Populate L and R for downflowing fluid heat balance for first element (which has the injection temperature specified)
                    L[0,0] = 1 / Deltat + u_down / Deltaz[0]*2  + 1/R_cp/(A_flow_annulus*rho_f*cp_f)
                    L[0,2] = -1 / (A_flow_annulus*rho_f*cp_f)
                    L[0,3] = -1 / R_cp / (A_flow_annulus*rho_f*cp_f)
                    R[0,0] = 1 / Deltat*Tw_down_previous[0] + u_down/Deltaz[0]*Tinj*2
                    
                    #Populate L and R for rock temperature equation for first element
                    L[1,0] = 1
                    L[1,1] = -1
                    L[1,2] = Rt
                    R[1,0] = 0
                    
                    #Populate L and R for SBT algorithm for first element
                    L[2,np.arange(2,4*N,4)] = NPCP[0,0:N]
                    L[2,1] = 1
                    R[2,0] =  - BBCPOP[0] - BB[0] + BBinitial[0]
                    
                    #Populate L and R for upflowing fluid heat balance for first element
                    L[3,3] = 1 / Deltat + u_up/Deltaz[0] + 1/R_cp/(A_flow_centerpipe*rho_f*cp_f);
                    L[3,0] = -1/R_cp/(A_flow_centerpipe*rho_f*cp_f);
                    L[3,7] = -u_up/Deltaz[0];
                    R[3,0] = 1/Deltat*Tw_up_previous[0];            
                    
                    for iiii in range(2, N+1):  #Populate L and R for remaining elements
                        #Heat balance equation for downflowing fluid
                        L[(iiii-1)*4,(iiii-1)*4] = 1/Deltat + u_down/Deltaz[iiii-1] + 1/R_cp/(A_flow_annulus*rho_f*cp_f)
                        L[(iiii-1)*4,2+(iiii-1)*4] = -1/(A_flow_annulus*rho_f*cp_f)
                        L[(iiii-1)*4,3+(iiii-1)*4] = -1/R_cp/(A_flow_annulus*rho_f*cp_f)
                        L[(iiii-1)*4,(iiii-2)*4] = -u_down/Deltaz[iiii-1]
                        R[(iiii-1)*4,0] = 1/Deltat*Tw_down_previous[iiii-1]
                        
                        #Rock temperature equation
                        L[1 + (iiii - 1) * 4,  (iiii - 1) * 4] = 1
                        L[1 + (iiii - 1) * 4, 1 + (iiii - 1) * 4] = -1
                        L[1 + (iiii - 1) * 4, 2 + (iiii - 1) * 4] = Rt
                        R[1 + (iiii - 1) * 4, 0] = 0
                        
                        #SBT equation              
                        L[2 + (iiii - 1) * 4, np.arange(2,4*N,4)] = NPCP[iiii-1, :N]
                        L[2 + (iiii - 1) * 4, 1 + (iiii - 1) * 4] = 1
                        R[2 + (iiii - 1) * 4, 0] = -BBCPOP[iiii-1] - BB[iiii-1] + BBinitial[iiii-1]                
                        
                        #Heat balance for upflowing fluid
                        L[3+(iiii-1)*4,3+(iiii-1)*4] = 1/Deltat + u_up/Deltaz[iiii-1] + 1/R_cp/(A_flow_centerpipe*rho_f*cp_f)
                        if iiii<N:
                            L[3+(iiii-1)*4,(iiii-1)*4] = -1/R_cp/(A_flow_centerpipe*rho_f*cp_f)
                            L[3+(iiii-1)*4,3+(iiii)*4] = -u_up/Deltaz[iiii-1]
                        else: #The bottom element has the downflowing fluid becoming the upflowing fluid
                            L[3+(iiii-1)*4,(iiii-1)*4] = -1/R_cp/(A_flow_centerpipe*rho_f*cp_f) - u_up/Deltaz[iiii-1];
                        
                        R[3+(iiii-1)*4,0] = 1/Deltat*Tw_up_previous[iiii-1];
                
                
                elif coaxialflowtype == 2: #CXC
                    #Populate L and R for heat balance upflowing fluid for first element
                    L[0,0] = 1/Deltat + u_up/Deltaz[0]  + 1/R_cp/(A_flow_annulus*rho_f*cp_f);
                    L[0,2] = -1/(A_flow_annulus*rho_f*cp_f);
                    L[0,3] = -1/R_cp/(A_flow_annulus*rho_f*cp_f);
                    L[0,4] = -u_up/Deltaz[0];
                    R[0,0] = 1/Deltat*Tw_up_previous[0];
                    
                    #Populate L and R for rock temperature equation for first element
                    L[1,0] = 1
                    L[1,1] = -1
                    L[1,2] = Rt
                    R[1,0] = 0            
                
                    #Populate L and R for SBT algorithm for first element
                    L[2,np.arange(2,4*N,4)] = NPCP[0,0:N]
                    L[2,1] = 1
                    R[2,0] =  - BBCPOP[0] - BB[0] + BBinitial[0]
                    
                    #Populate L and R for heat balance fluid down for first element
                    L[3,3] = 1/Deltat + u_down/Deltaz[0]*2 + 1/R_cp/(A_flow_centerpipe*rho_f*cp_f);
                    L[3,0] = -1/R_cp/(A_flow_centerpipe*rho_f*cp_f);
                    R[3,0] = 1/Deltat*Tw_down_previous[0] + u_down/Deltaz[0]*Tin*2;   
                    
                    for iiii in range(2, N+1):  #Populate L and R for remaining elements
                        #Heat balance upflowing fluid
                        L[(iiii-1)*4,(iiii-1)*4] = 1/Deltat + u_up/Deltaz[iiii-1] + 1/R_cp/(A_flow_annulus*rho_f*cp_f);
                        L[(iiii-1)*4,2+(iiii-1)*4] = -1/(A_flow_annulus*rho_f*cp_f);
                        R[(iiii-1)*4,0] = 1/Deltat*Tw_up_previous[iiii-1];
                        if iiii<N:
                            L[(iiii-1)*4,(iiii)*4] = -u_up/Deltaz[iiii-1];
                            L[(iiii-1)*4,3+(iiii-1)*4] = -1/R_cp/(A_flow_annulus*rho_f*cp_f);
                        else: #iiii==N is the bottom element where the downflowing fluid becomes the upflowing fluid
                            L[(iiii-1)*4,3+(iiii-1)*4] = -1/R_cp/(A_flow_annulus*rho_f*cp_f) - u_up/Deltaz[iiii-1];                
        
                        #Rock temperature equation
                        L[1 + (iiii - 1) * 4,  (iiii - 1) * 4] = 1
                        L[1 + (iiii - 1) * 4, 1 + (iiii - 1) * 4] = -1
                        L[1 + (iiii - 1) * 4, 2 + (iiii - 1) * 4] = Rt
                        R[1 + (iiii - 1) * 4, 0] = 0
                        
                        #SBT equation              
                        L[2 + (iiii - 1) * 4, np.arange(2,4*N,4)] = NPCP[iiii-1, :N]
                        L[2 + (iiii - 1) * 4, 1 + (iiii - 1) * 4] = 1
                        R[2 + (iiii - 1) * 4, 0] = -BBCPOP[iiii-1] - BB[iiii-1] + BBinitial[iiii-1]  
        
                        #Heat balance downflowing fluid
                        L[3+(iiii-1)*4,3+(iiii-1)*4] = 1/Deltat + u_down/Deltaz[iiii-1] + 1/R_cp/(A_flow_centerpipe*rho_f*cp_f);
                        L[3+(iiii-1)*4,(iiii-1)*4] = -1/R_cp/(A_flow_centerpipe*rho_f*cp_f);
                        L[3+(iiii-1)*4,3+(iiii-2)*4] = -u_down/Deltaz[iiii-1];
                        R[3+(iiii-1)*4,0] = 1/Deltat*Tw_down_previous[iiii-1];                
            
            elif clg_configuration == 2: #U-loop geometry
                #Populate L and R for fluid heat balance for first element (which has the injection temperature specified)
                L[0, 0] = 1 / Deltat + uvector[0] / Deltaz[0] * (fullyimplicit) * 2
                L[0, 2] = -4 / np.pi / Dvector[0]**2 / rho_f / cp_f
                R[0, 0] = 1 / Deltat * Twprevious[0] + uvector[0] / Deltaz[0] * Tinj * 2 - uvector[0] / Deltaz[0] * Twprevious[0] * (1 - fullyimplicit) * 2
            
                #Populate L and R for rock temperature equation for first element   
                L[1, 0] = 1
                L[1, 1] = -1
                L[1, 2] = Rtvector[0]
                R[1, 0] = 0
                # Populate L and R for SBT algorithm for first element
                L[2, np.arange(2,3*N,3)] = NPCP[0,0:N]
                L[2,1] = 1
                R[2, 0] = -BBCPOP[0] - BB[0] + BBinitial[0]
            
                for iiii in range(2, N+1):  
                    # Heat balance equation
                    L[0+(iiii - 1) * 3,  (iiii - 1) * 3] = 1 / Deltat + uvector[iiii-1] / Deltaz[iiii-1] / 2 * (fullyimplicit) * 2
                    L[0+(iiii - 1) * 3, 2 + (iiii - 1) * 3] = -4 / np.pi / Dvector[iiii-1] ** 2 / rho_f / cp_f
                
                    if iiii == len(xinj):  # Upcoming pipe has first element temperature sum of all incoming water temperatures
                        for j in range(len(lateralendpoints)):
                            L[0+ (iiii - 1) * 3, 0 + (lateralendpoints[j]) * 3] = -ulateral[j] / Deltaz[iiii-1] / 2 / lateralflowmultiplier * (fullyimplicit) * 2 * (radiuslateral/radiusvertical)**2
                            R[0+(iiii - 1) * 3, 0] = 1 / Deltat * Twprevious[iiii-1] + uvector[iiii-1] / Deltaz[iiii-1] * (
                                    -Twprevious[iiii-1] + (radiuslateral/radiusvertical)**2*np.sum(lateralflowallocation_norm[j] * Twprevious[lateralendpoints[j]])) / 2 * (
                                                            1 - fullyimplicit) * 2
                    else:
                        L[0+(iiii-1) * 3, 0 + (int(previouswaterelements[iiii-1])) * 3] = -uvector[iiii-1] / Deltaz[iiii-1] / 2 * (
                                fullyimplicit) * 2
                        R[0+(iiii-1) * 3, 0] = 1 / Deltat * Twprevious[iiii-1] + uvector[iiii-1] / Deltaz[iiii-1] * (
                                -Twprevious[iiii-1] + Twprevious[int(previouswaterelements[iiii-1])]) / 2 * (1 - fullyimplicit) * 2
            
                    # Rock temperature equation
                    L[1 + (iiii - 1) * 3,  (iiii - 1) * 3] = 1
                    L[1 + (iiii - 1) * 3, 1 + (iiii - 1) * 3] = -1
                    L[1 + (iiii - 1) * 3, 2 + (iiii - 1) * 3] = Rtvector[iiii-1]
                    R[1 + (iiii - 1) * 3, 0] = 0
                
                    # SBT equation 
                    L[2 + (iiii - 1) * 3, np.arange(2,3*N,3)] = NPCP[iiii-1, :N]
                    L[2 + (iiii - 1) * 3, 1 + (iiii - 1) * 3] = 1
                    R[2 + (iiii - 1) * 3, 0] = -BBCPOP[iiii-1] - BB[iiii-1] + BBinitial[iiii-1]
            
            
            # Solving the linear system of equations
            Sol = np.linalg.solve(L, R)    
        
        elif sbt_version == 2: #we need to perform iterative convergence for pressure, temperature, and fluid properties
            kk = 1
            maxrelativechange = 1
            TfluidupnodesOld = np.zeros(N+1)
            TfluiddownnodesOld = np.zeros(N+1)
            while kk <= maxnumberofiterations and maxrelativechange > reltolerance: #While loop iterates till either maximum number of iterations is reached or maximum relative change is less than user-defined target relative tolerance
                # Calculate frictional pressure drop (only turbulent flow is considered) [Pa]
                if np.any(Refluidupmidpoints < 2300):
                    print('Error: laminar flow in pipes; only turbulent flow models built-in for frictional pressure drop calculation')
                    print('Simulation terminated')
                    exit()
        
                if np.any(Refluiddownmidpoints < 2300):
                    print('Error: laminar flow in pipes; only turbulent flow models built-in for frictional pressure drop calculation')
                    print('Simulation terminated')
                    exit()
            
                fup = 1E-5 * np.ones(len(Refluidupmidpoints))  # Initial guess for upflowing turbulent flow friction factor
                fdown = 1E-5 * np.ones(len(Refluiddownmidpoints))  # Initial guess for downflowing turbulent flow friction factor
                for dd in range(1, 6):  # We assume 5 iterations are sufficient to converge turbulent friction factor
                    if coaxialflowtype == 1:  # CXA (injection in annulus; production from center pipe)
                        fup = 1 / (-2 * np.log10(eps_centerpipe / 3.7 / (2 * radiuscenterpipe) + 2.51 / Refluidupmidpoints / np.sqrt(fup))) ** 2
                        fdown = 1 / (-2 * np.log10(eps_annulus / 3.7 / Dh_annulus + 2.51 / Refluiddownmidpoints / np.sqrt(fdown))) ** 2
                    else:  # CXC (injection in center pipe; production from annulus)
                        fup = 1 / (-2 * np.log10(eps_annulus / 3.7 / Dh_annulus + 2.51 / Refluidupmidpoints / np.sqrt(fup))) ** 2
                        fdown = 1 / (-2 * np.log10(eps_centerpipe / 3.7 / (2 * radiuscenterpipe) + 2.51 / Refluiddownmidpoints / np.sqrt(fdown))) ** 2
        
                if coaxialflowtype == 1: #CXA (injection in annulus; production from center pipe)
                    DeltaP_frictionpipeup = fup*1/2*densityfluidupmidpoints*velocityfluidupmidpoints**2/(2*radiuscenterpipe)*Deltaz #Upflowing frictional pressure drop in pipe segments [Pa]
                    DeltaP_frictionpipedown = fdown*1/2*densityfluiddownmidpoints*velocityfluiddownmidpoints**2/(Dh_annulus)*Deltaz #Downflowing frictional pressure drop in pipe segments [Pa]
                elif coaxialflowtype == 2: #CXC (injection in center pipe; production from annulus)
                    DeltaP_frictionpipeup = fup*1/2*densityfluidupmidpoints*velocityfluidupmidpoints**2/(Dh_annulus)*Deltaz #Upflowing frictional pressure drop in pipe segments [Pa]
                    DeltaP_frictionpipedown = fdown*1/2*densityfluiddownmidpoints*velocityfluiddownmidpoints**2/(2*radiuscenterpipe)*Deltaz #Downflowing frictional pressure drop in pipe segments [Pa]
        
                #Calculate acceleration pressure change [Pa]
                DeltaP_accelerationup = densityfluidupnodes[1:]*velocityfluidupnodes[1:]**2 - densityfluidupnodes[:-1]*velocityfluidupnodes[:-1]**2
                DeltaP_accelerationdown = densityfluiddownnodes[1:]*velocityfluiddownnodes[1:]**2 - densityfluiddownnodes[:-1]*velocityfluiddownnodes[:-1]**2
                
                #Calculate nodal and midpoint fluid pressures [Pa]
                Pfluiddownnodes = (Pin * 1e5 
                    - np.cumsum(np.concatenate(([0], g * verticalchange * densityfluiddownmidpoints))) 
                    - np.cumsum(np.concatenate(([0], DeltaP_frictionpipedown))) 
                    - np.cumsum(np.concatenate(([0], DeltaP_accelerationdown))))
                
                CumSumPress = (-np.cumsum(np.concatenate(([0], g * verticalchange * densityfluidupmidpoints))) 
                + np.cumsum(np.concatenate(([0], DeltaP_frictionpipeup))) 
                + np.cumsum(np.concatenate(([0], DeltaP_accelerationup))))
                
                Pfluidupnodes = CumSumPress + (Pfluiddownnodes[-1] - CumSumPress[-1])
                Pfluidupmidpoints = 0.5 * (Pfluidupnodes[1:] + Pfluidupnodes[:-1])          #Midpoints pressures calculated as average of neighboring nodes [Pa]
                Pfluiddownmidpoints = 0.5 * (Pfluiddownnodes[1:] + Pfluiddownnodes[:-1])    #Midpoints pressures calculated as average of neighboring nodes [Pa]
        
                # Calculate thermal resistance in annulus and center pipe (assuming turbulent flow)
                Numidpointsup = 0.023 * (Refluidupmidpoints ** (4 / 5)) * (Prandtlfluidupmidpoints ** 0.4)  # Nusselt Number [-]
                Numidpointsdown = 0.023 * (Refluiddownmidpoints ** (4 / 5)) * (Prandtlfluiddownmidpoints ** 0.4)  # Nusselt Number [-]

                if coaxialflowtype == 1:  # CXA (injection in annulus; production from center pipe)
                    # Thermal resistance in annulus (downflowing)
                    Nu_down_o = Numidpointsdown  # %Based on Section 8.6 in Bergman (2011), for annulus turbulent flow, the Nusselt numbers for the inner and outer wall can be assumed the same
                    Nu_down_i = Numidpointsdown  
                    hmidpointsdown_o = Nu_down_o * thermalconductivityfluiddownmidpoints / Dh_annulus
                    hmidpointsdown_i = Nu_down_i * thermalconductivityfluiddownmidpoints / Dh_annulus
                    Rt = 1 / (np.pi * hmidpointsdown_o * radius * 2)  #Thermal resistance between annulus flow and surrounding rock (open-hole assumed)
                
                    # Thermal resistance in center pipe (upflowing)
                    Nu_up = Numidpointsup
                    hmidpointsup = Nu_up * thermalconductivityfluidupmidpoints / (2 * radiuscenterpipe)
                    R_cp = (
                        1 / (np.pi * hmidpointsup * 2 * radiuscenterpipe) +
                        np.log((2 * outerradiuscenterpipe) / (2 * radiuscenterpipe)) / (2 * np.pi * k_center_pipe) +
                        1 / (np.pi * hmidpointsdown_i * 2 * outerradiuscenterpipe)
                    )  # thermal resistance between annulus flow and center pipe flow
            
                elif coaxialflowtype == 2:  # CXC (injection in center pipe; production from annulus)
                    # Thermal resistance in annulus (upflowing)
                    Nu_up_o = Numidpointsup  # Based on Section 8.6 in Bergman (2011), for annulus turbulent flow, the Nusselt numbers for the inner and outer wall can be assumed the same
                    Nu_up_i = Numidpointsup  
                    hmidpointsup_o = Nu_up_o * thermalconductivityfluidupmidpoints / Dh_annulus
                    hmidpointsup_i = Nu_up_i * thermalconductivityfluidupmidpoints / Dh_annulus
                    Rt = 1 / (np.pi * hmidpointsup_o * radius * 2)  # Thermal resistance between annulus flow and surrounding rock
                
                    # Thermal resistance in center pipe (downflowing)
                    Nu_down = Numidpointsdown
                    hmidpointsdown = Nu_down * thermalconductivityfluiddownmidpoints / (2 * radiuscenterpipe)
                    R_cp = (
                        1 / (np.pi * hmidpointsdown * 2 * radiuscenterpipe) +
                        np.log((2 * outerradiuscenterpipe) / (2 * radiuscenterpipe)) / (2 * np.pi * k_center_pipe) +
                        1 / (np.pi * hmidpointsup_i * 2 * outerradiuscenterpipe)
                    )  # Thermal resistance between annulus flow and center pipe flow

                #Deltahstar is used in the fluid energy balance equation and specifies the different in enthalpy due to a difference in pressure
                if variablefluidproperties == 1: #The most accurate method uses the enthalpy property tables and is used when no constant fluid properties are specified.
                
                    # Calculate Deltahstar for downward flow
                    Deltahstardown = (
                        interpolator_enthalpy(np.array([[x, y] for x, y in zip(Pfluiddownnodes[1:], Tfluiddownnodes[:-1] + 273.15)]))
                        - interpolator_enthalpy(np.array([[x, y] for x, y in zip(Pfluiddownnodes[:-1], Tfluiddownnodes[:-1] + 273.15)]))
                    )
                    # Calculate Deltahstar for upward flow
                    Deltahstarup = (
                        interpolator_enthalpy(np.array([[x, y] for x, y in zip(Pfluidupnodes[:-1], Tfluidupnodes[1:] + 273.15)]))
                        - interpolator_enthalpy(np.array([[x, y] for x, y in zip(Pfluidupnodes[1:], Tfluidupnodes[1:] + 273.15)]))
                    )
                else: #If constant fluid properties are specified, then Deltahstar simplifies to 1/rho*(deltaP) (because the thermal expansion coefficient is zero) (the equations below show the full equation including the thermal expansion coefficient, so that it can also be used when fluid properties are not constant and the fluid is compressible)
                    Deltahstardown = (
                        1.0 / densityfluiddownnodes[:-1]
                        * (Pfluiddownnodes[1:] - Pfluiddownnodes[:-1])
                        * (1 - thermalexpansionfluiddownmidpoints * (Tfluiddownnodes[:-1] + 273.15))
                    )
                
                    Deltahstarup = (
                        1.0 / densityfluidupnodes[1:]
                        * (Pfluidupnodes[:-1] - Pfluidupnodes[1:])
                        * (1 - thermalexpansionfluidupmidpoints * (Tfluidupnodes[1:] + 273.15))
                    )
                
                #populate L and R
                if coaxialflowtype == 1: #CXA (injection in annulus, production from center pipe) (1:Tdown; 2:Tr; 3:Q; 4:Tup)
                    #Populate L and R for downflowing fluid energy balance for first element (which has the injection temperature specified)
                    L[0,0] = m*heatcapacityfluiddownmidpoints[0] + 1/R_cp[0]*Deltaz[0]/2
                    L[0,2] = -1*Deltaz[0]
                    L[0,3] = -1/R_cp[0]*Deltaz[0]/2
                    L[0,3+4] = -1/R_cp[0]*Deltaz[0]/2
                    R[0,0] = m*heatcapacityfluiddownmidpoints[0]*Tinj - 1/R_cp[0]*Deltaz[0]/2*Tinj - m*0.5*(velocityfluiddownnodes[1]**2-velocityfluiddownnodes[0]**2) - m*g*verticalchange[0] - m*Deltahstardown[0]
                    
                    #Populate L and R for rock temperature equation for first element
                    L[1,0] = 1/2
                    L[1,1] = -1
                    L[1,2] = Rt[0]
                    R[1,0] = -1/2*Tinj
                    
                    #Populate L and R for SBT algorithm for first element
                    L[2,np.arange(2,4*N,4)] = NPCP[0,0:N]
                    L[2,1] = 1
                    R[2,0] =  - BBCPOP[0] - BB[0] + BBinitial[0]                
                    
                    #Populate L and R for upflowing energy heat balance for first element
                    L[3,3] = m*heatcapacityfluidupmidpoints[0] + 1/R_cp[0]*Deltaz[0]/2
                    L[3,0] = -1/R_cp[0]*Deltaz[0]/2
                    L[3,7] = -m*heatcapacityfluidupmidpoints[0] + 1/R_cp[0]*Deltaz[0]/2
                    R[3,0] = -m*0.5*(velocityfluidupnodes[0]**2-velocityfluidupnodes[1]**2) + m*g*verticalchange[0] - m*Deltahstarup[0] + 1/R_cp[0]*Deltaz[0]/2*Tin
                
                    for iiii in range(2, N+1):  #Populate L and R for remaining elements (1:Tdown; 2:Tr; 3:Q; 4:Tup)
                        #Energy balance equation for downflowing fluid
                        L[(iiii-1)*4,2+(iiii-1)*4] = -1*Deltaz[iiii-1]
                        L[(iiii-1)*4,3+(iiii-1)*4] = -1/R_cp[iiii-1]*Deltaz[iiii-1]/2
                        if iiii<N:
                            L[(iiii-1)*4,3+(iiii)*4] = -1/R_cp[iiii-1]*Deltaz[iiii-1]/2
                            L[(iiii-1)*4,(iiii-1)*4] = m*heatcapacityfluiddownmidpoints[iiii-1] + 1/R_cp[iiii-1]*Deltaz[iiii-1]/2
                        else:
                            L[(iiii-1)*4,(iiii-1)*4] = m*heatcapacityfluiddownmidpoints[iiii-1]
                        
                        L[(iiii-1)*4,(iiii-2)*4] =  -m*heatcapacityfluiddownmidpoints[iiii-1] + 1/R_cp[iiii-1]*Deltaz[iiii-1]/2
                        R[(iiii-1)*4,0] = -m*0.5*(velocityfluiddownnodes[iiii]**2-velocityfluiddownnodes[iiii-1]**2) - m*g*verticalchange[iiii-1] - m*Deltahstardown[iiii-1]
                            
                        #Rock temperature equation
                        L[1+(iiii-1)*4,(iiii-2)*4] = 1/2
                        L[1+(iiii-1)*4,(iiii-1)*4] = 1/2
                        L[1+(iiii-1)*4,1+(iiii-1)*4] = -1
                        L[1+(iiii-1)*4,2+(iiii-1)*4] = Rt[iiii-1]
                        R[1+(iiii-1)*4,0] = 0
                            
                        #SBT equation              
                        L[2 + (iiii - 1) * 4, np.arange(2,4*N,4)] = NPCP[iiii-1, :N]
                        L[2 + (iiii - 1) * 4, 1 + (iiii - 1) * 4] = 1
                        R[2 + (iiii - 1) * 4, 0] = -BBCPOP[iiii-1] - BB[iiii-1] + BBinitial[iiii-1]  
                        
                        #Energy balance for upflowing fluid
                        L[3+(iiii-1)*4,3+(iiii-1)*4] = m*heatcapacityfluidupmidpoints[iiii-1]+1/R_cp[iiii-1]*Deltaz[iiii-1]/2
                        L[3+(iiii-1)*4,(iiii-2)*4] = -1/R_cp[iiii-1]*Deltaz[iiii-1]/2
                        if iiii<N:
                            L[3+(iiii-1)*4,(iiii-1)*4] = -1/R_cp[iiii-1]*Deltaz[iiii-1]/2
                            L[3+(iiii-1)*4,3+(iiii)*4] = -m*heatcapacityfluidupmidpoints[iiii-1] + 1/R_cp[iiii-1]*Deltaz[iiii-1]/2
                        else: #The bottom element has the downflowing fluid becoming the upflowing fluid
                            L[3+(iiii-1)*4,(iiii-1)*4] = -m*heatcapacityfluidupmidpoints[iiii-1]
                        R[3+(iiii-1)*4,0] = -m*0.5*(velocityfluidupnodes[iiii-1]**2-velocityfluidupnodes[iiii]**2) + m*g*verticalchange[iiii-1] - m*Deltahstarup[iiii-1]
                    
                elif coaxialflowtype == 2: #CXC (injection in center pipe, production from annulus) (1:Tup; 2:Tr; 3:Q; 4:Tdown)
                    #Populate L and R for upflowing fluid energy balance for first element
                    L[0,0] = m*heatcapacityfluidupmidpoints[0] + 1/R_cp[0]*Deltaz[0]/2
                    L[0,4] = -m*heatcapacityfluidupmidpoints[0] + 1/R_cp[0]*Deltaz[0]/2
                    L[0,2] = -1*Deltaz[0]
                    L[0,3] = -1/R_cp[0]*Deltaz[0]/2
                    R[0,0] = 1/R_cp[0]*Deltaz[0]/2*Tinj - m*0.5*(velocityfluidupnodes[0]**2-velocityfluidupnodes[1]**2) + m*g*verticalchange[0] - m*Deltahstarup[0]
                    
                    #Populate L and R for rock temperature equation for first element
                    L[1,0] = 1/2
                    L[1,1] = -1
                    L[1,2] = Rt[0]
                    L[1,4] = 1/2
                    
                    #Populate L and R for SBT algorithm for first element               
                    L[2,np.arange(2,4*N,4)] = NPCP[0,0:N]
                    L[2,1] = 1
                    R[2,0] =  - BBCPOP[0] - BB[0] + BBinitial[0]  

                    #Populate L and R for downflowing energy balance for first element (which has the injection temperature specified) 
                    L[3,3] = m*heatcapacityfluiddownmidpoints[0] + 1/R_cp[0]*Deltaz[0]/2
                    L[3,0] = -1/R_cp[0]*Deltaz[0]/2
                    L[3,4] = -1/R_cp[0]*Deltaz[0]/2
                    R[3,0] = m*heatcapacityfluiddownmidpoints[0]*Tinj - 1/R_cp[0]*Deltaz[0]/2*Tinj - m*0.5*(velocityfluiddownnodes[1]**2-velocityfluiddownnodes[0]**2) - m*g*verticalchange[0]-m*Deltahstardown[0]           
                    
                    for iiii in range(2, N+1): #Populate L and R for remaining elements (1:Tup; 2:Tr; 3:Q; 4:Tdown)
                        #Energy balance equation for upflowing fluid
                        L[(iiii-1)*4,(iiii-1)*4] = m*heatcapacityfluidupmidpoints[iiii-1]+1/R_cp[iiii-1]*Deltaz[iiii-1]/2
                        L[(iiii-1)*4,2+(iiii-1)*4] = -1*Deltaz[iiii-1]
                        L[(iiii-1)*4,3+(iiii-1)*4-4] = -1/R_cp[iiii-1]*Deltaz[iiii-1]/2
                        if iiii<N:
                            L[(iiii-1)*4,(iiii-1)*4+4] = -m*heatcapacityfluidupmidpoints[iiii-1]+1/R_cp[iiii-1]*Deltaz[iiii-1]/2
                            L[(iiii-1)*4,3+(iiii-1)*4] = -1/R_cp[iiii-1]*Deltaz[iiii-1]/2
                        else: #The bottom element has the downflowing fluid becoming the upflowing fluid
                            L[(iiii-1)*4,3+(iiii-1)*4] = -m*heatcapacityfluidupmidpoints[iiii-1]
                        R[(iiii-1)*4,0] = -m*0.5*(velocityfluidupnodes[iiii-1]**2-velocityfluidupnodes[iiii]**2)+m*g*verticalchange[iiii-1]-m*Deltahstarup[iiii-1]
                        
                        #Rock temperature equation
                        if iiii == N:
                            L[1+(iiii-1)*4,3+(iiii-1)*4] = 1/2 #At bottom element, the last downflowing fluid temperature becomes the first upflowing fluid temperature
                        else:
                            L[1+(iiii-1)*4,(iiii)*4] = 1/2
                        L[1+(iiii-1)*4,(iiii-1)*4] = 1/2
                        L[1+(iiii-1)*4,1+(iiii-1)*4] = -1
                        L[1+(iiii-1)*4,2+(iiii-1)*4] = Rt[iiii-1]
                        R[1+(iiii-1)*4,0] = 0
                        
                        #SBT equation              
                        L[2 + (iiii - 1) * 4, np.arange(2,4*N,4)] = NPCP[iiii-1, :N]
                        L[2 + (iiii - 1) * 4, 1 + (iiii - 1) * 4] = 1
                        R[2 + (iiii - 1) * 4, 0] = -BBCPOP[iiii-1] - BB[iiii-1] + BBinitial[iiii-1]   
                        
                        #Energy balance for downflowing fluid
                        L[3+(iiii-1)*4,3+(iiii-2)*4] = -m*heatcapacityfluiddownmidpoints[iiii-1]+1/R_cp[iiii-1]*Deltaz[iiii-1]/2
                        L[3+(iiii-1)*4,(iiii-1)*4] = -1/R_cp[iiii-1]*Deltaz[iiii-1]/2
                        if iiii<N:
                            L[3+(iiii-1)*4,3+(iiii-1)*4] = m*heatcapacityfluiddownmidpoints[iiii-1]+1/R_cp[iiii-1]*Deltaz[iiii-1]/2
                            L[3+(iiii-1)*4,(iiii)*4] = -1/R_cp[iiii-1]*Deltaz[iiii-1]/2
                        else: #The bottom element has the downflowing fluid becoming the upflowing fluid
                            L[3+(iiii-1)*4,3+(iiii-1)*4] = m*heatcapacityfluiddownmidpoints[iiii-1]
                        R[3+(iiii-1)*4,0] = - m*0.5*(velocityfluiddownnodes[iiii]**2-velocityfluiddownnodes[iiii-1]**2)-m*g*verticalchange[iiii-1]-m*Deltahstardown[iiii-1]           
                        
                        
                
                # Solving the linear system of equations
                Sol = np.linalg.solve(L, R)  
                if coaxialflowtype == 1: #CXA
                    Tfluiddownnodes = np.concatenate(([Tinj], Sol.ravel()[0::4]))
                    Tfluidupnodes =  np.concatenate((Sol.ravel()[3::4],[Tfluiddownnodes[-1]]))
                elif coaxialflowtype == 2: #CXC
                    Tfluiddownnodes = np.concatenate(([Tinj], Sol.ravel()[3::4]))
                    Tfluidupnodes = np.concatenate((Sol.ravel()[0::4],[Tfluiddownnodes[-1]]))
                
                Tfluiddownmidpoints = 0.5*Tfluiddownnodes[1:]+0.5*Tfluiddownnodes[:-1]
                Tfluidupmidpoints = 0.5*Tfluidupnodes[1:]+0.5*Tfluidupnodes[:-1]
                        
                #Update fluid properties
                densityfluiddownnodes = interpolator_density(np.array([[x, y] for x, y in zip(Pfluiddownnodes, Tfluiddownnodes + 273.15)]))
                densityfluidupnodes = interpolator_density(np.array([[x, y] for x, y in zip(Pfluidupnodes, Tfluidupnodes + 273.15)]))
                densityfluiddownmidpoints = interpolator_density(np.array([[x, y] for x, y in zip(Pfluiddownmidpoints, Tfluiddownmidpoints + 273.15)]))
                densityfluidupmidpoints = interpolator_density(np.array([[x, y] for x, y in zip(Pfluidupmidpoints, Tfluidupmidpoints + 273.15)]))
                viscosityfluiddownmidpoints = interpolator_viscosity(np.array([[x, y] for x, y in zip(Pfluiddownmidpoints, Tfluiddownmidpoints + 273.15)]))
                viscosityfluidupmidpoints = interpolator_viscosity(np.array([[x, y] for x, y in zip(Pfluidupmidpoints, Tfluidupmidpoints + 273.15)]))
                heatcapacityfluiddownmidpoints = interpolator_heatcapacity(np.array([[x, y] for x, y in zip(Pfluiddownmidpoints, Tfluiddownmidpoints + 273.15)]))
                heatcapacityfluidupmidpoints = interpolator_heatcapacity(np.array([[x, y] for x, y in zip(Pfluidupmidpoints, Tfluidupmidpoints + 273.15)]))
                thermalconductivityfluiddownmidpoints = interpolator_thermalconductivity(np.array([[x, y] for x, y in zip(Pfluiddownmidpoints, Tfluiddownmidpoints + 273.15)]))
                thermalconductivityfluidupmidpoints = interpolator_thermalconductivity(np.array([[x, y] for x, y in zip(Pfluidupmidpoints, Tfluidupmidpoints + 273.15)]))
                alphafluiddownmidpoints = thermalconductivityfluiddownmidpoints/densityfluiddownmidpoints/heatcapacityfluiddownmidpoints
                alphafluidupmidpoints = thermalconductivityfluidupmidpoints/densityfluidupmidpoints/heatcapacityfluidupmidpoints
                thermalexpansionfluiddownmidpoints = interpolator_thermalexpansion(np.array([[x, y] for x, y in zip(Pfluiddownmidpoints, Tfluiddownmidpoints + 273.15)]))
                thermalexpansionfluidupmidpoints = interpolator_thermalexpansion(np.array([[x, y] for x, y in zip(Pfluidupmidpoints, Tfluidupmidpoints + 273.15)]))
                Prandtlfluiddownmidpoints = viscosityfluiddownmidpoints/densityfluiddownmidpoints/alphafluiddownmidpoints
                Prandtlfluidupmidpoints = viscosityfluidupmidpoints/densityfluidupmidpoints/alphafluidupmidpoints

                if coaxialflowtype == 1: #CXA
                    Refluiddownmidpoints = densityfluiddownmidpoints*velocityfluiddownmidpoints*(Dh_annulus)/viscosityfluiddownmidpoints
                    Refluidupmidpoints = densityfluidupmidpoints*velocityfluidupmidpoints*(2*radiuscenterpipe)/viscosityfluidupmidpoints
                elif coaxialflowtype == 2: #CXC
                    Refluiddownmidpoints = densityfluiddownmidpoints*velocityfluiddownmidpoints*(2*radiuscenterpipe)/viscosityfluiddownmidpoints
                    Refluidupmidpoints = densityfluidupmidpoints*velocityfluidupmidpoints*(Dh_annulus)/viscosityfluidupmidpoints
                
                #Update fluid velocity
                if coaxialflowtype == 1: #CXA
                    velocityfluiddownmidpoints = m/A_flow_annulus/densityfluiddownmidpoints
                    velocityfluidupmidpoints = m/A_flow_centerpipe/densityfluidupmidpoints
                    velocityfluiddownnodes = m/A_flow_annulus/densityfluiddownnodes
                    velocityfluidupnodes = m/A_flow_centerpipe/densityfluidupnodes
                elif coaxialflowtype == 2: #CXC
                    velocityfluiddownmidpoints = m/A_flow_centerpipe/densityfluiddownmidpoints
                    velocityfluidupmidpoints = m/A_flow_annulus/densityfluidupmidpoints
                    velocityfluiddownnodes = m/A_flow_centerpipe/densityfluiddownnodes
                    velocityfluidupnodes = m/A_flow_annulus/densityfluidupnodes

                maxrelativechange = np.max(np.concatenate([np.abs(Tfluidupnodes - TfluidupnodesOld) / Tfluidupnodes, np.abs(Tfluiddownnodes - TfluiddownnodesOld) / Tfluiddownnodes]))
                linetoprint = f"Step = {i} | Iteration = {kk} | Max Rel Change = {maxrelativechange}"
                print(linetoprint)
                TfluidupnodesOld = Tfluidupnodes
                TfluiddownnodesOld = Tfluiddownnodes
                kk = kk+1
                            
                
        
        # Extracting Q array for current heat pulses, fluid temperature, and store fluid temperature for the next time step
        if sbt_version == 1:
            if clg_configuration == 1: #co-axial geometry 
                Q[:, i] = Sol.ravel()[2::4]
                if coaxialflowtype == 1: #CXA
                    Tw_down_Matrix[i,:] = Sol.ravel()[np.arange(0,4*N,4)] #Store the temperature for the downflowing water
                    Tw_down_previous = Sol.ravel()[np.arange(0,4*N,4)]    #Store the current downflowing water temperature to be used as previous water temperature in the next time step
                    Tw_up_Matrix[i,:] = Sol.ravel()[np.arange(3,4*N,4)]   #Store the temperature for the upflowing water
                    Tw_up_previous = Sol.ravel()[np.arange(3,4*N,4)]      #Store the current upflowing water temperature to be used as previous water temperature in the next time step
                elif  coaxialflowtype == 2: #CXC
                    Tw_down_Matrix[i,:] = Sol.ravel()[np.arange(3,4*N,4)] #Store the temperature for the downflowing water
                    Tw_down_previous = Sol.ravel()[np.arange(3,4*N,4)]   #Store the current downflowing water temperature to be used as previous water temperature in the next time step
                    Tw_up_Matrix[i,:] = Sol.ravel()[np.arange(0,4*N,4)]   #Store the temperature for the upflowing water
                    Tw_up_previous = Sol.ravel()[np.arange(0,4*N,4)]     #Store the current upflowing water temperature to be used as previous water temperature in the next time step
            
                #Calculate the fluid outlet temperature at the top of the first element based on the water temperature at the midpoint of the top element and the one below
                Toutput[i] = Tw_up_previous[0]+(Tw_up_previous[0]-Tw_up_previous[1])*0.5 
            
            elif clg_configuration == 2: #U-loop geometry
                Q[:, i] = Sol.ravel()[2::3]
                TwMatrix[i, :] = Sol.ravel()[np.arange(0,3*N,3)]
                Twprevious = Sol.ravel()[np.arange(0,3*N,3)]
                
                #Calculate the fluid outlet temperature at the top of the first element based on the water temperature at the midpoint of the top element and the one below
                top_element_index = len(xinj) + len(xprod) - 3
                Toutput[i] = Twprevious[top_element_index] + (Twprevious[top_element_index] - Twprevious[top_element_index - 1]) * 0.5
                
        elif sbt_version == 2:
            if clg_configuration == 1: #co-axial geometry 
                Q[:, i] = Sol.ravel()[2::4]
                Toutput[i] = Tfluidupnodes[0]                        #Store the fluid outlet temperature
                Poutput[i] = Pfluidupnodes[0]                        #Store the fluid outlet pressure
                
                Tfluidupnodesstore[:, i] = Tfluidupnodes             #Store upflowing nodal fluid temperatures
                Tfluidupmidpointsstore[:,i] = Tfluidupmidpoints      #Store upflowing midpoint fluid temperatures
                Tfluiddownnodesstore[:,i] = Tfluiddownnodes          #Store downflowing nodal fluid temperatures
                Tfluiddownmidpointsstore[:,i] = Tfluiddownmidpoints  #Store downflowing midpoint fluid temperatures
                
                Pfluidupnodesstore[:,i] = Pfluidupnodes              #Store upflowing nodal fluid pressures
                Pfluidupmidpointsstore[:,i] = Pfluidupmidpoints      #Store upflowing midpoint fluid pressures
                Pfluiddownnodesstore[:,i] = Pfluiddownnodes          #Store downflowing nodal fluid pressures
                Pfluiddownmidpointsstore[:,i] = Pfluiddownmidpoints  #Store downflowing midpoint fluid pressures
                
            elif clg_configuration == 2: #U-loop geometry
                tbd = 1
        
        #FMM algorithm for combining heat pulses
        #---------------------------------------
        if (FMM == 1 and i>50 and times[i] > FMMtriggertime):
            remainingtimes = times[(times >= combinedtimes[-1]) & (times < times[i])] 
            currentendtimespassed = times[i] - remainingtimes
            
            if len(remainingtimes)>40:
                if currentendtimespassed[24]>FMMtriggertime:
                    combinedtimes = np.append(combinedtimes,remainingtimes[24])
                    startindex = np.where(times == combinedtimes[-2])[0][0]
                    endindex = np.where(times == remainingtimes[24])[0][0]
                    newcombinedQ = np.sum(Q[:,startindex+1:endindex+1]*(times[startindex+1:endindex+1]-times[startindex:endindex])/(combinedtimes[-1]-combinedtimes[-2]),axis = 1) #weighted average
                    newcombinedQ = newcombinedQ.reshape(-1, 1)
                    combinedQ = np.hstack((combinedQ, newcombinedQ))
            
            
            #combine very old time pulses
            startindexforsecondlevel = np.where(combinedtimes == combinedtimes2ndlevel[-1])[0][0]
            if combinedtimes.size>30 and combinedtimes.size-startindexforsecondlevel>30:
                if (times[i] - combinedtimes[startindexforsecondlevel+20])>FMMtriggertime*5:
                    indicestodrop = np.arange(startindexforsecondlevel+1, startindexforsecondlevel+20)
                    weightedQ = np.sum(combinedQ[:,startindexforsecondlevel+1:startindexforsecondlevel+20+1]*\
                                        (combinedtimes[startindexforsecondlevel+1:startindexforsecondlevel+20+1]-combinedtimes[startindexforsecondlevel:startindexforsecondlevel+20])\
                                            /(combinedtimes[startindexforsecondlevel+20]-combinedtimes[startindexforsecondlevel]),axis = 1)
                    combinedtimes2ndlevel = np.append(combinedtimes2ndlevel,combinedtimes[startindexforsecondlevel+20])
                    combinedtimes = np.delete(combinedtimes, indicestodrop)
                    combinedQ[:,startindexforsecondlevel+20] = weightedQ
                    combinedQ = np.delete(combinedQ, indicestodrop, axis=1)


            #combine very very old time pulses
            startindexforthirdlevel = np.where(combinedtimes == combinedtimes3rdlevel[-1])[0][0]
            if combinedtimes.size>50 and combinedtimes.size-startindexforthirdlevel>50:
                if (times[i] - combinedtimes[startindexforthirdlevel+20])>FMMtriggertime*10:
                    indicestodrop = np.arange(startindexforthirdlevel+1, startindexforthirdlevel+20)
                    weightedQ = np.sum(combinedQ[:,startindexforthirdlevel+1:startindexforthirdlevel+20+1]*\
                                        (combinedtimes[startindexforthirdlevel+1:startindexforthirdlevel+20+1]-combinedtimes[startindexforthirdlevel:startindexforthirdlevel+20])\
                                            /(combinedtimes[startindexforthirdlevel+20]-combinedtimes[startindexforthirdlevel]),axis = 1)
                    combinedtimes3rdlevel = np.append(combinedtimes3rdlevel,combinedtimes[startindexforthirdlevel+20])
                    combinedtimes = np.delete(combinedtimes, indicestodrop)
                    combinedQ[:,startindexforthirdlevel+20] = weightedQ
                    combinedQ = np.delete(combinedQ, indicestodrop, axis=1)
            
            
            endindex = np.where(times == combinedtimes[-1])[0][0]
            remainingQ = Q[:,endindex+1:i+1]
            QFMM = np.hstack((combinedQ, remainingQ))
            timesFMM = np.append(combinedtimes,times[endindex+1:i+1])  
            
            
        else:
            QFMM = Q[:,0:i+1]
            timesFMM = times[0:i+1]
            

    
        # print('Time = ' + str(round(times[i]/3600/24/365*100)/100) + ' years')
        # print('Outlet Temperature = ' + str(round(Toutput[i],1)) + ' °C')
        # filename = 'python'+str(i)+'.mat'
        # if i>1:
        #     scipy.io.savemat(filename, dict(timestep=i,LPython=L,RPython=R,CPCPPython=CPCP,CPOPPython=CPOP,NPCPPython=NPCP,NPOPPython=NPOP,SolPython = Sol))
        
    #%% -----------------
    # 4. Post-Processing
    # The user can modify this section depending on the desired figures and simulation results
    #-------------------
    if sbt_version == 1:    
        HeatProduction = mstore * cp_f * (Toutput - Tinstore) / 1e6  # Calculates the heat production [MW]
    elif sbt_version == 2:
        if variablefluidproperties == 1: #Calculates the heat production as produced enthalpy minus injected enthalpy [MWth]
            HeatProduction = m*(interpolator_enthalpy(np.array([[x, y] for x, y in zip(Pfluidupnodesstore[0,:], Tfluidupnodesstore[0,:] + 273.15)])) - 
                                interpolator_enthalpy(np.array([[x, y] for x, y in zip(Pfluiddownnodesstore[0,:], Tfluiddownnodesstore[0,:] + 273.15)])))/1e6
        else: #For constant fluid properties, calculates the heat production as m*cp*DeltaT [MWth]
            HeatProduction = m*cp_f*(Toutput-Tin)/1e6
            
    AverageProductionTemperature = np.sum((times[1:] - times[:-1]) * Toutput[1:]) / times[-1]  # Calculates the weighted-average production temperature [deg.C]
    AverageHeatProduction = np.sum((times[1:] - times[:-1]) * HeatProduction[1:]) / times[-1]  # Calculates the weighted-average heat production [MW]
    line_to_print = f'Average production temperature = {AverageProductionTemperature:.2f} °C\n'
    # print(line_to_print, end='')

    line_to_print = f'Average heat production = {AverageHeatProduction:.2f} MWt\n'
    # print(line_to_print, end='')

    # line_to_print = f'Calculation time = {passedtime:.2f} s\n'
    # print(line_to_print)
    if is_plot:
        plot_final_fluid_temp_profile_v1(sbt_version=sbt_version, clg_configuration=clg_configuration, 
                                        Tw_up_previous=Tw_up_previous,
                                        Tw_down_previous=Tw_down_previous, Tfluiddownnodes=Tfluiddownnodes, 
                                        Deltaz=Deltaz, TwMatrix=TwMatrix, 
                                        numberoflaterals=numberoflaterals, coaxialflowtype=coaxialflowtype,
                                        interconnections=interconnections_new,
                                        lateralflowallocation=lateralflowallocation,
                                        xinj=xinj, xlat=xlat, xprod=xprod)
        
        plot_final_fluid_temp_profile_v2(sbt_version=sbt_version, 
                                        coaxialflowtype=coaxialflowtype, 
                                        Pfluiddownnodes=Pfluiddownnodes, Pfluidupnodes=Pfluidupnodes, 
                                        Deltaz=Deltaz)
        
        plot_heat_production(HeatProduction=HeatProduction, times=times)
        plot_production_temperature_linear(Toutput=Toutput, Tinstore=Tinstore, times=times)
        plot_production_tempterature_log(Toutput=Toutput, Tinstore=Tinstore, times=times)

    # print(Toutput + 273.15)
    # print(times[14:].shape)
    # print(Toutput[14:].shape)
    return times/365/24/3600, Toutput + 273.15 # return in Kelvin
    # return times[13:]/365/24/3600, Toutput[13:] + 273.15 # return in Kelvin # already done in clgs.py


