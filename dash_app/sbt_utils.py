# -------------------------------------------
# SBT inputs and computing functions
# -------------------------------------------

import math
import pandas as pd
import numpy as np
import scipy.io
from scipy.special import erfc, jv, yv
from scipy.interpolate import RegularGridInterpolator
from paths import inpath_dict

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
        times = np.array(times)
    elif mesh_fineness == 2:
        #Note 2: To capture the start-up effects, several small time steps are taken during the first 10,000 seconds in the time vector considered. To speed up the simulation, this can be avoided with limited impact on the long-term results. For example, an alternative time vector would be:
        times = [0] + list(range(100, 1000, 100)) + list(range(1000, 10000, 1000)) + list(np.logspace(np.log10(100*100), np.log10(20*365*24*3600), 75))
        times = np.array(times)

    fullyimplicit = None
    if clg_configuration == 2:
        fullyimplicit = 1            # Should be between 0 and 1. Only required when clg_configuration is 2. Most stable is setting it to 1 which results in a fully implicit Euler scheme when calculting the fluid temperature at each time step. With a value of 0, the convective term is modelled using explicit Euler. A value of 0.5 would model the convective term 50% explicit and 50% implicit, which may be slightly more accurate than fully implicit.

    variablefluidproperties = piperoughness = variableinjectiontemperature = variableflowrate = flowratefilename = None

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
        
        # AB: should these be editable to the user?
        # converge parameters
        reltolerance = HYPERPARAM4 # 1e-5                   # Target maximum acceptable relative tolerance each time step [-]. Lower tolerance will result in more accurate results but requires longer computational time
        maxnumberofiterations = HYPERPARAM5 # 15            # Maximum number of iterations each time step [-]. Each time step, solution converges until relative tolerance criteria is met or maximum number of time steps is reached.

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
        # NOTE: option for speeding up, change step size here
        verticaldepthsection = np.arange(0, -DrillingDepth_L1-1, -100) #removed *1000 on DrillingDepth_L1
        horizontalextentsection = np.arange(100, HorizontalExtent_L2+1, 100) #removed *1000 on HorizontalExtent_L2
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

        for i in range(int(numberoflaterals)):
            x = np.concatenate((x, xlat[:, i].reshape(-1,1)))
            y = np.concatenate((y, ylat[:, i].reshape(-1,1)))
            z = np.concatenate((z, zlat[:, i].reshape(-1,1)))

        # this could be slightly faster but not tested for accuracy
        # x = np.vstack([xinj, xprod, xlat[:, :numberoflaterals]])
        # y = np.vstack([yinj, yprod, ylat[:, :numberoflaterals]])
        # z = np.vstack([zinj, zprod, zlat[:, :numberoflaterals]])

    return locals()
    

def set_tube_geometry(sbt_version, clg_configuration, Diameter1, Diameter2, PipeParam3, PipeParam4, PipeParam5):
   
    """ 
        ## Tube Geometry
        Diameter1                        # If coaxial: radius, radiuscenterpipe and If U-loop: radiusvertical and radiuscenterpipe
        Diameter2                        #
    """

    radius = radiuscenterpipe = thicknesscenterpipe = coaxialflowtype = autoadjustlateralflowrates = None # coaxial
    numberoflaterals = radiusvertical = radiuslateral = lateralflowmultiplier = lateralflowallocation =  None # uloop 
    
    if clg_configuration == 1:  # co-axial geometry (1)
        radius = Diameter1/2  # 0.2286/2                # Wellbore radius [m] (everything is assumed open-hole)
        radiuscenterpipe = Diameter2/2  # 0.127/2       # Inner radius of inner pipe [m]
        thicknesscenterpipe = PipeParam3 # 0.0127       # Thickness of inner pipe [m]
        k_center_pipe = PipeParam4 # 0.006                                           # Thermal conductivity of insulation of center pipe wall [W/m/K]
        coaxialflowtype = PipeParam5 # 1                                             # 1 = CXA (fluid injection in annulus); 2 = CXC (fluid injection in center pipe)

    elif clg_configuration == 2: # U-loop geometry (2)

        radiusvertical = Diameter1/2 # 0.15                         # Radius of "vertical" injection and production well (open hole assumed for heat transfer) [m] (it is labeled here as vertical but it is allowed to be deviated)
        radiuslateral = Diameter2/2 # 0.30                          # Radius of laterals (open hole assumed for heat transfer) [m]
        
        numberoflaterals = int(PipeParam3) if PipeParam3 is not None else 1  # Number of laterals (must be integer) [-]
        numberoflaterals = max(1, numberoflaterals)  # Ensure at least 1 lateral
        lateralflowallocation = PipeParam4 # [1/3, 1/3, 1/3]                         # Distribution of flow accross laterals, must add up to 1 (it will get normalized below if sum does not equal to 1). Length of array must match number of laterals.
        # PipeParam5 is only for coaxial geometry (flow direction). For U-tube, lateralflowmultiplier should come from HyperParam5 or default to 1.
        # If PipeParam5 is a string (like "Inject in Annulus"), it's from coaxial UI and should be ignored for U-tube.
        if PipeParam5 is not None and isinstance(PipeParam5, (int, float)):
            try:
                lateralflowmultiplier = float(PipeParam5)
            except (ValueError, TypeError):
                lateralflowmultiplier = 1
        else:
            lateralflowmultiplier = 1  # Default for U-tube when PipeParam5 is None or a string 

        # TODO: ADD NEW PARAMETER (required if clg_configuration is 2)
        # USER EDITABLE? (v27)
        if sbt_version == 2:
            autoadjustlateralflowrates = 0               # Only used in SBT v2 U-loop. Must be 0 or 1. "0" means the flow rate in each lateral remains constant with a distribution as specified in lateralflowallocation. "1" means that the flow rate in each lateral is adjusted over time to ensure matching fluid pressures at the exit of each lateral. Sometimes this can cause convergence issues.

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
    k_f=0.68 # k_f = 0.667 
    mu_f=600*10**-6

    g = 9.81                       # Gravitational acceleration [m/s^2]
    gamma = 0.577215665            # Euler's constant
    alpha_f = k_f / rho_f / cp_f   # Fluid thermal diffusivity [m2/s]
    Pr_f = mu_f / rho_f / alpha_f  # Fluid Prandtl number [-]

    return locals()


def compute_tube_geometry(sbt_version, clg_configuration, fluid,
                                radius, radiuscenterpipe, thicknesscenterpipe,
                                x, y, z, xinj, yinj, zinj, xprod, yprod, zprod, xlat, ylat, zlat,
                                mdot,
                                piperoughness,
                                numberoflaterals, radiuslateral, radiusvertical, lateralflowmultiplier, lateralflowallocation):

    interconnections = radiusvector = None
    mlateral = mlateralold = None
    Deltaz = np.sqrt((x[1:] - x[:-1]) ** 2 + (y[1:] - y[:-1]) ** 2 + (z[1:] - z[:-1]) ** 2)  # Length of each segment [m]
    Deltaz = Deltaz.reshape(-1)

    if clg_configuration == 1: # COAXIAL
        # Ensure radiuscenterpipe and thicknesscenterpipe are not None for coaxial
        if radiuscenterpipe is None:
            raise ValueError(f"radiuscenterpipe is None for coaxial configuration. This should be set from Diameter2.")
        if thicknesscenterpipe is None:
            raise ValueError(f"thicknesscenterpipe is None for coaxial configuration. This should be set from PipeParam3.")
        
        outerradiuscenterpipe = radiuscenterpipe+thicknesscenterpipe    # Outer radius of inner pipe [m]
        A_flow_annulus = math.pi*(radius**2-outerradiuscenterpipe**2)   # Flow area of annulus pipe [m^2]
        A_flow_centerpipe = math.pi*radiuscenterpipe**2                 # Flow area of center pipe [m^2]
        Dh_annulus = 2*(radius-outerradiuscenterpipe)                   # Hydraulic diameter of annulus [m]
        
        # Validate geometry to prevent numerical instability
        if A_flow_annulus <= 0:
            raise ValueError(f"Invalid annulus flow area: {A_flow_annulus:.6f} m². "
                           f"radius={radius:.6f} m, outerradiuscenterpipe={outerradiuscenterpipe:.6f} m. "
                           f"This indicates invalid coaxial geometry (center pipe outer radius >= wellbore radius).")
        if Dh_annulus <= 0:
            raise ValueError(f"Invalid annulus hydraulic diameter: {Dh_annulus:.6f} m. "
                           f"radius={radius:.6f} m, outerradiuscenterpipe={outerradiuscenterpipe:.6f} m. "
                           f"This indicates invalid coaxial geometry (center pipe outer radius >= wellbore radius).")
        if A_flow_centerpipe <= 0:
            raise ValueError(f"Invalid center pipe flow area: {A_flow_centerpipe:.6f} m². "
                           f"radiuscenterpipe={radiuscenterpipe:.6f} m. "
                           f"This indicates invalid center pipe geometry.")
        
        # Warn if annulus is very thin (may cause numerical instability)
        min_annulus_area_ratio = 0.01  # 1% of wellbore area
        min_annulus_area = min_annulus_area_ratio * math.pi * radius**2
        if A_flow_annulus < min_annulus_area:
            print(f"[WARNING] Annulus flow area is very small ({A_flow_annulus:.6f} m² < {min_annulus_area:.6f} m²). "
                  f"This may cause numerical instability in the solver. "
                  f"radius={radius:.6f} m, outerradiuscenterpipe={outerradiuscenterpipe:.6f} m", flush=True)
        
        if sbt_version == 2:
            eps_annulus = Dh_annulus*piperoughness                      # Relative roughness annulus [-]
            eps_centerpipe = 2*radiuscenterpipe*piperoughness           # Relative roughness inner pipe [-]
        LoverR = Deltaz / radius                                        # Ratio of pipe segment length to radius along the wellbore [-]
        RelativeLengthChanges = (Deltaz[1:] - Deltaz[:-1]) / Deltaz[:-1]
        
    elif clg_configuration == 2: # U-LOOP

        # Ensure numberoflaterals is valid for interconnections calculation
        num_lat_int = int(numberoflaterals) if numberoflaterals is not None else 1
        num_lat_int = max(1, num_lat_int)  # Ensure at least 1
        if num_lat_int > 1:
            interconnections = np.concatenate((np.array([len(xinj)],dtype=int), np.array([len(xprod)],dtype=int), (np.ones(num_lat_int - 1, dtype=int) * len(xlat))))
        else:
            interconnections = np.concatenate((np.array([len(xinj)],dtype=int), np.array([len(xprod)],dtype=int)))
        interconnections = np.cumsum(interconnections)  # lists the indices of interconnections between inj, prod, and laterals (this will used to take care of the duplicate coordinates of the start and end points of the laterals)
        # Use original calculation that matches Deltaz structure after deletion
        # Use len(xlat) to match interconnections calculation (for 2D array, len() returns first dimension)
        radiusvector = np.concatenate([np.ones(len(xinj) + len(xprod) - 2) * radiusvertical, np.ones(numberoflaterals * len(xlat) - numberoflaterals) * radiuslateral])  # Stores radius of each element in a vector [m]
        Dvector = radiusvector * 2  # Diameter of each element [m]
        lateralflowallocation_norm = lateralflowallocation / np.sum(lateralflowallocation)  # Ensure the sum equals 1   
        
        if sbt_version == 2:
            radiusvectornodes = np.concatenate([np.full(len(xinj) + len(xprod), radiusvertical),np.full(numberoflaterals * xlat.shape[0] - 2 * numberoflaterals, radiuslateral)]) #Stores the element radius at each node [m] (The bottom nodes of the injector and producer are assume to have the radius of the injector and producer)
            Dvectornodes = radiusvectornodes*2 #Vector with element diameter at each node [m]
            FlowDistributionMidPoints = np.ones(len(xinj) + len(xprod) - 2)
            FlowDistributionNodes = np.ones(len(xinj)+len(xprod))
            for dd in range(numberoflaterals):
                # Ensure lateralflowallocation_norm[dd] is a scalar, not an array
                allocation_val = lateralflowallocation_norm[dd]
                if hasattr(allocation_val, '__len__') and not isinstance(allocation_val, str):
                    allocation_scalar = float(np.asarray(allocation_val).item())
                else:
                    allocation_scalar = float(allocation_val)
                multiplier_val = float(lateralflowmultiplier) if lateralflowmultiplier is not None else 1.0
                FlowDistributionMidPoints = np.concatenate([FlowDistributionMidPoints, multiplier_val * allocation_scalar * np.ones(xlat.shape[0] - 1)])
                FlowDistributionNodes = np.concatenate([FlowDistributionNodes, multiplier_val * allocation_scalar * np.ones(xlat.shape[0] - 2)])
            Area = math.pi*Dvector**2/4              #Vector with cross sectional flow area of each element at the element midpoints (all elements are assumed open hole) [m2]
            AreaNodes = math.pi*Dvectornodes**2/4  #Vector with cross sectional flow area of each element at the nodes (all elements are assumed open hole) [m2]
            eps = Dvector*piperoughness       #Vector with relative roughness at midpoint of each element [-]
            signofcorrection = 1              #Parameter used in the script for updating the lateral flow rates to equalize the fluid lateral outlet pressures
            secondaverageabslateralflowcorrection = 1 #Parameter used in the script for updating the lateral flow rates to equalize the fluid lateral outlet pressures
            mvector = mdot*FlowDistributionMidPoints #Array that stores the flow rate in each element (initially assumes uniform distribution of flow among laterals but this will be adjusted below (if autoadjustlateralflowrates = 1) to ensure identical pressure change accross all laterals) [kg/s]
            mnodalvector = mdot*FlowDistributionNodes #Array that stores the flow rate at each node (initialy assumes uniform distribution of flow among laterals but this will be adjusted below (if autoadjustlateralflowrates = 1) to ensure identical pressure change accross all laterals) [kg/s]
            mlateral = lateralflowallocation_norm*mdot*lateralflowmultiplier #%Array that stores the flow rate through each lateral (initially assumes uniform flow distribution accross the laterals)
            mlateralold = mlateral.copy()

        Deltaz = np.delete(Deltaz, interconnections - 1)  # Removes the phantom elements due to duplicate coordinates
        # radiusvector is already created with the correct length (matching Deltaz after deletion), so no deletion needed
        Dvector = radiusvector * 2  # Diameter of each element [m]
        
        # Verify Dvector and Deltaz have matching lengths
        if len(Dvector) != len(Deltaz):
            raise ValueError(f"Dvector length ({len(Dvector)}) does not match Deltaz length ({len(Deltaz)}) after deletion. radiusvector length: {len(radiusvector)}")
        
        # Update Area and eps for SBT v2 after Dvector is modified
        if sbt_version == 2:
            Area = math.pi*Dvector**2/4              #Vector with cross sectional flow area of each element at the element midpoints (all elements are assumed open hole) [m2]
            eps = Dvector*piperoughness       #Vector with relative roughness at midpoint of each element [-]
        
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
    
    # QUALITY CONTROL # AB: bring back later on UI
    # if smallestLoverR < 10:
    #     print('Warning: smallest ratio of segment length over radius is less than 10. Good practice is to keep this ratio larger than 10.')
    # if max(abs(RelativeLengthChanges)) > 0.5:
    #     print('Warning: abrupt change(s) in segment length detected, which may cause numerical instabilities. Good practice is to avoid abrupt length changes to obtain smooth results.')
        
    return locals()


def calc_tube_min_time_steps(clg_configuration, radius, radiusvector, alpha_m, Deltaz, LimitPointSourceModel, LimitCylinderModelRequired, LimitInfiniteModel):

    timeforpointssource = max(Deltaz)**2 / alpha_m * LimitPointSourceModel  # Calculates minimum time step size when point source model becomes applicable [s]

    if clg_configuration == 1: # co-axial geometry (1)
        timeforlinesource = radius**2 / alpha_m * LimitCylinderModelRequired  # Calculates minimum time step size when line source model becomes applicable [s]

    elif clg_configuration == 2: # U-loop geometry (2)
        timeforlinesource = max(radiusvector)**2 / alpha_m * LimitCylinderModelRequired  # Calculates minimum time step size when line source model becomes applicable [s]

    timeforfinitelinesource = max(Deltaz)**2 / alpha_m * LimitInfiniteModel  # Calculates minimum time step size when finite line source model should be considered [s]

    return timeforpointssource, timeforlinesource, timeforfinitelinesource,


def prepare_interpolators(sbt_version, variablefluidproperties, fluid, rho_f, cp_f, k_f, mu_f, variableflowrate=None):

    interpolator_density = interpolator_enthalpy = interpolator_entropy = None
    interpolator_heatcapacity = interpolator_heatcapacity = interpolator_phase = None 
    interpolator_thermalconductivity = interpolator_thermalexpansion = interpolator_viscosity = None

    if sbt_version == 1:
        if fluid == 2:  # CO2
            if variableflowrate == 1:  # Variable Mass Flow Rate Mode
                try:
                    mat1 = scipy.io.loadmat(inpath_dict["properties_CO2_pathname_sbtv2"])
                    mat2 = scipy.io.loadmat(inpath_dict["properties_CO2v2_pathname"])
                    Pvector = mat1['Pvector']
                    Tvector = mat1['Tvector']
                    density = mat1['density']
                    enthalpy = mat1['enthalpy']
                    entropy = mat1['entropy']
                    heatcapacity = mat1['heatcapacity']
                    phase = mat1['phase']
                    thermalconductivity = mat1['thermalconductivity']
                    thermalexpansion = mat1['thermalexpansion']
                    viscosity = mat1['viscosity']
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
                except Exception as e:
                    print(f"Error loading CO2 properties files for SBT v1 (Variable Mass Flow Rate Mode): {e}")
                    raise
            else:  # Constant Mass Flow Rate Mode
                try:
                    mat = scipy.io.loadmat(inpath_dict["properties_CO2_pathname_sbtv2"])
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
                except Exception as e:
                    print(f"Error loading CO2 properties file for SBT v1 (Constant Mass Flow Rate Mode): {e}")
                    raise
        elif fluid == 1:  # H2O
            try:
                mat = scipy.io.loadmat(inpath_dict["properties_H2O_pathname_sbtv2"])
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
            except Exception as e:
                print(f"Error loading H2O properties file for SBT v1: {e}")
                raise

    if sbt_version == 2:

        if variablefluidproperties == 0:  # For computational purposes, use constant fluid property tables
            # Define vectors for pressure and temperature (convert to numpy arrays)
            Pvector = np.array([1, 1e9])
            Tvector = np.array([1, 1e4])
        
            # Create 2x2 arrays with constant fluid properties
            density = np.array([[rho_f] * 2] * 2)
            enthalpy = np.array([[0] * 2] * 2)  # Enthalpy difference not used in constant mode
            entropy = np.array([[0] * 2] * 2)  # Entropy difference not used in constant mode
            heatcapacity = np.array([[cp_f] * 2] * 2)
            phase = np.array([[1] * 2] * 2)  # Assume liquid phase (1) for constant properties
            thermalconductivity = np.array([[k_f] * 2] * 2)
            viscosity = np.array([[mu_f] * 2] * 2)
            thermalexpansion = np.array([[0] * 2] * 2)  # Incompressible fluid has zero thermal expansion coefficient
        else:  # If variable fluid properties, import pre-generated tables
            # print('Loading fluid properties...')
            if fluid == 1:  # H2O
                try:
                    # mat = scipy.io.loadmat('properties_H2O.mat') 
                    mat = scipy.io.loadmat(inpath_dict["properties_H2O_pathname_sbtv2"])
                    Pvector = mat['Pvector']
                    # print("pVector", Pvector)
                    Tvector = mat['Tvector']
                    density = mat['density']
                    enthalpy = mat['enthalpy']
                    entropy = mat['entropy']
                    heatcapacity = mat['heatcapacity']
                    phase = mat['phase']
                    thermalconductivity = mat['thermalconductivity']
                    thermalexpansion = mat['thermalexpansion']
                    viscosity = mat['viscosity']
                    # print('Fluid properties for water loaded successfully')
                except Exception as e:
                    print(f"Error loading properties for water: {e}")
                    raise
            elif fluid == 2:  #CO2
                try:
                    # mat = scipy.io.loadmat('properties_CO2.mat') 
                    mat = scipy.io.loadmat(inpath_dict["properties_CO2_pathname_sbtv2"]) 
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


def get_profiles(sbt_version, times, variableinjectiontemperature, variableflowrate, flowratefilename, Tinj, mdot):

    # Read injection temperature profile if provided
    Tinstore = np.zeros(len(times))
    # print(variableinjectiontemperature)
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

def precalculations(clg_configuration, Deltaz, alpha_m, k_m, times, 
                    NoArgumentsFinitePipeCorrection, NoDiscrFinitePipeCorrection, NoDiscrInfCylIntegration,
                    NoArgumentsInfCylIntegration,
                    x,y,z,
                    timeforlinesource, radius, radiusvector, interconnections):

    interconnections_new = None
    delta_times = [t2 - t1 for t1, t2 in zip(times[:-1], times[1:])]
    fpcmaxarg = max(Deltaz)**2 / (4 * alpha_m * min(delta_times))

    fpcminarg = min(Deltaz)**2 / (4 * alpha_m * times[-1])
    # fpcmaxarg = max(Deltaz)**2 / (4 * alpha_m * (min(times[1:] - times[:-1])))
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
        besselminarg = alpha_m * min(delta_times) / radius**2
        # besselminarg = alpha_m * (min(times[1:] - times[:-1])) / radius**2
        besselmaxarg = alpha_m * timeforlinesource / radius**2

    
    elif clg_configuration == 2: # U-loop geometry (2)
        besselminarg = alpha_m * min(delta_times) / max(radiusvector)**2
        # besselminarg = alpha_m * (min(times[1:] - times[:-1])) / max(radiusvector)**2
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
    
    # print("Elements in the tube: ", elementcenters.shape) # for coaxial it's (30000, 3) and for uloop it's (89, 3) at the moment
    # print(elementcenters)

    SMatrix = np.zeros((N, N))  # Initializes the spacing matrix, which holds the distance between center points of each element [m]
    SoverL = np.zeros((N, N))  # Initializes the ratio of spacing to element length matrix
    for i in range(N):
        SMatrix[i, :] = np.sqrt((elementcenters[i, 0] - elementcenters[:, 0])**2 + (elementcenters[i, 1] - elementcenters[:, 1])**2 + (elementcenters[i, 2] - elementcenters[:, 2])**2)
        SoverL[i, :] = SMatrix[i, :] / Deltaz[i] #Calculates the ratio of spacing between two elements and element length
    
    return Amin1vector, argumentbesselvec, besselcylinderresult, elementcenters, SMatrix, SoverL, N, interconnections_new, finitecorrectiony

