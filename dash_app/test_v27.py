if clg_configuration == 1: #co-axial geometry
    outerradiuscenterpipe = radiuscenterpipe+thicknesscenterpipe #Outer radius of inner pipe [m]
    A_flow_annulus = math.pi*(radius**2-outerradiuscenterpipe**2)       #Flow area of annulus pipe [m^2]
    A_flow_centerpipe = math.pi*radiuscenterpipe**2         #Flow area of center pipe [m^2]
    Dh_annulus = 2*(radius-outerradiuscenterpipe) #Hydraulic diameter of annulus [m]    
    if sbt_version == 2:
        eps_annulus = Dh_annulus*piperoughness            #Relative roughness annulus [-]
        eps_centerpipe = 2*radiuscenterpipe*piperoughness #Relative roughness inner pipe [-]
    
elif clg_configuration == 2: #U-loop geometry
    interconnections = np.concatenate((np.array([len(xinj)],dtype=int), np.array([len(xprod)],dtype=int), (np.ones(numberoflaterals - 1, dtype=int) * len(xlat))))
    interconnections = np.cumsum(interconnections)  # lists the indices of interconnections between inj, prod, and laterals (this will used to take care of the duplicate coordinates of the start and end points of the laterals)
    radiusvector = np.concatenate([np.ones(len(xinj) + len(xprod) - 2) * radiusvertical, np.ones(numberoflaterals * len(xlat) - numberoflaterals) * radiuslateral])  # Stores radius of each element in a vector [m]
    Dvector = radiusvector * 2  # Diameter of each element [m]
    lateralflowallocation = lateralflowallocation / np.sum(lateralflowallocation)  # Ensure the sum equals 1
    if sbt_version == 2:
        radiusvectornodes = np.concatenate([np.full(len(xinj) + len(xprod), radiusvertical),np.full(numberoflaterals * xlat.shape[0] - 2 * numberoflaterals, radiuslateral)]) #Stores the element radius at each node [m] (The bottom nodes of the injector and producer are assume to have the radius of the injector and producer)
        Dvectornodes = radiusvectornodes*2 #Vector with element diameter at each node [m]
        FlowDistributionMidPoints = np.ones(len(xinj) + len(xprod) - 2)
        FlowDistributionNodes = np.ones(len(xinj)+len(xprod))
        for dd in range(numberoflaterals):
            FlowDistributionMidPoints = np.concatenate([FlowDistributionMidPoints,lateralflowmultiplier * lateralflowallocation[dd] * np.ones(xlat.shape[0] - 1)])
            FlowDistributionNodes = np.concatenate([FlowDistributionNodes,lateralflowmultiplier * lateralflowallocation[dd] * np.ones(xlat.shape[0] - 2)])
        Area = math.pi*Dvector**2/4              #Vector with cross sectional flow area of each element at the element midpoints (all elements are assumed open hole) [m2]
        AreaNodes = math.pi*Dvectornodes**2/4  #Vector with cross sectional flow area of each element at the nodes (all elements are assumed open hole) [m2]
        eps = Dvector*piperoughness       #Vector with relative roughness at midpoint of each element [-]
        signofcorrection = 1              #Parameter used in the script for updating the lateral flow rates to equalize the fluid lateral outlet pressures
        secondaverageabslateralflowcorrection = 1 #Parameter used in the script for updating the lateral flow rates to equalize the fluid lateral outlet pressures
        mvector = m*FlowDistributionMidPoints #Array that stores the flow rate in each element (initially assumes uniform distribution of flow among laterals but this will be adjusted below (if autoadjustlateralflowrates = 1) to ensure identical pressure change accross all laterals) [kg/s]
        mnodalvector = m*FlowDistributionNodes #Array that stores the flow rate at each node (initialy assumes uniform distribution of flow among laterals but this will be adjusted below (if autoadjustlateralflowrates = 1) to ensure identical pressure change accross all laterals) [kg/s]
        mlateral = lateralflowallocation*m*lateralflowmultiplier #%Array that stores the flow rate through each lateral (initially assumes uniform flow distribution accross the laterals)
        mlateralold = mlateral.copy()
        
if sbt_version == 2:
    if variablefluidproperties == 0:  # For computational purposes, use constant fluid property tables
        # Define vectors for pressure and temperature
        Pvector = np.array([[1, 1e9]])
        Tvector = np.array([[1, 1e4]])
    
        # Create 2x2 arrays with constant fluid properties
        density = np.array([[rho_f] * 2] * 2)
        heatcapacity = np.array([[cp_f] * 2] * 2)
        thermalconductivity = np.array([[k_f] * 2] * 2)
        viscosity = np.array([[mu_f] * 2] * 2)
        thermalexpansion = np.array([[0] * 2] * 2)  # Incompressible fluid has zero thermal expansion coefficient
    else:  # If variable fluid properties, import pre-generated tables
        print('Loading fluid properties...')
        if fluid == 1:  # H2O
            try:
                mat = scipy.io.loadmat('properties_H2O.mat') 
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
    #Prepare interpolators
    Pvector_1d = Pvector.ravel()
    Tvector_1d = Tvector.ravel()
    interpolator_density = RegularGridInterpolator((Pvector_1d, Tvector_1d), density)
    interpolator_heatcapacity = RegularGridInterpolator((Pvector_1d, Tvector_1d), heatcapacity)
    interpolator_thermalconductivity = RegularGridInterpolator((Pvector_1d, Tvector_1d), thermalconductivity)
    interpolator_thermalexpansion = RegularGridInterpolator((Pvector_1d, Tvector_1d), thermalexpansion)
    interpolator_viscosity = RegularGridInterpolator((Pvector_1d, Tvector_1d), viscosity)
    if variablefluidproperties == 1:
        interpolator_enthalpy = RegularGridInterpolator((Pvector_1d, Tvector_1d), enthalpy)
        interpolator_entropy = RegularGridInterpolator((Pvector_1d, Tvector_1d), entropy)
        interpolator_phase = RegularGridInterpolator((Pvector_1d, Tvector_1d), phase)

Deltaz = np.sqrt((x[1:] - x[:-1]) ** 2 + (y[1:] - y[:-1]) ** 2 + (z[1:] - z[:-1]) ** 2)  # Length of each segment [m]
Deltaz = Deltaz.reshape(-1)
if clg_configuration == 2:
    Deltaz = np.delete(Deltaz, interconnections - 1)  # Removes the phantom elements due to duplicate coordinates
TotalLength = np.sum(Deltaz)  # Total length of all elements (for informational purposes only) [m]

# Quality Control
if clg_configuration == 1: #co-axial geometry
    LoverR = Deltaz / radius  # Ratio of pipe segment length to radius along the wellbore [-]
elif clg_configuration == 2: #U-loop geometry
    LoverR = Deltaz / radiusvector  # Ratio of pipe segment length to radius along the wellbore [-]
smallestLoverR = np.min(LoverR)  # Smallest ratio of pipe segment length to pipe radius. This ratio should be larger than 10. [-]

if smallestLoverR < 10:
    print('Warning: smallest ratio of segment length over radius is less than 10. Good practice is to keep this ratio larger than 10.')

if clg_configuration == 1: #co-axial geometry
    RelativeLengthChanges = (Deltaz[1:] - Deltaz[:-1]) / Deltaz[:-1]
elif clg_configuration == 2: #U-loop geometry
    if numberoflaterals > 1:
        DeltazOrdered = np.concatenate((Deltaz[0:(interconnections[0]-1)], Deltaz[(interconnections[1]-2):(interconnections[2]-3)], Deltaz[(interconnections[0]-1):(interconnections[1]-2)]))
    else:
        DeltazOrdered = np.concatenate((Deltaz[0:interconnections[0] - 1], Deltaz[interconnections[1] - 1:-1], Deltaz[interconnections[0]:interconnections[1] - 2]))
    RelativeLengthChanges = (DeltazOrdered[1:] - DeltazOrdered[:-1]) / DeltazOrdered[:-1]

if max(abs(RelativeLengthChanges)) > 0.5:
    print('Warning: abrupt change(s) in segment length detected, which may cause numerical instabilities. Good practice is to avoid abrupt length changes to obtain smooth results.')

if clg_configuration == 2: #additional quality control for U-loop geometry
    for dd in range(1, numberoflaterals + 1):
        if abs(xinj[-1] - xlat[0][dd - 1]) > 1e-12 or abs(yinj[-1] - ylat[0][dd - 1]) > 1e-12 or abs(zinj[-1] - zlat[0][dd - 1]) > 1e-12:
            print(f'Error: Coordinate mismatch between bottom of injection well and start of lateral #{dd}')
        if abs(xprod[0] - xlat[-1][dd - 1]) > 1e-12 or abs(yprod[0] - ylat[-1][dd - 1]) > 1e-12 or abs(zprod[0] - zlat[-1][dd - 1]) > 1e-12:
            print(f'Error: Coordinate mismatch between bottom of production well and end of lateral #{dd}')
    if len(lateralflowallocation) != numberoflaterals:
        print('Error: Length of array "lateralflowallocation" does not match the number of laterals')
  
