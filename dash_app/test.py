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

    
def compute_tube_geometry(sbt_version, clg_configuration, radiuscenterpipe, thicknesscenterpipe, 
                                xinj, xprod, xlat, numberoflaterals, radiuslateral, lateralflowallocation):

    interconnections = radiusvector = None
    Deltaz = np.sqrt((x[1:] - x[:-1]) ** 2 + (y[1:] - y[:-1]) ** 2 + (z[1:] - z[:-1]) ** 2)  # Length of each segment [m]
    # Deltaz = Deltaz.reshape(-1) # AB update
    Deltaz = Deltaz.ravel() 

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
