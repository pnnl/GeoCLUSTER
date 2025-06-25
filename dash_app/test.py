def prepare_interpolators(sbt_version, variablefluidproperties, fluid, rho_f, cp_f, k_f, mu_f):

    interpolator_density = interpolator_enthalpy = interpolator_entropy = None
    interpolator_heatcapacity = interpolator_heatcapacity = interpolator_phase = None 
    interpolator_thermalconductivity = interpolator_thermalexpansion = interpolator_viscosity = None

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
        
        # Prepare interpolators
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
    
    return interpolator_density, interpolator_enthalpy, \
                interpolator_entropy, interpolator_heatcapacity, interpolator_heatcapacity, \
                    interpolator_phase, interpolator_thermalconductivity, interpolator_thermalexpansion, interpolator_viscosity
