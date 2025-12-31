#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import time
from sbt_v27 import run_sbt as run_sbt_final

print('\n') # look for ***!!!!!!****
print("scenario 1")
print("utube - H2O")
start = time.time()
times, Tout, Pout = run_sbt_final(
        ## Model Specifications 
        sbt_version=1, 
        
        mesh_fineness=0, # coarse -- done 

        HYPERPARAM1=0, # variableflowrate (0 not variable, 1 is) -- done
        HYPERPARAM2='MassFlowRate.xlsx', # 'MassFlowRate.xlsx'  --done     
        HYPERPARAM3=0, # variableinjectiontemperature -- done
        HYPERPARAM4='InjectionTemperatures.xlsx', # injectiontemperaturefilename -- done
        HYPERPARAM5=None, # for sbt v1 it's None -- done

        accuracy=1, # coarse - done

        ## Geologic Properties
        Tsurf=25, # done
        GeoGradient=65/1000, # done 
        k_m=3.00, # done
        c_m=790.0, # done
        rho_m=2800, # done

        ## Operations
        clg_configuration=2, # 1 coaxial and 2 is utube - done
        mdot=30, # mass flow rate -- done
        Tinj=55, # injection temperature -- done
        fluid=1, # h2o is 1 sco2 is 2 -- done 
        DrillingDepth_L1=3.5, # km - done
        HorizontalExtent_L2=10.0, # km - done
        Diameter1=0.35, ## diametervertical - done
        Diameter2=0.35, ## diameterlateral - done
        PipeParam3=1, ## numberoflaterals - done - 3 works
        PipeParam4=[1], ## lateralflowallocation - done - [1/3, 1/3, 1/3] works
        PipeParam5=1, ## lateralflowmultiplier - done
        )
print("Time Steps: ", len(times)) # 175 time steps (in seconds)
print("Tout: ", Tout[-1]-272.15) # last value (in Kelvin converted to Celcius)
print("Pout: ", Pout) # should be None, later in clgs.py set to 20 MPa constant in Pa
end = time.time()
print("TOTAL TIME: ", end - start) # runs in 0.23 seconds 
print('\n') 

print('\n') # look for ***!!!!!!**** -- seems to be less elevated than the database (db get 128 Tout but here get 107 Tout)
print("scenario 1")
print("utube - sCO2")
start = time.time()
times, Tout, Pout = run_sbt_final(
        ## Model Specifications 
        sbt_version=2, 
        
        mesh_fineness=0, # coarse -- done 

        HYPERPARAM1=20*10, # Fluid input pressure [bar] -- done
        HYPERPARAM2=1e-6, # Pipe/borehole roughness --done     
        HYPERPARAM3=0, # variablefluidproperties -- constant fluid properties or 1 variable -- done -- ONLY mode that can be varied and it's across both fluids!!!
        HYPERPARAM4=1e-5 , # reltolerance | Target maximum acceptable relative tolerance -- done
        HYPERPARAM5=15, # maxnumberofiterations |Maximum number of iterations each time step [-]
        
        accuracy=1, # coarse - done

        ## Geologic Properties
        Tsurf=25, # done
        GeoGradient=65/1000, # done 
        k_m=3.00, # done
        c_m=790.0, # done
        rho_m=2800, # done

        ## Operations
        clg_configuration=2, # 1 coaxial and 2 is utube - done
        mdot=30, # mass flow rate -- done
        Tinj=55, # injection temperature -- done
        fluid=2, # h2o is 1 sco2 is 2 -- done 
        DrillingDepth_L1=3.5, # km - done
        HorizontalExtent_L2=10.0, # km - done
        Diameter1=0.35, ## diametervertical - done
        Diameter2=0.35, ## diameterlateral - done
        PipeParam3=1, ## numberoflaterals - done - 3 works
        PipeParam4=[1], ## lateralflowallocation - done - [1/3, 1/3, 1/3] works
        PipeParam5=1, ## lateralflowmultiplier - done
        )
print("Time Steps: ", len(times)) # 175 time steps (in seconds)
print("Tout: ", Tout[-1]-272.15) # last value (in Kelvin converted to Celcius)
print("Pout: ", Pout[-1]/1e6) # should be PRESENT, later in clgs.py set to 20 MPa constant in Pa
end = time.time()
print("TOTAL TIME: ", end - start) # runs in 0.62 seconds 
print('\n') 


print('\n') # look for ***!!!!!!****
print("scenario 1")
print("coaxial - H2O")
start = time.time() # TODO!
times, Tout, Pout = run_sbt_final(
        ## Model Specifications 
        sbt_version=1, 
        
        mesh_fineness=0, # coarse -- done 

        HYPERPARAM1=0, # variableflowrate (0 not variable, 1 is) -- done
        HYPERPARAM2='MassFlowRate.xlsx', # 'MassFlowRate.xlsx'  --done     
        HYPERPARAM3=0, # variableinjectiontemperature -- done
        HYPERPARAM4='InjectionTemperatures.xlsx', # injectiontemperaturefilename -- done
        HYPERPARAM5=None, # for sbt v1 it's None -- done

        accuracy=1, # coarse - done

        ## Geologic Properties
        Tsurf=25, # done
        GeoGradient=65/1000, # done 
        k_m=3.00, # done
        c_m=790.0, # done
        rho_m=2800, # done

        ## Operations
        clg_configuration=1, # 1 coaxial and 2 is utube - done
        mdot=30, # mass flow rate -- done
        Tinj=55, # injection temperature -- done
        fluid=1, # h2o is 1 sco2 is 2 -- done 
        
        DrillingDepth_L1=3.5, # km - done
        HorizontalExtent_L2=10.0, # km - done
        Diameter1=0.445, ## diameter wellbore (annulus)- done
        Diameter2=0.201, ## diameter center pipe - done

        PipeParam3=0.013, ## thicknesscenterpipe - done 
        PipeParam4=0.025, ##  insulation thermal conductivity | k_center_pipe -- done
        PipeParam5="Inject in Annulus", ## coaxialflowtype where 1 is in inject in annulus and 2 center -- done
        )
print("Time Steps: ", len(times)) # 175 time steps (in seconds)
print("Tout: ", Tout[-1]-272.15) # last value (in Kelvin converted to Celcius)
print("Pout: ", Pout) # should be None, later in clgs.py set to 20 MPa constant in Pa
end = time.time()
print("TOTAL TIME: ", end - start) # runs in 0.30 seconds 
print('\n') 


print('\n') # look for ***!!!!!!****
print("scenario 1")
print("coaxial - sCO2")
start = time.time()
times, Tout, Pout = run_sbt_final(
        ## Model Specifications 
        sbt_version=2, 
        
        mesh_fineness=0, # coarse -- done 

        HYPERPARAM1=20*10, # Fluid input pressure [bar] -- done
        HYPERPARAM2=1e-6, # Pipe/borehole roughness --done     
        HYPERPARAM3=0, # variablefluidproperties -- constant fluid properties or 1 variable -- done -- ONLY mode that can be varied and it's across both fluids!!!
        HYPERPARAM4=1e-5 , # reltolerance | Target maximum acceptable relative tolerance -- done
        HYPERPARAM5=15, # maxnumberofiterations |Maximum number of iterations each time step [-]
        
        accuracy=1, # coarse - done

        ## Geologic Properties
        Tsurf=25, # done
        GeoGradient=65/1000, # done 
        k_m=3.00, # done
        c_m=790.0, # done
        rho_m=2800, # done

        ## Operations
        clg_configuration=1, # 1 coaxial and 2 is utube - done
        mdot=30, # mass flow rate -- done
        Tinj=55, # injection temperature -- done
        fluid=2, # h2o is 1 sco2 is 2 -- done 
        
        DrillingDepth_L1=3.5, # km - done
        HorizontalExtent_L2=10.0, # km - done
        Diameter1=0.445, ## diameter wellbore (annulus)- done
        Diameter2=0.201, ## diameter center pipe - done

        PipeParam3=0.013, ## thicknesscenterpipe - done 
        PipeParam4=0.025, ##  insulation thermal conductivity | k_center_pipe -- done
        PipeParam5="Inject in Annulus", ## coaxialflowtype where 1 is in inject in annulus and 2 center -- done
        )
print("Time Steps: ", len(times)) # 175 time steps (in seconds)
print("Tout: ", Tout[-1]-272.15) # last value (in Kelvin converted to Celcius)
print("Pout: ", Pout[-1]/1e6) # should be PRESENT, later in clgs.py set to 20 MPa constant in Pa
end = time.time()
print("TOTAL TIME: ", end - start) # runs in 0.62 seconds 
print('\n') 



