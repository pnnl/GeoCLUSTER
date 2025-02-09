#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from sbt import run_sbt

# example
# Tout = run_sbt(sbt_version=1, clg_configuration=1, modt=24, Tinj=30, GeoGradient=0.05, k_m=3)
# Tout = run_sbt(clg_configuration=2, sbt_version=1, 
#                mdot=20, Tinj=20, HorizontalExtent_L2=1, DrillingDepth_L1=2, 
#                Tsurf=20, GeoGradient=90/1000, k_m=2.83, c_m=825, rho_m=2875)

run_sbt(
            ## Model Specifications 
            sbt_version=1, mesh_fineness=0, HYPERPARAM1=0, HYPERPARAM2="MassFlowRate.xlsx", 
            HYPERPARAM3=0, HYPERPARAM4="InjectionTemperatures.xlsx", HYPERPARAM5=None, 
            accuracy=1,

             ## Operations
            clg_configuration=1, mdot=20, Tinj=20, fluid=1, ## Operations
            DrillingDepth_L1=3.5, HorizontalExtent_L2=10.0, #BoreholeDiameter=1, ## Wellbore Geometry
            Diameter1=0.2286, Diameter2=0.127, PipeParam3=0.0127, PipeParam4=0.006,
            PipeParam5=1, ## Tube Geometry

            ## Geologic Properties
            Tsurf=25, GeoGradient=50/1000, k_m=3.00, c_m=790.0, rho_m=2750, # k_m=2.83, c_m=825, rho_m=2875, 
            )

# run_sbt(
#             ## Model Specifications 
#             sbt_version=1, mesh_fineness=0, HYPERPARAM1=0, HYPERPARAM2="MassFlowRate.xlsx", 
#             HYPERPARAM3=0, HYPERPARAM4="InjectionTemperatures.xlsx", HYPERPARAM5=None, 
#             accuracy=1,

#              ## Operations
#             clg_configuration=2, # coaxial = 1 and uloop = 2
#             mdot=24, Tinj=30, fluid=1, ## Operations
#             DrillingDepth_L1=3.5, HorizontalExtent_L2=10.0, #BoreholeDiameter=1, ## Wellbore Geometry
#             # DrillingDepth_L1=3500, HorizontalExtent_L2=10000.0, #BoreholeDiameter=1, ## Wellbore Geometry
#             Diameter1=0.3500, Diameter2=0.3500, PipeParam3=1, PipeParam4=[1],  #0.006
#             PipeParam5=1, ## Tube Geometry

#             ## Geologic Properties
#             Tsurf=25, GeoGradient=50/1000, k_m=3.00, c_m=790.0, rho_m=2750, 
#             )

#             # NOTE: HDF5: fluid dependent properties assumed vs. sbt assumption is constant

# start_vals_d = {"mdot": 24.0, "L2": 10000, "L1": 3500 , "grad": 0.05, "D": 0.3500, "Tinj": 30.0, "k": 3.0}

# case = water, sbt1,
# run_sbt(
#             ## Model Specifications 
#             sbt_version=sbt_version,

#              ## Operations
#             clg_configuration=case, mdot=mdot, Tinj=Tinj, fluid=fluid, ## Operations
#             DrillingDepth_L1=L1, HorizontalExtent_L2=L2,
#             Diameter1=D, Diameter2=D,  
#             )