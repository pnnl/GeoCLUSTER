#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from dash import dcc, html, Input, Output, State, ctx, ALL
import dash_bootstrap_components as dbc
from dash.exceptions import PreventUpdate

PARAMETER_INFO = {
    # Heat Transfer and System Configuration
    "Heat Transfer Mode": {
        "definition": "Specify the way thermal energy moves due to a temperature difference. The main modes are conduction (through solids or contact) or convection (through moving fluids).",
        "recommended_range": "Conduction, Convection",
        "typical_value": "Conduction",
        "unit": "mode",
        "description": "Specify the way thermal energy moves due to a temperature difference. The main modes are conduction (through solids or contact) or convection (through moving fluids)."
    },
    
    "Heat Exchanger": {
        "definition": "Specify the geometry of the system. \"U-Tube\" uses two separate wells or laterals; \"Coaxial\" uses concentric pipes for flow-in and flow-out. Coaxial is often the default in simple single-well simulations, but U-Tube systems often have better long-term heat extraction in closed-loop systems.",
        "recommended_range": "U-Tube, Coaxial",
        "typical_value": "U-Tube",
        "unit": "mode",
        "description": "Specify the geometry of the system. \"U-Tube\" uses two separate wells or laterals; \"Coaxial\" uses concentric pipes for flow-in and flow-out. Coaxial is often the default in simple single-well simulations, but U-Tube systems often have better long-term heat extraction in closed-loop systems."
    },
    
    "Working Fluid": {
        "definition": "Specify a liquid or gas that carries heat through the system. H₂O (water) is typical for geothermal due to availability; sCO₂ (supercritical carbon dioxide) is included for advanced systems with enhanced heat extraction potential.",
        "recommended_range": "H2O, sCO2",
        "typical_value": "H2O",
        "unit": "fluid",
        "description": "Specify a liquid or gas that carries heat through the system. H₂O (water) is typical for geothermal due to availability; sCO₂ (supercritical carbon dioxide) is included for advanced systems with enhanced heat extraction potential."
    },
    
    "End-Use": {
        "definition": "Specify the intended application of the geothermal system output. Options include Heating (direct thermal use), Electricity (power generation), or All (both applications).",
        "recommended_range": "Heating, Electricity, All",
        "typical_value": "All",
        "unit": "mode",
        "description": "The end-use selection determines the economic calculations and system design parameters."
    },
    
    "Model Version": {
        "definition": "Select the computational model for simulating closed-loop geothermal systems.",
        "recommended_range": "Database, Simulator",
        "typical_value": "Database",
        "unit": "model",
        "description": "Database: Use pre-calculated database model for fast results from pre-computed scenarios. Simulator: Live model for simulating scenarios not in the database, including depths deeper than 5 km, geothermal gradients larger than 70°C/km, and number of laterals greater than 1."
    },
    
    # Geologic Properties
    "Surface Temperature (˚C)": {
        "definition": "Set the ground-level or ambient temperature. A value of 25°C is a typical average surface temperature in temperate regions during geothermal operation.",
        "recommended_range": "0-40°C",
        "typical_value": "25°C",
        "unit": "°C",
        "description": "Surface temperature affects the initial conditions for geothermal calculations and heat transfer modeling."
    },
    
    "Surface Temperature (˚F)": {
        "definition": "Set the ground-level or ambient temperature. A value of 77°F is a typical average surface temperature in temperate regions during geothermal operation.",
        "recommended_range": "32-104°F",
        "typical_value": "77°F",
        "unit": "°F",
        "description": "Surface temperature affects the initial conditions for geothermal calculations and heat transfer modeling."
    },
    
    "Geothermal Gradient (°C/m)": {
        "definition": "Set the rate at which temperature increases with depth. The western United States typically exhibits a higher geothermal gradient (~0.035 °C/m) than the eastern U.S. (~0.025 °C/m).",
        "recommended_range": "0.015-0.050°C/m",
        "typical_value": "0.03°C/m",
        "unit": "°C/m",
        "description": "Higher gradients indicate more favorable geothermal conditions for energy extraction."
    },
    
    "Geothermal Gradient (˚F/ft)": {
        "definition": "Set the rate at which temperature increases with depth. A value of 0.027°F/ft means that the temperature increases by 27°F for every 1000 feet of depth. This represents average conditions in continental crust and it is hot enough to run a small power plant or provide heating for buildings.",
        "recommended_range": "0.008-0.109°F/ft",
        "typical_value": "0.027°F/ft",
        "unit": "°F/ft",
        "description": "Higher gradients indicate more favorable geothermal conditions for energy extraction."
    },
    
    "Rock Thermal Conductivity (W/m-°C)": {
        "definition": "Set how quickly heat moves through rock. A value of 3 W/m-K represents moderately conductive rock, such as granite.",
        "recommended_range": "0.4-5.0 W/m-K",
        "typical_value": "3 W/m-K",
        "unit": "W/m-K",
        "description": "Higher conductivity improves heat transfer from the rock to the working fluid."
    },
    
    "Rock Thermal Conductivity (Btu/ft-h-˚F)": {
        "definition": "Set how quickly heat moves through rock. A value of 1.73 Btu/ft-h-°F represents moderately conductive rock, such as granite.",
        "recommended_range": "0.23-2.89 Btu/ft-h-°F",
        "typical_value": "1.73 Btu/ft-h-°F",
        "unit": "Btu/ft-h-°F",
        "description": "Higher conductivity improves heat transfer from the rock to the working fluid."
    },
    
    "Rock Specific Heat Capacity (J/kg-°C)": {
        "definition": "Set the amount of energy the rock can absorb or release when its temperature changes by 1°C, which determines how quickly the rock heats up or cools down in response to fluid circulation. A value of 0.051 J/kg-K represents an average for various dry rocks.",
        "recommended_range": "500-2000 J/kg-K",
        "typical_value": "0.051 J/kg-K",
        "unit": "J/kg-K",
        "description": "Affects the thermal storage capacity of the rock formation."
    },
    
    "Rock Specific Heat Capacity (Btu/lb-˚F)": {
        "definition": "Set the amount of energy the rock can absorb or release when its temperature changes by 1°F, which determines how quickly the rock heats up or cools down in response to fluid circulation. A value of 0.000012 Btu/lb-°F represents an average for various dry rocks.",
        "recommended_range": "0.119-0.477 Btu/lb-°F",
        "typical_value": "0.000012 Btu/lb-°F",
        "unit": "Btu/lb-°F",
        "description": "Affects the thermal storage capacity of the rock formation."
    },
    
    "Rock Density (kg/m3)": {
        "definition": "Set the mass of rock per unit volume, which affects heat storage and fluid flow behavior in a geothermal system. For example, denser rock holds more heat and changes temperature more slowly. A value of 790 kg/m³ implies the rock is highly porous, fractured, or contains gas-filled voids. Typical rock densities usually range from 2,000 to 3,000 kg/m³.",
        "recommended_range": "1000-3500 kg/m³",
        "typical_value": "2800 kg/m³",
        "unit": "kg/m³",
        "description": "Density influences the thermal mass and heat storage capacity of the formation."
    },
    
    "Rock Density (lb/ft3)": {
        "definition": "Set the mass of rock per unit volume, which affects heat storage and fluid flow behavior in a geothermal system. For example, denser rock holds more heat and changes temperature more slowly. A value of 49.3 lb/ft³ implies the rock is highly porous, fractured, or contains gas-filled voids. Typical rock densities usually range from 124.9 to 218.6 lb/ft³.",
        "recommended_range": "62.4-218.6 lb/ft³",
        "typical_value": "49.3 lb/ft³",
        "unit": "lb/ft³",
        "description": "Density influences the thermal mass and heat storage capacity of the formation."
    },
    
    # Wellbore Operations
    "Injection Temperature (˚C)": {
        "definition": "Set the temperature of the fluid entering the subsurface. A value of 30°C is a common injection temperature for low-enthalpy systems.",
        "recommended_range": "30-60°C",
        "typical_value": "30°C",
        "unit": "°C",
        "description": "Lower injection temperatures generally improve heat extraction efficiency."
    },
    
    "Injection Temperature (˚F)": {
        "definition": "Set the temperature of the fluid entering the subsurface. A value of 86°F is a common injection temperature for low-enthalpy systems.",
        "recommended_range": "86-140°F",
        "typical_value": "86°F",
        "unit": "°F",
        "description": "Lower injection temperatures generally improve heat extraction efficiency."
    },
    
    "Mass Flow Rate (kg/s)": {
        "definition": "Set the total mass of working fluid that circulates through the geothermal system every second. A value of 24 kg/s moves enough fluid to extract significant heat but keeps pumping requirements and pressure losses manageable.",
        "recommended_range": "5-300 kg/s",
        "typical_value": "24 kg/s",
        "unit": "kg/s",
        "description": "Higher flow rates increase heat extraction but may require more pumping power."
    },
    
    "Mass Flow Rate (lb/s)": {
        "definition": "Set the total mass of working fluid that circulates through the geothermal system every second. A value of 52.9 lb/s moves enough fluid to extract significant heat but keeps pumping requirements and pressure losses manageable.",
        "recommended_range": "11.0-661.4 lb/s",
        "typical_value": "52.9 lb/s",
        "unit": "lb/s",
        "description": "Higher flow rates increase heat extraction but may require more pumping power."
    },
    
    # Tube Geometry
    "Borehole Diameter (m)": {
        "definition": "Set the width of the hole drilled into the ground to access the geothermal reservoir. A value of 0.35 m can manage frictional pressure losses where lower widths can increase pressure drops and reduce heat transfer.",
        "recommended_range": "0.2159-0.4445 m",
        "typical_value": "0.35 m",
        "unit": "m",
        "description": "Larger diameters allow for higher flow rates but increase drilling costs."
    },
    
    "Borehole Diameter (ft)": {
        "definition": "Set the width of the hole drilled into the ground to access the geothermal reservoir. A value of 1.15 ft can manage frictional pressure losses where lower widths can increase pressure drops and reduce heat transfer.",
        "recommended_range": "0.71-1.46 ft",
        "typical_value": "1.15 ft",
        "unit": "ft",
        "description": "Larger diameters allow for higher flow rates but increase drilling costs."
    },
    
    "Wellbore Diameter Vertical (m)": {
        "definition": "Set the diameter of the vertical injection and production well of the U-tube design. A value of 0.444 m is a relatively large open-hole diameter for maximizing heat transfer surface area.",
        "recommended_range": "0.2159-0.4445 m",
        "typical_value": "0.444 m",
        "unit": "m",
        "description": "Affects the heat transfer area and flow resistance in the vertical section."
    },
    
    "Wellbore Diameter Vertical (ft)": {
        "definition": "Set the diameter of the vertical injection and production well of the U-tube design. A value of 1.46 ft is a relatively large open-hole diameter for maximizing heat transfer surface area.",
        "recommended_range": "0.71-1.46 ft",
        "typical_value": "1.46 ft",
        "unit": "ft",
        "description": "Affects the heat transfer area and flow resistance in the vertical section."
    },
    
    "Wellbore Diameter Lateral (m)": {
        "definition": "Set the diameter of the lateral branches of the U-tube design. A value of 0.444 m is a relatively large open-hole diameter for maximizing heat transfer surface area.",
        "recommended_range": "0.2159-0.4445 m",
        "typical_value": "0.444 m",
        "unit": "m",
        "description": "Influences heat transfer and flow characteristics in the lateral sections."
    },
    
    "Wellbore Diameter Lateral (ft)": {
        "definition": "Set the diameter of the lateral branches of the U-tube design. A value of 1.46 ft is a relatively large open-hole diameter for maximizing heat transfer surface area.",
        "recommended_range": "0.71-1.46 ft",
        "typical_value": "1.46 ft",
        "unit": "ft",
        "description": "Influences heat transfer and flow characteristics in the lateral sections."
    },
    
    # Base parameter names (without units) for pattern matching
    "Wellbore Diameter Vertical": {
        "definition": "Set the diameter of the vertical injection and production well of the U-tube design. A value of 0.444 m is a relatively large open-hole diameter for maximizing heat transfer surface area.",
        "recommended_range": "0.2159-0.4445 m",
        "typical_value": "0.444 m",
        "unit": "m",
        "description": "Affects the heat transfer area and flow resistance in the vertical section."
    },
    
    "Wellbore Diameter Lateral": {
        "definition": "Set the diameter of the lateral branches of the U-tube design. A value of 0.444 m is a relatively large open-hole diameter for maximizing heat transfer surface area.",
        "recommended_range": "0.2159-0.4445 m",
        "typical_value": "0.444 m",
        "unit": "m",
        "description": "Influences heat transfer and flow characteristics in the lateral sections."
    },
    
    "Horizontal Extent (m)": {
        "definition": "Set the horizontal length of the well. A value of 10 km represents long multi-lateral systems. A value of 50 km far exceeds directional drilling and would require massive pressure support and well integrity.",
        "recommended_range": "1000-50000 m",
        "typical_value": "10000 m",
        "unit": "m",
        "description": "Longer horizontal extents increase heat extraction area but require more drilling."
    },
    
    "Horizontal Extent (ft)": {
        "definition": "Set the horizontal length of the well. A value of 32,808 ft represents long multi-lateral systems. A value of 164,042 ft far exceeds directional drilling and would require massive pressure support and well integrity.",
        "recommended_range": "3,281-164,042 ft",
        "typical_value": "32,808 ft",
        "unit": "ft",
        "description": "Longer horizontal extents increase heat extraction area but require more drilling."
    },
    
    "Drilling Depth (m)": {
        "definition": "Set the depth of the hole drilling into the ground to access the geothermal reservoir. A value of 3.5 km targets mid-to-high enthalpy zones. The deeper the drill, the hotter the rock and higher the drilling cost.",
        "recommended_range": "1000-10000 m",
        "typical_value": "3500 m",
        "unit": "m",
        "description": "Deeper drilling accesses higher temperatures but increases costs significantly."
    },
    
    "Drilling Depth (ft)": {
        "definition": "Set the depth of the hole drilling into the ground to access the geothermal reservoir. A value of 11,483 ft targets mid-to-high enthalpy zones. The deeper the drill, the hotter the rock and higher the drilling cost.",
        "recommended_range": "3,281-32,808 ft",
        "typical_value": "11,483 ft",
        "unit": "ft",
        "description": "Deeper drilling accesses higher temperatures but increases costs significantly."
    },
    
    "Number of Laterals": {
        "definition": "Set the number of horizontal branches from the main well, with more laterals increasing reservoir contact. A value of 1 is the minimum number of laterals for a U-tube design.",
        "recommended_range": "1-20",
        "typical_value": "1",
        "unit": "count",
        "description": "Set the number of horizontal branches from the main well, with more laterals increasing reservoir contact. A value of 1 is the minimum number of laterals for a U-tube design."
    },
    
    "Lateral Flow Multiplier": {
        "definition": "Set the scaling factor for flow distribution in lateral sections.",
        "recommended_range": "0-1",
        "typical_value": "1",
        "unit": "multiplier",
        "description": "Scales the flow rate allocated to each lateral. A value of 1 means full flow allocation. The multiplier is automatically calculated as 1 divided by the number of laterals "
    },
    
    # Economic Parameters
    "Drilling Cost ($/m)": {
        "definition": "Set the cost per meter drilled. A value of $1,000/m represents an average cost.",
        "recommended_range": "0-4000 $/m",
        "typical_value": "1000 $/m",
        "unit": "$/m",
        "description": "A major component of geothermal project costs, varies with depth and geology."
    },
    
    "Drilling Cost ($/ft)": {
        "definition": "Set the cost per foot drilled. A value of $305/ft represents an average cost.",
        "recommended_range": "0-1,220 $/ft",
        "typical_value": "305 $/ft",
        "unit": "$/ft",
        "description": "A major component of geothermal project costs, varies with depth and geology."
    },
    
    # Base parameter name (without units) for pattern matching
    "Drilling Cost": {
        "definition": "Set the cost per meter drilled. A value of $1,000/m represents an average cost.",
        "recommended_range": "0-4000 $/m",
        "typical_value": "1000 $/m",
        "unit": "$/m",
        "description": "A major component of geothermal project costs, varies with depth and geology."
    },
    
  "Discount Rate (%)": {
        "definition": "Set the percentage that reflects how much future money is worth today, accounting for both the time and value of money and project risk. A 7% rate reflects moderate risk and cost of capital.",
        "recommended_range": "0-20%",
        "typical_value": "7%",
        "unit": "%",
        "description": "Higher rates make long-term projects less attractive economically."
    },
    
    "Lifetime (years)": {
        "definition": "Set the economic life of the project. A value of 40 years is common for geothermal projects with long-term resource stability.",
        "recommended_range": "10-40 years",
        "typical_value": "40 years",
        "unit": "years",
        "description": "Longer lifetimes improve project economics but increase uncertainty."
    },
    
    "Plant CAPEX ($/kWt)": {
        "definition": "Set the capital expenditure (CAPEX) per thermal kilowatt, which reflects the capital cost to build the surface plant for heat production, like heat exchangers, pumps, control systems, surface piping, and construction. A value of $100/kWt is a typical baseline for thermal-only systems, like district heating or industrial processes.",
        "recommended_range": "0-1000 $/kWt",
        "typical_value": "100 $/kWt",
        "unit": "$/kWt",
        "description": "Costs for heat exchangers, pumps, and other surface equipment."
    },
    
    "Plant CAPEX ($/kWe)": {
        "definition": "Set the capital expenditure (CAPEX) per electric kilowatt, which reflects the capital cost to build the surface plant for electricity production, like cooling towers, electrical systems, binary power cycle components, plant equipment, and construction. A value of $3,000/kWe is a typical baseline for flash or binary geothermal plants, making electricity for the grid or on-site industrial power.",
        "recommended_range": "0-10000 $/kWe",
        "typical_value": "3000 $/kWe",
        "unit": "$/kWe",
        "description": "Costs for turbines, generators, and power conversion equipment."
    },
    
    "Pre-cooling (˚C)": {
        "definition": "Set the temperature to which the working fluid is cooled before it's injected back underground. This temperature reflects the lowest temperature that can be consistently and economically achieved to help maximize the thermal gradient between the rock and the fluid. A value of 13°C is a realistic baseline for a system that includes ambient air cooling or mechanical chillers.",
        "recommended_range": "0-40°C",
        "typical_value": "13°C",
        "unit": "°C",
        "description": "Optimizes turbine efficiency and power output for sCO2 cycles."
    },
    
    "Pre-cooling (˚F)": {
        "definition": "Set the temperature to which the working fluid is cooled before it's injected back underground. This temperature reflects the lowest temperature that can be consistently and economically achieved to help maximize the thermal gradient between the rock and the fluid. A value of 55.4°F is a realistic baseline for a system that includes ambient air cooling or mechanical chillers.",
        "recommended_range": "32-104°F",
        "typical_value": "55.4°F",
        "unit": "°F",
        "description": "Optimizes turbine efficiency and power output for sCO2 cycles."
    },
    
    "Turbine Outlet Pressure (bar)": {
        "definition": "Set the pressure of the working fluid after it exits the turbine, determining how much energy can be extracted in the turbine and what condition (phase) the fluid is in for cooling and reinjection. A value of 80 bar keeps the working fluid in a dense supercritical or subcooled state.",
        "recommended_range": "75-200 bar",
        "typical_value": "80 bar",
        "unit": "bar",
        "description": "Critical parameter for sCO2 cycle efficiency and power output."
    },
    
    "Turbine Outlet Pressure (psi)": {
        "definition": "Set the pressure of the working fluid after it exits the turbine, determining how much energy can be extracted in the turbine and what condition (phase) the fluid is in for cooling and reinjection. A value of 1,160 psi keeps the working fluid in a dense supercritical or subcooled state.",
        "recommended_range": "1,088-2,901 psi",
        "typical_value": "1,160 psi",
        "unit": "psi",
        "description": "Critical parameter for sCO2 cycle efficiency and power output."
    },
    
    "Turbine Isentropic Efficiency (sCO2 electricity)": {
        "definition": "Efficiency of the turbine in converting thermal energy to mechanical work during the expansion process. The value of 90% represents a typical isentropic efficiency for sCO2 turbines, accounting for losses due to friction, heat transfer, and other irreversibilities.",
        "unit": "%",
        "description": "Higher efficiency increases net power output from the sCO2 cycle."
    },
    
    "Generator Efficiency (sCO2 electricity)": {
        "definition": "Efficiency of converting mechanical work from the turbine into electrical power. The value of 98% represents a typical generator efficiency, accounting for electrical losses in the conversion process.",
        "unit": "%",
        "description": "Higher efficiency increases net electrical power output from the system."
    },
    
    "Compressor Isentropic Efficiency (sCO2 electricity)": {
        "definition": "Efficiency of the compressor in pressurizing the working fluid with minimal entropy increase. The value of 90% represents a typical isentropic efficiency for sCO2 compressors, accounting for losses due to friction, heat transfer, and other irreversibilities.",
        "unit": "%",
        "description": "Higher efficiency reduces the work required to compress the fluid, increasing net power output."
    },
    
    # Model Fine-tuning Parameters
    "Mesh Fineness": {
        "definition": "Set the spatial resolution of the borehole based on how finely it is broken into discrete segments for simulation. A value of 0 (coarse) compared to 5 (research-grade accuracy), is the fastest option but least precise geometry.",
        "recommended_range": "0-5",
        "typical_value": "0",
        "unit": "index",
        "description": "Finer meshes provide more accurate results but require more computation time."
    },
    
    "Accuracy": {
        "definition": "Set the numerical accuracy level, from 1 (fastest) to 5 (most precise temperature and pressure).",
        "recommended_range": "1-5",
        "typical_value": "1",
        "unit": "index",
        "description": "Higher accuracy settings provide more precise results but increase computation time."
    },
    
    # Mode Parameters
    "Mass Flow Rate Mode": {
        "definition": "Set the profile of the fluid's mass that circulates through the geothermal system as either constant (steady operation) or variable (creates a time-dependent flow profile). Constant simplifies modeling.",
        "recommended_range": "Constant, Variable",
        "typical_value": "Constant",
        "unit": "mode",
        "description": "Constant mode provides steady-state operation for simplified modeling."
    },
    
    "Injection Temperature Mode": {
        "definition": "Set the profile of the fluid's temperature at the entrance of the subsurface as either constant (fixed temperature) or variable (creates a time-varying injection temperatures). Constant simplifies modeling.",
        "recommended_range": "Constant, Variable",
        "typical_value": "Constant",
        "unit": "mode",
        "description": "Constant mode provides fixed injection temperature for simplified modeling."
    },
    
    "Fluid Properties Mode": {
        "definition": "Set whether fluid properties are calculated at every time step as both a function of temperature and pressure (variable) or kept fixed (constant). Fluid properties include pressure values, temperature values, density, enthalpy, entropy, heat capacity, phase, thermal conductivity, thermal expansion, viscosity, and phase (liquid or vapor). Variable is recommended for high-fidelity modeling.",
        "recommended_range": "Constant, Variable",
        "typical_value": "Variable",
        "unit": "mode",
        "description": "Variable mode provides more accurate fluid property calculations."
    },
    
    # Additional Parameters
    "Pipe Roughness (m)": {
        "definition": "Set a measure of how smooth the inside surface of the pipe or borehole is, which affects how much friction the fluid encounters as it flows. A value of 1e-6 is very smooth – like new steel casting.",
        "recommended_range": "1e-6 to 3e-6 m",
        "typical_value": "0.000001 m",
        "unit": "m",
        "description": "Lower roughness values reduce friction losses and improve flow efficiency."
    },
    "Pipe Roughness (µm)": {
        "definition": "Set a measure of how smooth the inside surface of the pipe or borehole is, which affects how much friction the fluid encounters as it flows. A value of 1 µm (1e-6 m) is very smooth – like new steel casting.",
        "recommended_range": "1 to 3 µm (1e-6 to 3e-6 m)",
        "typical_value": "1 µm (0.000001 m)",
        "unit": "µm",
        "description": "Lower roughness values reduce friction losses and improve flow efficiency."
    },
    
    "Inlet Pressure (MPa)": {
        "definition": "Set the pressure of the fluid entering the system. A value of 10 MPa (or 100 bar) ensures flow in deep or high-resistant wells.",
        "recommended_range": "5-20 MPa",
        "typical_value": "10 MPa",
        "unit": "MPa",
        "description": "Higher inlet pressures can improve flow rates but increase pumping requirements."
    },
    
    "Wellbore Diameter (m)": {
        "definition": "Set the outer diameter of the coaxial borehole. A value of 0.444 m is a relatively large open-hole diameter for maximizing heat transfer surface area.",
        "recommended_range": "0.2159-0.4445 m",
        "typical_value": "0.458 m",
        "unit": "m",
        "description": "Affects the heat transfer area and flow resistance in coaxial systems."
    },
    
    "Center Pipe Diameter (m)": {
        "definition": "Set the inner diameter of the injection pipe. A value of 0.2 m (20 cm) is standard for high-velocity flow.",
        "recommended_range": "0.127-0.348 m",
        "typical_value": "0.2 m",
        "unit": "m",
        "description": "Affects the annular flow area and heat transfer characteristics in coaxial systems."
    },
    
    "Center Pipe Thickness (m)": {
        "definition": "Set the wall thickness of the inner pipe. A value of 0.013 m (13 mm) provides structural integrity and thermal isolation.",
        "recommended_range": "0.005-0.025 m",
        "typical_value": "0.013 m",
        "unit": "m",
        "description": "Thicker walls provide more structural strength but reduce flow area."
    },
    
    "Insulation Thermal Conductivity (W/m-°C)": {
        "definition": "Set the insulation performance of the inner pipe. A value of 0.025 W/m-K is a good approximation of a well-insulated pipe layer.",
        "recommended_range": "0.025-0.5 W/m-K",
        "typical_value": "0.025 W/m-K",
        "unit": "W/m-K",
        "description": "Lower conductivity values provide better thermal insulation."
    },
    
    "Coaxial Flow Type": {
        "definition": "Set the flow direction as either sending cold fluid down the outer pipe (\"Inject in Annulus\") or down the inner pipe (\"Inject in Center Pipe\").",
        "recommended_range": "Inject in Annulus, Inject in Center Pipe",
        "typical_value": "Inject in Annulus",
        "unit": "mode",
        "description": "Determines whether fluid is injected through the annular space or the center pipe."
    },
    
    "Lateral Flow Allocation": {
        "definition": "Set the flow allocation for each lateral branch of the well system.",
        "recommended_range": "0-1 per lateral",
        "typical_value": "Equal distribution (1/n for n laterals)",
        "unit": "fraction",
        "description": "Controls how flow is distributed across multiple lateral branches. For n laterals, the allocation values should sum to 1.0 (e.g., for 3 laterals: 1/3, 1/3, 1/3)."
    },
    
    "Mass Flow Rate Profile": {
        "definition": "User-provided Excel file containing mass flow rate profile over time.",
        "recommended_range": "Excel file with time and flow rate columns",
        "typical_value": "MassFlowRate.xlsx",
        "unit": "file",
        "description": "Required when Mass Flow Rate Mode is set to Variable. First column stores time in seconds, second column stores mass flow rate in kg/s."
    },
    
    "Injection Temperature Profile": {
        "definition": "User-provided Excel file containing injection temperature profile over time.",
        "recommended_range": "Excel file with time and temperature columns",
        "typical_value": "InjectionTemperatures.xlsx",
        "unit": "file",
        "description": "Required when Injection Temperature Mode is set to Variable. First column stores time in seconds, second column stores injection temperature in °C."
    },
    
    # CovHDF5 specific parameters
    "Permeability (HWR)": {
        "definition": "Set the horizontal-to-vertical permeability ratio (HWR) for the convection model. This dimensionless parameter controls how easily fluid flows horizontally versus vertically through the rock formation.",
        "recommended_range": "0.1-1.0",
        "typical_value": "0.5",
        "unit": "dimensionless",
        "description": "Higher values indicate more horizontal flow (anisotropic), while lower values indicate more vertical flow. A value of 1.0 represents isotropic permeability (equal horizontal and vertical flow)."
    },
    
    "CovHDF5 Mass Flow Rate (kg/s)": {
        "definition": "Set the total mass of working fluid that circulates through the geothermal system every second. A value of 6 kg/s moves enough fluid to extract significant heat but keeps pumping requirements and pressure losses manageable.",
        "recommended_range": "2-10 kg/s",
        "typical_value": "6 kg/s",
        "unit": "kg/s",
        "description": "Higher flow rates increase heat extraction but may require more pumping power."
    },
    
    "CovHDF5 Mass Flow Rate (lb/s)": {
        "definition": "Set the total mass of working fluid that circulates through the geothermal system every second. A value of 13.2 lb/s moves enough fluid to extract significant heat but keeps pumping requirements and pressure losses manageable.",
        "recommended_range": "4.4-22.0 lb/s",
        "typical_value": "13.2 lb/s",
        "unit": "lb/s",
        "description": "Higher flow rates increase heat extraction but may require more pumping power."
    },
    
    "CovHDF5 Horizontal Extent (m)": {
        "definition": "Set the horizontal distance of the well system. For CovHDF5, this represents the extent of the horizontal well section.",
        "recommended_range": "1,000-5,000 m",
        "typical_value": "2,500 m",
        "unit": "m",
        "description": "Longer horizontal sections increase heat extraction area but require more drilling."
    },
    
    "CovHDF5 Horizontal Extent (ft)": {
        "definition": "Set the horizontal distance of the well system. For CovHDF5, this represents the extent of the horizontal well section.",
        "recommended_range": "3,281-16,404 ft",
        "typical_value": "8,202 ft",
        "unit": "ft",
        "description": "Longer horizontal sections increase heat extraction area but require more drilling."
    },
    
    "CovHDF5 Drilling Depth (m)": {
        "definition": "Set the depth of the hole drilling into the ground to access the geothermal reservoir. A value of 3 km targets mid-to-high enthalpy zones. The deeper the drill, the hotter the rock and higher the drilling cost.",
        "recommended_range": "1,000-5,000 m",
        "typical_value": "3,000 m",
        "unit": "m",
        "description": "Deeper drilling accesses higher temperatures but increases costs significantly."
    },
    
    "CovHDF5 Drilling Depth (ft)": {
        "definition": "Set the depth of the hole drilling into the ground to access the geothermal reservoir. A value of 9,843 ft targets mid-to-high enthalpy zones. The deeper the drill, the hotter the rock and higher the drilling cost.",
        "recommended_range": "3,281-16,404 ft",
        "typical_value": "9,843 ft",
        "unit": "ft",
        "description": "Deeper drilling accesses higher temperatures but increases costs significantly."
    },
    
    "CovHDF5 Geothermal Gradient (°C/m)": {
        "definition": "Set the rate of temperature increase with depth. This represents how quickly the earth's temperature increases as you go deeper.",
        "recommended_range": "0.03-0.06°C/m",
        "typical_value": "0.045°C/m",
        "unit": "°C/m",
        "description": "Higher gradients provide access to higher temperatures at shallower depths, improving system efficiency."
    },
    
    "CovHDF5 Geothermal Gradient (°F/ft)": {
        "definition": "Set the rate of temperature increase with depth. This represents how quickly the earth's temperature increases as you go deeper.",
        "recommended_range": "0.016-0.033 °F/ft",
        "typical_value": "0.025 °F/ft",
        "unit": "°F/ft",
        "description": "Higher gradients provide access to higher temperatures at shallower depths, improving system efficiency."
    },
    
    "CovHDF5 Injection Temperature (°C)": {
        "definition": "Set the temperature of the working fluid when it enters the system. This is the temperature at which fluid is injected into the well.",
        "recommended_range": "30-60 °C",
        "typical_value": "45 °C",
        "unit": "°C",
        "description": "Lower injection temperatures generally improve heat extraction efficiency."
    },
    
    "CovHDF5 Injection Temperature (°F)": {
        "definition": "Set the temperature of the working fluid when it enters the system. This is the temperature at which fluid is injected into the well.",
        "recommended_range": "86-140 °F",
        "typical_value": "113 °F",
        "unit": "°F",
        "description": "Lower injection temperatures generally improve heat extraction efficiency."
    }
}

# Model-specific descriptions for Model Version popup
MODEL_DESCRIPTIONS = {
    "HDF5": {
        "title": "Database",
        "description": "Use pre-calculated database model for fast results from pre-computed scenarios. The database includes results for both water and supercritical CO2 working fluids."
    },
    "SBT V1.0": {
        "title": "Simulator - Slender-Body Theory V1.0",
        "description": "Fastest live model for simulating closed-loop geothermal scenarios and configurations that were not originally included in the pre-calculated database (Beckers et al., 2023). For example, with the Slender-Body Theory model, designs can be simulated for depths deeper than 5 km, geothermal gradients larger than 70°C/km, and number of laterals greater than 1, which were originally upper limits considered for the respective parameters when generating the database. This version supports water (H2O) only."
    },
    "SBT V2.0": {
        "title": "Simulator - Slender-Body Theory V2.0",
        "description": "Most comprehensive live model with support for both water and supercritical CO2 working fluids. This version extends the capabilities of V1.0 to include supercritical CO2 simulations."
    }
}

MODEL_LABELS = {
    "HDF5": "Database",
    "SBT V1.0": "Simulator",
    "SBT V2.0": "Simulator"
}

# Parameters whose titles naturally wrap to two lines and need extra icon offset
MULTILINE_INFO_PARAMS = {
    "Rock Thermal Conductivity (W/m-°C)",
    "Rock Thermal Conductivity (Btu/ft-h-˚F)",
    "Rock Specific Heat Capacity (J/kg-°C)",
    "Rock Specific Heat Capacity (Btu/lb-˚F)",
    "Insulation Thermal Conductivity (W/m-°C)",
}

MULTILINE_INFO_CUSTOM_OFFSETS = {
    "Rock Thermal Conductivity (W/m-°C)": "-10px",
    "Rock Thermal Conductivity (Btu/ft-h-˚F)": "-10px",
    "Rock Specific Heat Capacity (J/kg-°C)": "-10px",
    "Rock Specific Heat Capacity (Btu/lb-˚F)": "-10px",
    "Insulation Thermal Conductivity (W/m-°C)": "-15px",
}

# Parameters related to lateral configuration needing custom icon alignment
LATERAL_INFO_PARAMS = {
    "Number of Laterals",
    "Lateral Flow Multiplier",
    "Lateral Flow Allocation",
}

LATERAL_CUSTOM_TRANSFORMS = {
    "Lateral Flow Allocation": "translateX(-6px)",
}

def param_name_to_id_suffix(name: str) -> str:
    """
    Turn a PARAMETER_INFO key into a consistent id suffix.
    Example: "Drilling Depth (m)" -> "drilling-depth-m"
    """
    import re
    s = name.lower()
    s = s.replace('˚', 'deg').replace('°', 'deg')
    s = re.sub(r'[()$]', '', s)
    s = s.replace('/', '-')
    s = re.sub(r'\s+', '-', s)
    s = re.sub(r'-+', '-', s)
    s = s.strip('-')
    return s

def create_info_button(parameter_name, button_id=None):
    """Create an information button for a parameter."""
    if button_id is None:
        button_id = {
            "type": "info-btn",
            "param": param_name_to_id_suffix(parameter_name),
        }
    
    top_offset = "-3px"
    horizontal_transform = "translateX(-2px)"

    if parameter_name in MULTILINE_INFO_PARAMS:
        top_offset = MULTILINE_INFO_CUSTOM_OFFSETS.get(parameter_name, "-6px")
        # Only apply horizontal transform for Insulation Thermal Conductivity (needs more left adjustment)
        if parameter_name == "Insulation Thermal Conductivity (W/m-°C)":
            horizontal_transform = "translateX(-16px)"

    if parameter_name in LATERAL_INFO_PARAMS:
        top_offset = "1px"
        horizontal_transform = LATERAL_CUSTOM_TRANSFORMS.get(parameter_name, "translateX(-8px)")

    if parameter_name is None:
        top_offset = "-3px"

    return html.Div([
        dbc.Button(
            html.Img(src="/assets/info.svg", style={"width": "16px", "height": "16px"}),
            id=button_id,
            color="link",
            size="sm",
            className="ms-1 info-button",
            style={
                "textDecoration": "none",
                "padding": "0",
                "backgroundColor": "transparent",
                "border": "none",
                "display": "inline-flex",
                "alignItems": "center",
                "justifyContent": "center",
                "verticalAlign": "middle",
                "position": "relative",
                "top": top_offset,
                "transform": horizontal_transform
            }
        )
    ])



def create_enhanced_slider(DivID, ID, ptitle, min_v, max_v, mark_dict, start_v, div_style, parameter_name=None, step_i=None):
    """Create a slider with an information button."""
    info_button = create_info_button(parameter_name) if parameter_name else html.Div()
    
    slider_props = {
        "id": ID,
        "min": min_v,
        "max": max_v,
        "marks": mark_dict,
        "value": start_v,
        "tooltip": {"placement": "bottom", "always_visible": True}
    }
    
    if step_i is not None:
        slider_props["step"] = step_i
    
    return html.Div(id=DivID,
                    style=div_style,
                    children=[
                       html.Div(className="title-button-container", style={"display": "flex", "justifyContent": "flex-start", "alignItems": "center"}, children=[
                           html.P(ptitle, style={"fontWeight": "bold", "margin": 0}),
                           info_button
                       ]),
                       dcc.Slider(**slider_props),
                       ]
                    )

def create_enhanced_dropdown(DivID, ID, ptitle, options, disabled, div_style, parameter_name=None):
    """Create a dropdown with an information button."""
    info_button = create_info_button(parameter_name) if parameter_name else html.Div()
    
    if options and isinstance(options[0], dict):
        default_value = options[0]["value"]
    else:
        default_value = options[0] if options else None
    
    value = None if not options else default_value
    
    return html.Div(
            id=DivID,
            className="name-input-container-dd",
            style=div_style,
            children=[
                    html.Div(className="title-button-container", style={"display": "flex", "justifyContent": "flex-start", "alignItems": "center"}, children=[
                        html.P(ptitle, className="input-title", style={"margin": 0}),
                        info_button
                    ]),
                    dcc.Dropdown(
                            id=ID,
                            options=options,
                            value=value,
                            clearable=False,
                            searchable=False,
                            disabled=disabled,
                            className="select-dropdown"
                    ),
            ])

def create_enhanced_input_box(DivID, ID, ptitle, min_v, max_v, start_v, step_i, div_style, parameter_name=None):
    """Create an input box with an information button."""
    info_button = create_info_button(parameter_name) if parameter_name else html.Div()
    
    return html.Div(
            id=DivID,
            style=div_style,
            className="name-input-container",
            children=[
                html.Div(className="title-button-container", style={"display": "flex", "justifyContent": "flex-start", "alignItems": "center"}, children=[
                    html.P(ptitle, className="input-title", style={"margin": 0}),
                    info_button
                ]),
                dcc.Input(id=ID, disabled=True,
                            value=start_v, type='number', min=min_v, max=max_v, step=step_i, className="input-box"),
        ])

def create_info_modal():
    """Create the info modal component with click count store to prevent automatic opening when switching models or units"""
    return html.Div([
        dcc.Store(id="info-btn-last-clicks", data={}),
        dbc.Modal([
            dbc.ModalHeader(dbc.ModalTitle(id="info-modal-title")),
            dbc.ModalBody(id="info-modal-body"),
            dbc.ModalFooter(
                dbc.Button("Close", id="close-info-modal", className="ms-auto", n_clicks=0, 
                          style={"cursor": "pointer", "fontWeight": "bold"})
            ),
        ], id="info-modal", is_open=False, size="lg")
    ])

def register_info_modal_callbacks(app):
    """Register info modal callbacks using pattern-matching IDs."""
    suffix_to_param = {
        param_name_to_id_suffix(p): p for p in PARAMETER_INFO.keys()
    }

    @app.callback(
        [Output("info-modal", "is_open"),
         Output("info-modal-title", "children"),
         Output("info-modal-body", "children"),
         Output("info-btn-last-clicks", "data")],
        [
            Input({"type": "info-btn", "param": ALL}, "n_clicks"),
        ],
        [
            State("info-btn-last-clicks", "data"),
            State("info-modal", "is_open"),
            State("model-select", "value"),
        ],
        prevent_initial_call=True,
    )
    def toggle_info_modal(clicks_list, last_clicks, is_open, selected_model):
        """Handle info popup clicks for all parameters dynamically using pattern matching."""
        if not clicks_list:
            raise PreventUpdate

        last_clicks = last_clicks or {}

        triggered = ctx.triggered_id
        if not triggered or not isinstance(triggered, dict) or triggered.get("type") != "info-btn":
            raise PreventUpdate

        suffix = triggered.get("param")
        param = suffix_to_param.get(suffix)
        if not param:
            raise PreventUpdate

        max_clicks = max((c or 0) for c in clicks_list) if clicks_list else 0
        
        if max_clicks == 0:
            raise PreventUpdate
        
        last_clicks = last_clicks.copy()
        last_clicks["_max"] = max_clicks
        last_clicks[suffix] = max_clicks

        info = PARAMETER_INFO.get(param)
        if not info:
            raise PreventUpdate

        if param == "Model Version":
            header_style = {"fontSize": "16px", "fontWeight": "bold", "marginTop": "15px", "marginBottom": "8px"}
            
            hdf5_info = MODEL_DESCRIPTIONS.get("HDF5", {})
            sbt1_info = MODEL_DESCRIPTIONS.get("SBT V1.0", {})
            sbt2_info = MODEL_DESCRIPTIONS.get("SBT V2.0", {})
            
            if selected_model == "HDF5":
                hdf5_style = {"fontWeight": "bold"}
                sbt1_style = {"fontWeight": "normal"}
                sbt2_style = {"fontWeight": "normal"}
            elif selected_model == "SBT V1.0":
                hdf5_style = {"fontWeight": "normal"}
                sbt1_style = {"fontWeight": "bold"}
                sbt2_style = {"fontWeight": "normal"}
            elif selected_model == "SBT V2.0":
                hdf5_style = {"fontWeight": "normal"}
                sbt1_style = {"fontWeight": "normal"}
                sbt2_style = {"fontWeight": "bold"}
            else:
                hdf5_style = {"fontWeight": "normal"}
                sbt1_style = {"fontWeight": "normal"}
                sbt2_style = {"fontWeight": "normal"}
            
            modal_content = [
                html.H6(hdf5_info.get("title", "Database"), className="text-primary", style=header_style),
                html.P(hdf5_info.get("description", ""), className="mb-3", style=hdf5_style),
                
                html.H6(sbt1_info.get("title", "Simulator - Slender-Body Theory V1.0"), className="text-primary", style=header_style),
                html.P(sbt1_info.get("description", ""), className="mb-3", style=sbt1_style),
                
                html.H6(sbt2_info.get("title", "Simulator - Slender-Body Theory V2.0"), className="text-primary", style=header_style),
                html.P(sbt2_info.get("description", ""), className="mb-3", style=sbt2_style),
                
                html.H6("Simulator Model Selection", className="text-primary", style=header_style),
                html.P([
                    "When \"Simulator\" is selected, the Slender-Body Theory version is automatically determined based on the working fluid selection:",
                    html.Br(),
                    html.Br(),
                    "• ",
                    html.Strong("H2O selected:"),
                    " Slender-Body Theory V1.0 is used",
                    html.Br(),
                    "• ",
                    html.Strong("sCO2 selected or both H2O and sCO2 selected:"),
                    " Slender-Body Theory V2.0 is used",
                ], className="mb-3"),
            ]
            
            model_label = MODEL_LABELS.get(selected_model, "Model Version") if selected_model else "Model Version"
            bold_title = html.Strong(model_label)
            return True, bold_title, modal_content, last_clicks

        # Standard handling for other parameters
        modal_content = []
        info_sections = [
            ("Definition", "definition"),
            ("Recommended Range", "recommended_range"),
            ("Typical Value", "typical_value"),
            ("Description", "description"),
        ]

        for label, key in info_sections:
            value = info.get(key)
            if value:
                modal_content.append(html.H6(f"{label}:", className="text-primary"))
                modal_content.append(html.P(value, className="mb-3"))

        if not modal_content:
            modal_content.append(html.P("No additional information available yet.", className="mb-3"))

        return True, f"Information: {param}", modal_content, last_clicks

    @app.callback(
        Output("info-modal", "is_open", allow_duplicate=True),
        Input("close-info-modal", "n_clicks"),
        prevent_initial_call=True
    )
    def close_modal(n_clicks):
        if n_clicks:
            return False
        raise PreventUpdate 