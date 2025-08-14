#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# web app and interactive graphics libraries 
from dash import dcc, html
import dash_bootstrap_components as dbc

# ---------------------------
# Parameter Information Dictionary
# ---------------------------

PARAMETER_INFO = {
    # Heat Transfer and System Configuration
    "Heat Transfer Mode": {
        "definition": "Specify the way thermal energy moves due to a temperature difference. The main modes are conduction (through solids or contact) or convection (through moving fluids).",
        "recommended_range": "Conduction, Convection",
        "typical_value": "Conduction",
        "unit": "mode",
        "description": "Conduction is the primary mode for geothermal heat transfer through rock formations."
    },
    
    "Heat Exchanger": {
        "definition": "Specify the geometry of the system. \"U-Tube\" uses two separate wells or laterals; \"Coaxial\" uses concentric pipes for flow-in and flow-out. Coaxial is often the default in simple single-well simulations, but U-Tube systems often have better long-term heat extraction in closed-loop systems.",
        "recommended_range": "U-Tube, Coaxial",
        "typical_value": "U-Tube",
        "unit": "mode",
        "description": "U-Tube systems provide better long-term heat extraction in closed-loop geothermal systems."
    },
    
    "Working Fluid": {
        "definition": "Specify a liquid or gas that carries heat through the system. H₂O (water) is typical for geothermal due to availability; sCO₂ (supercritical carbon dioxide) is included for advanced systems with enhanced heat extraction potential.",
        "recommended_range": "H2O, sCO2",
        "typical_value": "H2O",
        "unit": "fluid",
        "description": "Water is the most common working fluid due to availability and favorable thermal properties."
    },
    
    # Geologic Properties
    "Surface Temperature (˚C)": {
        "definition": "Set the ground-level or ambient temperature. A value of 25°C is a typical average surface temperature in temperate regions during geothermal operation.",
        "recommended_range": "0-40°C",
        "typical_value": "25°C",
        "unit": "°C",
        "description": "Surface temperature affects the initial conditions for geothermal calculations and heat transfer modeling."
    },
    
    "Geothermal Gradient (K/m)": {
        "definition": "Set the rate at which temperature increases with depth. A value of 0.05 K/m means that the temperature increases by 50°C for every kilometer of depth. 50°C/km represents average conditions in continental crust and it is hot enough to run a small power plant or provide heating for buildings.",
        "recommended_range": "0.015-0.200 K/m",
        "typical_value": "0.05 K/m",
        "unit": "K/m",
        "description": "Higher gradients indicate more favorable geothermal conditions for energy extraction."
    },
    
    "Rock Thermal Conductivity (W/m-K)": {
        "definition": "Set how quickly heat moves through rock. A value of 3 W/m-K represents moderately conductive rock, such as granite.",
        "recommended_range": "0.4-5.0 W/m-K",
        "typical_value": "3 W/m-K",
        "unit": "W/m-K",
        "description": "Higher conductivity improves heat transfer from the rock to the working fluid."
    },
    
    "Rock Specific Heat Capacity (J/kg-K)": {
        "definition": "Set the amount of energy the rock can absorb or release when its temperature changes by 1°C, which determines how quickly the rock heats up or cools down in response to fluid circulation. A value of 0.051 J/kg-K represents an average for various dry rocks.",
        "recommended_range": "500-2000 J/kg-K",
        "typical_value": "0.051 J/kg-K",
        "unit": "J/kg-K",
        "description": "Affects the thermal storage capacity of the rock formation."
    },
    
    "Rock Density (kg/m3)": {
        "definition": "Set the mass of rock per unit volume, which affects heat storage and fluid flow behavior in a geothermal system. For example, denser rock holds more heat and changes temperature more slowly. A value of 790 kg/m³ implies the rock is highly porous, fractured, or contains gas-filled voids. Typical rock densities usually range from 2,000 to 3,000 kg/m³.",
        "recommended_range": "1000-3500 kg/m³",
        "typical_value": "790 kg/m³",
        "unit": "kg/m³",
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
    
    "Mass Flow Rate (kg/s)": {
        "definition": "Set the total mass of working fluid that circulates through the geothermal system every second. A value of 24 kg/s moves enough fluid to extract significant heat but keeps pumping requirements and pressure losses manageable.",
        "recommended_range": "5-300 kg/s",
        "typical_value": "24 kg/s",
        "unit": "kg/s",
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
    
    "Wellbore Radius Vertical (m)": {
        "definition": "Set the radius of the vertical injection and production well of the U-tube design. A value of 0.222 m is a relatively large open-hole radius for maximizing heat transfer surface area.",
        "recommended_range": "0.10795-0.22225 m",
        "typical_value": "0.222 m",
        "unit": "m",
        "description": "Affects the heat transfer area and flow resistance in the vertical section."
    },
    
    "Wellbore Radius Lateral (m)": {
        "definition": "Set the radius of the lateral branches of the U-tube design. A value of 0.222 m is a relatively large open-hole radius for maximizing heat transfer surface area.",
        "recommended_range": "0.10795-0.22225 m",
        "typical_value": "0.222 m",
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
    
    "Drilling Depth (m)": {
        "definition": "Set the depth of the hole drilling into the ground to access the geothermal reservoir. A value of 3.5 km targets mid-to-high enthalpy zones. The deeper the drill, the hotter the rock and higher the drilling cost.",
        "recommended_range": "1000-10000 m",
        "typical_value": "3500 m",
        "unit": "m",
        "description": "Deeper drilling accesses higher temperatures but increases costs significantly."
    },
    
    "Number of Laterals": {
        "definition": "Set the number of horizontal branches from the main well, with more laterals increasing reservoir contact. A value of 1 is the minimum number of laterals for a U-tube design.",
        "recommended_range": "1-20",
        "typical_value": "1",
        "unit": "count",
        "description": "More laterals increase heat extraction area but add complexity and cost."
    },
    
    "Lateral Flow Multiplier": {
        "definition": "Set the scaling factor for flow distribution in lateral sections.",
        "recommended_range": "0-1",
        "typical_value": "1",
        "unit": "multiplier",
        "description": "Used to optimize flow distribution for maximum heat extraction."
    },
    
    # Economic Parameters
    "Drilling Cost ($/m)": {
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
    
    "Turbine Outlet Pressure (bar)": {
        "definition": "Set the pressure of the working fluid after it exits the turbine, determining how much energy can be extracted in the turbine and what condition (phase) the fluid is in for cooling and reinjection. A value of 80 bar keeps the working fluid in a dense supercritical or subcooled state.",
        "recommended_range": "75-200 bar",
        "typical_value": "80 bar",
        "unit": "bar",
        "description": "Critical parameter for sCO2 cycle efficiency and power output."
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
    
    "Inlet Pressure (MPa)": {
        "definition": "Set the pressure of the fluid entering the system. A value of 10 MPa (or 100 bar) ensures flow in deep or high-resistant wells.",
        "recommended_range": "5-20 MPa",
        "typical_value": "10 MPa",
        "unit": "MPa",
        "description": "Higher inlet pressures can improve flow rates but increase pumping requirements."
    },
    
    "Wellbore Radius (m)": {
        "definition": "Set the outer radius of the coaxial borehole. A value of 0.222 m is a relatively large open-hole radius for maximizing heat transfer surface area.",
        "recommended_range": "0.10795-0.22225 m",
        "typical_value": "0.229 m",
        "unit": "m",
        "description": "Affects the heat transfer area and flow resistance in coaxial systems."
    },
    
    "Center Pipe Radius (m)": {
        "definition": "Set the inner radius of the injection pipe. A value of 0.010 m (1 cm) is standard for high-velocity flow.",
        "recommended_range": "0.0635-0.174 m",
        "typical_value": "0.010 m",
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
    
    "Insulation Thermal Conductivity (W/m-K)": {
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
    
    # Legacy parameters (keeping for backward compatibility)
    "Lateral Flow Allocation": {
        "definition": "The fraction of total flow allocated to lateral sections.",
        "recommended_range": "0-1",
        "typical_value": "0.33",
        "unit": "fraction",
        "description": "Controls how flow is distributed between vertical and lateral sections."
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
    }
}

# ---------------------------
# Helper Functions
# ---------------------------

def create_info_button(parameter_name, button_id=None):
    """
    Create an information button for a parameter.
    
    Args:
        parameter_name (str): The name of the parameter
        button_id (str): Optional custom button ID
        
    Returns:
        html.Div: A button component with info icon
    """
    if button_id is None:
        # Create a standardized button ID based on parameter name
        # First replace special characters, then spaces, then clean up multiple dashes
        button_id = parameter_name.lower()
        button_id = button_id.replace('˚', 'deg').replace('°', 'deg')
        button_id = button_id.replace('(', '').replace(')', '')
        button_id = button_id.replace('/', '-').replace('$', '')
        button_id = button_id.replace(' ', '-')
        # Clean up multiple consecutive dashes
        import re
        button_id = re.sub(r'-+', '-', button_id)
        button_id = f"info-btn-{button_id}"
    
    
    return html.Div([
        dbc.Button(
            "ℹ",
            id=button_id,
            color="link",
            size="sm",
            className="ms-1 info-button",
            style={
                "fontSize": "10px",
                "textDecoration": "none",
                "padding": "0",
                "color": "#17a2b8",
                "borderRadius": "50%",
                "border": "1px solid #17a2b8",
                "backgroundColor": "transparent",
                "width": "16px",
                "height": "16px",
                "display": "inline-flex",
                "alignItems": "center",
                "justifyContent": "center",
                "lineHeight": "1",
                "textAlign": "center",
                "transform": "translateX(-1px) translateY(-2px)", # Adjusted for centering and vertical position
                "position": "relative",
                "top": "-2px" # Adjusted for vertical position
            }
        )
    ])



def create_enhanced_slider(DivID, ID, ptitle, min_v, max_v, mark_dict, start_v, div_style, parameter_name=None, step_i=None):
    """
    Create a slider with an information button.
    
    Args:
        DivID (str): The div ID
        ID (str): The slider ID
        ptitle (str): The parameter title
        min_v (float): Minimum value
        max_v (float): Maximum value
        mark_dict (dict): Marks dictionary
        start_v (float): Starting value
        div_style (dict): Style dictionary
        parameter_name (str): Parameter name for info popup
        step_i (float): Step increment (for slider1 type)
        
    Returns:
        html.Div: Enhanced slider component with info button
    """
    info_button = create_info_button(parameter_name) if parameter_name else html.Div()
    
    slider_props = {
        "id": ID,
        "min": min_v,
        "max": max_v,
        "marks": mark_dict,
        "value": start_v,
        "tooltip": {"placement": "bottom", "always_visible": True}
    }
    
    # Add step if provided (for slider1 type)
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
    """
    Create a dropdown with an information button.
    
    Args:
        DivID (str): The div ID
        ID (str): The dropdown ID
        ptitle (str): The parameter title
        options (list): Dropdown options
        disabled (bool): Whether dropdown is disabled
        div_style (dict): Style dictionary
        parameter_name (str): Parameter name for info popup
        
    Returns:
        html.Div: Enhanced dropdown component with info button
    """
    info_button = create_info_button(parameter_name) if parameter_name else html.Div()
    
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
                            value=options[0],
                            clearable=False,
                            searchable=False,
                            disabled=disabled,
                            className="select-dropdown"
                    ),
            ])

def create_enhanced_input_box(DivID, ID, ptitle, min_v, max_v, start_v, step_i, div_style, parameter_name=None):
    """
    Create an input box with an information button.
    
    Args:
        DivID (str): The ID for the container div
        ID (str): The ID for the input component
        ptitle (str): The title/label for the input
        min_v (float): Minimum value
        max_v (float): Maximum value
        start_v (float): Starting value
        step_i (float): Step increment
        div_style (dict): Style for the container div
        parameter_name (str): The name of the parameter for info popup
        
    Returns:
        html.Div: An input box component with info button
    """
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