#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# web app and interactive graphics libraries 
from dash import dcc, html
import dash_bootstrap_components as dbc

# ---------------------------
# Parameter Information Dictionary
# ---------------------------

PARAMETER_INFO = {
    # Geologic Properties
    "Surface Temperature (˚C)": {
        "definition": "The temperature at the Earth's surface in the geothermal field.",
        "recommended_range": "0-40°C",
        "typical_value": "20°C",
        "unit": "°C",
        "description": "Surface temperature affects the initial conditions for geothermal calculations and heat transfer modeling."
    },
    
    "Geothermal Gradient (K/m)": {
        "definition": "The rate of temperature increase with depth in the Earth's crust.",
        "recommended_range": "0.015-0.200 K/m",
        "typical_value": "0.065 K/m",
        "unit": "K/m",
        "description": "Higher gradients indicate more favorable geothermal conditions for energy extraction."
    },
    
    "Rock Thermal Conductivity (W/m-K)": {
        "definition": "A measure of how well rock conducts heat.",
        "recommended_range": "0.4-5.0 W/m-K",
        "typical_value": "2.83 W/m-K",
        "unit": "W/m-K",
        "description": "Higher conductivity improves heat transfer from the rock to the working fluid."
    },
    
    "Rock Specific Heat Capacity (J/kg-K)": {
        "definition": "The amount of heat required to raise the temperature of 1 kg of rock by 1°C.",
        "recommended_range": "500-2000 J/kg-K",
        "typical_value": "825 J/kg-K",
        "unit": "J/kg-K",
        "description": "Affects the thermal storage capacity of the rock formation."
    },
    
    "Rock Density (kg/m3)": {
        "definition": "The mass per unit volume of the rock formation.",
        "recommended_range": "1000-3500 kg/m³",
        "typical_value": "2875 kg/m³",
        "unit": "kg/m³",
        "description": "Density influences the thermal mass and heat storage capacity of the formation."
    },
    
    # Wellbore Operations
    "Injection Temperature (˚C)": {
        "definition": "The temperature of the working fluid when injected into the wellbore.",
        "recommended_range": "30-60°C",
        "typical_value": "30°C",
        "unit": "°C",
        "description": "Lower injection temperatures generally improve heat extraction efficiency."
    },
    
    "Mass Flow Rate (kg/s)": {
        "definition": "The rate at which working fluid flows through the system.",
        "recommended_range": "5-300 kg/s",
        "typical_value": "20 kg/s",
        "unit": "kg/s",
        "description": "Higher flow rates increase heat extraction but may require more pumping power."
    },
    
    # Tube Geometry
    "Borehole Diameter (m)": {
        "definition": "The diameter of the drilled borehole.",
        "recommended_range": "0.2159-0.4445 m",
        "typical_value": "0.2286 m",
        "unit": "m",
        "description": "Larger diameters allow for higher flow rates but increase drilling costs."
    },
    
    "Wellbore Radius Vertical (m)": {
        "definition": "The radius of the vertical section of the wellbore.",
        "recommended_range": "0.10795-0.22225 m",
        "typical_value": "0.10795 m",
        "unit": "m",
        "description": "Affects the heat transfer area and flow resistance in the vertical section."
    },
    
    "Wellbore Radius Lateral (m)": {
        "definition": "The radius of the horizontal lateral sections of the wellbore.",
        "recommended_range": "0.10795-0.22225 m",
        "typical_value": "0.10795 m",
        "unit": "m",
        "description": "Influences heat transfer and flow characteristics in the lateral sections."
    },
    
    "Horizontal Extent (m)": {
        "definition": "The total horizontal distance covered by the wellbore laterals.",
        "recommended_range": "1000-50000 m",
        "typical_value": "1000 m",
        "unit": "m",
        "description": "Longer horizontal extents increase heat extraction area but require more drilling."
    },
    
    "Drilling Depth (m)": {
        "definition": "The vertical depth to which the wellbore is drilled.",
        "recommended_range": "1000-10000 m",
        "typical_value": "2000 m",
        "unit": "m",
        "description": "Deeper drilling accesses higher temperatures but increases costs significantly."
    },
    
    "Number of Laterals": {
        "definition": "The number of horizontal lateral branches from the main wellbore.",
        "recommended_range": "1-20",
        "typical_value": "3",
        "unit": "count",
        "description": "More laterals increase heat extraction area but add complexity and cost."
    },
    
    "Lateral Flow Allocation": {
        "definition": "The fraction of total flow allocated to lateral sections.",
        "recommended_range": "0-1",
        "typical_value": "0.33",
        "unit": "fraction",
        "description": "Controls how flow is distributed between vertical and lateral sections."
    },
    
    "Lateral Flow Multiplier": {
        "definition": "A factor that modifies the flow rate in lateral sections.",
        "recommended_range": "0-1",
        "typical_value": "1",
        "unit": "multiplier",
        "description": "Used to optimize flow distribution for maximum heat extraction."
    },
    
    # Economic Parameters
    "Drilling Cost ($/m)": {
        "definition": "The cost per meter of drilling the wellbore.",
        "recommended_range": "0-4000 $/m",
        "typical_value": "1000 $/m",
        "unit": "$/m",
        "description": "A major component of geothermal project costs, varies with depth and geology."
    },
    
    "Discount Rate (%)": {
        "definition": "The rate used to discount future cash flows to present value.",
        "recommended_range": "0-20%",
        "typical_value": "7%",
        "unit": "%",
        "description": "Higher rates make long-term projects less attractive economically."
    },
    
    "Lifetime (years)": {
        "definition": "The expected operational lifetime of the geothermal system.",
        "recommended_range": "10-40 years",
        "typical_value": "30 years",
        "unit": "years",
        "description": "Longer lifetimes improve project economics but increase uncertainty."
    },
    
    "Plant CAPEX ($/kWt)": {
        "definition": "Capital expenditure per kilowatt thermal for the surface plant.",
        "recommended_range": "0-1000 $/kWt",
        "typical_value": "100 $/kWt",
        "unit": "$/kWt",
        "description": "Costs for heat exchangers, pumps, and other surface equipment."
    },
    
    "Plant CAPEX ($/kWe)": {
        "definition": "Capital expenditure per kilowatt electric for the power plant.",
        "recommended_range": "0-10000 $/kWe",
        "typical_value": "3000 $/kWe",
        "unit": "$/kWe",
        "description": "Costs for turbines, generators, and power conversion equipment."
    },
    
    "Pre-cooling (˚C)": {
        "definition": "Temperature reduction before the working fluid enters the turbine.",
        "recommended_range": "0-40°C",
        "typical_value": "13°C",
        "unit": "°C",
        "description": "Optimizes turbine efficiency and power output for sCO2 cycles."
    },
    
    "Turbine Outlet Pressure (bar)": {
        "definition": "The pressure at the turbine outlet in the power cycle.",
        "recommended_range": "75-200 bar",
        "typical_value": "80 bar",
        "unit": "bar",
        "description": "Critical parameter for sCO2 cycle efficiency and power output."
    },
    
    # Model Fine-tuning Parameters
    "Mesh Fineness": {
        "definition": "The resolution of the numerical mesh used in simulations.",
        "recommended_range": "0-2",
        "typical_value": "0",
        "unit": "index",
        "description": "Finer meshes provide more accurate results but require more computation time."
    },
    
    "Accuracy": {
        "definition": "The numerical accuracy setting for the simulation solver.",
        "recommended_range": "1-5",
        "typical_value": "1",
        "unit": "index",
        "description": "Higher accuracy settings provide more precise results but increase computation time."
    },
    
    # Mode Parameters
    "Mass Flow Rate Mode": {
        "definition": "Whether the mass flow rate is constant or variable throughout the system.",
        "recommended_range": "Constant, Variable",
        "typical_value": "Constant",
        "unit": "mode",
        "description": "Variable mode allows for more complex flow patterns and optimization."
    },
    
    "Injection Temperature Mode": {
        "definition": "Whether the injection temperature is constant or variable.",
        "recommended_range": "Constant, Variable",
        "typical_value": "Constant",
        "unit": "mode",
        "description": "Variable temperature mode can optimize heat transfer efficiency."
    },
    
    "Fluid Properties Mode": {
        "definition": "Whether fluid properties are calculated as constant or variable.",
        "recommended_range": "Variable, Constant",
        "typical_value": "Variable",
        "unit": "mode",
        "description": "Variable mode provides more accurate fluid property calculations."
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
                       html.Div([
                           html.P(ptitle, style={"fontWeight": "bold", "display": "inline"}),
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
                    html.Div([
                        html.P(ptitle, className="input-title", style={"display": "inline"}),
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
                html.Div([
                    html.P(ptitle, className="input-title", style={"display": "inline"}),
                    info_button
                ]),
                dcc.Input(id=ID, disabled=True,
                            value=start_v, type='number', min=min_v, max=max_v, step=step_i, className="input-box"),
        ]) 