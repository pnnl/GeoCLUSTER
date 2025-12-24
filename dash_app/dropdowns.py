#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# web app and interactive graphics libraries 
from dash import dcc, html
from info_popups import create_info_button

# ---------------------------
# Define dropdown options.
# ---------------------------

# User-facing model options - these are the labels shown to users
# The values are internal model names, but we'll update them dynamically based on fluid
def get_model_options(fluid="All"):
    """
    Get model dropdown options based on current fluid selection.
    For Simulator, the value depends on fluid: H2O -> SBT V1.0, All/sCO2 -> SBT V2.0
    """
    simulator_value = "SBT V1.0" if fluid == "H2O" else "SBT V2.0"
    return [
        {"label": "Database", "value": "HDF5"},
        {"label": "Simulator", "value": simulator_value}
    ]

# Helper function to get internal model name based on user selection and fluid
def get_internal_model(user_model, fluid):
    """
    Maps user-facing model selection to internal model name.
    
    Args:
        user_model: "HDF5" (Database) or "Simulator"
        fluid: "All", "H2O", or "sCO2"
    
    Returns:
        Internal model name: "HDF5", "SBT V1.0", or "SBT V2.0"
    """
    if user_model == "HDF5":
        return "HDF5"
    elif user_model == "Simulator":
        # For Simulator, determine version based on fluid
        if fluid == "H2O":
            return "SBT V1.0"
        else:  # "All" or "sCO2"
            return "SBT V2.0"
    else:
        return "HDF5"  # Default fallback
interp_list = ["True", "False"]
case_list = ["utube", "coaxial"]
fluid_list = ["All", "H2O", "sCO2"]
end_use_list = ["All", "Heating", "Electricity"] 

param_list = ["Horizontal Extent (m)", "Vertical Extent (m)", "Geothermal Gradient (°C/m)", "Borehole Diameter (m)", 
                "Injection Temperature (˚C)", "Rock Thermal Conductivity (W/m-°C)"]

dropdown_list = ["Run Interpolation", "Model Version", "Heat-Exchanger", "Working Fluid", "End-Use"]

def dropdown_card():

    # -----------------------------------------------------------------------
    # A Div containing controls for graphs. Controls include the following:
    #
    #   Dropdowns for categories:
    #           interpolation, case, fluid
    #
    # -----------------------------------------------------------------------

    return html.Div(id="dropdown-container",
                    children=[

                        html.Div(id="dropdown-card0",
                            # style={'display': 'none'},
                            children=[
                                html.Div([
                                    html.P("Model Version", className="dropdown-text", style={"display": "inline-block", "margin": "0", "verticalAlign": "middle"}),
                                    create_info_button("Model Version")
                                ], style={"display": "flex", "alignItems": "center", "gap": "5px", "marginBottom": "0"}),
                                dcc.Dropdown(
                                    id="model-select",
                                    options=get_model_options(),  # Will be updated dynamically
                                    value="HDF5",  # Start with Database
                                    clearable=False,
                                    searchable=False
                                ),
                            ]
                        ),

                        html.Div(id="dropdown-card1",
                            style={'display': 'none'},
                            children=[
                                html.P("Run Interpolation", className="dropdown-text"),
                                dcc.Dropdown(
                                    id="interpolation-select",
                                    options=[{"label": i, "value": i} for i in interp_list],
                                    value=interp_list[0],
                                    clearable=False,
                                    searchable=False
                                ),
                            ]
                        ),
                        
                        html.Div(id="dropdown-card2",
                            children=[
                                html.Div([
                                    html.P("Heat-Exchanger", className="dropdown-text", style={"display": "inline-block", "margin": "0", "verticalAlign": "middle"}),
                                    create_info_button("Heat Exchanger")
                                ], style={"display": "flex", "alignItems": "center", "gap": "5px", "marginBottom": "0"}),
                                dcc.Dropdown(
                                        id="case-select",
                                        options=[{"label": i, "value": i} for i in case_list],
                                        value=case_list[0],
                                        clearable=False,
                                        searchable=False
                                    ),
                                ]
                            ),
                        html.Div(id="dropdown-card3",
                            children=[
                                html.Div([
                                    html.P("Working Fluid", className="dropdown-text", style={"display": "inline-block", "margin": "0", "verticalAlign": "middle"}),
                                    create_info_button("Working Fluid")
                                ], style={"display": "flex", "alignItems": "center", "gap": "5px", "marginBottom": "0"}),
                                dcc.Dropdown(
                                    id="fluid-select",
                                    options=[{"label": i, "value": i} for i in fluid_list],
                                    value=fluid_list[0],
                                    clearable=False,
                                    searchable=False
                                    )
                                ]
                            ),
                        html.Div(id="dropdown-card4",
                            children=[
                                html.Div([
                                    html.P("End-Use", className="dropdown-text", style={"display": "inline-block", "margin": "0", "verticalAlign": "middle"}),
                                    create_info_button("End-Use")
                                ], style={"display": "flex", "alignItems": "center", "gap": "5px", "marginBottom": "0"}),
                                dcc.Dropdown(
                                    id="end-use-select",
                                    options=[{"label": i, "value": i} for i in end_use_list],
                                    value=end_use_list[0],
                                    clearable=False,
                                    searchable=False
                                ),
                            ]
                        ),
                ]
            )


