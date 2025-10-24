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
        "definition": "Set the mass of rock per unit volume, which affects heat storage and fluid flow behavior in a geothermal system. For example, denser rock holds more heat and changes temperature more slowly. A value of 2750 kg/m³ represents typical rock density for geothermal formations.",
        "recommended_range": "1000-3500 kg/m³",
        "typical_value": "2750 kg/m³",
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
    "Lateral Flow Allocation": {
        "definition": "Set the distribution of flow among lateral branches in a multi-lateral well system.",
        "recommended_range": "0-1",
        "typical_value": "0.5",
        "unit": "fraction",
        "description": "Controls how flow is distributed between different lateral branches."
    },
    
    "Wellbore Radius (m)": {
        "definition": "Set the radius of the wellbore for coaxial systems.",
        "recommended_range": "0.1-0.5 m",
        "typical_value": "0.2 m",
        "unit": "m",
        "description": "Affects the flow area and heat transfer characteristics in coaxial systems."
    },
    
    "Center Pipe Radius (m)": {
        "definition": "Set the radius of the center pipe in coaxial well systems.",
        "recommended_range": "0.05-0.2 m",
        "typical_value": "0.1 m",
        "unit": "m",
        "description": "Determines the inner flow area and affects pressure losses in coaxial systems."
    },
    
    "Center Pipe Thickness (m)": {
        "definition": "Set the wall thickness of the center pipe in coaxial systems.",
        "recommended_range": "0.005-0.02 m",
        "typical_value": "0.01 m",
        "unit": "m",
        "description": "Affects structural integrity and thermal resistance in coaxial systems."
    },
    
    "Insulation Thermal Conductivity (W/m-K)": {
        "definition": "Set the thermal conductivity of insulation material used in coaxial systems.",
        "recommended_range": "0.02-0.1 W/m-K",
        "typical_value": "0.04 W/m-K",
        "unit": "W/m-K",
        "description": "Lower values provide better thermal insulation and reduce heat losses."
    },
    
    "Coaxial Flow Type": {
        "definition": "Set the flow direction in coaxial systems - whether fluid is injected through the annulus or center pipe.",
        "recommended_range": "Inject in Annulus, Inject in Center Pipe",
        "typical_value": "Inject in Annulus",
        "unit": "mode",
        "description": "Affects flow patterns and heat transfer characteristics in coaxial systems."
    },
    
    "Inlet Pressure (MPa)": {
        "definition": "Set the pressure at the inlet of the geothermal system.",
        "recommended_range": "5-50 MPa",
        "typical_value": "10 MPa",
        "unit": "MPa",
        "description": "Higher inlet pressures can improve flow rates but increase pumping requirements."
    },
    
    "Pipe Roughness (m)": {
        "definition": "Set the surface roughness of the pipe walls, which affects frictional pressure losses.",
        "recommended_range": "1e-6 to 1e-3 m",
        "typical_value": "1e-5 m",
        "unit": "m",
        "description": "Smoother pipes reduce pressure losses but may be more expensive to manufacture."
    },
    
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
                "transform": "translateX(-1px) translateY(-4px)",
                "position": "relative",
                "top": "-2px"
            }
        )
    ])

def create_info_modal():
    """Create the info modal component with timestamp store to prevent automatic opening when switching models or units"""
    return html.Div([
        dcc.Store(id="info-btn-last-ts", data=0),
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
         Output("info-btn-last-ts", "data")],
        [
            Input({"type": "info-btn", "param": ALL}, "n_clicks_timestamp"),
        ],
        [
            State("info-btn-last-ts", "data"),
            State("info-modal", "is_open"),
        ],
        prevent_initial_call=True,
    )
    def toggle_info_modal(ts_list, last_ts, is_open):
        """Handle info popup clicks for all parameters dynamically using pattern matching with timestamps."""
        # No buttons mounted or nothing ever clicked
        if not ts_list:
            raise PreventUpdate

        # Current max click time among mounted buttons
        current_max = max((t or 0) for t in ts_list)

        # If no new click since last time, do nothing
        if current_max <= (last_ts or 0):
            raise PreventUpdate

        # A *new* click happened; figure out which button
        triggered = ctx.triggered_id
        if not triggered or not isinstance(triggered, dict) or triggered.get("type") != "info-btn":
            raise PreventUpdate

        suffix = triggered.get("param")
        param = suffix_to_param.get(suffix)
        if not param:
            raise PreventUpdate

        info = PARAMETER_INFO.get(param)
        if not info:
            raise PreventUpdate

        modal_content = [
            html.H6("Definition:", className="text-primary"),
            html.P(info["definition"], className="mb-3"),
            html.H6("Recommended Range:", className="text-primary"),
            html.P(info["recommended_range"], className="mb-3"),
            html.H6("Typical Value:", className="text-primary"),
            html.P(f"{info['typical_value']}", className="mb-3"),
            html.H6("Description:", className="text-primary"),
            html.P(info["description"], className="mb-3"),
        ]
        return True, f"Information: {param}", modal_content, current_max

    @app.callback(
        Output("info-modal", "is_open", allow_duplicate=True),
        Input("close-info-modal", "n_clicks"),
        prevent_initial_call=True
    )
    def close_modal(n_clicks):
        if n_clicks:
            return False
        raise PreventUpdate
