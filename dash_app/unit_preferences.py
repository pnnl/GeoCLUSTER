#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Unit Preferences UI Components for GeoCLUSTER
Provides interactive components for users to select their preferred units
"""

from dash import dcc, html
import dash_bootstrap_components as dbc
from unit_conversions import get_available_units, get_unit_symbol, unit_converter

def create_unit_preferences_card():
    """Create the main unit preferences card"""
    
    return dbc.Card([
        dbc.CardHeader([
            html.H5("Unit Preferences", className="card-title"),
            html.P("Select your preferred units for different parameter types", className="card-subtitle text-muted")
        ]),
        dbc.CardBody([
            dbc.Row([
                dbc.Col([
                    create_temperature_unit_selector(),
                    create_length_unit_selector(),
                    create_mass_flow_unit_selector()
                ], width=6),
                dbc.Col([
                    create_thermal_conductivity_unit_selector(),
                    create_heat_capacity_unit_selector(),
                    create_density_unit_selector()
                ], width=6)
            ]),
            dbc.Row([
                dbc.Col([
                    create_pressure_unit_selector(),
                    create_geothermal_gradient_unit_selector()
                ], width=12)
            ], className="mt-3"),
            dbc.Row([
                dbc.Col([
                    html.Hr(),
                    html.P("Changes will apply to all sliders and displays", className="text-muted small"),
                    dbc.Button("Apply Preferences", id="apply-unit-preferences", 
                              color="primary", size="sm", className="mt-2")
                ], width=12)
            ])
        ])
    ], className="mb-3")

def create_temperature_unit_selector():
    """Create temperature unit selector"""
    available_units = get_available_units('temperature')
    
    return html.Div([
        html.Label("Temperature Units:", className="form-label fw-bold"),
        dcc.Dropdown(
            id="temperature-unit-selector",
            options=[{"label": f"{unit} ({get_unit_symbol(unit)})", "value": unit} 
                    for unit in available_units],
            value="C",
            clearable=False,
            className="mb-2"
        )
    ])

def create_length_unit_selector():
    """Create length unit selector"""
    available_units = get_available_units('length')
    
    return html.Div([
        html.Label("Length Units:", className="form-label fw-bold"),
        dcc.Dropdown(
            id="length-unit-selector",
            options=[{"label": f"{unit} ({get_unit_symbol(unit)})", "value": unit} 
                    for unit in available_units],
            value="m",
            clearable=False,
            className="mb-2"
        )
    ])

def create_mass_flow_unit_selector():
    """Create mass flow unit selector"""
    available_units = get_available_units('mass_flow')
    
    return html.Div([
        html.Label("Mass Flow Units:", className="form-label fw-bold"),
        dcc.Dropdown(
            id="mass-flow-unit-selector",
            options=[{"label": f"{unit} ({get_unit_symbol(unit)})", "value": unit} 
                    for unit in available_units],
            value="kg/s",
            clearable=False,
            className="mb-2"
        )
    ])

def create_thermal_conductivity_unit_selector():
    """Create thermal conductivity unit selector"""
    available_units = get_available_units('thermal_conductivity')
    
    return html.Div([
        html.Label("Thermal Conductivity Units:", className="form-label fw-bold"),
        dcc.Dropdown(
            id="thermal-conductivity-unit-selector",
            options=[{"label": f"{unit} ({get_unit_symbol(unit)})", "value": unit} 
                    for unit in available_units],
            value="W/m-K",
            clearable=False,
            className="mb-2"
        )
    ])

def create_heat_capacity_unit_selector():
    """Create heat capacity unit selector"""
    available_units = get_available_units('heat_capacity')
    
    return html.Div([
        html.Label("Heat Capacity Units:", className="form-label fw-bold"),
        dcc.Dropdown(
            id="heat-capacity-unit-selector",
            options=[{"label": f"{unit} ({get_unit_symbol(unit)})", "value": unit} 
                    for unit in available_units],
            value="J/kg-K",
            clearable=False,
            className="mb-2"
        )
    ])

def create_density_unit_selector():
    """Create density unit selector"""
    available_units = get_available_units('density')
    
    return html.Div([
        html.Label("Density Units:", className="form-label fw-bold"),
        dcc.Dropdown(
            id="density-unit-selector",
            options=[{"label": f"{unit} ({get_unit_symbol(unit)})", "value": unit} 
                    for unit in available_units],
            value="kg/m3",
            clearable=False,
            className="mb-2"
        )
    ])

def create_pressure_unit_selector():
    """Create pressure unit selector"""
    available_units = get_available_units('pressure')
    
    return html.Div([
        html.Label("Pressure Units:", className="form-label fw-bold"),
        dcc.Dropdown(
            id="pressure-unit-selector",
            options=[{"label": f"{unit} ({get_unit_symbol(unit)})", "value": unit} 
                    for unit in available_units],
            value="Pa",
            clearable=False,
            className="mb-2"
        )
    ])

def create_geothermal_gradient_unit_selector():
    """Create geothermal gradient unit selector"""
    available_units = get_available_units('geothermal_gradient')
    
    return html.Div([
        html.Label("Geothermal Gradient Units:", className="form-label fw-bold"),
        dcc.Dropdown(
            id="geothermal-gradient-unit-selector",
            options=[{"label": f"{unit} ({get_unit_symbol(unit)})", "value": unit} 
                    for unit in available_units],
            value="K/m",
            clearable=False,
            className="mb-2"
        )
    ])

def create_quick_unit_selector():
    """Create a compact unit selector for the main interface"""
    
    return html.Div(id="unit-selector-container",
                    children=[
                        html.Div(id="unit-selector-card0",
                            children=[
                                html.P("Units", className="dropdown-text"),
                                dcc.Dropdown(
                                    id="quick-unit-selector",
                                    options=[
                                        {"label": "Metric (SI)", "value": "metric"},
                                        {"label": "Imperial (US)", "value": "imperial"}
                                    ],
                                    value="metric",
                                    clearable=False,
                                    searchable=False
                                ),
                            ]
                        ),

                    ])

def get_unit_preferences_from_inputs(temperature_unit, length_unit, mass_flow_unit, 
                                   thermal_conductivity_unit, heat_capacity_unit, 
                                   density_unit, pressure_unit, geothermal_gradient_unit):
    """Extract unit preferences from input values"""
    return {
        'temperature': temperature_unit,
        'length': length_unit,
        'mass_flow': mass_flow_unit,
        'thermal_conductivity': thermal_conductivity_unit,
        'heat_capacity': heat_capacity_unit,
        'density': density_unit,
        'pressure': pressure_unit,
        'geothermal_gradient': geothermal_gradient_unit
    }

def apply_metric_units():
    """Apply metric (SI) units"""
    return {
        'temperature': 'C',
        'length': 'm',
        'mass_flow': 'kg/s',
        'thermal_conductivity': 'W/m-K',
        'heat_capacity': 'J/kg-K',
        'density': 'kg/m3',
        'pressure': 'Pa',
        'geothermal_gradient': 'K/m'
    }

def apply_imperial_units():
    """Apply imperial (US) units"""
    return {
        'temperature': 'F',
        'length': 'yd',
        'mass_flow': 'lb/s',
        'thermal_conductivity': 'Btu/yd-h-F',
        'heat_capacity': 'Btu/lb-F',
        'density': 'lb/yd3',
        'pressure': 'psi',
        'geothermal_gradient': 'F/yd'
    }

def update_unit_converter_preferences(preferences):
    """Update the global unit converter with new preferences"""
    unit_converter.update_preferences(preferences)
