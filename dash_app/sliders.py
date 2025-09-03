#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# web app and interactive graphics libraries 
from dash import dcc, html
import dash_bootstrap_components as dbc

# sourced scripts
from plots import * # u_sCO2, u_H2O, c_sCO2, c_H2O
from info_popups import create_enhanced_slider, create_enhanced_dropdown, create_enhanced_input_box
from unit_conversions import unit_converter, get_unit_symbol

# Conversion helper functions
def get_temperature_converted_values(base_value_celsius, target_unit):
    """Convert temperature values from Celsius to target unit"""
    if target_unit == 'C':
        return base_value_celsius
    elif target_unit == 'F':
        converted = unit_converter.convert_temperature(base_value_celsius, 'C', 'F')
        print(f"DEBUG: Converting {base_value_celsius}°C to {converted}°F")
        return converted
    else:
        return base_value_celsius

def get_length_converted_values(base_value_meters, target_unit):
    """Convert length values from meters to target unit"""
    if target_unit == 'm':
        return base_value_meters
    elif target_unit == 'ft':
        converted = unit_converter.convert_length(base_value_meters, 'm', 'ft')
        print(f"DEBUG: Converting {base_value_meters}m to {converted}ft")
        return converted
    elif target_unit == 'yd':
        converted = unit_converter.convert_length(base_value_meters, 'm', 'yd')
        print(f"DEBUG: Converting {base_value_meters}m to {converted}yd")
        return converted
    else:
        return base_value_meters

def get_mass_flow_converted_values(base_value_kg_s, target_unit):
    """Convert mass flow from kg/s to target unit"""
    if target_unit == 'kg/s':
        return base_value_kg_s
    elif target_unit == 'lb/s':
        converted = unit_converter.convert_mass_flow(base_value_kg_s, 'kg/s', 'lb/s')
        print(f"DEBUG: Converting {base_value_kg_s} kg/s to {converted} lb/s")
        return converted
    else:
        return base_value_kg_s

def get_thermal_conductivity_converted_values(base_value_w_m_k, target_unit):
    """Convert thermal conductivity from W/m-K to target unit"""
    if target_unit == 'W/m-K':
        return base_value_w_m_k
    elif target_unit == 'Btu/ft-h-F':
        converted = unit_converter.convert_thermal_conductivity(base_value_w_m_k, 'W/m-K', 'Btu/ft-h-F')
        print(f"DEBUG: Converting {base_value_w_m_k} W/m-K to {converted} Btu/ft-h-F")
        return converted
    elif target_unit == 'Btu/yd-h-F':
        converted = unit_converter.convert_thermal_conductivity(base_value_w_m_k, 'W/m-K', 'Btu/yd-h-F')
        print(f"DEBUG: Converting {base_value_w_m_k} W/m-K to {converted} Btu/yd-h-F")
        return converted
    else:
        return base_value_w_m_k

def get_heat_capacity_converted_values(base_value_j_kg_k, target_unit):
    """Convert heat capacity from J/kg-K to target unit"""
    if target_unit == 'J/kg-K':
        return base_value_j_kg_k
    elif target_unit == 'Btu/lb-F':
        converted = unit_converter.convert_heat_capacity(base_value_j_kg_k, 'J/kg-K', 'Btu/lb-F')
        print(f"DEBUG: Converting {base_value_j_kg_k} J/kg-K to {converted} Btu/lb-F")
        return converted
    else:
        return base_value_j_kg_k

def get_density_converted_values(base_value_kg_m3, target_unit):
    """Convert density from kg/m³ to target unit"""
    if target_unit == 'kg/m3':
        return base_value_kg_m3
    elif target_unit == 'lb/ft3':
        converted = unit_converter.convert_density(base_value_kg_m3, 'kg/m3', 'lb/ft3')
        print(f"DEBUG: Converting {base_value_kg_m3} kg/m³ to {converted} lb/ft³")
        return converted
    elif target_unit == 'lb/yd3':
        converted = unit_converter.convert_density(base_value_kg_m3, 'kg/m3', 'lb/yd3')
        print(f"DEBUG: Converting {base_value_kg_m3} kg/m³ to {converted} lb/yd³")
        return converted
    else:
        return base_value_kg_m3

def get_pressure_converted_values(base_value_pa, target_unit):
    """Convert pressure from Pa to target unit"""
    if target_unit == 'Pa':
        return base_value_pa
    elif target_unit == 'psi':
        converted = unit_converter.convert_pressure(base_value_pa, 'Pa', 'psi')
        print(f"DEBUG: Converting {base_value_pa} Pa to {converted} psi")
        return converted
    else:
        return base_value_pa

def get_geothermal_gradient_converted_values(base_value_k_m, target_unit):
    """Convert geothermal gradient from K/m to target unit"""
    if target_unit == 'K/m':
        return base_value_k_m
    elif target_unit == 'F/ft':
        converted = unit_converter.convert_geothermal_gradient(base_value_k_m, 'K/m', 'F/ft')
        print(f"DEBUG: Converting {base_value_k_m} K/m to {converted} F/ft")
        return converted
    elif target_unit == 'F/yd':
        converted = unit_converter.convert_geothermal_gradient(base_value_k_m, 'K/m', 'F/yd')
        print(f"DEBUG: Converting {base_value_k_m} K/m to {converted} F/yd")
        return converted
    else:
        return base_value_k_m

def get_drilling_cost_converted_values(base_value_per_meter, target_unit):
    """Convert drilling cost from $/m to $/ft or $/yd"""
    if target_unit == 'm':
        return base_value_per_meter
    elif target_unit == 'ft':
        converted = base_value_per_meter / 3.28084
        return converted
    elif target_unit == 'yd':
        converted = base_value_per_meter / 1.09361
        return converted
    else:
        return base_value_per_meter

def create_imperial_marks(min_val, max_val, unit):
    """Create min/max marks for imperial units"""
    
    if unit == 'ft':
        # For feet, show min and max with ft suffix
        # Round to 2 decimal places for feet to avoid showing too many decimals
        min_rounded = round(min_val, 2)
        max_rounded = round(max_val, 2)
        result = {min_rounded: f'{min_rounded:.2f} ft', max_rounded: f'{max_rounded:.2f} ft'}
        return result
    elif unit == 'yd':
        min_rounded = round(min_val, 2)
        max_rounded = round(max_val, 2)
        result = {min_rounded: f'{min_rounded:.2f} yd', max_rounded: f'{max_rounded:.2f} yd'}
        return result
    elif unit == 'F':
        # For Fahrenheit, show min and max with °F suffix
        return {min_val: f'{min_val:.0f}°F', max_val: f'{max_val:.0f}°F'}
    elif unit == 'psi':
        # For psi, show min and max with psi suffix
        return {min_val: f'{min_val:.0f} psi', max_val: f'{max_val:.0f} psi'}
    elif unit == 'lb/s':
        # For lb/s, show min and max with lb/s suffix
        return {min_val: f'{min_val:.1f} lb/s', max_val: f'{max_val:.1f} lb/s'}
    elif unit == 'Btu/ft-h-F':
        # For thermal conductivity, show min and max
        return {min_val: f'{min_val:.1f}', max_val: f'{max_val:.1f}'}
    elif unit == 'Btu/yd-h-F':
        return {min_val: f'{min_val:.1f}', max_val: f'{max_val:.1f}'}
    elif unit == 'Btu/lb-F':
        # For heat capacity, show min and max
        return {min_val: f'{min_val:.0f}', max_val: f'{max_val:.0f}'}
    elif unit == 'lb/ft3':
        # For density, show min and max
        return {min_val: f'{min_val:.0f}', max_val: f'{max_val:.0f}'}
    elif unit == 'lb/yd3':
        return {min_val: f'{min_val:.0f}', max_val: f'{max_val:.0f}'}
    elif unit == 'F/ft':
        # For geothermal gradient, show min and max
        return {min_val: f'{min_val:.2f}', max_val: f'{max_val:.2f}'}
    elif unit == 'F/yd':
        return {min_val: f'{min_val:.2f}', max_val: f'{max_val:.2f}'}
    elif unit == 'drillcost_ft':  # Special case for drilling cost in $/ft
        # For drilling cost in $/ft, show min and max with $/ft suffix
        min_rounded = round(min_val, 0)
        max_rounded = round(max_val, 0)
        return {min_rounded: f'${min_rounded:.0f}/ft', max_rounded: f'${max_rounded:.0f}/ft'}
    elif unit == 'drillcost_yd':  # Special case for drilling cost in $/yd
        min_rounded = round(min_val, 0)
        max_rounded = round(max_val, 0)
        return {min_rounded: f'${min_rounded:.0f}/yd', max_rounded: f'${max_rounded:.0f}/yd'}
    elif unit == 'ft':  # Regular length unit in feet
        # For feet, show min and max with ft suffix
        # Round to 2 decimal places for feet to avoid showing too many decimals
        min_rounded = round(min_val, 2)
        max_rounded = round(max_val, 2)
        result = {min_rounded: f'{min_rounded:.2f} ft', max_rounded: f'{max_rounded:.2f} ft'}
        return result
    elif unit == 'yd':  # Regular length unit in yards
        min_rounded = round(min_val, 2)
        max_rounded = round(max_val, 2)
        result = {min_rounded: f'{min_rounded:.2f} yd', max_rounded: f'{max_rounded:.2f} yd'}
        return result
    else:
        # Default case - just show min and max
        return {min_val: str(min_val), max_val: str(max_val)}

# -------------------------------------------------------------------------------------------
# Create dictionaries to assign the values and boundaries of the parameter siders.
#
#   mdot - mass flow rates [kg/s]
#   L2 - horizontal extent [m]
#   L1 - drilling depth [m] 
#   grad - thermal gradient [K/m] 
#   D - borehole diameter [m]
#   Tinj - injection temperature [K]
#   k - rock thermal conductivity [W/m-K]
#
# -------------------------------------------------------------------------------------------


def create_steps(arg_arr, str_round_place, val_round_place):
    """
    OUTPUT: Returns a newly annnotated dictionary for a slider component.
    """
    
    value_list = [round(elem, val_round_place) for elem in arg_arr.tolist()]

    # print(value_list)

    if value_list[-1] >= 1000:
        arg_arr = arg_arr / 1000
        value_list = [int(item) for item in value_list]

     # print(value_list)

    # exception
    if value_list[0] > 2 and value_list[0] < 100:
        value_list[0] = int(value_list[0])
        value_list[-1] = int(value_list[-1])

    list_str = [str_round_place.format(x) for x in arg_arr.tolist()]

    for i in range(len(list_str)):
        if value_list[i] >= 1000:
            list_str[i] = list_str[i] + "k"
        if i >= 1 and i < len(list_str) - 1:
            list_str[i] = ""

    slider_dict = dict(zip(value_list, list_str))

    return slider_dict

def get_dynamic_slider_labels():
    """Get slider labels with current unit preferences"""
    temp_unit = get_unit_symbol(unit_converter.user_preferences.get('temperature', 'C'))
    length_unit = get_unit_symbol(unit_converter.user_preferences.get('length', 'm'))
    mass_flow_unit = get_unit_symbol(unit_converter.user_preferences.get('mass_flow', 'kg/s'))
    thermal_cond_unit = get_unit_symbol(unit_converter.user_preferences.get('thermal_conductivity', 'W/m-K'))
    heat_cap_unit = get_unit_symbol(unit_converter.user_preferences.get('heat_capacity', 'J/kg-K'))
    density_unit = get_unit_symbol(unit_converter.user_preferences.get('density', 'kg/m3'))
    pressure_unit = get_unit_symbol(unit_converter.user_preferences.get('pressure', 'Pa'))
    grad_unit = get_unit_symbol(unit_converter.user_preferences.get('geothermal_gradient', 'K/m'))
    
    return {
        'wellbore_operations': [
            f"Injection Temperature ({temp_unit})", 
            f"Mass Flow Rate ({mass_flow_unit})", 
            f"Horizontal Extent ({length_unit})", 
            f"Drilling Depth ({length_unit})"
        ],
        'tube_geometry': [
            f"Wellbore Radius Vertical ({length_unit})", 
            f"Wellbore Radius Lateral ({length_unit})",
            "Number of Laterals", 
            "Lateral Flow Allocation", 
            "Lateral Flow Multiplier"
        ],
        'economic_params': [
            f"Drilling Cost ($/{get_unit_symbol(unit_converter.user_preferences.get('length', 'm'))})", 
            "Discount Rate (%)", 
            "Lifetime (years)", 
            "Plant CAPEX ($/kWt)", 
            "Plant CAPEX ($/kWe)", 
            f"Pre-cooling ({temp_unit})", 
            f"Turbine Outlet Pressure ({pressure_unit})"
        ],
        'geologic_properties': [
            f"Surface Temperature ({temp_unit})", 
            f"Geothermal Gradient ({grad_unit})", 
            f"Rock Thermal Conductivity ({thermal_cond_unit})", 
            f"Rock Specific Heat Capacity ({heat_cap_unit})", 
            f"Rock Density ({density_unit})"
        ]
    }

# Initialize with default units
wellbore_operations_l = ["Injection Temperature (˚C)", "Mass Flow Rate (kg/s)", "Horizontal Extent (m)", "Drilling Depth (m)"]
tube_geometry_l = ["Wellbore Radius Vertical (m)", "Wellbore Radius Lateral (m)", "Number of Laterals", "Lateral Flow Allocation", "Lateral Flow Multiplier"]
economic_params_l = ["Drilling Cost ($/m)", "Discount Rate (%)", "Lifetime (years)", "Plant CAPEX ($/kWt)", "Plant CAPEX ($/kWe)", "Pre-cooling (˚C)", "Turbine Outlet Pressure (bar)"]
geologic_properties_l = ["Surface Temperature (˚C)", "Geothermal Gradient (K/m)", "Rock Thermal Conductivity (W/m-K)", "Rock Specific Heat Capacity (J/kg-K)", "Rock Density (kg/m3)"]

def update_slider_labels():
    """Update slider labels when units change"""
    global wellbore_operations_l, tube_geometry_l, economic_params_l, geologic_properties_l
    slider_labels = get_dynamic_slider_labels()
    wellbore_operations_l = slider_labels['wellbore_operations']
    tube_geometry_l = slider_labels['tube_geometry']
    economic_params_l = slider_labels['economic_params']
    geologic_properties_l = slider_labels['geologic_properties']
    return slider_labels
                            
model_finetuning_l = ["Mesh Fineness", "Accuracy", "Mass Flow Rate Mode", "Mass Flow Rate Profile", 
                                            "Injection Temperature Mode", "Injection Temperature Profile"]

slider_list = wellbore_operations_l + tube_geometry_l + economic_params_l + geologic_properties_l + model_finetuning_l

mdot_dict = create_steps(arg_arr=u_sCO2.mdot, str_round_place='{:.1f}', val_round_place=1)
L2_dict = create_steps(arg_arr=u_sCO2.L2, str_round_place='{:.0f}', val_round_place=0)
L1_dict = create_steps(arg_arr=u_sCO2.L1, str_round_place='{:.0f}', val_round_place=0)
grad_dict = create_steps(arg_arr=u_sCO2.grad, str_round_place='{:.2f}', val_round_place=2)
D_dict = create_steps(arg_arr=u_sCO2.D, str_round_place='{:.4f}', val_round_place=4)
Tinj_dict = create_steps(arg_arr=u_sCO2.Tinj - 273.15, str_round_place='{:.1f}', val_round_place=2)
k_dict = create_steps(arg_arr=u_sCO2.k, str_round_place='{:.1f}', val_round_place=1)

# slider labels for economic paramters
drillcost_dict = create_steps(arg_arr=np.arange(start=0, stop=4100, step=200), str_round_place='{:.0f}', val_round_place=0)
discount_dict = create_steps(arg_arr=np.arange(start=0, stop=20.1, step=1), str_round_place='{:.1f}', val_round_place=1)
lifetime_dict = create_steps(arg_arr=np.arange(start=10, stop=41, step=5), str_round_place='{:.0f}', val_round_place=0)
kwt_dict = create_steps(arg_arr=np.arange(start=0, stop=1010, step=100), str_round_place='{:.0f}', val_round_place=0)
kwe_dict = create_steps(arg_arr=np.arange(start=0, stop=10100, step=1000), str_round_place='{:.0f}', val_round_place=0)
precool_dict = create_steps(arg_arr=np.arange(start=0, stop=40.1, step=5), str_round_place='{:.0f}', val_round_place=0)
turb_pout_dict = create_steps(arg_arr=np.arange(start=75, stop=201, step=10), str_round_place='{:.0f}', val_round_place=0)

discount_dict = {0: '0', 1: '', 2: '', 3: '', 4: '', 5: '', 6: '', 7: '', 8: '', 9: '', 10: '', 
                    11: '', 12: '', 13: '', 14: '', 15: '', 16: '', 17: '', 18: '', 19: '', 20: '20'}
precool_dict = {0: '0', 5: '', 10: '', 15: '', 20: '', 25: '', 30: '', 35: '', 40: '40'}

fineness_dict = {0: '0 (coarse)', 1: '', 2: '2 (fine)'} #, 3: '', 4: '', 5: '5 (fine)'}
accuracy_dict = {1: '1 (coarse)', 2: '', 3: '', 4: '', 5: '5 (fine)'}

# new slider labels
# grad_dict = {0.01: '0.01', 0.1: '0.1'}
# k_dict = {0.1: '0.1', 7.0: '7.0'} # 1.5 4.5
Tsurf_dict = {0: '0', 40: '40'}
c_dict = {500: '500', 2000: '2000'}
rho_dict = {1000: '1000', 3500: '3500'}

# Tinj_dict =  {20: '20', 200: '200'}
# mdot_dict = {1: '1', 2000: '2000'}
# L2_dict = {1000: '1k', 50000: '50k'}
# L1_dict = {1000: '1k', 10000: '10k'}

radius_vertical_dict = {0.10795: '0.10795', 0.22225: '0.22225'}
radius_lateral_dict = {0.10795: '0.10795', 0.22225: '0.22225'}
radius_centerpipe_dict = {0.0635: '0.0635', 0.174: '0.174'}
thickness_centerpipe_dict = {0.005: '0.005', 0.025: '0.025'}

inlet_pressure_dict = {5: '5', 20: '20'}
pipe_roughness_dict =  {0.000001: '0.000001', 0.000003: '0.000003'}

# TODO: need to make it general across parameters 
start_vals_hdf5 = {"Tsurf": 25, "c": 790.0, "rho": 2750, "n-laterals": 1, "lateral-flow": 1, "lateral-multiplier": 1}
start_vals_d = {"mdot": 24.0, "L2": 10000, "L1": 3500 , "Tinj": 30.0, "grad": 0.05, "D": 0.3500, "k": 3.0}
# start_vals_d = {"mdot": 20.0, "L2": 1000, "L1": 2000 , "grad": 0.05, "D": 0.2280, "Tinj": 20.0, "k": 2.83} #
start_vals_sbt = {"mesh": 0, "accuracy": 1, "mass-mode": 0, "temp-mode": 0,
                    "radius-vertical": 0.3500, "radius-lateral": 0.3500,
                    "radius": 0.2286,
                    "radiuscenterpipe": 0.127, 
                    "thicknesscenterpipe": 0.0127,
                    "k_center_pipe": 0.006,
                    "coaxialflowtype": 1,
                    "inletpressure": 10,
                    "piperoughness": 0.000001 # 1e-6
                    } 
start_vals_econ = {"drillcost": 1000, "discount-rate": 7.0, "lifetime": 40, "kwt": 100,
                    "kwe": 3000, "precool": 13, "turb-pout": 80
                    }

# ---------------------------
# Global style element(s).
# ---------------------------
div_block_style = {'display': 'block'}
div_none_style = {'display': 'none'}
p_bold_style = {"fontWeight": "bold"}


def slider1(DivID, ID, ptitle, min_v, max_v, mark_dict, step_i, start_v, div_style):

    # ---------------------------------------------------------------------------
    # Create a Div with the name and the slider stacked **with** the option to 
    # define steps.
    # ---------------------------------------------------------------------------

    return html.Div(id=DivID,
                    style=div_style,
                    children=[
                       html.P(ptitle, style=p_bold_style),
                       dcc.Slider(id=ID,
                       min=min_v, max=max_v,
                       marks=mark_dict, 
                       step=step_i,
                       value=start_v,
                       tooltip={"placement": "bottom", "always_visible": True}),
                       ]
                    )


def slider2(DivID, ID, ptitle, min_v, max_v, mark_dict, start_v, div_style):

    # ---------------------------------------------------------------------------
    # Create a Div with the name and the slider stacked **without** the option to 
    # define steps.
    # ---------------------------------------------------------------------------

    return html.Div(id=DivID,
                    style=div_style,
                    children=[
                       html.P(ptitle, style=p_bold_style),
                       dcc.Slider(id=ID,
                       min=min_v, max=max_v, #step=None,
                       marks=mark_dict,
                       value=start_v,
                       tooltip={"placement": "bottom", "always_visible": True}
                       ),
                       ]
                    )


def input_box(DivID, ID, ptitle, min_v, max_v, start_v, step_i, div_style):

    return html.Div(
            id=DivID,
            style=div_style,
            className="name-input-container",
            children=[
                html.P(ptitle, className="input-title"),
                dcc.Input(id=ID, disabled=True,
                            value=start_v, type='number', min=min_v, max=max_v, step=step_i, className="input-box"),
        ])


def dropdown_box(DivID, ID, ptitle, options, disabled, div_style):

        return html.Div(
                id=DivID,
                className="name-input-container-dd",
                style=div_style,
                children=[
                        html.P(ptitle, className="input-title"),
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

def slider_card():

    # -----------------------------------------------------------------------
    # A Div containing controls for graphs. Controls include the following:
    #
    #   Sliders for thermal parameters:
    #           mdot, L2, L1, grad, D, Tinj, and k parameters.
    #
    #   Sliders for economic parameters;
    #       drillcost, discount-rate, lifetime, kwt, kwe, precool, turb-pout
    #
    # -----------------------------------------------------------------------
    
    # Get current unit preferences dynamically
    from unit_conversions import unit_converter

    return html.Div(
        children=[
            dbc.Button(
                    "ⓘ Finetune System Values",
                    id="collapse-button3",
                    className="mb-33",
                    color="primary",
                    n_clicks=0,
                ),
            dbc.Collapse(
                dbc.Card(
                    dbc.CardBody(
                        children=[

                            ##
                            html.Div(id="slider-container",
                                    children=[
                                    html.Div(id="therm-sliders",
                                        children=[

                                            html.Div(
                                                id="geologic-params-div",
                                                className="params-div",
                                                children=[
                                                    html.P("GEOLOGIC PROPERTIES", className="param-class-name"),
                                                    create_enhanced_slider(DivID="Tsurf-select-div", ID="Tsurf-select", ptitle=f"Surface Temperature ({get_unit_symbol(unit_converter.user_preferences.get('temperature', 'C'))})", min_v=0, max_v=40.0, 
                                                            mark_dict=create_imperial_marks(0, 40.0, unit_converter.user_preferences.get('temperature', 'C')) if unit_converter.user_preferences.get('temperature', 'C') == 'F' else Tsurf_dict, start_v=start_vals_hdf5["Tsurf"], div_style=div_none_style, parameter_name="Surface Temperature (˚C)"),
                                                   
                                                    html.Div(
                                                            id="grad-container",
                                                            children=[
                                                                    create_enhanced_slider(DivID="grad-select-div", ID="grad-select", ptitle=f"Geothermal Gradient ({get_unit_symbol(unit_converter.user_preferences.get('geothermal_gradient', 'K/m'))})", 
                                                                            min_v=get_geothermal_gradient_converted_values(0.03, unit_converter.user_preferences.get('geothermal_gradient', 'K/m')), 
                                                                            max_v=get_geothermal_gradient_converted_values(0.07, unit_converter.user_preferences.get('geothermal_gradient', 'K/m')), 
                                                                            mark_dict=create_imperial_marks(
                                                                                get_geothermal_gradient_converted_values(0.03, unit_converter.user_preferences.get('geothermal_gradient', 'K/m')),
                                                                                get_geothermal_gradient_converted_values(0.07, unit_converter.user_preferences.get('geothermal_gradient', 'K/m')),
                                                                                unit_converter.user_preferences.get('geothermal_gradient', 'K/m')
                                                                            ) if unit_converter.user_preferences.get('geothermal_gradient', 'K/m') in ['F/ft', 'F/yd'] else grad_dict, 
                                                                            start_v=get_geothermal_gradient_converted_values(start_vals_d["grad"], unit_converter.user_preferences.get('geothermal_gradient', 'K/m')), 
                                                                            div_style=div_block_style, parameter_name="Geothermal Gradient (K/m)")
                                                            ]),

                                                    html.Div(
                                                            id="k-container",
                                                            children=[
                                                                    create_enhanced_slider(DivID="k-select-div", ID="k-select", ptitle=f"Rock Thermal Conductivity ({get_unit_symbol(unit_converter.user_preferences.get('thermal_conductivity', 'W/m-K'))})", 
                                                                            min_v=get_thermal_conductivity_converted_values(1.5, unit_converter.user_preferences.get('thermal_conductivity', 'W/m-K')), 
                                                                            max_v=get_thermal_conductivity_converted_values(4.5, unit_converter.user_preferences.get('thermal_conductivity', 'W/m-K')), 
                                                                            mark_dict=create_imperial_marks(
                                                                                get_thermal_conductivity_converted_values(1.5, unit_converter.user_preferences.get('thermal_conductivity', 'W/m-K')),
                                                                                get_thermal_conductivity_converted_values(4.5, unit_converter.user_preferences.get('thermal_conductivity', 'W/m-K')),
                                                                                unit_converter.user_preferences.get('thermal_conductivity', 'W/m-K')
                                                                            ) if unit_converter.user_preferences.get('thermal_conductivity', 'W/m-K') in ['Btu/ft-h-F', 'Btu/yd-h-F'] else k_dict, 
                                                                            start_v=get_thermal_conductivity_converted_values(start_vals_d["k"], unit_converter.user_preferences.get('thermal_conductivity', 'W/m-K')), 
                                                                            div_style=div_block_style, parameter_name="Rock Thermal Conductivity (W/m-K)")
                                                                    
                                                            ]),
                                                    create_enhanced_slider(DivID="c-select-div", ID="c-select", ptitle=f"Rock Specific Heat Capacity ({get_unit_symbol(unit_converter.user_preferences.get('heat_capacity', 'J/kg-K'))})", 
                                                            min_v=get_heat_capacity_converted_values(500, unit_converter.user_preferences.get('heat_capacity', 'J/kg-K')), 
                                                            max_v=get_heat_capacity_converted_values(2000, unit_converter.user_preferences.get('heat_capacity', 'J/kg-K')), 
                                                            mark_dict=create_imperial_marks(
                                                                get_heat_capacity_converted_values(500, unit_converter.user_preferences.get('heat_capacity', 'J/kg-K')),
                                                                get_heat_capacity_converted_values(2000, unit_converter.user_preferences.get('heat_capacity', 'J/kg-K')),
                                                                unit_converter.user_preferences.get('heat_capacity', 'J/kg-K')
                                                            ) if unit_converter.user_preferences.get('heat_capacity', 'J/kg-K') == 'Btu/lb-F' else c_dict, step_i=1, 
                                                            start_v=get_heat_capacity_converted_values(start_vals_hdf5["c"], unit_converter.user_preferences.get('heat_capacity', 'J/kg-K')), 
                                                            div_style=div_none_style, parameter_name="Rock Specific Heat Capacity (J/kg-K)"),
                                                    create_enhanced_slider(DivID="rho-select-div", ID="rho-select", ptitle=f"Rock Density ({get_unit_symbol(unit_converter.user_preferences.get('density', 'kg/m3'))})", 
                                                            min_v=get_density_converted_values(1000, unit_converter.user_preferences.get('density', 'kg/m3')), 
                                                            max_v=get_density_converted_values(3500, unit_converter.user_preferences.get('density', 'kg/m3')), 
                                                            mark_dict=create_imperial_marks(
                                                                get_density_converted_values(1000, unit_converter.user_preferences.get('density', 'kg/m3')),
                                                                get_density_converted_values(3500, unit_converter.user_preferences.get('density', 'kg/m3')),
                                                                unit_converter.user_preferences.get('density', 'kg/m3')
                                                            ) if unit_converter.user_preferences.get('density', 'kg/m3') in ['lb/ft3', 'lb/yd3'] else rho_dict, step_i=1, 
                                                            start_v=get_density_converted_values(start_vals_hdf5["rho"], unit_converter.user_preferences.get('density', 'kg/m3')), 
                                                            div_style=div_none_style, parameter_name="Rock Density (kg/m3)"),
                                                ]
                                            ),
                                        
                                            html.Div(
                                                id="wellbore-params-div", # wellbore operations
                                                className="params-div",
                                                children=[
                                                    html.P("WELLBORE OPERATIONS", className="param-class-name"),

                                                    html.Div(
                                                            id="Tinj-container",
                                                            children=[
                                                                create_enhanced_slider(DivID="Tinj-select-div", ID="Tinj-select", ptitle=f"Injection Temperature ({get_unit_symbol(unit_converter.user_preferences.get('temperature', 'C'))})", 
                                                                        min_v=get_temperature_converted_values(30.0, unit_converter.user_preferences.get('temperature', 'C')), 
                                                                        max_v=get_temperature_converted_values(60.0, unit_converter.user_preferences.get('temperature', 'C')), 
                                                                        mark_dict=create_imperial_marks(
                                                                            get_temperature_converted_values(30.0, unit_converter.user_preferences.get('temperature', 'C')),
                                                                            get_temperature_converted_values(60.0, unit_converter.user_preferences.get('temperature', 'C')),
                                                                            unit_converter.user_preferences.get('temperature', 'C')
                                                                        ) if unit_converter.user_preferences.get('temperature', 'C') == 'F' else Tinj_dict, 
                                                                        start_v=get_temperature_converted_values(30.0, unit_converter.user_preferences.get('temperature', 'C')), 
                                                                        div_style=div_block_style, parameter_name="Injection Temperature (˚C)")
                                                            ]),
                                                    html.Div(
                                                            id="mdot-container",
                                                            children=[        
                                                                create_enhanced_slider(DivID="mdot-select-div", ID="mdot-select", ptitle=f"Mass Flow Rate ({get_unit_symbol(unit_converter.user_preferences.get('mass_flow', 'kg/s'))})", 
                                                                        min_v=get_mass_flow_converted_values(u_sCO2.mdot[0], unit_converter.user_preferences.get('mass_flow', 'kg/s')), 
                                                                        max_v=get_mass_flow_converted_values(u_sCO2.mdot[-1], unit_converter.user_preferences.get('mass_flow', 'kg/s')), 
                                                                        mark_dict=create_imperial_marks(
                                                                            get_mass_flow_converted_values(u_sCO2.mdot[0], unit_converter.user_preferences.get('mass_flow', 'kg/s')),
                                                                            get_mass_flow_converted_values(u_sCO2.mdot[-1], unit_converter.user_preferences.get('mass_flow', 'kg/s')),
                                                                            unit_converter.user_preferences.get('mass_flow', 'kg/s')
                                                                        ) if unit_converter.user_preferences.get('mass_flow', 'kg/s') == 'lb/s' else mdot_dict, 
                                                                        start_v=get_mass_flow_converted_values(start_vals_d["mdot"], unit_converter.user_preferences.get('mass_flow', 'kg/s')), 
                                                                        div_style=div_block_style, parameter_name="Mass Flow Rate (kg/s)")
                                                            ]),
                                                ]

                                            ),

                                            html.Div(
                                                id="tube-params-div", # tube geometry 
                                                className="params-div",
                                                children=[
                                                    html.P("TUBE GEOMETRY", className="param-class-name"),
                                                    html.Div(
                                                            id="diameter-container",
                                                            children=[ 
                                                                create_enhanced_slider(DivID="diameter-select-div", ID="diameter-select", ptitle=f"Borehole Diameter ({get_unit_symbol(unit_converter.user_preferences.get('length', 'm'))})", 
                                                                        min_v=get_length_converted_values(0.2159, unit_converter.user_preferences.get('length', 'm')), 
                                                                        max_v=get_length_converted_values(0.4445, unit_converter.user_preferences.get('length', 'm')), 
                                                                        mark_dict=create_imperial_marks(
                                                                            get_length_converted_values(0.2159, unit_converter.user_preferences.get('length', 'm')),
                                                                            get_length_converted_values(0.4445, unit_converter.user_preferences.get('length', 'm')),
                                                                            unit_converter.user_preferences.get('length', 'm')
                                                                        ) if unit_converter.user_preferences.get('length', 'm') in ['ft', 'yd'] else D_dict, 
                                                                        start_v=get_length_converted_values(start_vals_d["D"], unit_converter.user_preferences.get('length', 'm')), 
                                                                        div_style=div_block_style, parameter_name="Borehole Diameter (m)")
                                                            ]),
                                                    html.Div(
                                                            id="Diameter1-container",
                                                            children=[
                                                                create_enhanced_slider(DivID="radius-vertical-select-div", ID="radius-vertical-select", ptitle=f"Wellbore Radius Vertical ({get_unit_symbol(unit_converter.user_preferences.get('length', 'm'))})", min_v=0.10795, max_v=0.22225,
                                                                mark_dict=create_imperial_marks(
                                                                    get_length_converted_values(0.10795, unit_converter.user_preferences.get('length', 'm')),
                                                                    get_length_converted_values(0.22225, unit_converter.user_preferences.get('length', 'm')),
                                                                    unit_converter.user_preferences.get('length', 'm')
                                                                ) if unit_converter.user_preferences.get('length', 'm') in ['ft', 'yd'] else radius_vertical_dict, step_i=0.001, start_v=start_vals_sbt["radius-vertical"], div_style=div_none_style, parameter_name="Wellbore Radius Vertical (m)")
                                                            ]),
                                                    html.Div(
                                                            id="Diameter2-container",
                                                            children=[
                                                                create_enhanced_slider(DivID="radius-lateral-select-div", ID="radius-lateral-select", ptitle=f"Wellbore Radius Lateral ({get_unit_symbol(unit_converter.user_preferences.get('length', 'm'))})", min_v=0.10795, max_v=0.22225,
                                                                        mark_dict=create_imperial_marks(
                                                                            get_length_converted_values(0.10795, unit_converter.user_preferences.get('length', 'm')),
                                                                            get_length_converted_values(0.22225, unit_converter.user_preferences.get('length', 'm')),
                                                                            unit_converter.user_preferences.get('length', 'm')
                                                                        ) if unit_converter.user_preferences.get('length', 'm') in ['ft', 'yd'] else radius_lateral_dict, step_i=0.001, start_v=start_vals_sbt["radius-lateral"], div_style=div_none_style, parameter_name="Wellbore Radius Lateral (m)")
                                                            ]),
                                                    html.Div(
                                                            id="L2-container",
                                                            children=[ 
                                                                create_enhanced_slider(DivID="L2-select-div", ID="L2-select", ptitle=f"Horizontal Extent ({get_unit_symbol(unit_converter.user_preferences.get('length', 'm'))})", 
                                                                        min_v=get_length_converted_values(u_sCO2.L2[0], unit_converter.user_preferences.get('length', 'm')), 
                                                                        max_v=get_length_converted_values(u_sCO2.L2[-1], unit_converter.user_preferences.get('length', 'm')), 
                                                                        mark_dict=create_imperial_marks(
                                                                            get_length_converted_values(u_sCO2.L2[0], unit_converter.user_preferences.get('length', 'm')),
                                                                            get_length_converted_values(u_sCO2.L2[-1], unit_converter.user_preferences.get('length', 'm')),
                                                                            unit_converter.user_preferences.get('length', 'm')
                                                                        ) if unit_converter.user_preferences.get('length', 'm') in ['ft', 'yd'] else L2_dict, 
                                                                        start_v=get_length_converted_values(start_vals_d["L2"], unit_converter.user_preferences.get('length', 'm')), 
                                                                        div_style=div_block_style, parameter_name="Horizontal Extent (m)")
                                                            ]),
                                                    html.Div(
                                                            id="L1-container",
                                                            children=[ 
                                                                create_enhanced_slider(DivID="L1-select-div", ID="L1-select", ptitle=f"Drilling Depth ({get_unit_symbol(unit_converter.user_preferences.get('length', 'm'))})", 
                                                                        min_v=get_length_converted_values(u_sCO2.L1[0], unit_converter.user_preferences.get('length', 'm')), 
                                                                        max_v=get_length_converted_values(u_sCO2.L1[-1], unit_converter.user_preferences.get('length', 'm')), 
                                                                        mark_dict=create_imperial_marks(
                                                                            get_length_converted_values(u_sCO2.L1[0], unit_converter.user_preferences.get('length', 'm')),
                                                                            get_length_converted_values(u_sCO2.L1[-1], unit_converter.user_preferences.get('length', 'm')),
                                                                            unit_converter.user_preferences.get('length', 'm')
                                                                        ) if unit_converter.user_preferences.get('length', 'm') in ['ft', 'yd'] else L1_dict, 
                                                                        start_v=get_length_converted_values(start_vals_d["L1"], unit_converter.user_preferences.get('length', 'm')), 
                                                                        div_style=div_block_style, parameter_name="Drilling Depth (m)")
                                                            ]),
                                                    html.Div(
                                                            id="num-lat-container",
                                                            children=[ 
                                                                create_enhanced_input_box(DivID="num-lat-div", ID="n-laterals-select", ptitle="Number of Laterals", 
                                                                            min_v=0, max_v=20, start_v=start_vals_hdf5["n-laterals"], step_i=1, div_style=div_none_style, parameter_name="Number of Laterals")
                                                            ]),
                                                    html.Div(
                                                            id="lat-allo-container",
                                                            children=[ 
                                                                create_enhanced_input_box(DivID="lat-allocation-div", ID="lateral-flow-select", ptitle="Lateral Flow Allocation", 
                                                                            min_v=0, max_v=1, start_v=start_vals_hdf5["lateral-flow"], step_i=0.01, div_style=div_none_style, parameter_name="Lateral Flow Allocation")
                                                            ]),
                                                    html.Div(
                                                            id="lat-flow-container",
                                                            children=[
                                                                create_enhanced_input_box(DivID="lat-flow-mul-div", ID="lateral-multiplier-select", ptitle="Lateral Flow Multiplier", 
                                                                                        min_v=0, max_v=1, start_v=start_vals_hdf5["lateral-multiplier"], step_i=0.05, div_style=div_none_style, parameter_name="Lateral Flow Multiplier")
                                                            ]),

                                                ]
                                            ),

                                    ]),

                                    html.Div(id="econ-sliders",
                                        children=[
                                        
                                            html.Div(
                                                    id="economics-params-div",
                                                    className="params-div",
                                                    children=[
                                                        html.P("ECONOMIC PARAMETERS", className="param-class-name"),
                                                                                                                 create_enhanced_slider(DivID="drillcost-div", ID="drillcost-select", ptitle=f"Drilling Cost ($/{get_unit_symbol(unit_converter.user_preferences.get('length', 'm'))})", 
                                                                 min_v=get_drilling_cost_converted_values(0, unit_converter.user_preferences.get('length', 'm')), 
                                                                 max_v=get_drilling_cost_converted_values(4000, unit_converter.user_preferences.get('length', 'm')), 
                                                                 mark_dict=create_imperial_marks(
                                                                     get_drilling_cost_converted_values(0, unit_converter.user_preferences.get('length', 'm')),
                                                                     get_drilling_cost_converted_values(4000, unit_converter.user_preferences.get('length', 'm')),
                                                                     f"drillcost_{unit_converter.user_preferences.get('length', 'm')}"
                                                                 ) if unit_converter.user_preferences.get('length', 'm') in ['ft', 'yd'] else drillcost_dict, 
                                                                 start_v=get_drilling_cost_converted_values(start_vals_econ["drillcost"], unit_converter.user_preferences.get('length', 'm')), 
                                                                 div_style=div_block_style, parameter_name="Drilling Cost ($/m)"),
                                                        create_enhanced_slider(DivID="discount-rate-div", ID="discount-rate-select", ptitle="Discount Rate (%)", min_v=0, max_v=20, 
                                                                mark_dict=discount_dict, start_v=start_vals_econ["discount-rate"], div_style=div_block_style, parameter_name="Discount Rate (%)"),
                                                        # slider2(DivID="lifetime-div", ID="lifetime-select", ptitle="Lifetime (years)", min_v=10, max_v=40, 
                                                        #         mark_dict=lifetime_dict, start_v=30),
                                                        create_enhanced_slider(DivID="lifetime-div", ID="lifetime-select", ptitle="Lifetime (years)", min_v=10, max_v=40, 
                                                                mark_dict=lifetime_dict, step_i=1, start_v=start_vals_econ["lifetime"], div_style=div_block_style, parameter_name="Lifetime (years)"),
                                                        create_enhanced_slider(DivID="kwt-div", ID="kwt-select", ptitle="Plant CAPEX ($/kWt)", min_v=0, max_v=1000, 
                                                                mark_dict=kwt_dict, start_v=start_vals_econ["kwt"], div_style=div_block_style, parameter_name="Plant CAPEX ($/kWt)"),
                                                        create_enhanced_slider(DivID="kwe-div", ID="kwe-select", ptitle="Plant CAPEX ($/kWe)", min_v=0, max_v=10000, 
                                                                mark_dict=kwe_dict, start_v=start_vals_econ["kwe"], div_style=div_block_style, parameter_name="Plant CAPEX ($/kWe)"),
                                                    ]
                                                    ),
                                            html.Div(id="sCO2-card",
                                                    className="params-div",
                                                    children=[
                                                        html.P("ⓘ Multiple LCOE minima exist. Dial here to explore:", id="sCO2-text"),  # Run the optimizer 
                                                        create_enhanced_slider(DivID="precool-div", ID="precool-select", ptitle=f"Pre-cooling ({get_unit_symbol(unit_converter.user_preferences.get('temperature', 'C'))})", min_v=0, max_v=40, 
                                                                mark_dict=precool_dict, start_v=start_vals_econ["precool"], div_style=div_block_style, parameter_name="Pre-cooling (˚C)"),
                                                        create_enhanced_slider(DivID="turb-pout-div", ID="turb-pout-select", ptitle=f"Turbine Outlet Pressure ({get_unit_symbol(unit_converter.user_preferences.get('pressure', 'bar'))})", 
                                                                min_v=get_pressure_converted_values(75, unit_converter.user_preferences.get('pressure', 'bar')), 
                                                                max_v=get_pressure_converted_values(200, unit_converter.user_preferences.get('pressure', 'bar')), 
                                                                mark_dict=turb_pout_dict, 
                                                                start_v=get_pressure_converted_values(start_vals_econ["turb-pout"], unit_converter.user_preferences.get('pressure', 'bar')), 
                                                                div_style=div_block_style, parameter_name="Turbine Outlet Pressure (bar)"),
                                                        html.Div(id="check-visual-card",
                                                                children=[
                                                                        dcc.Checklist(id="checklist",
                                                                            options=[' '],
                                                                            value=[' '],
                                                                            inline=True
                                                                        ),
                                                                        # html.P("Display heat added or removed from the CLG system to tune for LCOE minima.", id="TS-diagram-text"),
                                                                        html.P("In CLGS, working fluids drive turbines to generate electricity. \
                                                                            Display a cycle for sCO2 that depends on turbine cooling and turbine outlet pressure. \
                                                                            Use the cycle to tune for LCOE minima and maximize net electric power by minimizing excess heating.", id="TS-diagram-text"),
                                                            ])
                                                    ]),

                                            html.Div(
                                                id="model-params-div",
                                                className="params-div",
                                                style=div_none_style,
                                                children=[
                                                    html.P("MODEL FINE-TUNING", className="param-class-name"),
                                                    create_enhanced_slider(DivID="mesh-div", ID="mesh-select", ptitle="Mesh Fineness", min_v=0, max_v=2, 
                                                                mark_dict=fineness_dict, step_i=1, start_v=start_vals_sbt["mesh"], div_style=div_block_style, parameter_name="Mesh Fineness"),
                                                    create_enhanced_slider(DivID="accuracy-div", ID="accuracy-select", ptitle="Accuracy", min_v=1, max_v=5, 
                                                                mark_dict=accuracy_dict, step_i=1,start_v=start_vals_sbt["accuracy"], div_style=div_block_style, parameter_name="Accuracy"),
                                                
                                                    html.Div(
                                                            id="hyperparam1-container",
                                                            children=[
                                                                create_enhanced_dropdown(DivID="mass-flow-mode-div", ID="mass-mode-select", ptitle="Mass Flow Rate Mode", 
                                                                                                options=["Constant", "Variable"], disabled=True, div_style=div_block_style, parameter_name="Mass Flow Rate Mode")
                                                        ]),
                                                    html.Div(
                                                        id="hyperparam3-container",
                                                        children=[
                                                                create_enhanced_dropdown(DivID="temp-flow-mode-div", ID="temp-mode-select", ptitle="Injection Temperature Mode", 
                                                                                        options=["Constant", "Variable"], disabled=True, div_style=div_block_style, parameter_name="Injection Temperature Mode")
                                                        ]),
                                                    html.Div(
                                                        id="hyperparam5-container",
                                                        children=[
                                                                create_enhanced_dropdown(DivID="fluid-mode-div", ID="fluid-mode-select", ptitle="Fluid Properties Mode", 
                                                                                                options=["Variable", "Constant"], disabled=True, div_style=div_none_style, parameter_name="Fluid Properties Mode")
                                                    ]),
                                                ]

                                            ),
                                    ]),
                                ]),



                        #     html.Div(id="slider-container",
                        #             children=[
                        #             html.Div(id="therm-sliders",
                        #                 children=[

                        #                 slider2(DivID="mdot-select-div", ID="mdot-select", ptitle="Mass Flow Rate (kg/s)", min_v=u_sCO2.mdot[0], max_v=u_sCO2.mdot[-1], 
                        #                         mark_dict=mdot_dict, start_v=start_vals_d["mdot"]),
                        #                 slider2(DivID="L2-select-div", ID="L2-select", ptitle="Horizontal Extent (m)", min_v=u_sCO2.L2[0], max_v=u_sCO2.L2[-1], 
                        #                         mark_dict=L2_dict, start_v=start_vals_d["L2"]),
                        #                 slider2(DivID="L1-select-div", ID="L1-select", ptitle="Drilling Depth (m)", min_v=u_sCO2.L1[0], max_v=u_sCO2.L1[-1], 
                        #                         mark_dict=L1_dict, start_v=start_vals_d["L1"]),
                        #                 slider2(DivID="grad-select-div", ID="grad-select", ptitle="Geothermal Gradient (K/m)", min_v=u_sCO2.grad[0], max_v=u_sCO2.grad[-1], 
                        #                         mark_dict=grad_dict, start_v=start_vals_d["grad"]),
                        #                 # slider2(DivID="diameter-select-div", ID="diameter-select", ptitle="Borehole Diameter (m)", min_v=u_sCO2.D[0], max_v=u_sCO2.D[-1], 
                        #                 #         mark_dict=D_dict, start_v=0.3500),
                        #                 slider1(DivID="diameter-select-div", ID="diameter-select", ptitle="Borehole Diameter (m)", min_v=0.2159, max_v=0.4445, 
                        #                         mark_dict=D_dict, step_i=0.002, start_v=start_vals_d["D"]),
                        #                 # slider2(DivID="Tinj-select-div", ID="Tinj-select", ptitle="Injection Temperature (˚C)", min_v=u_sCO2.Tinj[0] - 273.15, max_v=u_sCO2.Tinj[-1] - 273.15, 
                        #                 #         mark_dict=Tinj_dict, start_v=303.15-273.15),
                        #                 slider2(DivID="Tinj-select-div", ID="Tinj-select", ptitle="Injection Temperature (˚C)", min_v=30.0, max_v=60.0, 
                        #                         mark_dict=Tinj_dict, start_v=start_vals_d["Tinj"]),
                        #                 # slider1(DivID="Tinj-select-div", ID="Tinj-select", ptitle="Injection Temperature (˚C)", min_v=30, max_v=60, 
                        #                 #     mark_dict=Tinj_dict, step_i=0.2,start_v=30),
                        #                 slider2(DivID="k-select-div", ID="k-select", ptitle="Rock Thermal Conductivity (W/m-K)", min_v=u_sCO2.k[0], max_v=u_sCO2.k[-1], 
                        #                         mark_dict=k_dict, start_v=start_vals_d["k"])
                        #             ]),

                        #             html.Div(id="econ-sliders",
                        #             children=[

                        #             slider2(DivID="drillcost-div", ID="drillcost-select", ptitle="Drilling Cost ($/m)", min_v=0, max_v=4000, 
                        #                     mark_dict=drillcost_dict, start_v=1000),
                        #             slider2(DivID="discount-rate-div", ID="discount-rate-select", ptitle="Discount Rate (%)", min_v=0, max_v=20, 
                        #                     mark_dict=discount_dict, start_v=7.0),
                        #             # slider2(DivID="lifetime-div", ID="lifetime-select", ptitle="Lifetime (years)", min_v=10, max_v=40, 
                        #             #         mark_dict=lifetime_dict, start_v=30),
                        #             slider1(DivID="lifetime-div", ID="lifetime-select", ptitle="Lifetime (years)", min_v=10, max_v=40, 
                        #                     mark_dict=lifetime_dict, step_i=1, start_v=40),
                        #             slider2(DivID="kwt-div", ID="kwt-select", ptitle="Plant CAPEX ($/kWt)", min_v=0, max_v=1000, 
                        #                     mark_dict=kwt_dict, start_v=100),
                        #             slider2(DivID="kwe-div", ID="kwe-select", ptitle="Plant CAPEX ($/kWe)", min_v=0, max_v=10000, 
                        #                     mark_dict=kwe_dict, start_v=3000),
                                    
                        #             html.Div(id="sCO2-card",
                        #                      children=[
                        #                         html.P("ⓘ Multiple LCOE minima exist. Dial here to explore:", id="sCO2-text"),  # Run the optimizer 
                        #                         slider2(DivID="precool-div", ID="precool-select", ptitle="Pre-cooling (˚C)", min_v=0, max_v=40, 
                        #                                 mark_dict=precool_dict, start_v=13),
                        #                         slider2(DivID="turb-pout-div", ID="turb-pout-select", ptitle="Turbine Outlet Pressure (bar)", min_v=75, max_v=200, 
                        #                                 mark_dict=turb_pout_dict, start_v=80),
                        #                         html.Div(id="check-visual-card",
                        #                                  children=[
                        #                                         dcc.Checklist(id="checklist",
                        #                                             options=[' '],
                        #                                             value=[' '],
                        #                                             inline=True
                        #                                         ),
                        #                                         # html.P("Display heat added or removed from the CLG system to tune for LCOE minima.", id="TS-diagram-text"),
                        #                                         html.P("In CLGS, working fluids drive turbines to generate electricity. \
                        #                                             Display a cycle for sCO2 that depends on turbine cooling and turbine outlet pressure. \
                        #                                             Use the cycle to tune for LCOE minima and maximize net electric power by minimizing excess heating.", id="TS-diagram-text"),
                        #                             ])
                        #                     ])

                        #             ]),
                        #         ])



                            ]
                    )),
                id="collapse3",
                is_open=True,
            ),
        ]
    )
