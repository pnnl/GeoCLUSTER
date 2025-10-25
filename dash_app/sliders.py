#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# web app and interactive graphics libraries 
from dash import dcc, html
import dash_bootstrap_components as dbc

# sourced scripts
from plots import * # u_sCO2, u_H2O, c_sCO2, c_H2O

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

wellbore_operations_l = ["Injection Temperature (˚C)", "Mass Flow Rate (kg/s)", "Horizontal Extent (m)", "Drilling Depth (m)"]

tube_geometry_l = ["Wellbore Radius Vertical (m)", "Wellbore Radius Lateral (m)", # "Borehole Diameter (m)",
                                "Number of Laterals", "Lateral Flow Allocation", "Lateral Flow Multiplier"]

economic_params_l = ["Drilling Cost ($/m)", "Discount Rate (%)", "Lifetime (years)", "Plant CAPEX ($/kWt)", 
                    "Plant CAPEX ($/kWe)", "Pre-cooling (˚C)", "Turbine Outlet Pressure (bar)"]

geologic_properties_l = ["Surface Temperature (˚C)", "Geothermal Gradient (°C/m)", "Rock Thermal Conductivity (W/m-K)", 
                            "Rock Specific Heat Capacity (J/kg-K)", "Rock Density (kg/m3)"]
                            
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
start_vals_hdf5 = {"Tsurf": 25, "c": 790.0, "rho": 2800, "n-laterals": 1, "lateral-flow": 1, "lateral-multiplier": 1}
start_vals_d = {"mdot": 24.0, "L2": 10000, "L1": 3500 , "Tinj": 50.0, "grad": 0.03, "D": 0.3500, "k": 3.0}
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
p_bold_style = {"fontWeight": "bold", "textAlign": "center", "width": "100%"}


def slider1(DivID, ID, ptitle, min_v, max_v, mark_dict, step_i, start_v, div_style, parameter_name=None, custom_title=False):

    # ---------------------------------------------------------------------------
    # Create a Div with the name and the slider stacked **with** the option to 
    # define steps.
    # ---------------------------------------------------------------------------
    
    from info_popups import create_info_button
    
    info_button = create_info_button(parameter_name) if parameter_name else html.Div()

    return html.Div(id=DivID,
                    style=div_style,
                    children=[
                       html.Div(className="title-button-container", style={"display": "flex", "justifyContent": "center", "alignItems": "center"}, children=[
                           html.Div([
                               html.Div([
                                   html.Div("Rock Specific Heat Capacity", style={"textAlign": "center", "fontWeight": "bold"}),
                                   html.Div("(J/kg-K)", style={"textAlign": "center", "fontWeight": "bold"})
                               ]) if custom_title and "Rock Specific Heat Capacity" in ptitle else html.P(ptitle, style=p_bold_style),
                               info_button
                           ], style={"display": "flex", "alignItems": "center", "gap": "5px"})
                       ]),
                       dcc.Slider(id=ID,
                       min=min_v, max=max_v,
                       marks=mark_dict, 
                       step=step_i,
                       value=start_v,
                       tooltip={"placement": "bottom", "always_visible": True}),
                       ]
                    )


def slider2(DivID, ID, ptitle, min_v, max_v, mark_dict, start_v, div_style, parameter_name=None, custom_title=False):

    # ---------------------------------------------------------------------------
    # Create a Div with the name and the slider stacked **without** the option to 
    # define steps.
    # ---------------------------------------------------------------------------
    
    from info_popups import create_info_button
    
    info_button = create_info_button(parameter_name) if parameter_name else html.Div()

    return html.Div(id=DivID,
                    style=div_style,
                    children=[
                       html.Div(className="title-button-container", style={"display": "flex", "justifyContent": "center", "alignItems": "center"}, children=[
                           html.Div([
                               html.Div([
                                   html.Div("Rock Thermal Conductivity", style={"textAlign": "center", "fontWeight": "bold"}),
                                   html.Div("(W/m-K)", style={"textAlign": "center", "fontWeight": "bold"})
                               ]) if custom_title and "Rock Thermal Conductivity" in ptitle else html.P(ptitle, style=p_bold_style),
                               info_button
                           ], style={"display": "flex", "alignItems": "center", "gap": "5px"})
                       ]),
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
                                                    slider2(DivID="Tsurf-select-div", ID="Tsurf-select", ptitle="Surface Temperature (˚C)", min_v=0, max_v=40.0, 
                                                            mark_dict=Tsurf_dict, start_v=start_vals_hdf5["Tsurf"], div_style=div_none_style, parameter_name="Surface Temperature (˚C)"),
                                                   
                                                    html.Div(
                                                            id="grad-container",
                                                            children=[
                                                    slider2(DivID="grad-select-div", ID="grad-select", ptitle="Geothermal Gradient (°C/m)", min_v=u_sCO2.grad[0], max_v=u_sCO2.grad[-1], 
                                                            mark_dict=grad_dict, start_v=start_vals_d["grad"], div_style=div_block_style, parameter_name="Geothermal Gradient (°C/m)")
                                                            ]),

                                                    html.Div(
                                                            id="k-container",
                                                            children=[
                                                    slider2(DivID="k-select-div", ID="k-select", ptitle="Rock Thermal Conductivity (W/m-K)", min_v=u_sCO2.k[0], max_v=u_sCO2.k[-1], 
                                                            mark_dict=k_dict, start_v=start_vals_d["k"], div_style=div_block_style, parameter_name="Rock Thermal Conductivity (W/m-K)", custom_title=True)
                                                                    
                                                            ]),
                                                    slider1(DivID="c-select-div", ID="c-select", ptitle="Rock Specific Heat Capacity (J/kg-K)", min_v=500, max_v=2000, 
                                                            mark_dict=c_dict, step_i=1, start_v=start_vals_hdf5["c"], div_style=div_none_style, parameter_name="Rock Specific Heat Capacity (J/kg-K)", custom_title=True),
                                                    slider1(DivID="rho-select-div", ID="rho-select", ptitle="Rock Density (kg/m3)", min_v=1000, max_v=3500, 
                                                            mark_dict=rho_dict, step_i=1, start_v=start_vals_hdf5["rho"], div_style=div_none_style, parameter_name="Rock Density (kg/m3)"),
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
                                                                slider2(DivID="Tinj-select-div", ID="Tinj-select", ptitle="Injection Temperature (˚C)", min_v=30.0, max_v=60.0, 
                                                                        mark_dict=Tinj_dict, start_v=30.0, div_style=div_block_style, parameter_name="Injection Temperature (˚C)")
                                                            ]),
                                                    html.Div(
                                                            id="mdot-container",
                                                            children=[        
                                                                slider2(DivID="mdot-select-div", ID="mdot-select", ptitle="Mass Flow Rate (kg/s)", min_v=u_sCO2.mdot[0], max_v=u_sCO2.mdot[-1], 
                                                                        mark_dict=mdot_dict, start_v=start_vals_d["mdot"], div_style=div_block_style, parameter_name="Mass Flow Rate (kg/s)")
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
                                                                slider1(DivID="diameter-select-div", ID="diameter-select", ptitle="Borehole Diameter (m)", min_v=0.2159, max_v=0.4445, 
                                                                        mark_dict=D_dict, step_i=0.002, start_v=start_vals_d["D"], div_style=div_block_style, parameter_name="Borehole Diameter (m)")
                                                            ]),
                                                    html.Div(
                                                            id="Diameter1-container",
                                                            children=[
                                                                slider1(DivID="radius-vertical-select-div", ID="radius-vertical-select", ptitle="Wellbore Radius Vertical (m)", min_v=0.10795, max_v=0.22225,
                                                                mark_dict=radius_vertical_dict, step_i=0.001, start_v=start_vals_sbt["radius-vertical"], div_style=div_none_style, parameter_name="Wellbore Radius Vertical (m)")
                                                            ]),
                                                    html.Div(
                                                            id="Diameter2-container",
                                                            children=[
                                                                slider1(DivID="radius-lateral-select-div", ID="radius-lateral-select", ptitle="Wellbore Radius Lateral (m)", min_v=0.10795, max_v=0.22225,
                                                                        mark_dict=radius_lateral_dict, step_i=0.001, start_v=start_vals_sbt["radius-lateral"], div_style=div_none_style, parameter_name="Wellbore Radius Lateral (m)")
                                                            ]),
                                                    html.Div(
                                                            id="L2-container",
                                                            children=[ 
                                                                slider2(DivID="L2-select-div", ID="L2-select", ptitle="Horizontal Extent (m)", min_v=u_sCO2.L2[0], max_v=u_sCO2.L2[-1], 
                                                                        mark_dict=L2_dict, start_v=start_vals_d["L2"], div_style=div_block_style, parameter_name="Horizontal Extent (m)")
                                                            ]),
                                                    html.Div(
                                                            id="L1-container",
                                                            children=[ 
                                                                slider2(DivID="L1-select-div", ID="L1-select", ptitle="Drilling Depth (m)", min_v=u_sCO2.L1[0], max_v=u_sCO2.L1[-1], 
                                                                        mark_dict=L1_dict, start_v=start_vals_d["L1"], div_style=div_block_style, parameter_name="Drilling Depth (m)")
                                                            ]),
                                                    html.Div(
                                                            id="num-lat-container",
                                                            children=[ 
                                                                input_box(DivID="num-lat-div", ID="n-laterals-select", ptitle="Number of Laterals", 
                                                                            min_v=0, max_v=20, start_v=start_vals_hdf5["n-laterals"], step_i=1, div_style=div_none_style)
                                                            ]),
                                                    html.Div(
                                                            id="lat-allo-container",
                                                            children=[ 
                                                                input_box(DivID="lat-allocation-div", ID="lateral-flow-select", ptitle="Lateral Flow Allocation", 
                                                                            min_v=0, max_v=1, start_v=start_vals_hdf5["lateral-flow"], step_i=0.01, div_style=div_none_style)
                                                            ]),
                                                    html.Div(
                                                            id="lat-flow-container",
                                                            children=[
                                                                input_box(DivID="lat-flow-mul-div", ID="lateral-multiplier-select", ptitle="Lateral Flow Multiplier", 
                                                                                        min_v=0, max_v=1, start_v=start_vals_hdf5["lateral-multiplier"], step_i=0.05, div_style=div_none_style)
                                                            ]),
                                                    # html.Div(
                                                    #     id="num-lat-div",
                                                    #     style=div_none_style,
                                                    #     className="name-input-container",
                                                    #     children=[
                                                    #         html.P("Number of Laterals", className="input-title"),
                                                    #         dcc.Input(id="n-laterals-select", 
                                                    #                     disabled=True,
                                                    #                     value=start_vals_hdf5["n-laterals"], type='number', min=0, max=20, step=1, className="input-box"),
                                                    # ]),
                                                    # html.Div(
                                                    #     id="lat-allocation-div",
                                                    #     style=div_none_style,
                                                    #     className="name-input-container",
                                                    #     children=[
                                                    #         html.P("Lateral Flow Allocation", className="input-title"),
                                                    #         dcc.Input(id="lateral-flow-select", 
                                                    #                     disabled=True,
                                                    #                     value=start_vals_hdf5["lateral-flow"], type='number', min=0, max=1, step=0.01, className="input-box"),
                                                    # ]),
                                                    # html.Div(
                                                    #     id="lat-flow-mul-div",
                                                    #     style=div_none_style,
                                                    #     className="name-input-container",
                                                    #     children=[
                                                    #         html.P("Lateral Flow Multiplier", className="input-title"),
                                                    #         dcc.Input(id="lateral-multiplier-select", 
                                                    #                     disabled=True, 
                                                    #                     value=start_vals_hdf5["lateral-multiplier"], type='number', min=0, max=1, step=0.05, className="input-box"),
                                                    # ]),

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
                                                        slider2(DivID="drillcost-div", ID="drillcost-select", ptitle="Drilling Cost ($/m)", min_v=0, max_v=4000, 
                                                                mark_dict=drillcost_dict, start_v=start_vals_econ["drillcost"], div_style=div_block_style, parameter_name="Drilling Cost ($/m)"),
                                                        slider2(DivID="discount-rate-div", ID="discount-rate-select", ptitle="Discount Rate (%)", min_v=0, max_v=20, 
                                                                mark_dict=discount_dict, start_v=start_vals_econ["discount-rate"], div_style=div_block_style, parameter_name="Discount Rate (%)"),
                                                        # slider2(DivID="lifetime-div", ID="lifetime-select", ptitle="Lifetime (years)", min_v=10, max_v=40, 
                                                        #         mark_dict=lifetime_dict, start_v=30),
                                                        slider1(DivID="lifetime-div", ID="lifetime-select", ptitle="Lifetime (years)", min_v=10, max_v=40, 
                                                                mark_dict=lifetime_dict, step_i=1, start_v=start_vals_econ["lifetime"], div_style=div_block_style, parameter_name="Lifetime (years)"),
                                                        slider2(DivID="kwt-div", ID="kwt-select", ptitle="Plant CAPEX ($/kWt)", min_v=0, max_v=1000, 
                                                                mark_dict=kwt_dict, start_v=start_vals_econ["kwt"], div_style=div_block_style, parameter_name="Plant CAPEX ($/kWt)"),
                                                        slider2(DivID="kwe-div", ID="kwe-select", ptitle="Plant CAPEX ($/kWe)", min_v=0, max_v=10000, 
                                                                mark_dict=kwe_dict, start_v=start_vals_econ["kwe"], div_style=div_block_style, parameter_name="Plant CAPEX ($/kWe)"),
                                                    ]
                                                    ),
                                            html.Div(id="sCO2-card",
                                                    className="params-div",
                                                    children=[
                                                        html.P("ⓘ Multiple LCOE minima exist. Dial here to explore:", id="sCO2-text"),  # Run the optimizer 
                                                        slider2(DivID="precool-div", ID="precool-select", ptitle="Pre-cooling (˚C)", min_v=0, max_v=40, 
                                                                mark_dict=precool_dict, start_v=start_vals_econ["precool"], div_style=div_block_style, parameter_name="Pre-cooling (˚C)"),
                                                        slider2(DivID="turb-pout-div", ID="turb-pout-select", ptitle="Turbine Outlet Pressure (bar)", min_v=75, max_v=200, 
                                                                mark_dict=turb_pout_dict, start_v=start_vals_econ["turb-pout"], div_style=div_block_style, parameter_name="Turbine Outlet Pressure (bar)"),
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
                                                    slider1(DivID="mesh-div", ID="mesh-select", ptitle="Mesh Fineness", min_v=0, max_v=2, 
                                                                mark_dict=fineness_dict, step_i=1, start_v=start_vals_sbt["mesh"], div_style=div_block_style, parameter_name="Mesh Fineness"),
                                                    slider1(DivID="accuracy-div", ID="accuracy-select", ptitle="Accuracy", min_v=1, max_v=5, 
                                                                mark_dict=accuracy_dict, step_i=1,start_v=start_vals_sbt["accuracy"], div_style=div_block_style, parameter_name="Accuracy"),
                                                
                                                    html.Div(
                                                            id="hyperparam1-container",
                                                            children=[
                                                                dropdown_box(DivID="mass-flow-mode-div", ID="mass-mode-select", ptitle="Mass Flow Rate Mode", 
                                                                                                options=["Constant", "Variable"], disabled=True, div_style=div_block_style)
                                                        ]),
                                                    html.Div(
                                                        id="hyperparam3-container",
                                                        children=[
                                                                dropdown_box(DivID="temp-flow-mode-div", ID="temp-mode-select", ptitle="Injection Temperature Mode", 
                                                                                        options=["Constant", "Variable"], disabled=True, div_style=div_block_style)
                                                        ]),
                                                    html.Div(
                                                        id="hyperparam5-container",
                                                        children=[
                                                                dropdown_box(DivID="fluid-mode-div", ID="fluid-mode-select", ptitle="Fluid Properties Mode", 
                                                                                                options=["Variable", "Constant"], disabled=True, div_style=div_none_style)
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
