#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# -------------------------------------------------------------------------------------
# GeoCLUSTER, data tool and app that runs visuals on Closed-Loop Geothermal Data
# -------------------------------------------------------------------------------------

# TODO: get_Ts_diagram debugging, tab checks, visible and invisible checks
# TODO: really it should only be running the plotting when fluid or end-use == "All"?
# running unncessary amount of times maybe? This means running things more times

# --------------------
# Libraries.
# --------------------

# import pkg_resources

# installed_packages = pkg_resources.working_set
# installed_packages_list = sorted(["%s==%s" % (i.key, i.version)
#    for i in installed_packages])
# print(installed_packages_list)
# print(help("modules"))

import time

# web app and interactive graphics libraries 
import dash
from dash.dependencies import Input, Output, State
from dash import Dash, dcc, html, ctx
import dash_daq as daq                  # Adds more data acquisition (DAQ) and controls to dash callbacks 
import dash_bootstrap_components as dbc # Adds bootstrap components for more web themes and templates
from dash.exceptions import PreventUpdate
from flask_talisman import Talisman
from flask import send_from_directory
from flask_compress import Compress

# sourced scripts
from paths import inpath_dict
from write2excel import write_excelsheet
from sliders import * # u_sCO2, u_H2O, c_sCO2, c_H2O, and imports functions from plots.py
from dropdowns import *
from text import *
from tables import generate_summary_table
from plots import generate_econ_lineplots, generate_subsurface_lineplots, generate_subsurface_contours
from info_popups import PARAMETER_INFO, create_info_button, create_enhanced_slider, create_enhanced_input_box, create_enhanced_dropdown, create_info_modal, register_info_modal_callbacks
from unit_conversions import unit_converter, convert_value, get_unit_symbol
from unit_preferences import create_unit_preferences_card, get_unit_preferences_from_inputs, apply_metric_units, apply_imperial_units


# -----------------------------------------------------------------------------
# Create dash app with CSS styles and HTML components.
#
#   initalizer:  dash.Dash() 
#   layout:      Describes the HTML layout of the app.
#   callbacks:   Interactivity of the app.
#
# -----------------------------------------------------------------------------

requests_pathname_prefix = inpath_dict["requests_pathname_prefix"]
url_base_pathname = inpath_dict["url_base_pathname"]
geoCLUSTER_results_pathname = inpath_dict["geoCLUSTER_results_pathname"]
properties_H2O_pathname = inpath_dict["properties_H2O_pathname"]
properties_CO2v2_pathname = inpath_dict["properties_CO2v2_pathname"]
additional_properties_CO2v2_pathname = inpath_dict["additional_properties_CO2v2_pathname"]
tmatrix_pathname = inpath_dict["tmatrix_pathname"]

# sha384 for bootstrap@5.3.1/dist/css/bootstrap.min.css (jsDelivr)
BOOTSTRAP_CSS = {
    "rel": "stylesheet",
    "href": "https://cdn.jsdelivr.net/npm/bootstrap@5.3.1/dist/css/bootstrap.min.css",
    "integrity": "sha384-4bw+/aepP/YC94hEpVNVgiZdgIC5+VKNBQNGCHeKRQN+PtmoHDEXuppvnDJzQIu9",
    "crossorigin": "anonymous",
}
# print(dbc.themes.BOOTSTRAP)


app = dash.Dash(__name__, assets_folder='assets', 
                # external_stylesheets=[dbc.themes.BOOTSTRAP], 
                external_stylesheets=[BOOTSTRAP_CSS],          # ← use the dict, not dbc.themes.BOOTSTRAP
                # url_base_pathname=url_base_pathname, # not needed
                requests_pathname_prefix=requests_pathname_prefix,
                suppress_callback_exceptions=True,
                meta_tags=[{'name': 'viewport', 'content': 'width=device-width'},
                           {'name': 'description', 
                            'content': """Welcome to the Geothermal Closed-Loop User Tool in Energy Research (GeoCLUSTER). 
                                        This research was funded by the Geothermal Technologies Office (GTO) within the 
                                        Office of Energy Efficiency and  Renewable Energy (EERE) at the U.S. Department of 
                                        Energy (DOE) to form a collaborative study of closed-loop geothermal systems."""}
                            ])
# Talisman(
#     app.server,                 # <-- the Flask instance
#     force_https=True,           # redirects HTTP → HTTPS (308)
#     strict_transport_security=True,
#     strict_transport_security_max_age=31536000,     # one year
#     strict_transport_security_include_subdomains=True,
#     strict_transport_security_preload=True           # for Chrome preload list
# )

app.title = "GeoCLUSTER | Geothermal Closed-Loop User Tool in Energy Research"

# -----------------------------------------------
# Globl styles, configurations, and colors.
# -----------------------------------------------

tabs = ["about-tab", "energy-time-tab", "energy-tab", "econ-time-tab", "summary-tab"] # there are values and ids


plotly_config = {'displaylogo': False,
                'modeBarButtonsToRemove': ['autoScale', 'resetScale'], # High-level: zoom, pan, select, zoomIn, zoomOut, autoScale, resetScale
                'toImageButtonOptions': {
                    'format': 'png', # one of png, svg, jpeg, webp
                    'filename': 'custom_image',
                    'height': None,
                    'width': None,
                    'scale': 6 # Multiply title/legend/axis/canvas sizes by this factor
                        }
                  }


darkergrey = "#383838"
lightbrown = '#ede6dd'

dropdown_guidance_style = {'position': 'relative',
                            'marginTop': '28px',
                            'height': '25px',
                            'width': '100%',
                            "paddingBottom": "7px",
                            "paddingLeft": "20px",
                            "paddingRight": "20px",
                            'border':'1.5px white solid',
                            'backgroundColor': lightbrown,
                            "color": darkergrey,
                            "fontWeight": "bold",
                            "fontSize": "11px"}

carousel_image_style = {'objectFit': 'contain'}


# -----------------------------------------------
# HTML components.
# -----------------------------------------------

def description_card():

    # -----------------------------------------------------------------------
    # A Div containing dashboard title & descriptions.
    # -----------------------------------------------------------------------

    return html.Div(
        id="description-card",
        children=[
            html.Img(id="logo", src=app.get_asset_url('logo3.png')),
            html.Hr(id="hr-break1"),
        ],
    )

def generate_control_card():

    # -----------------------------------------------------------------------
    # A Div containing controls for graphs. Controls include the following:
    #
    #   Scenario buttons, sliders, and dropdowns.
    # -----------------------------------------------------------------------

    return html.Div(
        id="control-card",
        children=[

            html.Div(id='scenario1-div',
                     children=[

                        html.Button(html.Div(children=[
                                            html.P('Run At Commercial Scale'),
                                            html.P(scenario1_text, id="scenario1_text")]),
                        id='btn-nclicks-1', n_clicks=0),
            ]),
            # html.Div(id='scenario2-div',
            #          children=[
            #             html.Button('Optimize Power Output', 
            #             id='btn-nclicks-2', n_clicks=0),
            #             html.P(scenario2_text, id="scenario2_text"), 
            # ]),
            html.Div(id='scenario3-div',
                     children=[
                        html.Button(html.Div(children=[
                                            html.P('Optimize Economic Competitiveness'),
                                            html.P(scenario3_text, id="scenario3_text")]),
                            # 'Optimize Economic Competitiveness', 
                        id='btn-nclicks-3', n_clicks=0),
                        # html.P(scenario3_text, id="scenario3_text"),
            ]),            
            html.Hr(id="hr-break2"), 
            dropdown_card(),
            html.Br(),
            html.Br(),
            html.Div(id="slider-card", children=slider_card()),
            # disclaimer
            html.P(disclaimer_text, id='disclaimer-text'),
        ], 
    )



def graph_guidance_card(btnID, cardID, dropdown_children):

    # -----------------------------------------------------------------------
    # A Div containing graph guidance dropdowns for extra user context.
    # -----------------------------------------------------------------------

    return html.Div(
        [
            dbc.Button(
                "ⓘ Graph Guidance",
                id=btnID,
                color="primary",
                n_clicks=0,
                style=dropdown_guidance_style
            ),
            dbc.Collapse(
                dbc.Card(
                    children=[
                        dbc.CardBody(
                            children=[
                                dropdown_children
                            ]
                        )]
                    ),
                id=cardID,
                is_open=False,
            ),
        ]
    )

# -------------------------------------------------------------------------------------------------
# Defines tab contents for 5 of the following tabs:
# 
#   About, Subsurface Results, Subsurface Contours, Economic Results, and Summary
#
# -------------------------------------------------------------------------------------------------

about_tab = dcc.Tab(label='About',
                        id="about-tab",
                        value='about-tab',
                        selected_className="active_tabs",
                        children=[
                                html.Hr(className="tab-hr"),
                                html.Div(
                                    id="demo",
                                    children=[
                                                html.Div(id='intro2', 
                                                        children=[html.P("About Our Research", id='ab-title1'),
                                                                    html.P(note, id='ab-note'), 
                                                                    html.P(note2, id='ab-note2'),
                                                                    html.P("Navigating the Results", id='ab-title2'),
                                                                    html.P(note3, id='ab-note3'), 
                                                                    html.P("Resources", id='ab-title3'),
                                                                    html.Label([html.P("Download the contributing", id="shorttext1"),
                                                                                html.A('papers', href='https://gdr.openei.org/submissions/1473', id='hyperlink1'),
                                                                                html.P("and", id="shorttext2"),
                                                                                html.A('code', href='https://github.com/pnnl/GeoCLUSTER', id='hyperlink2'),
                                                                                html.P(".", id="shorttext3"),
                                                                                ], id='ab-note4')
                                                ]),
                                                html.Div(id='image-container', 
                                                        children=[html.Img(id="cluster-img", src=app.get_asset_url('CLGWG3.png'),),
                                                                    dbc.Carousel(id="carousel-ride",
                                                                        items=[
                                                                            {"key": "1", 
                                                                                "src": app.get_asset_url('CLGS-1.3.png'), 
                                                                                "img_style": carousel_image_style
                                                                                },
                                                                            {"key": "2", 
                                                                                "src": app.get_asset_url('CLGS-3.2.png'),
                                                                                "img_style": carousel_image_style
                                                                            },
                                                                        ],
                                                                        variant="dark",
                                                                        ride="carousel",
                                                                        interval=3000,
                                                                        controls=True,
                                                                        indicators=False),
                                                                    ]
                                                        )
                                                ],
                                    ),
                ])


energy_time_tab = dcc.Tab(label='Subsurface Results',
                        id="subsurface-tab",
                        value='energy-time-tab',
                        selected_className="active_tabs",
                        children=[
                            html.Div(className="extra-space"),
                            graph_guidance_card(btnID="collapse-button2", cardID="collapse2", dropdown_children=html.P(dropdown_text1)),
                            html.Div(id="error_block_div1"), 
                            html.Br(),
                            html.Div(
                                    id="graphics-parent",
                                    children=[
                                            html.Div(id="graphics-container",
                                            children=[
                                                    dcc.Graph(id="geothermal_time_plots",
                                                            config=plotly_config),
                                                    dbc.RadioItems(
                                                            options=[
                                                                {"label": "Auto Scale", "value": 1},
                                                                {"label": "Full Scale", "value": 2},
                                                            ],
                                                            value=1,
                                                            id="radio-graphic-control3",
                                                        ),
                                            ]

                                        )
                                ])

                ])


energy_tab = dcc.Tab(label='Subsurface Contours',
                        id="contour-tab",
                        value='energy-tab',
                        selected_className="active_tabs",
                        children=[
                            html.Div(className="extra-space"),
                            # html.Hr(),
                            html.Div(id="dropdown-card6",
                                    children=[
                                        graph_guidance_card(btnID="collapse-button", cardID="collapse", dropdown_children=html.P(dropdown_text2)),
                                        html.Div(id="error_block_div2"), 
                                        html.P("Select Parameter", id='select-text'),
                                        dcc.Dropdown(
                                            id="param-select",
                                            options=[{"label": i, "value": i} for i in param_list],
                                            value="Horizontal Extent (m)",
                                            clearable=False,
                                            searchable=False,
                                        ),
                                    ]
                                ), 
                        dcc.Graph(id="geothermal_plots",
                                    config=plotly_config)

                ])

economics_time_tab = dcc.Tab(label='Economic Results',
                        id="econ-tab",
                        value='economics-time-tab',
                        selected_className="active_tabs",
                            children=[
                            html.Div(className="extra-space"),
                            # html.Hr(),
                            graph_guidance_card(btnID="collapse-button4", cardID="collapse4", 
                                                dropdown_children=
                                                    html.Div(
                                                        children=[
                                                            html.P(dropdown_text1.replace("5 options", "7 options")),
                                                            html.P(dropdown_econ_text1),
                                                            html.Img(id="formulas-img", src=app.get_asset_url('lcoe-lcoh-formulas.png')),
                                                            dcc.Markdown(dropdown_econ_markdown_text1, mathjax=True),
                                                            html.P(dropdown_econ_text2)
                                                            ]
                                                        )
                                                    ),
                            html.Div(id="warning_block_div3"),
                            html.Div(id="error_block_div3"), 
                            html.Br(),
                            html.Div(
                                    id="graphics-parent-econ",
                                    children=[
                                        html.Div(id="graphics-container-econ",
                                                children=[
                                                    dcc.Graph(id="econ_plots",
                                                                config=plotly_config),
                                                    dbc.RadioItems(
                                                                options=[
                                                                    {"label": "Auto Scale", "value": 1},
                                                                    {"label": "Full Scale", "value": 2},
                                                                ],
                                                                value=1,
                                                                id="radio-graphic-control4",
                                                            ),
                                                    html.P("Temperature-entropy (T-S) Diagram", id="ts-text"),
                                                    # html.P("sCO2 only", id="ts-subtext")
                                                    # dcc.Graph(id="ts_plot",
                                                    #             config=plotly_config),
                                                ]

                                            )
                                    ])

                ])

summary_tab = dcc.Tab(label='Summary',
                        id="summary-tab",
                        value='summary-tab',
                        selected_className="active_tabs",
                        children=[
                                html.Hr(className="tab-hr"),
                                html.Button("Download Results", id="btn_xlsx"),
                                dcc.Download(id="download-dataframe-xlsx"),
                                html.Br(),
                                html.Br(),
                                dcc.Graph(id="table",
                                        config=plotly_config)
                ])

# -----------------------------------------------------------------------------
# Define dash app layout here.
# -----------------------------------------------------------------------------

app.layout = html.Div(
    id="app-container",
    children=[
        dcc.Store(id='econ-memory'),
        dcc.Store(id='econ-results'),
        dcc.Store(id='econ-errors'),
        dcc.Store(id='thermal-memory'),
        dcc.Store(id='thermal-results-mass'),
        dcc.Store(id='thermal-results-time'),
        dcc.Store(id='thermal-results-errors'),
        dcc.Store(id='thermal-contours-errors'),
        dcc.Store(id='summary-memory'),
        dcc.Store(id='TandP-data'),

        # Left column
        html.Div(
            id="left-column",
            className="four columns",
            children=[description_card(), 
                      generate_control_card()
                      ]
        ),
        # Right column
        html.Div(
            id="right-column",
            className="eight columns",
            children=[
                html.Div(
                    id="geothermal_card",
                    children=#[dcc.Tabs(id='tabs', value='tab-1', children=[])]
                              [dcc.Tabs(id="tabs", 
                                        value='about-tab', 
                                        children=[about_tab,
                                                energy_time_tab,
                                                energy_tab,
                                                economics_time_tab,
                                                summary_tab,
                              ])                    ],
                ),
            ],
        ),
        
        # Information Modal
        create_info_modal(),
    ],
)

# Register info popup callbacks
register_info_modal_callbacks(app)

# -----------------------------------------------------------------------------
# Define dash app callbacks begin here.
# -----------------------------------------------------------------------------

@app.callback(
    Output(component_id="download-dataframe-xlsx", component_property="data"),
    [Input(component_id="btn_xlsx", component_property="n_clicks"),
     Input(component_id='summary-memory', component_property='data'),
     Input(component_id='thermal-results-mass', component_property='data'),
     Input(component_id='thermal-results-time', component_property='data'),
     Input(component_id='econ-results', component_property='data')],
    prevent_initial_call=True,
)

def generate_summary(n_clicks, df_tbl, df_mass_flow_rate, df_time, df_econ):

    if "btn_xlsx" == ctx.triggered_id:
        try:
            print("Download button clicked - starting Excel generation...")
            print(f"Summary data type: {type(df_tbl)}")
            print(f"Summary data: {df_tbl}")
            print(f"Thermal mass data: {df_mass_flow_rate}")
            print(f"Thermal time data: {df_time}")
            print(f"Economic data: {df_econ}")
            
            # Check if we have valid data
            if df_tbl is None:
                print("Warning: Summary data is None")
                df_tbl = {}
            
            if df_mass_flow_rate is None:
                print("Warning: Thermal mass data is None")
                df_mass_flow_rate = []
                
            if df_time is None:
                print("Warning: Thermal time data is None")
                df_time = []
                
            if df_econ is None:
                print("Warning: Economic data is None")
                df_econ = {}
            
            # Create the Excel file
            write_excelsheet(df_summary=df_tbl, df_subsurf_res_mass=df_mass_flow_rate, 
                                    df_subsurf_res_time=df_time, df_econ=df_econ, 
                                            geoCLUSTER_results_pathname=geoCLUSTER_results_pathname)
            
            print(f"Excel file created successfully at: {geoCLUSTER_results_pathname}")
            
            # Try to read the file and send as data instead of file path (firewall-friendly)
            try:
                import pandas as pd
                # Read the Excel file we just created
                excel_data = pd.read_excel(geoCLUSTER_results_pathname, sheet_name=None)
                
                # Get the main summary sheet for download
                if 'SBT Parameters & Results' in excel_data:
                    main_df = excel_data['SBT Parameters & Results']
                elif 'Summary' in excel_data:
                    main_df = excel_data['Summary']
                else:
                    # Fallback to first available sheet
                    main_df = list(excel_data.values())[0]
                
                print("Sending Excel data as DataFrame (firewall-friendly method)")
                # Use send_data_frame which is more reliable and firewall-friendly
                return dcc.send_data_frame(
                    main_df.to_excel, 
                    filename="geoCLUSTER_results.xlsx", 
                    index=False
                )
                
            except Exception as file_error:
                print(f"Error reading Excel file: {file_error}")
                print("Falling back to file path method...")
                # Fallback to original method
                return dcc.send_file(geoCLUSTER_results_pathname)
            
        except Exception as e:
            print(f"Error in download callback: {str(e)}")
            import traceback
            traceback.print_exc()
            # Return a simple error message instead of crashing
            return None
    else:
        raise PreventUpdate


@app.callback(
    Output(component_id="collapse", component_property="is_open"),
    [Input(component_id="collapse-button", component_property="n_clicks")],
    [State(component_id="collapse", component_property="is_open")],
)
def toggle_collapse(n, is_open):

    # ----------------------------------------------
    # Subsurface contours graph guidance.
    # ----------------------------------------------

    if n:
        return not is_open
    return is_open

@app.callback(
    Output(component_id="collapse2", component_property="is_open"),
    [Input(component_id="collapse-button2", component_property="n_clicks")],
    [State(component_id="collapse2", component_property="is_open")],
)
def toggle_collapse(n, is_open):

    # ----------------------------------------------
    # Subsurface results graph guidance.
    # ----------------------------------------------

    if n:
        return not is_open
    return is_open


@app.callback(
    Output(component_id="collapse3", component_property="is_open"),
    [Input(component_id="collapse-button3", component_property="n_clicks")],
    [State(component_id="collapse3", component_property="is_open")],
)
def toggle_collapse(n, is_open):

    # ----------------------------------------------
    # Open to finetune system values.
    # ----------------------------------------------

    if n:
        return not is_open
    return is_open


@app.callback(
    Output(component_id="collapse4", component_property="is_open"),
    [Input(component_id="collapse-button4", component_property="n_clicks")],
    [State(component_id="collapse4", component_property="is_open")],
)
def toggle_collapse(n, is_open):

    # ----------------------------------------------
    # Economic results graph guidance.
    # ----------------------------------------------

    if n:
        return not is_open
    return is_open

@app.callback(
   [Output(component_id="scenario1-div", component_property="style"),
   Output(component_id="scenario3-div", component_property="style"),
   Output(component_id="hr-break1", component_property="style"),
   ],
   [Input(component_id="model-select", component_property="value")
    ],
   prevent_initial_call=True
   )

def update_tabs(selected_model):

    if selected_model == "HDF5": 

        return {'display': 'block'}, {'display': 'block'}, {'display': 'block'}

    elif selected_model == "SBT V1.0" or selected_model == "SBT V2.0": 
        
        return {'display': 'block'}, {'display': 'block'}, {'display': 'block'}

    elif selected_model == "CovHDF5":
        
        return {'display': 'none'}, {'display': 'none'}, {'display': 'none'}

    else:
        raise PreventUpdate


@app.callback(
#    [
    Output(component_id="graphics-parent", component_property="children"),
    # Output(component_id="graphics-parent-econ", component_property="children"),
#    ],
   Input(component_id="model-select", component_property="value"),
   prevent_initial_call=True
   )

def update_loading(selected_model):

    if selected_model == "HDF5": 

        return html.Div(id="graphics-container",
                                            children=[
                                                    dcc.Graph(id="geothermal_time_plots",
                                                            config=plotly_config),
                                                    dbc.RadioItems(
                                                            options=[
                                                                {"label": "Auto Scale", "value": 1},
                                                                {"label": "Full Scale", "value": 2},
                                                            ],
                                                            value=1,
                                                            id="radio-graphic-control3",
                                                        ),
                                            ]
                                        )

    elif selected_model == "SBT V1.0" or selected_model == "SBT V2.0": 
        
        return dcc.Loading(
                parent_className="loader-wrapper",
                type='circle',
                children=[
                        html.Div(id="graphics-container",
                            children=[
                                    dcc.Graph(id="geothermal_time_plots",
                                            config=plotly_config),
                                    dbc.RadioItems(
                                            options=[
                                                {"label": "Auto Scale", "value": 1},
                                                {"label": "Full Scale", "value": 2},
                                            ],
                                            value=1,
                                            id="radio-graphic-control3",
                                        ),
                            ]

                    )
                ]
            )
    else:
        raise PreventUpdate


@app.callback(
   Output(component_id="contour-tab", component_property="style"),
   [Input(component_id="model-select", component_property="value")
    ],
   prevent_initial_call=True
   )

def update_tabs(selected_model):

    if selected_model == "HDF5": 

        return {'display': 'block'}

    elif selected_model == "SBT V1.0" or selected_model == "SBT V2.0" : 
            
        return {'display': 'none'}


@app.callback(
   Output(component_id="fluid-select", component_property="value", allow_duplicate=True),
   [Input(component_id="tabs", component_property="value"),
    Input(component_id='btn-nclicks-1', component_property='n_clicks'),
    Input(component_id='btn-nclicks-3', component_property='n_clicks'),
    Input(component_id="fluid-select", component_property="value")
    ],
    prevent_initial_call=True
    )

def retain_entry_between_tabs(tab, btn1, btn3, fluid):

    if tab != "energy-tab" and fluid == "All":

        if "btn-nclicks-1" == ctx.triggered_id:
            return "H2O"

        if "btn-nclicks-3" == ctx.triggered_id:
            return "H2O"
        else:
            raise PreventUpdate
    else:
        raise PreventUpdate

@app.callback(
    Output(component_id="end-use-select", component_property="value"),
   [Input(component_id="tabs", component_property="value"),
    Input(component_id='btn-nclicks-1', component_property='n_clicks'),
    Input(component_id='btn-nclicks-3', component_property='n_clicks'),
    Input(component_id="end-use-select", component_property="value")
    ],
    prevent_initial_call=True
    )

def flip_to_tab(tab, btn1, btn3, end_use):

    if tab != "energy-tab" and end_use == "All":

        if "btn-nclicks-1" == ctx.triggered_id:
            return "Heating"

        if "btn-nclicks-3" == ctx.triggered_id:
            return "Heating"
        else:
            raise PreventUpdate

    else:
        raise PreventUpdate


@app.callback(
    [Output(component_id="fluid-select", component_property="value", allow_duplicate=True),
    Output(component_id="fluid-select", component_property="options", allow_duplicate=True)
    ],
    [Input(component_id="model-select", component_property="value"),
    ],
    
    prevent_initial_call=True
    )

def update_working_fluid(model):

    if model == "SBT V2.0":
        fluid_list = ["All", "H2O", "sCO2"]
        if ctx.triggered_id == "model-select":
            return "All", [{"label": i, "value": i} for i in fluid_list]
        else:
            raise PreventUpdate
    elif model == "SBT V1.0":
        fluid_list = ["H2O"]
        if ctx.triggered_id == "model-select":
            return "H2O", [{"label": i, "value": i} for i in fluid_list]
        else:
            raise PreventUpdate
    elif model == "CovHDF5":
        fluid_list = ["H2O"]  # CovHDF5 only supports H2O
        if ctx.triggered_id == "model-select":
            return "H2O", [{"label": i, "value": i} for i in fluid_list]
        else:
            raise PreventUpdate
    else:
        raise PreventUpdate


@app.callback(
   [Output(component_id='fluid-select', component_property='options'),
    Output(component_id='fluid-select', component_property='value'),
    Output(component_id='interpolation-select', component_property='options'),
    Output(component_id='interpolation-select', component_property='value')
   ],
   [Input(component_id="tabs", component_property="value"),
    Input(component_id='fluid-select', component_property='value'),
    Input(component_id="model-select", component_property="value")
    ])

def change_dropdown(at, fluid, model):
    
    if model == "HDF5":

        if at == "energy-time-tab":
            fluid_list = ["All", "H2O", "sCO2"]
            interpol_list = ["True"]
            return [{"label": i, "value": i} for i in fluid_list], fluid, [{"label": i, "value": i} for i in interpol_list], interpol_list[0]

        elif at == "about-tab":
            fluid_list = ["All", "H2O", "sCO2"]
            interpol_list = ["True"]
            return [{"label": i, "value": i} for i in fluid_list], fluid, [{"label": i, "value": i} for i in interpol_list], interpol_list[0]
    
        elif at == "energy-tab":
            fluid_list = ["H2O", "sCO2"]
            interpol_list = ["True"]
            if fluid != "All":
                return [{"label": i, "value": i} for i in fluid_list], fluid, [{"label": i, "value": i} for i in interpol_list], interpol_list[0]
            else:
                return [{"label": i, "value": i} for i in fluid_list], fluid_list[0], [{"label": i, "value": i} for i in interpol_list], interpol_list[0]

        elif at == "economics-time-tab":
            fluid_list = ["All", "H2O", "sCO2"]
            interpol_list = ["True"]
            return [{"label": i, "value": i} for i in fluid_list], fluid, [{"label": i, "value": i} for i in interpol_list], interpol_list[0]

        elif at == "summary-tab":
            fluid_list = ["All", "H2O", "sCO2"]
            interpol_list = ["True"]
            return [{"label": i, "value": i} for i in fluid_list], fluid, [{"label": i, "value": i} for i in interpol_list], interpol_list[0]
            # raise PreventUpdate

    elif model == "CovHDF5":
        # CovHDF5 only supports H2O fluid
        fluid_list = ["H2O"]
        interpol_list = ["True"]  # CovHDF5 uses interpolation
        
        if at == "energy-time-tab":
            return [{"label": i, "value": i} for i in fluid_list], "H2O", [{"label": i, "value": i} for i in interpol_list], interpol_list[0]
        elif at == "about-tab":
            return [{"label": i, "value": i} for i in fluid_list], "H2O", [{"label": i, "value": i} for i in interpol_list], interpol_list[0]
        elif at == "energy-tab":
            return [{"label": i, "value": i} for i in fluid_list], "H2O", [{"label": i, "value": i} for i in interpol_list], interpol_list[0]
        elif at == "economics-time-tab":
            return [{"label": i, "value": i} for i in fluid_list], "H2O", [{"label": i, "value": i} for i in interpol_list], interpol_list[0]
        elif at == "summary-tab":
            return [{"label": i, "value": i} for i in fluid_list], "H2O", [{"label": i, "value": i} for i in interpol_list], interpol_list[0]
        else:
            raise PreventUpdate

    if model == "SBT V1.0" or model == "SBT V2.0":

        raise PreventUpdate


@app.callback(
   Output(component_id="tabs", component_property="value"),
   [Input(component_id="tabs", component_property="value"),
    Input(component_id='btn-nclicks-1', component_property='n_clicks'),
     # Input(component_id='btn-nclicks-2', component_property='n_clicks'),
     Input(component_id='btn-nclicks-3', component_property='n_clicks')
     ])

def flip_to_tab(tab, btn1, btn3):

    if tab == "about-tab":

        if "btn-nclicks-1" == ctx.triggered_id:
            return "economics-time-tab"

        # if "btn-nclicks-2" == ctx.triggered_id:
        #     return "energy-time-tab"

        if "btn-nclicks-3" == ctx.triggered_id:
            return "economics-time-tab"
        else:
            raise PreventUpdate

    else:
        raise PreventUpdate


@app.callback(
    [Output(component_id='mdot-select', component_property='value'),
     Output(component_id='L2-select', component_property='value'),
     Output(component_id='L1-select', component_property='value'),
     Output(component_id='grad-select', component_property='value'),
     Output(component_id='diameter-select', component_property='value'),
     Output(component_id='Tinj-select', component_property='value'),
     Output(component_id='k-select', component_property='value'),

     Output(component_id='drillcost-select', component_property='value'),
     Output(component_id='discount-rate-select', component_property='value'),
     Output(component_id='lifetime-select', component_property='value'),
     Output(component_id='kwt-select', component_property='value'),
     Output(component_id='kwe-select', component_property='value'),
     Output(component_id='precool-select', component_property='value'),
     Output(component_id='turb-pout-select', component_property='value'),
     
     Output(component_id='Tsurf-select', component_property='value'),
     Output(component_id='c-select', component_property='value'),
     Output(component_id='rho-select', component_property='value'),
    #  Output(component_id='radius-vertical-select', component_property='value'),
    #  Output(component_id='radius-lateral-select', component_property='value'),
     Output(component_id='n-laterals-select', component_property='value'),
     Output(component_id='lateral-flow-select', component_property='value'),
     Output(component_id='lateral-multiplier-select', component_property='value'),

   ],
    [Input(component_id='btn-nclicks-1', component_property='n_clicks'),
     # Input(component_id='btn-nclicks-2', component_property='n_clicks'),
     Input(component_id='btn-nclicks-3', component_property='n_clicks'),
     Input(component_id="tabs", component_property="value"),
     Input(component_id='case-select', component_property='value'),
     Input(component_id='fluid-select', component_property='value'),
     Input(component_id='end-use-select', component_property='value'),

     Input(component_id='model-select', component_property='value')
     ]
)

def update_slider_with_btn(btn1, btn3, at, case, fluid, end_use, model):

    # ----------------------------------------------------------------------------------------------
    # Defines scenario button values when clicked.  
    #
    # Commercial Scale:
    # For example, consider a u-shaped configuration installed in a geothermal reservoir, with a 
    # horizontal extent of 10 km, vertical depth of 2.5 km, geothermal gradient of 0.065 ˚C/m, borehole 
    # diameter of 0.3 m, inlet temperature of 40˚C, and rock thermal conductivity of 3.5 W/m K, 
    # parameters within the range of anticipated values for commercial scale geothermal applications
    #
    #
    # ----------------------------------------------------------------------------------------------

    default_output = (25, # Tsurf
                      790, # c
                      2750, # rho
                      1, # n-laterals
                      1, # lateral-flow
                      1 # lateral-multiplier
                      )

    if model == "HDF5":
        # output = ('utube', 24, 10000, 3500, 0.050, 0.35, 30, 3, 1000, 7.0, 40, 100, 3000, 13, 80) # return to default

        if "btn-nclicks-1" == ctx.triggered_id:
            output = (24, 10000, 2500, 0.050, 0.30, 40, 3.5,  1000, 7.0, 40, 100, 3000, 13, 80)
            return output + default_output
        
        # elif "btn-nclicks-2" == ctx.triggered_id: 
        #     if case == 'coaxial':
        #         output = (77, 20000, 5000, 0.070, 0.44, 30, 4.5,  1000, 7.0, 40, 100, 3000, 13, 80)
        #         return output
        #     if case == "utube":
        #         output = (100, 20000, 5000, 0.070, 0.44, 30, 4.5,  1000, 7.0, 40, 100, 3000, 13, 80)
        #         return output

        elif "btn-nclicks-3" == ctx.triggered_id:

            if case == 'coaxial' and fluid == "H2O" or fluid == "All":

                if end_use == "Electricity":
                    output = (39.2, 20000, 5000, 0.070, 0.444, 60, 4.5,  1000, 7.0, 40, 100, 3000, 13, 80)
                if end_use == "Heating" or end_use == "All":
                    output = (73.4, 13000, 5000, 0.070, 0.444, 30, 4.5,  1000, 7.0, 40, 100, 3000, 13, 80)
                return output + default_output

            if case == 'utube' and fluid == "H2O" or fluid == "All":

                if end_use == "Electricity":
                    output = (43, 20000, 5000, 0.070, 0.44, 60, 4.5,  1000, 7.0, 40, 100, 3000, 13, 80)
                if end_use == "Heating" or end_use == "All":
                    output = (100, 20000, 5000, 0.070, 0.44, 30, 4.5,  1000, 7.0, 40, 100, 3000, 13, 80)
                return output + default_output

            if case == 'coaxial' and fluid == "sCO2":
                
                if end_use == "Electricity":
                    output = (69.6, 13000, 5000, 0.070, 0.44, 60, 4.5,  1000, 7.0, 40, 100, 3000, 13, 80)
                if end_use == "Heating" or end_use == "All":
                    output = (69.6, 6000, 5000, 0.070, 0.44, 45, 4.5,  1000, 7.0, 40, 100, 3000, 13, 80)
                return output + default_output
                
            if case == 'utube' and fluid == "sCO2":
                
                if end_use == "Electricity":
                    output = (100, 20000, 5000, 0.070, 0.44, 60, 4.5,  1000, 7.0, 40, 100, 3000, 13, 80)
                if end_use == "Heating" or end_use == "All":
                    output = (100, 11000, 5000, 0.070, 0.44, 45, 4.5,  1000, 7.0, 40, 100, 3000, 13, 80)
                return output + default_output

        else:
            raise PreventUpdate

    elif model == "CovHDF5":
        if "btn-nclicks-1" == ctx.triggered_id:
            output = (6, 2500, 3000, 0.045, 0.5, 45)
            k_value = 3.05
            return output + (k_value,) + default_output
        
        elif "btn-nclicks-3" == ctx.triggered_id:
            output = (4, 1500, 2000, 0.035, 0.3, 35)
            k_value = 3.05
            return output + (k_value,) + default_output
        
        else:
            raise PreventUpdate

    elif model == "SBT V1.0" or model == "SBT V2.0":
        raise PreventUpdate

@app.callback(
    [
    Output(component_id='model-params-div', component_property='style'),
    Output(component_id="Tsurf-select-div", component_property="style"),
    Output(component_id="c-select-div", component_property="style"),
    Output(component_id="rho-select-div", component_property="style"),
    Output(component_id="radius-vertical-select-div", component_property="style"),
    Output(component_id="radius-lateral-select-div", component_property="style"),
    Output(component_id="diameter-select-div", component_property="style"), # show HDF5, hide SBT
    Output(component_id="fluid-mode-div", component_property="style"), # hide HDF5 and v1, show v2
    Output(component_id='num-lat-div', component_property='style'),
    Output(component_id='lat-allocation-div', component_property='style'),
    # Output(component_id='lateral-flow-select-div', component_property='style'),
    Output(component_id='lat-flow-mul-div', component_property='style'),
    ],
   [Input(component_id="model-select", component_property="value"),
    ],
    prevent_initial_call=True
    )

def show_model_params(model):

    # print("show_model_params: ", model)

    b = {'display': 'block'}
    n = {'display': 'none'}

    if model == "HDF5":
        return n, n, n, n, n, n, b, n, n, n, n

    if model == "SBT V1.0":
        return b, b, b, b, b, b, n, b, b, n, n
    
    if model == "SBT V2.0":
        return b, b, b, b, b, b, n, b, b, n, b
    
    if model == "CovHDF5":
        # CovHDF5: hide SBT params, show standard params (but different ones)
        return n, n, n, n, n, n, b, n, n, n, n
    

@app.callback(
    [
    Output(component_id='economics-params-div', component_property='style', allow_duplicate=True),
    Output(component_id='kwt-div', component_property='style'),
    Output(component_id='kwe-div', component_property='style'),
    ],
   [
    Input(component_id="tabs", component_property="value"),
    # Input(component_id="model-select", component_property="value"),
    Input(component_id="fluid-select", component_property="value"),
    Input(component_id="end-use-select", component_property="value")
    ],
    prevent_initial_call=True,
    )

def econ_sliders_visibility(tab, fluid, end_use):

    b = {'display': 'block'}
    n = {'display': 'none'}

    econ_parms_div_style = {
        'display': 'block',
        "border": "solid 3px #c4752f",
        "borderRadius": "10px",
        "marginBottom": "5px",
        "marginRight": "5px",
        "paddingBottom": "5px",
    }

    econ_parms_div_style_2 = {
        'display': 'block',
        "borderTop": "solid 3px #c4752f",
        "borderLeft": "solid 3px #c4752f",
        "borderRight": "solid 3px #c4752f",
        "borderRadius": "10px 10px 0px 0px",
        "marginBottom": "5px",
        "marginRight": "5px",
        "paddingBottom": "5px"
    }

    if tab == "energy-time-tab" or tab == "energy-tab":
        # print("bye econ!")
        return n, n, n
    else: #
        # print("hi econ!") # how should econ look under different conditions?
        if fluid == "All" and end_use == "All":
            return econ_parms_div_style_2, b, b
        if fluid == "All" and end_use == "Heating":
            return econ_parms_div_style, b, n
        if fluid == "All" and end_use == "Electricity":
            return econ_parms_div_style_2, n, b
        elif fluid == "H2O" and end_use == "All":
            return econ_parms_div_style, b, b
        elif fluid == "H2O" and end_use == "Heating":
            return econ_parms_div_style, b, n
        elif fluid == "H2O" and end_use == "Electricity":
            return econ_parms_div_style, n, b
        elif fluid == "sCO2" and end_use == "All":
            return econ_parms_div_style_2, b, b
        elif fluid == "sCO2" and end_use == "Heating":
            return econ_parms_div_style, b, n
        elif fluid == "sCO2" and end_use == "Electricity":
            return econ_parms_div_style_2, n, b
        else:
            return b, b, b


@app.callback(
   [Output(component_id='sCO2-card', component_property='style'),
    Output(component_id='sCO2-text', component_property='style'),
    Output(component_id='check-visual-card', component_property='style'),
   ],
   [
    Input(component_id="tabs", component_property="value"),
    Input(component_id="fluid-select", component_property="value"),
    Input(component_id="end-use-select", component_property="value")
    ])

def show_hide_detailed_card(tab, fluid, end_use):

    # print("show_hide_detailed_card, tab: ", tab, end_use)

    if tab == "energy-time-tab" or tab == "energy-tab":
        # print("white 1")
        return {'border': 'solid 0px white'}, {'display': 'none'}, {'display': 'none'}

    if tab == "economics-time-tab" or tab == "about-tab" or tab == "summary-tab":

        if fluid == "H2O" or end_use == "Heating":
            # print("white 2")
            return {'border': 'solid 0px white'}, {'display': 'none'}, {'display': 'none'}
        else:
            # print("orange")
            return {'border': 'solid 3px #c4752f', 'display': 'block'}, {'display': 'block'}, {'display': 'inline-block'}
        

@app.callback(
   [Output(component_id='mdot-select-div', component_property='style'),
    Output(component_id='L2-select-div', component_property='style'),
    Output(component_id='L1-select-div', component_property='style'),
    Output(component_id='grad-select-div', component_property='style'),
    Output(component_id='diameter-select-div', component_property='style', allow_duplicate=True),
    Output(component_id='Tinj-select-div', component_property='style'), # AB DEBUG HERE
    Output(component_id='k-select-div', component_property='style'),

    Output(component_id='drillcost-div', component_property='style'),
    Output(component_id='discount-rate-div', component_property='style'),
    Output(component_id='lifetime-div', component_property='style'),
    Output(component_id='precool-div', component_property='style'),
    Output(component_id='turb-pout-div', component_property='style'),

    Output(component_id='wellbore-params-div', component_property='style'),
   ],
   [Input(component_id='param-select', component_property='value'),
    Input(component_id="tabs", component_property="value"),
    Input(component_id="fluid-select", component_property="value"),
    Input(component_id="end-use-select", component_property="value"),
    Input(component_id="model-select", component_property="value")
    ],
    prevent_initial_call=True,
    )

def show_hide_element(visibility_state, tab, fluid, end_use, model):

    # ----------------------------------------------------------------------------------------------
    # Reveals or hides sliders depending on which tab selected and which dropdowns.
    # ----------------------------------------------------------------------------------------------

    # print("show_hide_element: ", model, tab, fluid, end_use, visibility_state)
    # fluid getting changed makes this run "twice" but it's working as expected.
    
    b = {'display': 'block'}
    n = {'display': 'none'}

    if model == "HDF5":
        if tab == "about-tab":
            if fluid == "H2O" or end_use == "Heating":
                return b, b, b, b, b, b, b,\
                        b, b, b, n, n, b

            else:
                return b, b, b, b, b, b, b,\
                        b, b, b, b, b, b
        
        elif tab == "energy-time-tab":
            return b, b, b, b, b, b, b, \
                    n, n, n, n, n, b
        
        elif tab == "energy-tab":
            
            # Handle parameter visibility based on parameter name, not index
            # Handle parameter visibility based on parameter name
            if visibility_state:
                # Check parameter name directly (works for both metric and imperial)
                if "Horizontal Extent" in visibility_state:
                    return n, n, b, b, b, b, b, \
                            n, n, n, n, n, b
                elif "Vertical Extent" in visibility_state:
                    return n, b, n, b, b, b, b, \
                            n, n, n, n, n, b
                elif "Geothermal Gradient" in visibility_state:
                    return n, b, b, n, b, b, b, \
                            n, n, n, n, n, b
                elif "Borehole Diameter" in visibility_state:
                    return n, b, b, b, n, b, b, \
                            n, n, n, n, n, b
                elif "Injection Temperature" in visibility_state:
                    return n, b, b, b, b, n, b, \
                            n, n, n, n, n, b
                elif "Thermal Conductivity" in visibility_state:
                    return n, b, b, b, b, b, n, \
                            n, n, n, n, n, b
            
            # Default case - show all sliders
            return b, b, b, b, b, b, b, \
                    n, n, n, n, n, b
        
        elif tab == "economics-time-tab":
            if fluid == "H2O":
                if end_use == "All":
                    return b, b, b, b, b, b, b, \
                            b, b, b, n, n, b
                if end_use == "Heating":
                    return b, b, b, b, b, b, b, \
                            b, b, b, n, n, b
                if end_use == "Electricity":
                    return b, b, b, b, b, b, b, \
                            b, b, b, n, n, b

            else:
                if end_use == "All":
                    return b, b, b, b, b, b, b, \
                            b, b, b, b, b, b
                if end_use == "Heating":
                    return b, b, b, b, b, b, b, \
                            b, b, b, n, n, b
                if end_use == "Electricity":
                    return b, b, b, b, b, b, b, \
                            b, b, b, b, b, b

        elif tab == "summary-tab":
            if fluid == "H2O":
                return b, b, b, b, b, b, b, \
                        b, b, b, n, n, b
            else:
                return b, b, b, b, b, b, b, \
                        b, b, b, b, b, b
                        
    elif model == "SBT V1.0":
        if tab == "about-tab":
            if fluid == "H2O" or end_use == "Heating":
                return b, b, b, b, b, b, b,\
                        b, b, b, n, n, b

            else:
                return b, b, b, b, b, b, b,\
                        b, b, b, b, b, b
        
        elif tab == "energy-time-tab":
            return b, b, b, b, b, b, b, \
                    n, n, n, n, n, b
        
        elif tab == "energy-tab":
            
            # Handle parameter visibility based on parameter name, not index
            # Handle parameter visibility based on parameter name
            if visibility_state:
                # Check parameter name directly (works for both metric and imperial)
                if "Horizontal Extent" in visibility_state:
                    return n, n, b, b, b, b, b, \
                            n, n, n, n, n, b
                elif "Vertical Extent" in visibility_state:
                    return n, b, n, b, b, b, b, \
                            n, n, n, n, n, b
                elif "Geothermal Gradient" in visibility_state:
                    return n, b, b, n, b, b, b, \
                            n, n, n, n, n, b
                elif "Borehole Diameter" in visibility_state:
                    return n, b, b, b, n, b, b, \
                            n, n, n, n, n, b
                elif "Injection Temperature" in visibility_state:
                    return n, b, b, b, b, n, b, \
                            n, n, n, n, n, b
                elif "Thermal Conductivity" in visibility_state:
                    return n, b, b, b, b, b, n, \
                            n, n, n, n, n, b
            
            # Default case - show all sliders
            return b, b, b, b, b, b, b, \
                    n, n, n, n, n, b
        
        elif tab == "economics-time-tab":
            if fluid == "H2O":
                if end_use == "All":
                    return b, b, b, b, b, b, b, \
                            b, b, b, n, n, b
                if end_use == "Heating":
                    return b, b, b, b, b, b, b, \
                            b, b, b, n, n, b
                if end_use == "Electricity":
                    return b, b, b, b, b, b, b, \
                            b, b, b, n, n, b

            else:
                if end_use == "All":
                    return b, b, b, b, b, b, b, \
                            b, b, b, b, b, b
                if end_use == "Heating":
                    return b, b, b, b, b, b, b, \
                            b, b, b, n, n, b
                if end_use == "Electricity":
                    return b, b, b, b, b, b, b, \
                            b, b, b, b, b, b

        elif tab == "summary-tab":
            if fluid == "H2O":
                return b, b, b, b, b, b, b, \
                        b, b, b, n, n, b
            else:
                return b, b, b, b, b, b, b, \
                        b, b, b, b, b, b
    else:
        raise PreventUpdate


@app.callback(
   [
    Output(component_id='Tsurf-select-div', component_property='children'),
    Output(component_id='grad-container', component_property='children'),
    Output(component_id='k-container', component_property='children'),
    Output(component_id='Tinj-container', component_property='children'),
    Output(component_id='mdot-container', component_property='children'),
    Output(component_id='diameter-container', component_property='children'),
    Output(component_id='L2-container', component_property='children'),
    Output(component_id='L1-container', component_property='children'),
   ],
   [Input(component_id="model-select", component_property="value"),
    Input(component_id="case-select", component_property="value"),
    Input(component_id="quick-unit-selector", component_property="value")],
   prevent_initial_call=True
    )

def update_slider_ranges(model, case, unit_system):


    # Define styles for showing/hiding sliders
    div_block_style = {"width": "98%", "margin": "auto", "margin-bottom": "10px", "display": "block"}
    div_none_style = {"width": "98%", "margin": "auto", "margin-bottom": "10px", "display": "none"}

    # Default case - return empty containers if no model selected
    if not model:
        empty_container = html.Div(style={'display': 'none'})
        return empty_container, empty_container, empty_container, empty_container, empty_container, empty_container, empty_container, empty_container

    if model == "HDF5":
        # For HDF5, we need to update the sliders with unit-aware values
        if unit_system == "imperial":
            Tsurf_container = create_enhanced_slider(DivID="Tsurf-select-div", ID="Tsurf-select", ptitle="Surface Temperature (°F)", 
                                                    min_v=32, max_v=104, 
                                                    mark_dict={32: '32', 104: '104'}, 
                                                    start_v=77, 
                                                    div_style=div_none_style, parameter_name="Surface Temperature (˚F)")
            grad_container = create_enhanced_slider(DivID="grad-select-div", ID="grad-select", ptitle="Geothermal Gradient (°F/ft)", 
                                                    min_v=0.0082, max_v=0.0546, 
                                                    mark_dict={0.0082: '0.008', 0.0546: '0.055'}, 
                                                    start_v=0.0328, 
                                                    div_style=div_block_style, parameter_name="Geothermal Gradient (˚F/ft)")
            k_container = create_enhanced_slider(DivID="k-select-div", ID="k-select", ptitle="Rock Thermal Conductivity (Btu/ft-h-°F)", 
                                                    min_v=0.231, max_v=2.89, 
                                                    mark_dict={0.231: '0.2', 2.89: '2.9'}, 
                                                    start_v=1.73, 
                                                    div_style=div_block_style, parameter_name="Rock Thermal Conductivity (Btu/ft-h-˚F)")
            Tinj_container = create_enhanced_slider(DivID="Tinj-select-div", ID="Tinj-select", ptitle="Injection Temperature (°F)", 
                                                    min_v=86, max_v=140, 
                                                    mark_dict={86: '86', 140: '140'}, 
                                                    start_v=86, 
                                                    div_style=div_block_style, parameter_name="Injection Temperature (˚F)")
            mdot_container = create_enhanced_slider(DivID="mdot-select-div", ID="mdot-select", ptitle="Mass Flow Rate (lb/s)", 
                                                    min_v=0.066, max_v=154, 
                                                    mark_dict={0.066: '0.07', 154: '154'}, 
                                                    start_v=53, 
                                                    div_style=div_block_style, parameter_name="Mass Flow Rate (lb/s)")
            diameter_container = create_enhanced_slider(DivID="diameter-select-div", ID="diameter-select", ptitle="Borehole Diameter (ft)", 
                                                    min_v=0.708, max_v=1.458, 
                                                    mark_dict={0.708: '0.7', 1.458: '1.5'}, step_i=0.002, 
                                                    start_v=1.15, 
                                                    div_style=div_block_style, parameter_name="Borehole Diameter (ft)")
            L2_container = create_enhanced_slider(DivID="L2-select-div", ID="L2-select", ptitle="Horizontal Extent (ft)", 
                                                    min_v=708, max_v=65617, 
                                                    mark_dict={708: '708', 65617: '65k'}, 
                                                    start_v=32808, 
                                                    div_style=div_block_style, parameter_name="Horizontal Extent (ft)")
            L1_container = create_enhanced_slider(DivID="L1-select-div", ID="L1-select", ptitle="Drilling Depth (ft)", 
                                                    min_v=3281, max_v=26247, 
                                                    mark_dict={3281: '3.3k', 26247: '26k'}, 
                                                    start_v=11483, 
                                                    div_style=div_block_style, parameter_name="Drilling Depth (ft)")
        else:
            # Metric HDF5 ranges
            Tsurf_container = create_enhanced_slider(DivID="Tsurf-select-div", ID="Tsurf-select", ptitle="Surface Temperature (°C)", 
                                                    min_v=0, max_v=40, 
                                                    mark_dict={0: '0', 40: '40'}, 
                                                    start_v=25, 
                                                    div_style=div_none_style, parameter_name="Surface Temperature (˚C)")
            grad_container = create_enhanced_slider(DivID="grad-select-div", ID="grad-select", ptitle="Geothermal Gradient (K/m)", 
                                                    min_v=0.015, max_v=0.1, 
                                                    mark_dict={0.015: '0.015', 0.1: '0.1'}, 
                                                    start_v=0.06, 
                                                    div_style=div_block_style, parameter_name="Geothermal Gradient (K/m)")
            k_container = create_enhanced_slider(DivID="k-select-div", ID="k-select", ptitle="Rock Thermal Conductivity (W/m-K)", 
                                                    min_v=0.4, max_v=5.0, 
                                                    mark_dict={0.4: '0.4', 5.0: '5.0'}, 
                                                    start_v=3.0, 
                                                    div_style=div_block_style, parameter_name="Rock Thermal Conductivity (W/m-K)")
            Tinj_container = create_enhanced_slider(DivID="Tinj-select-div", ID="Tinj-select", ptitle="Injection Temperature (°C)", 
                                                    min_v=30.0, max_v=60.0, 
                                                    mark_dict={30: '30', 60: '60'}, 
                                                    start_v=30.0, 
                                                    div_style=div_block_style, parameter_name="Injection Temperature (˚C)")
            mdot_container = create_enhanced_slider(DivID="mdot-select-div", ID="mdot-select", ptitle="Mass Flow Rate (kg/s)", 
                                                    min_v=0.03, max_v=70, 
                                                    mark_dict={0.03: '0.03', 70: '70'}, 
                                                    start_v=24, 
                                                    div_style=div_block_style, parameter_name="Mass Flow Rate (kg/s)")
            diameter_container = create_enhanced_slider(DivID="diameter-select-div", ID="diameter-select", ptitle="Borehole Diameter (m)", 
                                                    min_v=0.2159, max_v=0.4445, 
                                                    mark_dict={0.2159: '0.22', 0.4445: '0.44'}, step_i=0.002, 
                                                    start_v=0.35, 
                                                    div_style=div_block_style, parameter_name="Borehole Diameter (m)")
            L2_container = create_enhanced_slider(DivID="L2-select-div", ID="L2-select", ptitle="Horizontal Extent (m)", 
                                                    min_v=216, max_v=20000, 
                                                    mark_dict={216: '216', 20000: '20k'}, 
                                                    start_v=10000, 
                                                    div_style=div_block_style, parameter_name="Horizontal Extent (m)")
            L1_container = create_enhanced_slider(DivID="L1-select-div", ID="L1-select", ptitle="Drilling Depth (m)", 
                                                    min_v=1000, max_v=8000, 
                                                    mark_dict={1000: '1k', 8000: '8k'}, 
                                                    start_v=3500, 
                                                    div_style=div_block_style, parameter_name="Drilling Depth (m)")
        
        return Tsurf_container, grad_container, k_container, Tinj_container, mdot_container, diameter_container, L2_container, L1_container

    elif model == "SBT V1.0" or model == "SBT V2.0":
        # Hardcoded values for SBT models
        if unit_system == "imperial":
            Tsurf_container = create_enhanced_slider(DivID="Tsurf-select-div", ID="Tsurf-select", ptitle="Surface Temperature (°F)", 
                                                    min_v=32, max_v=104, 
                                                    mark_dict={32: '32', 104: '104'}, 
                                                    start_v=77, 
                                                    div_style=div_block_style, parameter_name="Surface Temperature (˚F)")
            grad_container = create_enhanced_slider(DivID="grad-select-div", ID="grad-select", ptitle="Geothermal Gradient (°F/ft)", 
                                                    min_v=0.0082, max_v=0.11, 
                                                    mark_dict={0.0082: '0.008', 0.11: '0.11'}, 
                                                    start_v=0.033, 
                                                    div_style=div_block_style, parameter_name="Geothermal Gradient (˚F/ft)")
            k_container = create_enhanced_slider(DivID="k-select-div", ID="k-select", ptitle="Rock Thermal Conductivity (Btu/ft-h-°F)", 
                                                    min_v=0.231, max_v=2.89, 
                                                    mark_dict={0.231: '0.2', 2.89: '2.9'}, 
                                                    start_v=1.73, 
                                                    div_style=div_block_style, parameter_name="Rock Thermal Conductivity (Btu/ft-h-˚F)")
            Tinj_container = create_enhanced_slider(DivID="Tinj-select-div", ID="Tinj-select", ptitle="Injection Temperature (°F)", 
                                                    min_v=86, max_v=212, 
                                                    mark_dict={86: '86', 212: '212'}, 
                                                    start_v=86, 
                                                    div_style=div_block_style, parameter_name="Injection Temperature (˚F)")
        else:
            # Metric SBT ranges
            Tsurf_container = create_enhanced_slider(DivID="Tsurf-select-div", ID="Tsurf-select", ptitle="Surface Temperature (°C)", 
                                                    min_v=0, max_v=40, 
                                                    mark_dict={0: '0', 40: '40'}, 
                                                    start_v=25, 
                                                    div_style=div_block_style, parameter_name="Surface Temperature (˚C)")
            grad_container = create_enhanced_slider(DivID="grad-select-div", ID="grad-select", ptitle="Geothermal Gradient (K/m)", 
                                                    min_v=0.015, max_v=0.200, 
                                                    mark_dict={0.015: '0.015', 0.200: '0.2'}, 
                                                    start_v=0.06, 
                                                    div_style=div_block_style, parameter_name="Geothermal Gradient (K/m)")
            k_container = create_enhanced_slider(DivID="k-select-div", ID="k-select", ptitle="Rock Thermal Conductivity (W/m-K)", 
                                                    min_v=0.4, max_v=5.0, 
                                                    mark_dict={0.4: '0.4', 5.0: '5.0'}, 
                                                    start_v=3.0, 
                                                    div_style=div_block_style, parameter_name="Rock Thermal Conductivity (W/m-K)")
            Tinj_container = create_enhanced_slider(DivID="Tinj-select-div", ID="Tinj-select", ptitle="Injection Temperature (°C)", 
                                                    min_v=30.0, max_v=100.0, 
                                                    mark_dict={30: '30', 100: '100'}, 
                                                    start_v=30.0, 
                                                    div_style=div_block_style, parameter_name="Injection Temperature (˚C)")
        
        # Tube geometry sliders for SBT models
        if unit_system == "imperial":
            mdot_container = create_enhanced_slider(DivID="mdot-select-div", ID="mdot-select", ptitle="Mass Flow Rate (lb/s)", 
                                                   min_v=22, max_v=441, mark_dict={22: '22', 441: '441'}, start_v=53, 
                                                   div_style=div_block_style, parameter_name="Mass Flow Rate (lb/s)")
            diameter_container = create_enhanced_slider(DivID="diameter-select-div", ID="diameter-select", ptitle="Borehole Diameter (ft)", 
                                                       min_v=0.708, max_v=1.458, mark_dict={0.708: '0.7', 1.458: '1.5'}, start_v=1.15, 
                                                       div_style=div_block_style, parameter_name="Borehole Diameter (ft)")
            L2_container = create_enhanced_slider(DivID="L2-select-div", ID="L2-select", ptitle="Horizontal Extent (ft)", 
                                                 min_v=1640, max_v=65617, mark_dict={1640: '1.6k', 65617: '65k'}, start_v=32808, 
                                                 div_style=div_block_style, parameter_name="Horizontal Extent (ft)")
            L1_container = create_enhanced_slider(DivID="L1-select-div", ID="L1-select", ptitle="Drilling Depth (ft)", 
                                                 min_v=3281, max_v=26247, mark_dict={3281: '3.3k', 26247: '26k'}, start_v=11483, 
                                                 div_style=div_block_style, parameter_name="Drilling Depth (ft)")
        else:
            mdot_container = create_enhanced_slider(DivID="mdot-select-div", ID="mdot-select", ptitle="Mass Flow Rate (kg/s)", 
                                                   min_v=10, max_v=200, mark_dict={10: '10', 200: '200'}, start_v=100, 
                                                   div_style=div_block_style, parameter_name="Mass Flow Rate (kg/s)")
            diameter_container = create_enhanced_slider(DivID="diameter-select-div", ID="diameter-select", ptitle="Borehole Diameter (m)", 
                                                       min_v=0.2159, max_v=0.4445, mark_dict={0.2159: '0.22', 0.4445: '0.44'}, start_v=0.3, 
                                                       div_style=div_block_style, parameter_name="Borehole Diameter (m)")
            L2_container = create_enhanced_slider(DivID="L2-select-div", ID="L2-select", ptitle="Horizontal Extent (m)", 
                                                 min_v=500, max_v=5000, mark_dict={500: '500', 5000: '5k'}, start_v=2500, 
                                                 div_style=div_block_style, parameter_name="Horizontal Extent (m)")
            L1_container = create_enhanced_slider(DivID="L1-select-div", ID="L1-select", ptitle="Drilling Depth (m)", 
                                                 min_v=1000, max_v=8000, mark_dict={1000: '1k', 8000: '8k'}, start_v=5000, 
                                                 div_style=div_block_style, parameter_name="Drilling Depth (m)")
                       
        return Tsurf_container, grad_container, k_container, Tinj_container, mdot_container, diameter_container, L2_container, L1_container

    elif model == "CovHDF5":
        # For CovHDF5, create sliders with permeability parameter instead of borehole diameter and thermal conductivity
        if unit_system == "imperial":
            # Imperial CovHDF5 ranges
            mdot_container = create_enhanced_slider(DivID="mdot-select-div", ID="mdot-select", ptitle="Mass Flow Rate (lb/s)",
                                                    min_v=4.4, max_v=22.0,
                                                    mark_dict={4.4: '4.4', 21.6: '22'},
                                                    start_v=13.2,
                                                    div_style=div_block_style, parameter_name="CovHDF5 Mass Flow Rate (lb/s)")
            L2_container = create_enhanced_slider(DivID="L2-select-div", ID="L2-select", ptitle="Horizontal Extent (ft)",
                                                    min_v=3281, max_v=16404,
                                                    mark_dict={3281: '3.3k', 8202: '8.2k', 16404: '16k'},
                                                    start_v=8202,
                                                    div_style=div_block_style, parameter_name="CovHDF5 Horizontal Extent (ft)")
            L1_container = create_enhanced_slider(DivID="L1-select-div", ID="L1-select", ptitle="Drilling Depth (ft)",
                                                    min_v=3281, max_v=16404,
                                                    mark_dict={3281: '3.3k', 16404: '16k'},
                                                    start_v=9843,
                                                    div_style=div_block_style, parameter_name="CovHDF5 Drilling Depth (ft)")
            grad_container = create_enhanced_slider(DivID="grad-select-div", ID="grad-select", ptitle="Geothermal Gradient (°F/ft)",
                                                    min_v=0.0164, max_v=0.0328,
                                                    mark_dict={0.0164: '0.016', 0.0328: '0.033'},
                                                    start_v=0.0246,
                                                    div_style=div_block_style, parameter_name="CovHDF5 Geothermal Gradient (°F/ft)")
            perm_container = create_enhanced_slider(DivID="diameter-select-div", ID="diameter-select", ptitle="Permeability (HWR)",
                                                    min_v=0.1, max_v=1.0,
                                                    mark_dict={0.1: '0.1', 0.5: '0.5', 0.999: '1.0'},
                                                    start_v=0.5,
                                                    div_style=div_block_style, parameter_name="Permeability (HWR)")
            Tinj_container = create_enhanced_slider(DivID="Tinj-select-div", ID="Tinj-select", ptitle="Injection Temperature (°F)",
                                                    min_v=86, max_v=140,
                                                    mark_dict={86: '86', 140: '140'},
                                                    start_v=113,
                                                    div_style=div_block_style, parameter_name="CovHDF5 Injection Temperature (°F)")
        else:
            # Metric CovHDF5 ranges
            mdot_container = create_enhanced_slider(DivID="mdot-select-div", ID="mdot-select", ptitle="Mass Flow Rate (kg/s)", 
                                                    min_v=2, max_v=10, 
                                                    mark_dict={2: '2', 10: '10'}, 
                                                    start_v=6, 
                                                    div_style=div_block_style, parameter_name="CovHDF5 Mass Flow Rate (kg/s)")
            L2_container = create_enhanced_slider(DivID="L2-select-div", ID="L2-select", ptitle="Horizontal Extent (m)", 
                                                    min_v=1000, max_v=5000, 
                                                    mark_dict={1000: '1k', 2500: '2.5k', 5000: '5k'}, 
                                                    start_v=2500, 
                                                    div_style=div_block_style, parameter_name="CovHDF5 Horizontal Extent (m)")
            L1_container = create_enhanced_slider(DivID="L1-select-div", ID="L1-select", ptitle="Drilling Depth (m)", 
                                                    min_v=1000, max_v=5000, 
                                                    mark_dict={1000: '1k', 5000: '5k'}, 
                                                    start_v=3000, 
                                                    div_style=div_block_style, parameter_name="CovHDF5 Drilling Depth (m)")
            grad_container = create_enhanced_slider(DivID="grad-select-div", ID="grad-select", ptitle="Geothermal Gradient (K/m)", 
                                                    min_v=0.03, max_v=0.06, 
                                                    mark_dict={0.03: '0.03', 0.06: '0.06'}, 
                                                    start_v=0.045, 
                                                    div_style=div_block_style, parameter_name="CovHDF5 Geothermal Gradient (K/m)")
            perm_container = create_enhanced_slider(DivID="diameter-select-div", ID="diameter-select", ptitle="Permeability (HWR)", 
                                                    min_v=0.1, max_v=1.0, 
                                                    mark_dict={0.1: '0.1', 0.5: '0.5', 0.999: '1.0'}, 
                                                    start_v=0.5, 
                                                    div_style=div_block_style, parameter_name="Permeability (HWR)")
            Tinj_container = create_enhanced_slider(DivID="Tinj-select-div", ID="Tinj-select", ptitle="Injection Temperature (°C)", 
                                                    min_v=30, max_v=60, 
                                                    mark_dict={30: '30', 60: '60'}, 
                                                    start_v=45, 
                                                    div_style=div_block_style, parameter_name="CovHDF5 Injection Temperature (°C)")
        
        # For CovHDF5, we don't need Tsurf, k (thermal conductivity), so we'll create hidden sliders
        Tsurf_container = create_enhanced_slider(DivID="Tsurf-select-div", ID="Tsurf-select", ptitle="Surface Temperature (°C)", 
                                                min_v=25, max_v=25, 
                                                mark_dict={25: '25'}, 
                                                start_v=25, 
                                                div_style={'display': 'none'}, parameter_name="Surface Temperature (°C)")
        # Create a hidden k-select slider with default value for CovHDF5
        k_container = create_enhanced_slider(DivID="k-select-div", ID="k-select", ptitle="Rock Thermal Conductivity (W/m-K)", 
                                            min_v=3.0, max_v=3.1, 
                                            mark_dict={3.0: '3.0', 3.1: '3.1'}, 
                                            start_v=3.05, 
                                            div_style={'display': 'none'}, parameter_name="Rock Thermal Conductivity (W/m-K)")
        
        return Tsurf_container, grad_container, k_container, Tinj_container, mdot_container, perm_container, L2_container, L1_container

    else:
        # Default case - return empty containers
        empty_container = html.Div(style={'display': 'none'})
        return empty_container, empty_container, empty_container, empty_container, empty_container, empty_container, empty_container, empty_container

@app.callback(
   [
     Output(component_id='hyperparam1-container', component_property='children'), 
     Output(component_id='hyperparam3-container', component_property='children'),
     Output(component_id='hyperparam5-container', component_property='children')
   ],
   [Input(component_id="model-select", component_property="value"),
    ],
   prevent_initial_call=True
    )
def update_sliders_hyperparms(model):

    if model == "SBT V1.0":

        hyperparam1 = create_enhanced_dropdown(DivID="mass-flow-mode-div", ID="mass-mode-select", ptitle="Mass Flow Rate Mode", 
                                                                                                options=["Constant", "Variable"], disabled=True, div_style=div_block_style, parameter_name=None)
        hyperparam3 = create_enhanced_dropdown(DivID="temp-flow-mode-div", ID="temp-mode-select", ptitle="Injection Temperature Mode", 
                                                                                        options=["Constant", "Variable"], disabled=True, div_style=div_block_style, parameter_name=None)
        hyperparam5 = create_enhanced_dropdown(DivID="fluid-mode-div", ID="fluid-mode-select", ptitle="Fluid Properties Mode", 
                                                                                                options=["Constant", "Variable"], disabled=True, div_style=div_none_style, parameter_name="Fluid Properties Mode")

        return hyperparam1, hyperparam3, hyperparam5
    
    elif model == "SBT V2.0":

        inlet_pressure_dict = {5: '5', 20: '20'}
        pipe_roughness_dict = {1e-6: '1e-6', 3e-6: '3e-6'}

        hyperparam1 = create_enhanced_slider(DivID="mass-flow-mode-div", ID="mass-mode-select", ptitle="Inlet Pressure (MPa)", min_v=5, max_v=20,
                                                            mark_dict=inlet_pressure_dict, step_i=0.1, start_v=start_vals_sbt["inletpressure"], div_style=div_block_style, parameter_name=None)
        hyperparam3 = create_enhanced_slider(DivID="temp-flow-mode-div", ID="temp-mode-select", ptitle="Pipe Roughness (m)", min_v=1e-6, max_v=3e-6,
                                                            mark_dict=pipe_roughness_dict, step_i=0.000001, start_v=start_vals_sbt["piperoughness"], div_style=div_block_style, parameter_name=None)
        hyperparam5 = create_enhanced_dropdown(DivID="fluid-mode-div", ID="fluid-mode-select", ptitle="Fluid Properties Mode", 
                                                                                                options=["Variable", "Constant"], disabled=True, div_style=div_block_style, parameter_name="Fluid Properties Mode")

        return hyperparam1, hyperparam3, hyperparam5

    else:
        raise PreventUpdate


# -----------------------------------------------------------------------------
# Define dash app plotting callbacks.
# -----------------------------------------------------------------------------

@app.callback(
    [Output(component_id="geothermal_time_plots", component_property="figure"),
     Output(component_id='thermal-memory', component_property='data'),
     Output(component_id='thermal-results-mass', component_property='data'),
     Output(component_id='thermal-results-time', component_property='data'),
     Output(component_id='thermal-results-errors', component_property='data'),
     Output(component_id="TandP-data", component_property="data"),
     ],
    [Input(component_id="interpolation-select", component_property="value"),
     Input(component_id="fluid-select", component_property="value"),
     Input(component_id="case-select", component_property="value"),
     Input(component_id="mdot-select", component_property="value"),
     Input(component_id="L2-select", component_property="value"),
     Input(component_id="L1-select", component_property="value"),
     Input(component_id="grad-select", component_property="value"),
     Input(component_id="diameter-select", component_property="value"),
     Input(component_id="Tinj-select", component_property="value"),
     Input(component_id="k-select", component_property="value"),
     Input(component_id='radio-graphic-control3', component_property='value'),
     Input(component_id='model-select', component_property='value'),
     
     # more variables 
     Input(component_id='Tsurf-select', component_property='value'),
     Input(component_id='c-select', component_property='value'),
     Input(component_id='rho-select', component_property='value'),

     Input(component_id='radius-vertical-select', component_property='value'), # diameter1
     Input(component_id='radius-lateral-select', component_property='value'), # diameter2
     Input(component_id='n-laterals-select', component_property='value'), # PipeParam3
     Input(component_id='lateral-flow-select', component_property='value'), # PipeParam4
     Input(component_id='lateral-multiplier-select', component_property='value'), # PipeParam5

     Input(component_id='mesh-select', component_property='value'),
     Input(component_id='accuracy-select', component_property='value'),
     Input(component_id='mass-mode-select', component_property='value'),
     Input(component_id='temp-mode-select', component_property='value'),
     Input(component_id='fluid-mode-select', component_property='value'),
     Input(component_id="quick-unit-selector", component_property="value"),

    ],
)

def update_subsurface_results_plots(interp_time, fluid, case, mdot, L2, L1, grad, D, Tinj, k_m, scale, model,
                        Tsurf, c_m, rho_m, 
                        # radius_vertical, radius_lateral, n_laterals, lateral_flow, lateral_multiplier,
                        Diameter1, Diameter2, PipeParam3, PipeParam4, PipeParam5,
                        mesh, accuracy, 
                        # mass_mode, temp_mode
                        HyperParam3, HyperParam4, HyperParam5, units
                        ):

    # -----------------------------------------------------------------------------
    # Creates and displays Plotly subplots of the subsurface results.
    # -----------------------------------------------------------------------------



    # Handle missing SBT parameters gracefully
    # Only set defaults if we're actually using SBT models
    if model in ["SBT V1.0", "SBT V2.0"]:
        # For SBT models, ensure we have valid parameters
        if not Diameter1 or Diameter1 is None:
            Diameter1 = 0.10795  # Default value
        if not Diameter2 or Diameter2 is None:
            Diameter2 = 0.10795  # Default value
        if not PipeParam3 or PipeParam3 is None:
            PipeParam3 = 1  # Default value
        if not PipeParam4 or PipeParam4 is None:
            PipeParam4 = 0.5  # Default value
        if not PipeParam5 or PipeParam5 is None:
            PipeParam5 = 0.5  # Default value
    else:
        # For HDF5, these parameters aren't used, so set to None
        Diameter1 = None
        Diameter2 = None
        PipeParam3 = None
        PipeParam4 = None
        PipeParam5 = None


    # Convert imperial values back to metric for backend calculations
    # Simple conversion function to convert imperial values back to metric
    def convert_to_metric(value, from_unit, to_unit, converter_func):
        """Convert a value from imperial to metric units"""
        if from_unit != to_unit:
            return converter_func(value, from_unit, to_unit)
        return value
    
    # Get current unit preferences
    from unit_conversions import unit_converter
    
    # Use the units parameter directly instead of inferring from preferences
    unit_system = units
    
    # Convert imperial slider values to metric for calculations (same approach as contour plots)
    if units == "imperial":
        # Convert imperial slider inputs to metric for internal calculations
        L2_m = L2 * 0.3048  # ft to m
        L1_m = L1 * 0.3048  # ft to m
        D_m = D * 0.3048  # ft to m
        mdot_kg_s = mdot * 0.45359237  # lb/s to kg/s
        Tinj_c = (Tinj - 32.0) * (5.0/9.0)  # °F to °C
        grad_k_m = grad * (5.0/9.0) / 0.3048  # °F/ft to K/m
        k_w_m_k = k_m * 1.730735  # BTU/(hr·ft·°F) to W/m·K
        c_j_kg_k = c_m * 4186.8  # BTU/(lb·°F) to J/(kg·K)
        rho_kg_m3 = rho_m * 16.0185  # lb/ft³ to kg/m³
        Tsurf_c = (Tsurf - 32.0) * (5.0/9.0)  # °F to °C
    else:
        # Metric values - pass through as-is
        L2_m = L2
        L1_m = L1
        D_m = D
        mdot_kg_s = mdot
        Tinj_c = Tinj
        grad_k_m = grad
        k_w_m_k = k_m
        c_j_kg_k = c_m
        rho_kg_m3 = rho_m
        Tsurf_c = Tsurf
    
    try:
        # print('subsurface')
        # if HDF5:
        # start = time.time()
        subplots, forty_yr_TPmeans_dict, df_mass_flow_rate, df_time, err_subres_dict, TandP_dict = generate_subsurface_lineplots(
            interp_time, fluid, case, mdot_kg_s, L2_m, L1_m, grad_k_m, D_m, Tinj_c, k_w_m_k, scale, model,
            Tsurf_c, c_j_kg_k, rho_kg_m3, 
            # radius_vertical, radius_lateral, 
            Diameter1, Diameter2, PipeParam3, PipeParam4, PipeParam5,
            mesh, accuracy, HyperParam3, HyperParam4, HyperParam5, unit_system
        )
        # if SBT:
        # end = time.time()
        # print("run generate_subsurface_lineplots:", end - start)
    except Exception as e:
        import traceback
        traceback.print_exc()
        # Return empty/default values on error
        from plotly.subplots import make_subplots
        subplots = make_subplots(rows=2, cols=3)
        forty_yr_TPmeans_dict = {}
        df_mass_flow_rate = {}
        df_time = {}
        err_subres_dict = {'Error': str(e)}
        TandP_dict = {}

    return subplots, forty_yr_TPmeans_dict, df_mass_flow_rate, df_time, err_subres_dict, TandP_dict



@app.callback(
    [Output(component_id="geothermal_plots", component_property="figure"),
     Output(component_id="thermal-contours-errors", component_property="data"),
    ],
    [Input(component_id="interpolation-select", component_property="value"),
     Input(component_id="fluid-select", component_property="value"),
     Input(component_id="case-select", component_property="value"),
     Input(component_id="param-select", component_property="value"),
     Input(component_id="mdot-select", component_property="value"),
     Input(component_id="L2-select", component_property="value"),
     Input(component_id="L1-select", component_property="value"),
     Input(component_id="grad-select", component_property="value"),
     Input(component_id="diameter-select", component_property="value"),
     Input(component_id="Tinj-select", component_property="value"),
     Input(component_id="k-select", component_property="value"),
     Input(component_id="quick-unit-selector", component_property="value"),
    ],
)

def update_subsurface_contours_plots(interp_time, fluid, case, param, mdot, L2, L1, grad, D, Tinj, k_m, units):

    # -----------------------------------------------------------------------------
    # Creates and displays Plotly subplots of the subsurface contours.
    # -----------------------------------------------------------------------------

    # Convert imperial slider values to metric for calculations
    if units == "imperial":
        # Convert imperial slider inputs to metric for internal calculations
        mdot_metric = mdot * 0.45359237  # lb/s to kg/s
        L2_metric = L2 * 0.3048  # ft to m
        L1_metric = L1 * 0.3048  # ft to m
        grad_metric = grad * (5.0/9.0) / 0.3048  # °F/ft to K/m
        D_metric = D * 0.3048  # ft to m
        Tinj_metric = (Tinj - 32.0) * (5.0/9.0)  # °F to °C
        k_metric = k_m * 1.730735  # BTU/(hr·ft·°F) to W/m·K
    else:
        # Metric values - pass through as-is
        mdot_metric = mdot
        L2_metric = L2
        L1_metric = L1
        grad_metric = grad
        D_metric = D
        Tinj_metric = Tinj
        k_metric = k_m

    # print('contours')
    subplots, err_subcontour_dict = generate_subsurface_contours(
        interp_time, fluid, case, param, mdot_metric, L2_metric, L1_metric, grad_metric, D_metric, Tinj_metric, k_metric, units
    )

    return subplots, err_subcontour_dict


@app.callback(
    Output("param-select", "options"),
    Output("param-select", "value"),
    Input("quick-unit-selector", "value"),
    prevent_initial_call=True,
)
def update_param_dropdown_options(units):
    """Update parameter dropdown options based on unit system"""
    from dropdowns import get_param_list
    
    if units == "imperial":
        options = [{"label": i, "value": i} for i in get_param_list("imperial")]
        value = "Horizontal Extent (ft)"
    else:
        options = [{"label": i, "value": i} for i in get_param_list("metric")]
        value = "Horizontal Extent (m)"
    
    return options, value


@app.callback(
    [Output(component_id="econ_plots", component_property="figure"),
     Output(component_id='econ-memory', component_property='data'),
     Output(component_id='econ-results', component_property='data'),
     Output(component_id='econ-errors', component_property='data'),
     # Output(component_id='ts_plot', component_property='figure')
     ],
    [
     Input(component_id="TandP-data", component_property="data"),

     Input(component_id="interpolation-select", component_property="value"),
     Input(component_id="fluid-select", component_property="value"),
     Input(component_id="case-select", component_property="value"),
     Input(component_id="end-use-select", component_property="value"),

     # Economic parameters only - these should trigger economic updates
     Input(component_id="drillcost-select", component_property="value"),
     Input(component_id="discount-rate-select", component_property="value"),
     Input(component_id="lifetime-select", component_property="value"),
     Input(component_id="kwt-select", component_property="value"),
     Input(component_id="kwe-select", component_property="value"),
     Input(component_id="precool-select", component_property="value"),
     Input(component_id="turb-pout-select", component_property="value"),
     Input(component_id='radio-graphic-control4', component_property='value'),
     Input(component_id="checklist", component_property="value"),
     Input(component_id='model-select', component_property='value'),
     Input(component_id="quick-unit-selector", component_property="value"),  # Add this input
     
     # Thermal parameters as State - these won't trigger the callback but are available for calculations
     State(component_id="mdot-select", component_property="value"),
     State(component_id="L2-select", component_property="value"),
     State(component_id="L1-select", component_property="value"),
     State(component_id="grad-select", component_property="value"),
     State(component_id="diameter-select", component_property="value"),
     State(component_id="Tinj-select", component_property="value"),
     State(component_id="k-select", component_property="value"),
    ],
)

def update_econ_plots(TandP_dict,
                     interp_time, fluid, case, end_use,
                     Drilling_cost_per_m, Discount_rate, Lifetime, 
                     Direct_use_heat_cost_per_kWth, Power_plant_cost_per_kWe, Pre_Cooling_Delta_T, Turbine_outlet_pressure,
                     scale, checklist, model, units,  # Add units parameter
                     mdot, L2, L1, grad, D, Tinj, k_m):

    # -----------------------------------------------------------------------------
    # Creates and displays Plotly subplots of the economic results.
    # -----------------------------------------------------------------------------

    # Check for None values that could cause errors during unit conversion
    required_params = [mdot, L2, L1, grad, D, Tinj, k_m, Drilling_cost_per_m, Discount_rate, Lifetime, 
                      Direct_use_heat_cost_per_kWth, Power_plant_cost_per_kWe, Pre_Cooling_Delta_T, Turbine_outlet_pressure]
    if any(param is None for param in required_params):
        from plotly.graph_objects import Figure
        empty_fig = Figure()
        return empty_fig, {}, {}, {}

    # Convert imperial values to metric for calculations (same pattern as subsurface plots)
    if units == "imperial":
        # Convert thermal parameters from imperial to metric
        mdot_metric = mdot * 0.45359237  # lb/s to kg/s
        L2_metric = L2 * 0.3048  # ft to m
        L1_metric = L1 * 0.3048  # ft to m
        grad_metric = grad * (5.0/9.0) / 0.3048  # °F/ft to K/m
        D_metric = D * 0.3048  # ft to m
        Tinj_metric = (Tinj - 32.0) * (5.0/9.0)  # °F to °C
        k_metric = k_m * 1.730735  # BTU/(hr·ft·°F) to W/m·K
        
        # Convert economic parameters from imperial to metric
        Drilling_cost_per_m_metric = Drilling_cost_per_m / 0.3048  # $/ft to $/m
        Pre_Cooling_Delta_T_metric = Pre_Cooling_Delta_T * (5.0/9.0)  # °F to °C
    else:
        # Metric values - pass through as-is
        mdot_metric = mdot
        L2_metric = L2
        L1_metric = L1
        grad_metric = grad
        D_metric = D
        Tinj_metric = Tinj
        k_metric = k_m
        Drilling_cost_per_m_metric = Drilling_cost_per_m
        Pre_Cooling_Delta_T_metric = Pre_Cooling_Delta_T

    # print('economics')

    if checklist == [' ']:
        is_plot_ts = True
    else:
        is_plot_ts = False

    try:
        economics_fig, econ_data_dict, econ_values_dict, err_econ_dict = generate_econ_lineplots(
            TandP_dict,
            interp_time, case, end_use, fluid, 
            mdot_metric, L2_metric, L1_metric, grad_metric, D_metric, Tinj_metric, k_metric,  # Use metric values
            Drilling_cost_per_m_metric, Discount_rate, Lifetime,  # Use metric drilling cost
            Direct_use_heat_cost_per_kWth, Power_plant_cost_per_kWe, Pre_Cooling_Delta_T_metric, Turbine_outlet_pressure,  # Use metric pre-cooling
            scale, 
            properties_H2O_pathname, 
            properties_CO2v2_pathname, 
            additional_properties_CO2v2_pathname,
            tmatrix_pathname,
            model,
            is_plot_ts,
            units
        )
        return economics_fig, econ_data_dict, econ_values_dict, err_econ_dict
    except Exception as e:
        from plotly.graph_objects import Figure
        empty_fig = Figure()
        error_dict = {'Error': str(e)}
        return empty_fig, {}, {}, error_dict


@app.callback(
   Output(component_id='ts-text', component_property='style'), 
   [
    Input(component_id="fluid-select", component_property="value"),
    Input(component_id="end-use-select", component_property="value"),
    Input(component_id="checklist", component_property="value"),
    ]
    )

def update_plot_title(fluid, end_use, checklist):

    if checklist == [' ']:
        is_title_show = True
    else:
        is_title_show = False

    if not is_title_show:

        return {'display': 'none'}

    if fluid == "H2O" or end_use == "Heating": 

        return {'display': 'none'}

    if end_use == "Electricity":

        return {'display': 'block', 'marginTop':'-620px'}

@app.callback(
    [Output(component_id="table", component_property="figure"),
     Output(component_id='summary-memory', component_property='data')
     ],
    [Input(component_id="interpolation-select", component_property="value"),
     Input(component_id="fluid-select", component_property="value"),
     Input(component_id="case-select", component_property="value"),
     # Input("end-use-select", "value"),

     Input(component_id="mdot-select", component_property="value"),
     Input(component_id="L2-select", component_property="value"),
     Input(component_id="L1-select", component_property="value"),
     Input(component_id="grad-select", component_property="value"),
     Input(component_id="diameter-select", component_property="value"),
     Input(component_id="Tinj-select", component_property="value"),
     Input(component_id="k-select", component_property="value"),

     Input(component_id="drillcost-select", component_property="value"),
     Input(component_id="discount-rate-select", component_property="value"),
     Input(component_id="lifetime-select", component_property="value"),
     Input(component_id="kwt-select", component_property="value"),
     Input(component_id="kwe-select", component_property="value"),
     Input(component_id="precool-select", component_property="value"),
     Input(component_id="turb-pout-select", component_property="value"),
     
     # SBT model parameters
     Input(component_id='Tsurf-select', component_property='value'),
     Input(component_id='c-select', component_property='value'),
     Input(component_id='rho-select', component_property='value'),
     Input(component_id='radius-vertical-select', component_property='value'),
     Input(component_id='radius-lateral-select', component_property='value'),
     Input(component_id='n-laterals-select', component_property='value'),
     Input(component_id='lateral-flow-select', component_property='value'),
     Input(component_id='lateral-multiplier-select', component_property='value'),
     Input(component_id='mesh-select', component_property='value'),
     Input(component_id='accuracy-select', component_property='value'),
     Input(component_id='mass-mode-select', component_property='value'),
     Input(component_id='temp-mode-select', component_property='value'),
     Input(component_id='fluid-mode-select', component_property='value'),
     
     Input(component_id='econ-memory', component_property='data'),
     Input(component_id='econ-results', component_property='data'),
     Input(component_id='thermal-memory', component_property='data'),
     Input(component_id='model-select', component_property='value'),
     Input(component_id='TandP-data', component_property='data'),
     Input(component_id="quick-unit-selector", component_property="value"),
    ],
)

def update_table(interp_time, fluid, case, mdot, L2, L1, grad, D, Tinj, k,
                 Drilling_cost_per_m, Discount_rate, Lifetime, 
                 Direct_use_heat_cost_per_kWth, Power_plant_cost_per_kWe, Pre_Cooling_Delta_T, Turbine_outlet_pressure,
                 Tsurf, c_m, rho_m, Diameter1, Diameter2, PipeParam3, PipeParam4, PipeParam5,
                 mesh, accuracy, HyperParam3, HyperParam4, HyperParam5,
                 econ_dict, econ_results, thermal_dict, model, tandp_data, units):

    thermal_params = [mdot, L2, L1, grad, D, Tinj, k]
    # Check for None values that could cause errors during unit conversion
    if any(param is None for param in thermal_params):
        from plotly.graph_objects import Table
        empty_table = Table()
        return empty_table, {}

    econ_params = [Drilling_cost_per_m, Discount_rate, Lifetime, Direct_use_heat_cost_per_kWth, Power_plant_cost_per_kWe, Pre_Cooling_Delta_T, Turbine_outlet_pressure]
    # Check for None economic parameters too
    if any(param is None for param in econ_params):
        from plotly.graph_objects import Table
        empty_table = Table()
        return empty_table, {}

    # Convert imperial slider values to metric for calculations (same approach as other callbacks)
    if units == "imperial":
        # Convert imperial slider inputs to metric for internal calculations
        L2_m = L2 * 0.3048  # ft to m
        L1_m = L1 * 0.3048  # ft to m
        D_m = D * 0.3048  # ft to m
        mdot_kg_s = mdot * 0.45359237  # lb/s to kg/s
        Tinj_c = (Tinj - 32.0) * (5.0/9.0)  # °F to °C
        grad_k_m = grad * (5.0/9.0) / 0.3048  # °F/ft to K/m
        k_w_m_k = k * 1.730735  # BTU/(hr·ft·°F) to W/m·K
        c_j_kg_k = c_m * 4186.8 if c_m is not None else 4186.8  # BTU/(lb·°F) to J/(kg·K)
        rho_kg_m3 = rho_m * 16.0185 if rho_m is not None else 40000  # lb/ft³ to kg/m³
        Tsurf_c = (Tsurf - 32.0) * (5.0/9.0) if Tsurf is not None else 15  # °F to °C
    else:
        # Metric values - pass through as-is
        L2_m = L2
        L1_m = L1  
        D_m = D
        mdot_kg_s = mdot
        Tinj_c = Tinj
        grad_k_m = grad
        k_w_m_k = k
        c_j_kg_k = c_m if c_m is not None else 1000  # Default heat capacity
        rho_kg_m3 = rho_m if rho_m is not None else 2500  # Default density
        Tsurf_c = Tsurf if Tsurf is not None else 15  # Default surface temperature
    

    # Add TandP data to thermal_dict for SBT models
    if model != "HDF5" and tandp_data:
        if thermal_dict is None:
            thermal_dict = {}
        thermal_dict['TandP-data'] = tandp_data

    # Check if we have the required data before proceeding (this is normal on app startup)
    if not thermal_dict or not econ_dict:
        from plotly.graph_objects import Table
        tbl = Table()
        summary_dict = {}
        return tbl, summary_dict
    
    try:
        tbl, summary_dict = generate_summary_table(
                    mdot_kg_s, L2_m, L1_m, grad_k_m, D_m, Tinj_c, k_w_m_k, Drilling_cost_per_m, Discount_rate, Lifetime, 
                    Direct_use_heat_cost_per_kWth, Power_plant_cost_per_kWe, Pre_Cooling_Delta_T, Turbine_outlet_pressure, 
                    interp_time, case, fluid, model,
                    thermal_dict, econ_dict,
                    Tsurf=Tsurf_c, c_m=c_j_kg_k, rho_m=rho_kg_m3, Diameter1=Diameter1, Diameter2=Diameter2, 
                    PipeParam3=PipeParam3, PipeParam4=PipeParam4, PipeParam5=PipeParam5,
                    mesh=mesh, accuracy=accuracy, HyperParam3=HyperParam3, HyperParam4=HyperParam4, HyperParam5=HyperParam5,
                    units=units
        )
    except Exception as e:
        import traceback
        traceback.print_exc()
        # Return empty/default values on error
        from plotly.graph_objects import Table
        tbl = Table()
        summary_dict = {}

    return tbl, summary_dict


@app.callback(
    [Output(component_id='error_block_div1', component_property='children'),
     Output(component_id='error_block_div2', component_property='children'),
     Output(component_id='error_block_div3', component_property='children'),
     ],
    [Input(component_id='thermal-results-errors', component_property='data'),
     Input(component_id='thermal-contours-errors', component_property='data'),
     Input(component_id='econ-errors', component_property='data'),
     Input(component_id='econ-results', component_property='data'),
    ]
)

def update_error_divs(err_sub_dict, err_contour_dict, err_econ_dict, econ_results_dict):


    try:
        err_div1 = html.Div(#id="error_block_div1",
                            style={'display': 'none'})

        err_div2 = html.Div(#id="error_block_div2",
                            style={'display': 'none'})

        err_div3 = html.Div(#id="error_block_div3",
                            style={'display': 'none'})

        error_style = {'display': 'block',
                        'width': '100%',
                        "paddingTop": "8px",
                        "paddingLeft": "20px",
                        "paddingRight": "20px",
                        "paddingBottom": "5px",
                        'backgroundColor': lightbrown,
                        'color': darkergrey,
                        }
        
        if err_sub_dict and err_sub_dict != {}:
            try:
                error_message = next(iter(err_sub_dict.values()))
                if "No outputs" in error_message:
                    error_message = "No outputs could be calculated because there is not enough data at these limits. Consider changing parameter value(s)."

                err_div1 = html.Div(#id="error_block_div1",
                                    style=error_style,
                                    children=[
                                        html.Img(id="error-img1", src=app.get_asset_url('error.png')), 
                                        dcc.Markdown("**Did not plot visual(s).**", style={'display': 'inline-block'}),
                                        html.P(error_message),
                                        ]
                                    )
            except Exception as e:
                err_div1 = html.Div(style={'display': 'none'})

        if err_contour_dict and err_contour_dict != {}:
            try:
                error_message = next(iter(err_contour_dict.values()))
                err_div2 = html.Div(#id="error_block_div2",
                                    style=error_style,
                                    children=[
                                        html.Img(id="error-img2", src=app.get_asset_url('error.png')), 
                                        dcc.Markdown("**Did not plot visual(s).**", style={'display': 'inline-block'}),
                                        html.P(error_message),
                                        ]
                                    )
            except Exception as e:
                err_div2 = html.Div(style={'display': 'none'})

        if err_econ_dict and err_econ_dict != {}:
            try:
                # Check if plots are actually being generated by looking at econ_results_dict
                plots_have_data = False
                if econ_results_dict and isinstance(econ_results_dict, dict):
                    # Check if there are meaningful economic values (not just "-")
                    lcoe_sco2 = econ_results_dict.get('LCOE sCO2', '')
                    lcoe_h2o = econ_results_dict.get('LCOE H2O', '')
                    lcoh_sco2 = econ_results_dict.get('LCOH sCO2', '')
                    lcoh_h2o = econ_results_dict.get('LCOH H2O', '')
                    
                    # If any of these have meaningful values (not "-"), then plots are being generated
                    plots_have_data = (lcoe_sco2 != "-" or lcoe_h2o != "-" or 
                                     lcoh_sco2 != "-" or lcoh_h2o != "-")
                
                # Only show warning if plots are NOT being generated
                if not plots_have_data:
                    error_message = next(iter(err_econ_dict.values()))
                    if error_message and error_message.strip():
                        if "object has no attribute" in error_message:
                            error_message = "No outputs could be calculated because there is not enough data at these limits. Consider changing parameter value(s)."
                        err_div3 = html.Div(#id="error_block_div3",
                                            style=error_style,
                                            children=[
                                                html.Img(id="error-img3", src=app.get_asset_url('error.png')),
                                                dcc.Markdown("**Did not plot visual(s).**", style={'display': 'inline-block'}),
                                                html.P("No outputs could be calculated because there is not enough data at these limits. Consider changing parameter value(s)."),
                                                ]
                                            )
                    else:
                        err_div3 = html.Div(style={'display': 'none'})
                else:
                    # Plots are being generated, don't show warning
                    err_div3 = html.Div(style={'display': 'none'})
                    
            except Exception as e:
                err_div3 = html.Div(style={'display': 'none'})

        return err_div1, err_div2, err_div3
        
    except Exception as e:
        import traceback
        traceback.print_exc()
        # Return safe default values
        return html.Div(style={'display': 'none'}), html.Div(style={'display': 'none'}), html.Div(style={'display': 'none'})


@app.callback(
    Output(component_id='warning_block_div3', component_property='children'),
    Input(component_id='econ-memory', component_property='data'),
)

def update_error_divs(levelized_cost_dict):
    
    warning_div3 = html.Div(style={'display': 'none'})

    if levelized_cost_dict is None:
        return warning_div3
        
    # Check if levelized_cost_dict exists and has the required keys
    if not isinstance(levelized_cost_dict, dict):
        return warning_div3

    error_style = {'display': 'block',
                    'width': '100%',
                    "paddingTop": "8px",
                    "paddingLeft": "20px",
                    "paddingRight": "20px",
                    "paddingBottom": "5px",
                    'backgroundColor': lightbrown,
                    'color': darkergrey,
                    }

    # Don't show economic warning - plots are working fine
    # The LCOE values being "-" is normal when calculations have issues but plots still display

    return warning_div3





# -----------------------------------------------------------------------------
# App runs here. Define configurations, proxies, etc.
# -----------------------------------------------------------------------------

# Info popup callbacks are now handled in info_popups.py

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

@app.callback(
    Output("quick-unit-selector", "value"),
    Input("quick-unit-selector", "value"),
    prevent_initial_call=True
)
def handle_quick_unit_changes(unit_system):
    """Handle quick unit system changes"""
    
    if unit_system == "metric":
        # Apply metric units
        preferences = apply_metric_units()
        unit_converter.update_preferences(preferences)
        print(f"Applied metric units: {preferences}")
        return "metric"
    
    elif unit_system == "imperial":
        # Apply imperial units
        preferences = apply_imperial_units()
        unit_converter.update_preferences(preferences)
        print(f"Applied imperial units: {preferences}")
        return "imperial"
    
    return unit_system

@app.callback(
   [
    Output(component_id='drillcost-div', component_property='children'),
    Output(component_id='precool-div', component_property='children'),
    Output(component_id='turb-pout-div', component_property='children'),
   ],
   [Input(component_id="quick-unit-selector", component_property="value")],
   prevent_initial_call=True
    )
def update_economic_sliders(unit_system):
    """Update economic sliders with hardcoded values based on unit system"""
    
    # Define div style
    div_block_style = {"width": "98%", "margin": "auto", "margin-bottom": "10px", "display": "block"}
    
    if unit_system == "imperial":
        drillcost_container = html.Div(
            key=f"drillcost-imperial",
            children=[create_enhanced_slider(
                DivID="drillcost-div", ID="drillcost-select", 
                ptitle="Drilling Cost ($/ft)", 
                min_v=0, max_v=1220, 
                mark_dict={0: '0', 1220: '1.2k'}, 
                start_v=610, 
                div_style=div_block_style, 
                parameter_name="Drilling Cost ($/ft)"
            )]
        )
        precool_container = html.Div(
            key=f"precool-imperial",
            children=[create_enhanced_slider(
                DivID="precool-div", ID="precool-select", 
                ptitle="Pre-cooling (°F)", 
                min_v=32, max_v=104, 
                mark_dict={32: '32', 104: '104'}, 
                start_v=32, 
                div_style=div_block_style, 
                parameter_name="Pre-cooling (˚F)"
            )]
        )
        turb_pout_container = html.Div(
            key=f"turb-pout-imperial",
            children=[create_enhanced_slider(
                DivID="turb-pout-div", ID="turb-pout-select", 
                ptitle="Turbine Outlet Pressure (psia)", 
                min_v=1088, max_v=2901, 
                mark_dict={1088: '1088', 2901: '2901'}, 
                start_v=1160, 
                div_style=div_block_style, 
                parameter_name="Turbine Outlet Pressure (psia)"
            )]
        )
    else:
        # Metric drilling cost ($/m)
        drillcost_container = html.Div(
            key=f"drillcost-metric",
            children=[create_enhanced_slider(
                DivID="drillcost-div", ID="drillcost-select", 
                ptitle="Drilling Cost ($/m)", 
                min_v=0, max_v=4000, 
                mark_dict={0: '0', 4000: '4k'}, 
                start_v=2000, 
                div_style=div_block_style, 
                parameter_name="Drilling Cost ($/m)"
            )]
        )
        # Metric pre-cooling (°C)
        precool_container = html.Div(
            key=f"precool-metric",
            children=[create_enhanced_slider(
                DivID="precool-div", ID="precool-select", 
                ptitle="Pre-cooling (°C)", 
                min_v=0, max_v=40, 
                mark_dict={0: '0', 40: '40'}, 
                start_v=0, 
                div_style=div_block_style, 
                parameter_name="Pre-cooling (˚C)"
            )]
        )
        turb_pout_container = html.Div(
            key=f"turb-pout-metric",
            children=[create_enhanced_slider(
                DivID="turb-pout-div", ID="turb-pout-select", 
                ptitle="Turbine Outlet Pressure (bar)", 
                min_v=75, max_v=200, 
                mark_dict={75: '75', 200: '200'}, 
                start_v=80, 
                div_style=div_block_style, 
                parameter_name="Turbine Outlet Pressure (bar)"
            )]
        )
    
    return drillcost_container, precool_container, turb_pout_container