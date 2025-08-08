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
from info_popups import PARAMETER_INFO, create_info_button, create_enhanced_slider, create_enhanced_input_box, create_enhanced_dropdown


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
            slider_card(),
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
                            html.Div(id="dropdown-card5",
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
        dbc.Modal([
            dbc.ModalHeader(dbc.ModalTitle(id="info-modal-title")),
            dbc.ModalBody(id="info-modal-body"),
            dbc.ModalFooter(
                dbc.Button("Close", id="close-info-modal", className="ms-auto", n_clicks=0, style={"cursor": "pointer", "fontWeight": "bold"})
            ),
        ], id="info-modal", is_open=False, size="lg"),
    ],
)


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
        write_excelsheet(df_summary=df_tbl, df_subsurf_res_mass=df_mass_flow_rate, 
                                    df_subsurf_res_time=df_time, df_econ=df_econ, 
                                            geoCLUSTER_results_pathname=geoCLUSTER_results_pathname)
        return dcc.send_file(geoCLUSTER_results_pathname)
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

        print(" ----------------------------- ")
        print("HDF5")

        return {'display': 'block'}, {'display': 'block'}, {'display': 'block'}

    elif selected_model == "SBT V1.0" or selected_model == "SBT V2.0" : 
        
        print(" ----------------------------- ")
        print("SBT")
        # TODO: update tabs styline
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
            
            if visibility_state == param_list[0]:
                return n, n, b, b, b, b, b, \
                        n, n, n, n, n, b
            if visibility_state == param_list[1]:
                return n, b, n, b, b, b, b, \
                        n, n, n, n, n, b
            if visibility_state == param_list[2]:
                return n, b, b, n, b, b, b, \
                        n, n, n, n, n, b
            if visibility_state == param_list[3]:
                return n, b, b, b, n, b, b, \
                        n, n, n, n, n, b
            if visibility_state == param_list[4]: # "Injection Temperature (˚C)"
                return n, b, b, b, b, n, b, \
                        n, n, n, n, n, n
            if visibility_state == param_list[5]:
                return n, b, b, b, b, b, n, \
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
        raise PreventUpdate
    else:
        raise PreventUpdate


@app.callback(
   [
    Output(component_id='grad-container', component_property='children', allow_duplicate=True),
    Output(component_id='k-container', component_property='children', allow_duplicate=True),
    Output(component_id='Tinj-container', component_property='children', allow_duplicate=True),
    Output(component_id='mdot-container', component_property='children', allow_duplicate=True),
    Output(component_id='diameter-container', component_property='children'),
    Output(component_id='L2-container', component_property='children'),
    Output(component_id='L1-container', component_property='children'),
   ],
   [Input(component_id="model-select", component_property="value")],
   prevent_initial_call=True
    )

def update_slider_ranges(model):

    grad_dict = {0.01: '0.01', 0.1: '0.1'}
    k_dict = {0.1: '0.1', 7.0: '7.0'} # 1.5 4.5
    Tinj_dict = {30: '30', 60: '60'}
    mdot_dict = {5: '5', 300: '300'}
    D_dict = {0.2159: '0.2159', 0.4445: '0.4445'}
    L2_dict = {1000: '1000', 50000: '50000'}
    L1_dict = {1000: '1000', 10000: '10000'}

    if model == "HDF5":

        grad_container = create_enhanced_slider(DivID="grad-select-div", ID="grad-select", ptitle="Geothermal Gradient (K/m)", min_v=u_sCO2.grad[0], max_v=u_sCO2.grad[-1], 
                                                mark_dict=grad_dict, start_v=start_vals_d["grad"], div_style=div_block_style, parameter_name="Geothermal Gradient (K/m)")
        k_container = create_enhanced_slider(DivID="k-select-div", ID="k-select", ptitle="Rock Thermal Conductivity (W/m-K)", min_v=u_sCO2.k[0], max_v=u_sCO2.k[-1], 
                                                mark_dict=k_dict, start_v=start_vals_d["k"], div_style=div_block_style, parameter_name="Rock Thermal Conductivity (W/m-K)")
        Tinj_container = create_enhanced_slider(DivID="Tinj-select-div", ID="Tinj-select", ptitle="Injection Temperature (˚C)", min_v=30.0, max_v=60.0, 
                                                mark_dict=Tinj_dict, start_v=30.0, div_style=div_block_style, parameter_name="Injection Temperature (˚C)")
        mdot_container = create_enhanced_slider(DivID="mdot-select-div", ID="mdot-select", ptitle="Mass Flow Rate (kg/s)", min_v=u_sCO2.mdot[0], max_v=u_sCO2.mdot[-1], 
                                                mark_dict=mdot_dict, start_v=start_vals_d["mdot"], div_style=div_block_style, parameter_name="Mass Flow Rate (kg/s)")
        diameter_container = create_enhanced_slider(DivID="diameter-select-div", ID="diameter-select", ptitle="Borehole Diameter (m)", min_v=0.2159, max_v=0.4445, 
                                                mark_dict=D_dict, step_i=0.002, start_v=start_vals_d["D"], div_style=div_block_style, parameter_name="Borehole Diameter (m)")
        L2_container = create_enhanced_slider(DivID="L2-select-div", ID="L2-select", ptitle="Horizontal Extent (m)", min_v=u_sCO2.L2[0], max_v=u_sCO2.L2[-1], 
                                                mark_dict=L2_dict, start_v=start_vals_d["L2"], div_style=div_block_style, parameter_name="Horizontal Extent (m)")
        L1_container = create_enhanced_slider(DivID="L1-select-div", ID="L1-select", ptitle="Drilling Depth (m)", min_v=u_sCO2.L1[0], max_v=u_sCO2.L1[-1], 
                                                mark_dict=L1_dict, start_v=start_vals_d["L1"], div_style=div_block_style, parameter_name="Drilling Depth (m)")                       
                                 
        return grad_container, k_container, Tinj_container, mdot_container, diameter_container, L2_container, L1_container

    elif model == "SBT V1.0" or model == "SBT V2.0":

        grad_dict = {0.015: '0.015', 0.200: '0.200'}
        k_dict = {0.4: '0.4', 5.0: '5.0'}
        Tinj_dict = {30: '30', 100: '100'}
        mdot_dict = {5: '5', 300: '300'}
        D_dict = {0.2159: '0.2159', 0.4445: '0.4445'}
        L2_dict = {1000: '1000', 50000: '50000'}
        L1_dict = {1000: '1000', 10000: '10000'}

        grad_container = create_enhanced_slider(DivID="grad-select-div", ID="grad-select", ptitle="Geothermal Gradient (K/m)", # min_v=0.01, max_v=0.1,
                                                                     min_v=0.015, max_v=0.200, 
                                                                mark_dict=grad_dict, start_v=start_vals_d["grad"], div_style=div_block_style, parameter_name="Geothermal Gradient (K/m)")
        k_container = create_enhanced_slider(DivID="k-select-div", ID="k-select", ptitle="Rock Thermal Conductivity (W/m-K)", #min_v=0.1, max_v=7.0, 
                                                                     min_v=0.4, max_v=5.0, 
                                                                        mark_dict=k_dict, start_v=start_vals_d["k"], div_style=div_block_style, parameter_name="Rock Thermal Conductivity (W/m-K)")
        Tinj_container = create_enhanced_slider(DivID="Tinj-select-div", ID="Tinj-select", ptitle="Injection Temperature (˚C)", min_v=30.0, max_v=100.0, 
                                                                    # min_v=20.0, max_v=200.0, 
                                                                        mark_dict=Tinj_dict, start_v=start_vals_d["Tinj"], div_style=div_block_style, parameter_name="Injection Temperature (˚C)")
        mdot_container = create_enhanced_slider(DivID="mdot-select-div", ID="mdot-select", ptitle="Mass Flow Rate (kg/s)", min_v=5, max_v=300,
                                                                        # min_v=u_sCO2.mdot[0], max_v=u_sCO2.mdot[-1], 
                                                                    mark_dict=mdot_dict, start_v=start_vals_d["mdot"], div_style=div_block_style, parameter_name="Mass Flow Rate (kg/s)")
        diameter_container = create_enhanced_slider(DivID="diameter-select-div", ID="diameter-select", ptitle="Borehole Diameter (m)", min_v=0.2159, max_v=0.4445, 
                                                                        mark_dict=D_dict, step_i=0.002, start_v=start_vals_d["D"], div_style=div_none_style, parameter_name="Borehole Diameter (m)")
        L2_container =  create_enhanced_slider(DivID="L2-select-div", ID="L2-select", ptitle="Horizontal Extent (m)", min_v=1000, max_v=50000,
                                                                        #min_v=u_sCO2.L2[0], max_v=u_sCO2.L2[-1], 
                                                                    mark_dict=L2_dict, start_v=start_vals_d["L2"], div_style=div_block_style, parameter_name="Horizontal Extent (m)")
        L1_container = create_enhanced_slider(DivID="L1-select-div", ID="L1-select", ptitle="Drilling Depth (m)", min_v=1000, max_v=10000,
                                                                        #min_v=u_sCO2.L1[0], max_v=u_sCO2.L1[-1], 
                                                                    mark_dict=L1_dict, start_v=start_vals_d["L1"], div_style=div_block_style, parameter_name="Drilling Depth (m)")

        return grad_container, k_container, Tinj_container, mdot_container, diameter_container, L2_container, L1_container
    else:
        raise PreventUpdate


@app.callback(
   [
     Output(component_id='Diameter1-container', component_property='children'), # Diameter1
     Output(component_id='Diameter2-container', component_property='children'), # Diameter2
     Output(component_id='num-lat-container', component_property='children'), # PipeParam3
     Output(component_id='lat-allo-container', component_property='children'), # PipeParam4
     Output(component_id='lat-flow-container', component_property='children'), # PipeParam5
   ],
   [Input(component_id="model-select", component_property="value"),
    Input(component_id="case-select", component_property="value")],
   prevent_initial_call=True
    )
def update_sliders_heat_exchanger(model, case):

    if model == "SBT V1.0" or model == "SBT V2.0": 
        
        # Define dictionaries for SBT models
        radius_vertical_dict = {0.10795: '0.10795', 0.22225: '0.22225'}
        radius_lateral_dict = {0.10795: '0.10795', 0.22225: '0.22225'}

        if case == "utube":

            radius_vertical = create_enhanced_slider(DivID="radius-vertical-select-div", ID="radius-vertical-select", ptitle="Wellbore Radius Vertical (m)", min_v=0.10795, max_v=0.22225,
                                                                mark_dict=radius_vertical_dict, step_i=0.001, start_v=start_vals_sbt["radius-vertical"], div_style=div_block_style, parameter_name="Wellbore Radius Vertical (m)")
            radius_lateral = create_enhanced_slider(DivID="radius-lateral-select-div", ID="radius-lateral-select", ptitle="Wellbore Radius Lateral (m)", min_v=0.10795, max_v=0.22225,
                                                                mark_dict=radius_lateral_dict, step_i=0.001, start_v=start_vals_sbt["radius-lateral"], div_style=div_block_style, parameter_name="Wellbore Radius Lateral (m)")
            n_laterals = create_enhanced_input_box(DivID="num-lat-div", ID="n-laterals-select", ptitle="Number of Laterals", 
                                    min_v=0, max_v=20, start_v=start_vals_hdf5["n-laterals"], step_i=1, div_style=div_block_style, parameter_name="Number of Laterals")
            lateral_flow = create_enhanced_input_box(DivID="lat-allocation-div", ID="lateral-flow-select", ptitle="Lateral Flow Allocation", 
                                    min_v=0, max_v=1, start_v=start_vals_hdf5["lateral-flow"], step_i=0.01, div_style=div_block_style, parameter_name="Lateral Flow Allocation")
            lateral_multiplier = create_enhanced_input_box(DivID="lat-flow-mul-div", ID="lateral-multiplier-select", ptitle="Lateral Flow Multiplier", 
                                    min_v=0, max_v=1, start_v=start_vals_hdf5["lateral-multiplier"], step_i=0.05, div_style=div_block_style, parameter_name="Lateral Flow Multiplier")

            return radius_vertical, radius_lateral, n_laterals, lateral_flow, lateral_multiplier
        
        elif case == "coaxial":

            radius_centerpipe_dict = {0.0635: '0.0635', 0.174: '0.174'}
            thickness_centerpipe_dict = {0.005: '0.005', 0.025: '0.025'}
            insulation_thermal_k_dict = {0.025: '0.025', 0.50: '0.5'}

            radius = create_enhanced_slider(DivID="radius-vertical-select-div", ID="radius-vertical-select", ptitle="Wellbore Radius (m)", min_v=0.10795, max_v=0.22225,
                                                                mark_dict=radius_vertical_dict, step_i=0.001, start_v=start_vals_sbt["radius"], div_style=div_block_style, parameter_name="Wellbore Radius (m)")
            
            radiuscenterpipe = create_enhanced_slider(DivID="radius-lateral-select-div", ID="radius-lateral-select", ptitle="Center Pipe Radius (m)", min_v=0.0635, max_v=0.174, #  # Center Pipe Radius (coaxial)	0.0635	0.174
                                                                mark_dict=radius_centerpipe_dict, step_i=0.001, start_v=start_vals_sbt["radiuscenterpipe"], div_style=div_block_style, parameter_name="Center Pipe Radius (m)")
            
            thicknesscenterpipe = create_enhanced_slider(DivID="num-lat-div", ID="n-laterals-select", ptitle="Center Pipe Thickness (m)", min_v=0.005, max_v=0.025,
                                                                mark_dict=thickness_centerpipe_dict, step_i=0.001, start_v=start_vals_sbt["thicknesscenterpipe"], div_style=div_block_style, parameter_name="Center Pipe Thickness (m)")
            
            k_center_pipe = create_enhanced_slider(DivID="lateral-flow-select-div", ID="lateral-flow-select", ptitle="Insulation Thermal Conductivity (W/m-K)", min_v=0.025, max_v=0.5,
                                                                mark_dict=insulation_thermal_k_dict, step_i=0.001, start_v=start_vals_sbt["k_center_pipe"], div_style=div_block_style, parameter_name="Insulation Thermal Conductivity (W/m-K)")
            coaxialflowtype = create_enhanced_dropdown(DivID="lat-flow-mul-div", ID="lateral-multiplier-select", ptitle="Coaxial Flow Type", 
                                                                    options=["Inject in Annulus", "Inject in Center Pipe"], disabled=False, div_style=div_block_style, parameter_name="Coaxial Flow Type")
            # slider1(DivID="lat-flow-mul-div", ID="lateral-multiplier-select", ptitle="Coaxial Flow Type", min_v=1, max_v=2,
                                                             # mark_dict=radius_vertical_dict, step_i=1, start_v=start_vals_sbt["coaxialflowtype"], div_style=div_block_style)
                                                     
            return radius, radiuscenterpipe, thicknesscenterpipe, k_center_pipe, coaxialflowtype

            # 1 = CXA (fluid injection in annulus); 2 = CXC (fluid injection in center pipe)
    
    elif model == "HDF5":

        # Define dictionaries for HDF5 models
        radius_vertical_dict = {0.2: '0.2', 0.6: '0.6'}
        radius_lateral_dict = {0.2: '0.2', 0.6: '0.6'}

        radius_vertical = create_enhanced_slider(DivID="radius-vertical-select-div", ID="radius-vertical-select", ptitle="Wellbore Radius Vertical (m)", min_v=0.2, max_v=0.6,
                                                            mark_dict=radius_vertical_dict, step_i=0.001, start_v=start_vals_sbt["radius-vertical"], div_style=div_none_style, parameter_name="Wellbore Radius Vertical (m)")
        radius_lateral = create_enhanced_slider(DivID="radius-lateral-select-div", ID="radius-lateral-select", ptitle="Wellbore Radius Lateral (m)", min_v=0.2, max_v=0.6,
                                                            mark_dict=radius_lateral_dict, step_i=0.001, start_v=start_vals_sbt["radius-lateral"], div_style=div_none_style, parameter_name="Wellbore Radius Lateral (m)")
        n_laterals = create_enhanced_input_box(DivID="num-lat-div", ID="n-laterals-select", ptitle="Number of Laterals", 
                                min_v=0, max_v=20, start_v=start_vals_hdf5["n-laterals"], step_i=1, div_style=div_none_style, parameter_name="Number of Laterals")
        lateral_flow = create_enhanced_input_box(DivID="lat-allocation-div", ID="lateral-flow-select", ptitle="Lateral Flow Allocation", 
                                min_v=0, max_v=1, start_v=start_vals_hdf5["lateral-flow"], step_i=0.01, div_style=div_none_style, parameter_name="Lateral Flow Allocation")
        lateral_multiplier = create_enhanced_input_box(DivID="lat-flow-mul-div", ID="lateral-multiplier-select", ptitle="Lateral Flow Multiplier", 
                                min_v=0, max_v=1, start_v=start_vals_hdf5["lateral-multiplier"], step_i=0.05, div_style=div_none_style, parameter_name="Lateral Flow Multiplier")

        return radius_vertical, radius_lateral, n_laterals, lateral_flow, lateral_multiplier

    else:
        raise PreventUpdate

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

    ],
)

def update_subsurface_results_plots(interp_time, fluid, case, mdot, L2, L1, grad, D, Tinj, k_m, scale, model,
                        Tsurf, c_m, rho_m, 
                        # radius_vertical, radius_lateral, n_laterals, lateral_flow, lateral_multiplier,
                        Diameter1, Diameter2, PipeParam3, PipeParam4, PipeParam5,
                        mesh, accuracy, 
                        # mass_mode, temp_mode
                        HyperParam3, HyperParam4, HyperParam5
                        ):

    # -----------------------------------------------------------------------------
    # Creates and displays Plotly subplots of the subsurface results.
    # -----------------------------------------------------------------------------

    # print('subsurface')
    # if HDF5:
    # start = time.time()
    subplots, forty_yr_TPmeans_dict, df_mass_flow_rate, df_time, err_subres_dict, TandP_dict = generate_subsurface_lineplots(
        interp_time, fluid, case, mdot, L2, L1, grad, D, Tinj, k_m, scale, model,
        Tsurf, c_m, rho_m, 
        # radius_vertical, radius_lateral, n_laterals, lateral_flow, lateral_multiplier,
        Diameter1, Diameter2, PipeParam3, PipeParam4, PipeParam5,
        mesh, accuracy, HyperParam3, HyperParam4, HyperParam5
    )
    # if SBT:
    # end = time.time()
    # print("run generate_subsurface_lineplots:", end - start)

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
    ],
)

def update_subsurface_contours_plots(interp_time, fluid, case, param, mdot, L2, L1, grad, D, Tinj, k_m):

    # -----------------------------------------------------------------------------
    # Creates and displays Plotly subplots of the subsurface contours.
    # -----------------------------------------------------------------------------

    # print('contours')
    subplots, err_subcontour_dict = generate_subsurface_contours(
        interp_time, fluid, case, param, mdot, L2, L1, grad, D, Tinj, k_m
    )

    return subplots, err_subcontour_dict


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
     Input(component_id='radio-graphic-control4', component_property='value'),
     Input(component_id="checklist", component_property="value"),
     Input(component_id='model-select', component_property='value'),
     
    ],
)

def update_econ_plots(TandP_dict,
                     interp_time, fluid, case, end_use,
                      mdot, L2, L1, grad, D, Tinj, k_m,
                      Drilling_cost_per_m, Discount_rate, Lifetime, 
                      Direct_use_heat_cost_per_kWth, Power_plant_cost_per_kWe, Pre_Cooling_Delta_T, Turbine_outlet_pressure,
                      scale, checklist, model
                      ):

    # -----------------------------------------------------------------------------
    # Creates and displays Plotly subplots of the economic results.
    # -----------------------------------------------------------------------------

    # print('economics')

    if checklist == [' ']:
        is_plot_ts = True
    else:
        is_plot_ts = False

    economics_fig, econ_data_dict, econ_values_dict, err_econ_dict = generate_econ_lineplots(
        TandP_dict,
        interp_time, case, end_use, fluid, 
        mdot, L2, L1, grad, D, Tinj, k_m,
        Drilling_cost_per_m, Discount_rate, Lifetime, 
        Direct_use_heat_cost_per_kWth, Power_plant_cost_per_kWe, Pre_Cooling_Delta_T, Turbine_outlet_pressure,
        scale, 
        properties_H2O_pathname, 
        properties_CO2v2_pathname, 
        additional_properties_CO2v2_pathname,
        tmatrix_pathname,
        model,
        is_plot_ts
    )

    return economics_fig, econ_data_dict, econ_values_dict, err_econ_dict


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
    # Output(component_id="table", component_property="figure"),
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
     Input(component_id='econ-memory', component_property='data'),
     Input(component_id='thermal-memory', component_property='data'),
     Input(component_id='model-select', component_property='value'),
     Input(component_id='TandP-data', component_property='data'),
    ],
)

def update_table(interp_time, fluid, case, mdot, L2, L1, grad, D, Tinj, k,
                 Drilling_cost_per_m, Discount_rate, Lifetime, 
                 Direct_use_heat_cost_per_kWth, Power_plant_cost_per_kWe, Pre_Cooling_Delta_T, Turbine_outlet_pressure,
                 econ_dict, thermal_dict, model, tandp_data):

    print(f"DEBUG SUMMARY: model={model}, tandp_data keys={list(tandp_data.keys()) if tandp_data else 'None'}")
    print(f"DEBUG SUMMARY: thermal_dict keys={list(thermal_dict.keys()) if thermal_dict else 'None'}")
    print(f"DEBUG SUMMARY: econ_dict keys={list(econ_dict.keys()) if econ_dict else 'None'}")

    # Add TandP data to thermal_dict for SBT models
    if model != "HDF5" and tandp_data:
        thermal_dict['TandP-data'] = tandp_data
        print(f"DEBUG SUMMARY: Added TandP-data to thermal_dict")

    tbl, summary_dict = generate_summary_table(
                mdot, L2, L1, grad, D, Tinj, k, Drilling_cost_per_m, Discount_rate, Lifetime, 
                Direct_use_heat_cost_per_kWth, Power_plant_cost_per_kWe, Pre_Cooling_Delta_T, Turbine_outlet_pressure, 
                interp_time, case, fluid, model,
                thermal_dict, econ_dict
    )

    print(f"DEBUG SUMMARY: Table generated successfully for {model}")
    print(f"DEBUG SUMMARY: Table figure type: {type(tbl)}")
    print(f"DEBUG SUMMARY: Summary dict keys: {list(summary_dict.keys()) if summary_dict else 'None'}")

    return tbl, summary_dict


@app.callback(
    [Output(component_id='error_block_div1', component_property='children'),
     Output(component_id='error_block_div2', component_property='children'),
     Output(component_id='error_block_div3', component_property='children'),
     ],
    [Input(component_id='thermal-results-errors', component_property='data'),
     Input(component_id='thermal-contours-errors', component_property='data'),
     Input(component_id='econ-errors', component_property='data'),
    ]
)

def update_error_divs(err_sub_dict, err_contour_dict, err_econ_dict):
    
    # print(err_sub_dict)
    # print(err_contour_dict)
    # print(err_econ_dict)
    # print('\n')

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
    if err_sub_dict != {}:


        error_message = next(iter(err_sub_dict.values()))

        if "No outputs" in error_message:
            error_message = "No outputs were able to be calculated because there are not enough data at these limits. Consider changing parameter value(s)."

        err_div1 = html.Div(#id="error_block_div1",
                            style=error_style,
                            children=[
                                html.Img(id="error-img1", src=app.get_asset_url('error.png')), 
                                dcc.Markdown("**Did not plot visual(s).**", style={'display': 'inline-block'}),
                                html.P(error_message),
                                ]
                        )


    if err_contour_dict != {}:

        error_message = next(iter(err_contour_dict.values()))

        err_div2 = html.Div(#id="error_block_div2",
                            style=error_style,
                            children=[
                                html.Img(id="error-img2", src=app.get_asset_url('error.png')), 
                                dcc.Markdown("**Did not plot visual(s).**", style={'display': 'inline-block'}),
                                html.P(error_message),
                                ]
                        )

    if err_econ_dict != {}:

        error_message = next(iter(err_econ_dict.values()))

        if "object has no attribute" in error_message:
            error_message = "No outputs were able to be calculated because there are not enough data at these limits. Consider changing parameter value(s)."

        err_div3 = html.Div(#id="error_block_div3",
                            style=error_style,
                            children=[
                                html.Img(id="error-img3", src=app.get_asset_url('error.png')),
                                dcc.Markdown("**Did not plot visual(s).**", style={'display': 'inline-block'}),
                                html.P(error_message),
                                ]
                        )

    return err_div1, err_div2, err_div3


@app.callback(
    Output(component_id='warning_block_div3', component_property='children'),
    Input(component_id='econ-memory', component_property='data'),
)

def update_error_divs(levelized_cost_dict):
    
    warning_div3 = html.Div(style={'display': 'none'})

    error_style = {'display': 'block',
                    'width': '100%',
                    "paddingTop": "8px",
                    "paddingLeft": "20px",
                    "paddingRight": "20px",
                    "paddingBottom": "5px",
                    'backgroundColor': lightbrown,
                    'color': darkergrey,
                    }

    if levelized_cost_dict['LCOE sCO2'] == "9999.00" or levelized_cost_dict['LCOE H2O'] == "9999.00":
        
        warning_div3 = html.Div(#id="error_block_div3",
                            style=error_style,
                            children=[
                                html.Img(id="warning-img", src=app.get_asset_url('warning.png')), 
                                dcc.Markdown("**LCOE is too high.**", style={'display': 'inline-block'}),
                                html.P("Outlet Temperature (°C) may be too low or the system is losing heat (negative kWe)."),
                                ]
                        )


    return warning_div3





# -----------------------------------------------------------------------------
# App runs here. Define configurations, proxies, etc.
# -----------------------------------------------------------------------------

@app.callback(
    [Output("info-modal", "is_open"),
     Output("info-modal-title", "children"),
     Output("info-modal-body", "children")],
    [
        Input("info-btn-surface-temperature-degc", "n_clicks"),
        Input("info-btn-geothermal-gradient-k-m", "n_clicks"),
        Input("info-btn-rock-thermal-conductivity-w-m-k", "n_clicks"),
        Input("info-btn-rock-specific-heat-capacity-j-kg-k", "n_clicks"),
        Input("info-btn-rock-density-kg-m3", "n_clicks"),
        Input("info-btn-injection-temperature-degc", "n_clicks"),
        Input("info-btn-mass-flow-rate-kg-s", "n_clicks"),
        Input("info-btn-borehole-diameter-m", "n_clicks"),
        Input("info-btn-horizontal-extent-m", "n_clicks"),
        Input("info-btn-drilling-depth-m", "n_clicks"),
        Input("info-btn-drilling-cost-m", "n_clicks"),
        Input("info-btn-discount-rate-%", "n_clicks"),
        Input("info-btn-lifetime-years", "n_clicks"),
        Input("info-btn-plant-capex-kwt", "n_clicks"),
        Input("info-btn-plant-capex-kwe", "n_clicks"),
        Input("info-btn-pre-cooling-degc", "n_clicks"),
        Input("info-btn-turbine-outlet-pressure-bar", "n_clicks"),
        Input("info-btn-mesh-fineness", "n_clicks"),
        Input("info-btn-accuracy", "n_clicks"),
        Input("info-btn-number-of-laterals", "n_clicks"),
        Input("info-btn-lateral-flow-allocation", "n_clicks"),
        Input("info-btn-lateral-flow-multiplier", "n_clicks"),
        Input("info-btn-fluid-properties-mode", "n_clicks")
    ],
    [State("info-modal", "is_open")],
    prevent_initial_call=True,
    suppress_callback_exceptions=True
)
def toggle_info_modal(*args):
    info_clicks = args[:-1]
    is_open = args[-1]
    
    # Get the triggered button ID
    from dash import ctx
    triggered_id = ctx.triggered_id if ctx.triggered_id else None
    
    # Debug: print the triggered ID and all click values
    print(f"Tooltip callback triggered. ID: {triggered_id}")
    print(f"All click values: {info_clicks}")
    
    # Check if any button was actually clicked (n_clicks > 0)
    button_clicked = False
    if triggered_id:
        # Find the index of the triggered button in the args
        button_ids = [
            "info-btn-surface-temperature-degc",
            "info-btn-geothermal-gradient-k-m",
            "info-btn-rock-thermal-conductivity-w-m-k",
            "info-btn-rock-specific-heat-capacity-j-kg-k",
            "info-btn-rock-density-kg-m3",
            "info-btn-injection-temperature-degc",
            "info-btn-mass-flow-rate-kg-s",
            "info-btn-borehole-diameter-m",
            "info-btn-horizontal-extent-m",
            "info-btn-drilling-depth-m",
            "info-btn-drilling-cost-m",
            "info-btn-discount-rate-%",
            "info-btn-lifetime-years",
            "info-btn-plant-capex-kwt",
            "info-btn-plant-capex-kwe",
            "info-btn-pre-cooling-degc",
            "info-btn-turbine-outlet-pressure-bar",
            "info-btn-mesh-fineness",
            "info-btn-accuracy",
            "info-btn-number-of-laterals",
            "info-btn-lateral-flow-allocation",
            "info-btn-lateral-flow-multiplier",
            "info-btn-fluid-properties-mode"
        ]
        
        try:
            button_index = button_ids.index(triggered_id)
            # Handle case where the component might not exist (None value)
            if button_index < len(info_clicks) and info_clicks[button_index] is not None and info_clicks[button_index] > 0:
                button_clicked = True
                print(f"Button {triggered_id} was clicked (n_clicks: {info_clicks[button_index]})")
            else:
                print(f"Button {triggered_id} was not actually clicked (n_clicks: {info_clicks[button_index] if button_index < len(info_clicks) else 'N/A'})")
        except (ValueError, IndexError):
            print(f"Button {triggered_id} not found in button list")
    
    parameter_names = [
        "Surface Temperature (˚C)",
        "Geothermal Gradient (K/m)",
        "Rock Thermal Conductivity (W/m-K)",
        "Rock Specific Heat Capacity (J/kg-K)",
        "Rock Density (kg/m3)",
        "Injection Temperature (˚C)",
        "Mass Flow Rate (kg/s)",
        "Borehole Diameter (m)",
        "Horizontal Extent (m)",
        "Drilling Depth (m)",
        "Drilling Cost ($/m)",
        "Discount Rate (%)",
        "Lifetime (years)",
        "Plant CAPEX ($/kWt)",
        "Plant CAPEX ($/kWe)",
        "Pre-cooling (˚C)",
        "Turbine Outlet Pressure (bar)",
        "Mesh Fineness",
        "Accuracy",
        "Number of Laterals",
        "Lateral Flow Allocation",
        "Lateral Flow Multiplier",
        "Fluid Properties Mode"
    ]
    
    # Map button IDs to parameter names
    button_to_param = {
        "info-btn-surface-temperature-degc": "Surface Temperature (˚C)",
        "info-btn-geothermal-gradient-k-m": "Geothermal Gradient (K/m)",
        "info-btn-rock-thermal-conductivity-w-m-k": "Rock Thermal Conductivity (W/m-K)",
        "info-btn-rock-specific-heat-capacity-j-kg-k": "Rock Specific Heat Capacity (J/kg-K)",
        "info-btn-rock-density-kg-m3": "Rock Density (kg/m3)",
        "info-btn-injection-temperature-degc": "Injection Temperature (˚C)",
        "info-btn-mass-flow-rate-kg-s": "Mass Flow Rate (kg/s)",
        "info-btn-borehole-diameter-m": "Borehole Diameter (m)",
        "info-btn-wellbore-radius-vertical-m": "Wellbore Radius Vertical (m)",
        "info-btn-wellbore-radius-lateral-m": "Wellbore Radius Lateral (m)",
        "info-btn-horizontal-extent-m": "Horizontal Extent (m)",
        "info-btn-drilling-depth-m": "Drilling Depth (m)",
        "info-btn-drilling-cost-m": "Drilling Cost ($/m)",
        "info-btn-discount-rate-%": "Discount Rate (%)",
        "info-btn-lifetime-years": "Lifetime (years)",
        "info-btn-plant-capex-kwt": "Plant CAPEX ($/kWt)",
        "info-btn-plant-capex-kwe": "Plant CAPEX ($/kWe)",
        "info-btn-pre-cooling-degc": "Pre-cooling (˚C)",
        "info-btn-turbine-outlet-pressure-bar": "Turbine Outlet Pressure (bar)",
        "info-btn-mesh-fineness": "Mesh Fineness",
        "info-btn-accuracy": "Accuracy",
        "info-btn-number-of-laterals": "Number of Laterals",
        "info-btn-lateral-flow-allocation": "Lateral Flow Allocation",
        "info-btn-lateral-flow-multiplier": "Lateral Flow Multiplier",
        "info-btn-fluid-properties-mode": "Fluid Properties Mode"
    }
    
    if button_clicked and triggered_id and triggered_id in button_to_param:
        param = button_to_param[triggered_id]
        info = PARAMETER_INFO.get(param, None)
        if info:
            print(f"Found info for parameter: {param}")
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
            return True, f"Information: {param}", modal_content
        else:
            print(f"No info found for parameter: {param}")
    else:
        if triggered_id:
            print(f"Button {triggered_id} was not actually clicked or not found in mapping")
        else:
            print("No button was triggered")
    
    # If no button was clicked, return current state
    return is_open, "", []





# Separate callback for close button to ensure it works
@app.callback(
    Output("info-modal", "is_open", allow_duplicate=True),
    Input("close-info-modal", "n_clicks"),
    prevent_initial_call=True
)
def close_modal(n_clicks):
    if n_clicks and n_clicks > 0:
        return False
    raise PreventUpdate



server = app.server 
# from app import server as application # in the wsgi.py file -- this targets the Flask server of Dash app

compress = Compress()
compress.init_app(app.server)     # gzip all static assets

app.server.config.update(
    SESSION_COOKIE_SECURE=True,
    SESSION_COOKIE_HTTPONLY=True,
    SESSION_COOKIE_SAMESITE="Lax"
)

"""" 
1. Dash/Flask code no longer issuing a Set-Cookie header.
2. No client-side script (such as the old Google-Analytics snippet) is writing a cookie anymore.
"""
    
@server.route("/.well-known/apple-app-site-association")
def aasa():
    return send_from_directory(
        "static",
        "apple-app-site-association",   # apple-app-site-association is purely a hint to iOS / iPadOS about when it may open a native app instead of Safari
        mimetype="application/json"
    )

if __name__ == '__main__':
    # app.run_server(port=8060, debug=True) 
    app.run_server(
        # host="127.0.0.1",
        port=8060,
        debug=False, # needs to be False in production
        ssl_context="adhoc" 
    )

    """
    ssl_context="adhoc"

    Flask, and more specifically Werkzeug, support the use of on-the-fly certificates, 
    which are useful to quickly serve an application over HTTPS without having to mess 
    with certificates. All you need to do, is add ssl_context='adhoc' to your app.run() call
    """
    # EXAMPLE: *************
    # app.run_server(port=8050, proxy="http://127.0.0.1:8059::https://<site>/<page_name>")



#  show_hide_element
    # Output(component_id='sCO2-card', component_property='style', allow_duplicate=True),
    #  Output(component_id='Tsurf-select-div', component_property='style'),
    #  Output(component_id='c-select-div', component_property='style'),
    #  Output(component_id='rho-select-div', component_property='style'),
    #  Output(component_id='radius-vertical-select-div', component_property='style'),
    #  Output(component_id='radius-lateral-select-div', component_property='style'),
    #  Output(component_id='n-laterals-select-div', component_property='style'),
    #  Output(component_id='lateral-flow-select-div', component_property='style'),
    #  Output(component_id='lateral-multiplier-select-div', component_property='style'),