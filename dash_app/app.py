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

import time
import traceback
from fractions import Fraction

# web app and interactive graphics libraries
import dash
import dash_bootstrap_components as dbc  # Adds bootstrap components for more web themes and templates
# import dash._callback_context as ctx
import dash_daq as daq  # Adds more data acquisition (DAQ) and controls to dash callbacks
from dash import Dash, ctx, dcc, html, no_update
from dash.dependencies import Input, Output, State, ALL, MATCH
from dash.exceptions import PreventUpdate
from dropdowns import *
# from dropdowns import get_model_options
from flask import send_from_directory
from flask_compress import Compress
from flask_talisman import Talisman
import plotly.graph_objects as go
from plotly.subplots import make_subplots

# sourced scripts
from paths import inpath_dict
from sliders import *  # u_sCO2, u_H2O, c_sCO2, c_H2O, and imports functions from plots.py
from tables import generate_summary_table
from text import *
from info_popups import (
    PARAMETER_INFO,
    create_info_button,
    create_info_modal,
    register_info_modal_callbacks,
)
from write2excel import write_excelsheet

is_print = False

if is_print:
    print(" ============= App start =============")
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
additional_properties_CO2v2_pathname = inpath_dict[
    "additional_properties_CO2v2_pathname"
]
tmatrix_pathname = inpath_dict["tmatrix_pathname"]

# sha384 for bootstrap@5.3.1/dist/css/bootstrap.min.css (jsDelivr)
BOOTSTRAP_CSS = {
    "rel": "stylesheet",
    "href": "https://cdn.jsdelivr.net/npm/bootstrap@5.3.1/dist/css/bootstrap.min.css",
    "integrity": "sha384-4bw+/aepP/YC94hEpVNVgiZdgIC5+VKNBQNGCHeKRQN+PtmoHDEXuppvnDJzQIu9",
    "crossorigin": "anonymous",
}
# print(dbc.themes.BOOTSTRAP)


app = dash.Dash(
    __name__,
    assets_folder="assets",
    # external_stylesheets=[dbc.themes.BOOTSTRAP],
    external_stylesheets=[BOOTSTRAP_CSS],  # ← use the dict, not dbc.themes.BOOTSTRAP
    # url_base_pathname=url_base_pathname, # not needed
    requests_pathname_prefix=requests_pathname_prefix,
    suppress_callback_exceptions=True,
    meta_tags=[
        {"name": "viewport", "content": "width=device-width"},
        {
            "name": "description",
            "content": """Welcome to the Geothermal Closed-Loop User Tool in Energy Research (GeoCLUSTER). 
                                        This research was funded by the Geothermal Technologies Office (GTO) within the 
                                        Office of Energy Efficiency and  Renewable Energy (EERE) at the U.S. Department of 
                                        Energy (DOE) to form a collaborative study of closed-loop geothermal systems.""",
        },
    ],
)
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

tabs = [
    "about-tab",
    "energy-time-tab",
    "energy-tab",
    "econ-time-tab",
    "summary-tab",
]  # there are values and ids


plotly_config = {
    "displaylogo": False,
    "modeBarButtonsToRemove": [
        "autoScale",
        "resetScale",
    ],  # High-level: zoom, pan, select, zoomIn, zoomOut, autoScale, resetScale
    "toImageButtonOptions": {
        "format": "png",  # one of png, svg, jpeg, webp
        "filename": "custom_image",
        "height": None,
        "width": None,
        "scale": 6,  # Multiply title/legend/axis/canvas sizes by this factor
    },
}



darkergrey = "#383838"
lightbrown = "#ede6dd"

dropdown_guidance_style = {
    "position": "relative",
    "marginTop": "28px",
    "height": "25px",
    "width": "100%",
    "paddingBottom": "7px",
    "paddingLeft": "20px",
    "paddingRight": "20px",
    "border": "1.5px white solid",
    "backgroundColor": lightbrown,
    "color": darkergrey,
    "fontWeight": "bold",
    "fontSize": "11px",
}

carousel_image_style = {"objectFit": "contain"}


HIDDEN_LATERAL_FLOW_STYLE = {"position": "absolute", "left": "-9999px", "visibility": "hidden"}


def hidden_lateral_flow_input(value=1.0):
    """Create the hidden lateral flow input so callbacks always have a target component."""

    return dcc.Input(
        id="lateral-flow-select",
        value=value,
        type="number",
        min=0,
        max=1,
        step=0.01,
        style=HIDDEN_LATERAL_FLOW_STYLE,
    )


def lateral_flow_placeholder(value=1.0, style=None):
    container_style = style if style is not None else div_none_style
    return html.Div(
        id="lat-allocation-div",
        style=container_style,
        children=[hidden_lateral_flow_input(value)],
    )


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
            html.Img(id="logo", src=app.get_asset_url("logo3.png")),
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
            html.Div(
                id="scenario1-div",
                children=[
                    html.Button(
                        html.Div(
                            children=[
                                html.P("Run At Commercial Scale"),
                                html.P(scenario1_text, id="scenario1_text"),
                            ]
                        ),
                        id="btn-nclicks-1",
                        n_clicks=0,
                    ),
                ],
            ),
            # html.Div(id='scenario2-div',
            #          children=[
            #             html.Button('Optimize Power Output',
            #             id='btn-nclicks-2', n_clicks=0),
            #             html.P(scenario2_text, id="scenario2_text"),
            # ]),
            html.Div(
                id="scenario3-div",
                children=[
                    html.Button(
                        html.Div(
                            children=[
                                html.P("Optimize Economic Competitiveness"),
                                html.P(scenario3_text, id="scenario3_text"),
                            ]
                        ),
                        # 'Optimize Economic Competitiveness',
                        id="btn-nclicks-3",
                        n_clicks=0,
                    ),
                    # html.P(scenario3_text, id="scenario3_text"),
                ],
            ),
            html.Hr(id="hr-break2"),
            dropdown_card(),
            html.Br(),
            html.Br(),
            slider_card(),
            # disclaimer
            html.P(disclaimer_text, id="disclaimer-text"),
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
                style=dropdown_guidance_style,
            ),
            dbc.Collapse(
                dbc.Card(children=[dbc.CardBody(children=[dropdown_children])]),
                id=cardID,
                is_open=False,
            ),
        ],
        style={"display": "block"}
    )


def create_resource_card(card_id, title, description, link=None, pdf_path=None, expanded_content=None):
    """Create an expandable resource card with title, description, and expand icon."""
    expand_icon_id = {"type": "resource-expand-icon", "index": card_id}
    collapse_id = {"type": "resource-collapse", "index": card_id}
    button_id = {"type": "resource-button", "index": card_id}
    
    # Determine if card should be clickable (has link or pdf) or expandable (has expanded_content)
    is_clickable = link is not None or pdf_path is not None
    has_expandable_content = expanded_content is not None
    
    # Build expanded content if provided
    if expanded_content is None and not is_clickable:
        expanded_content = html.Div()
    
    # Create click handler URL
    click_url = None
    if link:
        click_url = link
    elif pdf_path:
        click_url = app.get_asset_url(pdf_path)
    
    header_style = {
        "position": "relative",
        "cursor": "pointer",
        "padding": "15px",
        "border": "1px solid #C9C9C9",
        "borderRadius": "16px",
        "backgroundColor": "#fff",
    }
    
    # Wrap header in link if clickable
    header_children = [
        html.Div(
            style={"paddingRight": "35px"},
            children=[
                html.H6(title, style={"margin": "0", "fontWeight": "bold", "fontSize": "16px", "color": "rgb(50, 50, 50)"}),
                html.Div(description, style={"margin": "5px 0 0 0", "fontSize": "14px", "color": "rgb(50, 50, 50)", "lineHeight": "1.4"}),
            ]
        ),
        html.Img(
            id=expand_icon_id if has_expandable_content else f"{card_id}-icon",
            src=app.get_asset_url("expand_closed.svg"),
            className="resource-expand-icon",
            style={
                "width": "21px",
                "height": "21px",
                "position": "absolute",
                "top": "15px",
                "right": "15px",
            },
        ),
    ]
    
    # Check if description contains clipboard-wrapper (special case for API card)
    has_clipboard = False
    if isinstance(description, html.Div):
        # Check if any child has clipboard-wrapper id
        def check_for_clipboard(children):
            if isinstance(children, list):
                for child in children:
                    if check_for_clipboard(child):
                        return True
            elif hasattr(children, 'id') and children.id == 'clipboard-wrapper':
                return True
            elif hasattr(children, 'children'):
                return check_for_clipboard(children.children)
            return False
        has_clipboard = check_for_clipboard(description)
    
    if is_clickable and not has_expandable_content:
        if has_clipboard:
            # For cards with clipboard, use div with onClick to open link, but exclude clipboard clicks
            header_element = html.Div(
                className="resource-card-header",
                style=header_style,
                children=header_children,
                **{"data-href": click_url, "data-target": "_blank" if link else None}
            )
        else:
            # Card opens link directly - wrap in anchor tag
            header_element = html.A(
                className="resource-card-header",
                href=click_url,
                target="_blank",
                rel="noopener noreferrer",
                style=header_style,
                children=header_children,
            )
    else:
        # Card expands - use div with click handler
        div_props = {
            "className": "resource-card-header",
            "style": header_style,
            "children": header_children,
        }
        if has_expandable_content:
            div_props["id"] = button_id
            div_props["n_clicks"] = 0
        header_element = html.Div(**div_props)
    
    card_children = [header_element]
    if has_expandable_content:
        card_children.append(
            dbc.Collapse(
                id=collapse_id,
                is_open=False,
                children=expanded_content,
            )
        )
    
    return html.Div(
        id=card_id,
        className="resource-card",
        children=card_children,
    )


# -------------------------------------------------------------------------------------------------
# Defines tab contents for 5 of the following tabs:
#
#   About, Subsurface Results, Subsurface Contours, Economic Results, and Summary
#
# -------------------------------------------------------------------------------------------------

about_tab = dcc.Tab(
    label="About",
    id="about-tab",
    value="about-tab",
    selected_className="active_tabs",
    children=[
        html.Hr(className="tab-hr"),
        html.Div(
            id="demo",
            children=[
                html.Div(
                    id="intro2",
                    children=[
                        html.P("About Our Research", id="ab-title1"),
                        html.P(note, id="ab-note"),
                        html.P(note2, id="ab-note2"),
                        html.P("Closed-Loop Geothermal Working Group", id="ab-title4"),
                        html.P("This research was funded by the Geothermal Technologies Office (GTO) within the Office of Energy Efficiency and Renewable Energy (EERE) at the U.S. Department of Energy (DOE) to form a collaborative study of CLGSs involving four national laboratories and two universities.", id="ab-note5"),
                        html.P("Navigating the Results", id="ab-title2"),
                        html.P(note3, id="ab-note3"),
                    ],
                ),
                html.Div(
                    id="image-container",
                    children=[
                        html.P("Closed-Loop Geothermal System", id="ab-title5"),
                        dbc.Carousel(
                            id="carousel-ride",
                            items=[
                                {
                                    "key": "1",
                                    "src": app.get_asset_url("CLGS-1.3.png"),
                                    "img_style": carousel_image_style,
                                },
                                {
                                    "key": "2",
                                    "src": app.get_asset_url("CLGS-3.2.png"),
                                    "img_style": carousel_image_style,
                                },
                            ],
                            variant="dark",
                            interval=3000,
                            controls=True,
                            indicators=False,
                        ),
                        html.P("Acknowledgments", id="ab-title6", style={"marginTop": "30px"}),
                        html.P([
                            "This research was funded by the Geothermal Technologies Office (GTO) within the Office of Energy Efficiency and Renewable Energy (EERE) at the U.S. Department of Energy (DOE) to form a collaborative study of CLGSs involving four national laboratories and two universities: ",
                            "Idaho National Laboratory, ",
                            "National Laboratory of the Rockies, ",
                            "Pacific Northwest National Laboratory, ",
                            "Sandia National Laboratories, ",
                            "Stanford University, ",
                            "and PennState. ",
                            "Last updated December 2025.",
                        ], id="ab-note6"),
                    ],
                ),
            ],
        ),
        html.Div(
            id="resources-section",
            style={"width": "100%", "clear": "both", "padding": "0"},
            children=[
                html.P("Resources", id="ab-title3"),
                html.P("Explore our APIs, research papers, and open-source code.", id="ab-note4"),
                html.Div(
                    id="resource-cards-container",
                    children=[
                        html.Div(
                            id="resource-column-1",
                            children=[
                                create_resource_card(
                                    "resource-card-1",
                                    "GeoCLUSTER Code Repository",
                                    "Access the code behind the Python-based closed-loop geothermal techno-economic simulator with customizable reservoir and wellbore models. Last updated 12/2025.",
                                    link="https://github.com/pnnl/GeoCLUSTER",
                                ),
                                create_resource_card(
                                    "resource-card-2",
                                    "GeoCLUSTER APIs",
                                    "Access the closed-loop geothermal techno-economic model programmatically for scenario runs, individually or in batch. Last updated 12/2025.",
                                    link="https://colab.research.google.com/drive/1MDtSh6ymGeOTGAXI57BygN2D-PXFD-bX?usp=sharing",
                                ),
                                create_resource_card(
                                    "resource-card-6",
                                    "Geothermal Rising Conference, 2025",
                                    [
                                        "Raquel S.P. Hakes, Radoslav Bozinoski, Jassim Aljubran, Gabriela B. Anleu, Ryan P. Abernathey, Anastasia Bernat, and Aaron C. Buchko. ",
                                        html.Span("Multi-Site Techno-Economic Analysis of Closed-Loop Geothermal Systems.", style={"fontWeight": "bold"}),
                                        " Geothermal Rising Conference, 2025.",
                                    ],
                                    pdf_path="pdfs/Hakes et al - 2025 GRC - FINAL (1).pdf",
                                ),
                                create_resource_card(
                                    "resource-card-3",
                                    "Stanford Geothermal Workshop, 2025",
                                    [
                                        "Anastasia Bernat, Alexander Buchko, Koenraad Beckers, and Aaron Moreno. ",
                                        html.Span("GeoCLUSTER v2.0: A Closed-Loop, Techno-Economic Simulator Supporting New Case Studies.", style={"fontWeight": "bold"}),
                                        " Proceedings of the 50th Workshop on Geothermal Reservoir Engineering, Stanford University, February 10–12, 2025.",
                                    ],
                                    pdf_path="pdfs/GeoCLUSTER v2.0 A Closed-Loop.pdf",
                                ),
                                create_resource_card(
                                    "resource-card-5",
                                    "Stanford Geothermal Workshop, 2025",
                                    [
                                        "Gabriela A. Bran-Anleu, Raquel S.P. Hakes, Radoslav Bozinoski, and Koenraad Beckers. ",
                                        html.Span("A Parametric Study of L-Shape Coaxial Closed-Loop Geothermal Systems with Reservoir Convection.", style={"fontWeight": "bold"}),
                                        " Proceedings of the 50th Workshop on Geothermal Reservoir Engineering, Stanford University, February 12-14, 2025.",
                                    ],
                                    pdf_path="pdfs/A parametric study of Lshape.pdf",
                                ),
                                create_resource_card(
                                    "resource-card-14",
                                    "SBT Code",
                                    [
                                        "Access the slender-body theory (SBT) simulator - a Python-based tool for simulating the production temperatures, production pressures and heat extraction with closed-loop geothermal systems. Last updated 01/2025. ",
                                        html.Br(),
                                        "National Renewable Energy Laboratory.",
                                    ],
                                    link="https://github.com/NREL/SBT",
                                ),
                                create_resource_card(
                                    "resource-card-13",
                                    "DOE Report, 2024",
                                    [
                                        "U.S. Department of Energy. ",
                                        html.Span("Pathways to Commercial Liftoff: Next-Generation Geothermal Power.", style={"fontWeight": "bold"}),
                                        " March 2024.",
                                    ],
                                    pdf_path="pdfs/DOE Report 2024_Pathways to Commercial liftoff.pdf",
                                ),
                                create_resource_card(
                                    "resource-card-4",
                                    "Geothermics, 2024",
                                    [
                                        "Mark White, Yaroslav Vasyliv, Koenraad Beckers, Mario Martinez, Paolo Balestra, Carlo Parisi, Chad Augustine, Gabriela Bran-Anleu, Roland Horne, Laura Pauley, Giorgia Bettini, Theron Marshall, and Anastasia Bernat. ",
                                        html.Span("Numerical Investigation of Closed-Loop Geothermal Systems in Deep Geothermal Reservoirs.", style={"fontWeight": "bold"}),
                                        " Geothermics 116 (2024): 102852.",
                                    ],
                                    pdf_path="pdfs/Numerical investigation of .pdf",
                                ),
                                create_resource_card(
                                    "resource-card-15",
                                    "Geothermal Data Repository, 2023",
                                    [
                                        "Koenraad Beckers, Roland Horne, Chad Augustine, Laura Pauley, Doug Hollett, Andy Adams, Doug Blankenship, Zach Frone, Sean Porse, Seunghwan Baek, Paolo Balestra, Anastasia Bernat, Giorgia Bettini, Gabriela Bran-Anleu, Alec Kucala, Brian Kyanjo, Theron Marshall, Mario Martinez, Travis McLing, Carlo Parisi, Sam Subia, Yaroslav Vasyliv, and Mark White. ",
                                        html.Span("Closed Loop Geothermal Working Group: GeoCLUSTER App, Subsurface Simulation Results, and Publications.", style={"fontWeight": "bold"}),
                                        " Pacific Northwest National Laboratory, February, 3, 2023. Distributed by Geothermal Data Repository.",
                                    ],
                                    link="https://doi.org/10.15121/1972213",
                                ),
                            ],
                        ),
                        html.Div(
                            id="resource-column-2",
                            children=[
                                create_resource_card(
                                    "resource-card-7",
                                    "Stanford Geothermal Workshop, 2023",
                                    [
                                        "Carlo Parisi, Paolo Balestra, Brian Kyanjo, Theron D. Marshall, Travis L. Meling, and Mark D. White. ",
                                        html.Span("Closed Loop Geothermal Analysis Modeling and Simulation Using Idaho National Laboratory RELAP5-3D-FALCON Coupled Codes.", style={"fontWeight": "bold"}),
                                        " Proceedings of the 48th Workshop on Geothermal Reservoir Engineering, Stanford University, February 6-8, 2023.",
                                    ],
                                    pdf_path="pdfs/RELAP5-3D-Falcon.pdf",
                                ),
                                create_resource_card(
                                    "resource-card-8",
                                    "Stanford Geothermal Workshop, 2023",
                                    [
                                        "Koenraad Beckers, Yaroslav Vasyliv, Gabriela A. Bran-Anleu, Mario Martinez, Chad Augustine, Mark White, and the Closed-Loop Geothermal Working Group. ",
                                        html.Span("Tabulated Database of Closed-Loop Geothermal Systems Performance for Cloud-Based Technical and Economic Modeling of Heat Production and Electricity Generation.", style={"fontWeight": "bold"}),
                                        " Proceedings of the 48th Workshop on Geothermal Reservoir Engineering, Stanford University, February 6-8, 2023.",
                                    ],
                                    pdf_path="pdfs/Tabulated Database of Closed-Loop.pdf",
                                ),
                                create_resource_card(
                                    "resource-card-9",
                                    "Stanford Geothermal Workshop, 2023",
                                    [
                                        "Mark White, Mario Martinez, Yaroslav Vasyliv, Koenraad Beckers, Gabriela A. Bran-Anleu, Carlo Parisi, Paolo Balestra, Roland Horne, Chad Augustine, Laura Pauley, Giorgia Bettini, Theron Marshall, and the Closed Loop Geothermal Working Group. ",
                                        html.Span("Closed-Loop Geothermal Working Group Study – Understanding Thermal Performance and Economic Forecasts via Numerical Simulation.", style={"fontWeight": "bold"}),
                                        " Proceedings of the 48th Workshop on Geothermal Reservoir Engineering, Stanford University, February 6–8, 2023.",
                                    ],
                                    pdf_path="pdfs/Closed-loop Geothermal Working Group Study - Understanding Thermal Performance.pdf",
                                ),
                                create_resource_card(
                                    "resource-card-12",
                                    "Geothermal Rising Conference, 2021",
                                    [
                                        "Mark White, Mario Martinez, Yaroslav Vasyliv, Gabriela A. Bran-Anleu, Carlo Parisi, Paolo Balestra, Roland Horne, Chad Augustine, Laura Pauley, Doug Hollett, Giorgia Bettini, Theron Marshall, and the Closed Loop Geothermal Working Group. ",
                                        html.Span("Thermal and Mechanical Energy Performance Analysis of Closed-Loop Systems in Hot-Dry-Rock and Hot-Wet-Rock Reservoirs.", style={"fontWeight": "bold"}),
                                        " Proceedings of the Geothermal Rising Conference, Vol. 45, 2021.",
                                    ],
                                    pdf_path="pdfs/Thermal and Mechanical Energy Performance Analysis of.pdf",
                                ),
                                create_resource_card(
                                    "resource-card-11",
                                    "Geothermal Rising Conference, 2021",
                                    [
                                        "Yaroslav Vasyliv, Gabriela A. Bran-Anleu, Alec Kucala, Sam Subia, and Mario J. Martinez. ",
                                        html.Span("Analysis and Optimization of a Closed Loop Geothermal System in Hot Rock Reservoirs.", style={"fontWeight": "bold"}),
                                        " Proceedings of the Geothermal Rising Conference, Vol. 45, 2021.",
                                    ],
                                    pdf_path="pdfs/Analysis and Optimization.pdf",
                                ),
                                create_resource_card(
                                    "resource-card-10",
                                    "Geothermal Rising Conference, 2021",
                                    [
                                        "Carlo Parisi, Paolo Balestra, and Theron D. Marshall. ",
                                        html.Span("Geothermal Analysis Modeling and Simulation Using Idaho National Laboratory RELAP5-3D-PRONGHORN Coupled Codes.", style={"fontWeight": "bold"}),
                                        " Proceedings of the Geothermal Rising Conference, Vol. 45, 2021.",
                                    ],
                                    pdf_path="pdfs/Geothermal analysis modeling and simulation _ Geothermal rising 2021.pdf",
                                ),
                            ],
                        ),
                    ],
                ),
            ],
        ),
    ],
)


energy_time_tab = dcc.Tab(
    label="Subsurface Results",
    id="subsurface-tab",
    value="energy-time-tab",
    selected_className="active_tabs",
    children=[
        html.Div(className="extra-space"),
        graph_guidance_card(
            btnID="collapse-button2",
            cardID="collapse2",
            dropdown_children=html.P(dropdown_text1),
        ),
        html.Div(id="error_block_div1"),
        html.Br(),
        html.Div(
            id="graphics-parent",
            children=[
                html.Div(
                    id="graphics-container",
                    children=[
                        dcc.Graph(id="geothermal_time_plots", config=plotly_config),
                        dbc.RadioItems(
                            options=[
                                {"label": "Auto Scale", "value": 1},
                                {"label": "Full Scale", "value": 2},
                            ],
                            value=1,
                            id="radio-graphic-control3",
                        ),
                    ],
                )
            ],
        ),
    ],
)


energy_tab = dcc.Tab(
    label="Subsurface Contours",
    id="contour-tab",
    value="energy-tab",
    selected_className="active_tabs",
    children=[
        html.Div(className="extra-space"),
        # html.Hr(),
        html.Div(
            id="dropdown-card5",
            children=[
                graph_guidance_card(
                    btnID="collapse-button",
                    cardID="collapse",
                    dropdown_children=html.P(dropdown_text2),
                ),
                html.Div(id="error_block_div2"),
                html.P("Select Parameter", id="select-text"),
                dcc.Dropdown(
                    id="param-select",
                    options=[{"label": i, "value": i} for i in param_list],
                    value="Horizontal Extent (m)",
                    clearable=False,
                    searchable=False,
                ),
            ],
        ),
        dcc.Graph(id="geothermal_plots", config=plotly_config),
    ],
)

economics_time_tab = dcc.Tab(
    label="Economic Results",
    id="econ-tab",
    value="economics-time-tab",
    selected_className="active_tabs",
    children=[
        html.Div(className="extra-space"),
        # html.Hr(),
        graph_guidance_card(
            btnID="collapse-button4",
            cardID="collapse4",
            dropdown_children=html.Div(
                children=[
                    html.P(dropdown_text1.replace("5 options", "7 options")),
                    html.P(dropdown_econ_text1),
                    html.Img(
                        id="formulas-img",
                        src=app.get_asset_url("lcoe-lcoh-formulas.png"),
                    ),
                    dcc.Markdown(dropdown_econ_markdown_text1, mathjax=True),
                    html.P(dropdown_econ_text2),
                ]
            ),
        ),
        html.Div(id="warning_block_div3"),
        html.Div(id="error_block_div3"),
        html.Br(),
        html.Div(
            id="graphics-parent-econ",
            children=[
                html.Div(
                    id="graphics-container-econ",
                    children=[
                        dcc.Graph(id="econ_plots", config={**plotly_config, "responsive": True}),
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
                    ],
                )
            ],
        ),
    ],
)

summary_tab = dcc.Tab(
    label="Summary",
    id="summary-tab",
    value="summary-tab",
    selected_className="active_tabs",
    children=[
        html.Hr(className="tab-hr"),
        html.Button("Download Results", id="btn_xlsx"),
        dcc.Download(id="download-dataframe-xlsx"),
        html.Br(),
        html.Br(),
        dcc.Graph(id="table", config=plotly_config),
    ],
)

# -----------------------------------------------------------------------------
# Define dash app layout here.
# -----------------------------------------------------------------------------

app.layout = html.Div(
    id="app-container",
    children=[
        dcc.Store(id="econ-memory"),
        dcc.Store(id="econ-results"),
        dcc.Store(id="econ-errors"),
        dcc.Store(id="thermal-memory"),
        dcc.Store(id="thermal-results-mass"),
        dcc.Store(id="thermal-results-time"),
        dcc.Store(id="thermal-results-errors"),
        dcc.Store(id="thermal-contours-errors"),
        dcc.Store(id="summary-memory"),
        dcc.Store(id="TandP-data"),
        dcc.Store(id="slider-values-store", data={}),  # Store slider values per model
        dcc.Store(id="sbt-params-store", data={}),  # Store SBT-specific parameters for Excel export
        dcc.Store(
            id="fluid-selection-store",
            data={"preferred": "All", "last_specific": "H2O"},
        ),  # Remember user fluid preferences across tabs
        dcc.Store(id="plot-inputs-cache", data={}),
        dcc.Store(id="plot-params-store"),
        dcc.Store(id="plot-run-trigger"),
        dcc.Store(id="see-all-params-state", data={"expanded": False}),
        dcc.Store(id="calculation-request-id", data=0),
        dcc.Store(id="clipboard-init", data=0),
        dcc.Store(id="econ-resize-ping", data=0),
        dcc.Store(id="mass-mode-select", data="Constant"),
        dcc.Store(id="temp-mode-select", data="Constant"),
        # Left column
        html.Div(
            id="left-column",
            className="four columns",
            children=[description_card(), generate_control_card()],
        ),
        # Right column
        html.Div(
            id="right-column",
            className="eight columns",
            children=[
                dcc.Loading(
                    id="main-loader",
                    type="circle",
                    parent_className="loader-wrapper",
                    style={"minHeight": "100px"},
                    children=[
                        html.Div(
                            id="geothermal_card",
                            children=[
                                dcc.Tabs(
                                    id="tabs",
                                    value="about-tab",
                                    children=[
                                        about_tab,
                                        energy_time_tab,
                                        energy_tab,
                                        economics_time_tab,
                                        summary_tab,
                                    ],
                                )
                            ],
                        ),
                    ],
                ),
            ],
        ),
        # Information Modal
        create_info_modal(),
        # Clipboard toast notification
        dbc.Toast(
            id="clipboard-toast",
            header="Copied!",
            children="URL copied to clipboard",
            is_open=False,
            dismissable=True,
            duration=2000,
            style={"position": "fixed", "top": 66, "right": 10, "width": 350, "zIndex": 9999},
        ),
        # Hidden button to trigger toast
        html.Button(id="clipboard-toast-trigger-btn", style={"display": "none"}, n_clicks=0),
    ],
)

# Register info popup callbacks
register_info_modal_callbacks(app)


# -----------------------------------------------------------------------------
# Define dash app callbacks begin here.
# -----------------------------------------------------------------------------

# Resource cards expand/collapse callbacks
@app.callback(
    Output({"type": "resource-collapse", "index": MATCH}, "is_open"),
    [Input({"type": "resource-button", "index": MATCH}, "n_clicks")],
    [State({"type": "resource-collapse", "index": MATCH}, "is_open")],
    prevent_initial_call=True,
)
def toggle_resource_card(n_clicks, is_open):
    if is_print:
        print("toggle_resource_card")
    """Toggle resource card expand/collapse."""
    if n_clicks:
        return not is_open
    return is_open


# Combined client-side callback for clipboard and card click handling
app.clientside_callback(
    """
    function(data) {
        setTimeout(function() {
            try {
                // Show toast when clipboard is clicked
                var clipboardBtn = document.querySelector('#clipboard-wrapper button');
                if (clipboardBtn && !clipboardBtn.dataset.toastListenerAdded) {
                    clipboardBtn.dataset.toastListenerAdded = 'true';
                    clipboardBtn.addEventListener('click', function() {
                        setTimeout(function() {
                            var triggerBtn = document.getElementById('clipboard-toast-trigger-btn');
                            if (triggerBtn) {
                                triggerBtn.click();
                            }
                        }, 100);
                    });
                }
                
                // Handle div-based card headers (for cards with clipboard)
                var divHeaders = document.querySelectorAll('.resource-card-header[data-href]');
                divHeaders.forEach(function(header) {
                    if (!header.dataset.listenerAdded) {
                        header.dataset.listenerAdded = 'true';
                        header.addEventListener('click', function(e) {
                            var clickedElement = e.target;
                            var clipboardWrapper = document.getElementById('clipboard-wrapper');
                            if (clipboardWrapper && (clickedElement === clipboardWrapper || clipboardWrapper.contains(clickedElement))) {
                                return;
                            }
                            var href = header.getAttribute('data-href');
                            var target = header.getAttribute('data-target');
                            if (href) {
                                if (target === '_blank') {
                                    window.open(href, '_blank', 'noopener,noreferrer');
                                } else {
                                    window.location.href = href;
                                }
                            }
                        });
                    }
                });
                
                // Prevent clipboard clicks from triggering card links
                var wrapper = document.getElementById('clipboard-wrapper');
                if (wrapper && !wrapper.dataset.propagationListenerAdded) {
                    wrapper.dataset.propagationListenerAdded = 'true';
                    wrapper.addEventListener('click', function(e) {
                        e.stopPropagation();
                        e.stopImmediatePropagation();
                    }, true);
                    
                    var clipboardBtn2 = wrapper.querySelector('button');
                    if (clipboardBtn2) {
                        clipboardBtn2.addEventListener('click', function(e) {
                            e.stopPropagation();
                            e.stopImmediatePropagation();
                        }, true);
                    }
                }
            } catch(e) {
                console.log('Client-side callback error:', e);
            }
        }, 1000);
        return window.dash_clientside.no_update;
    }
    """,
    Output('clipboard-init', 'data'),
    Input('clipboard-init', 'data'),
)


@app.callback(
    Output("clipboard-toast", "is_open"),
    Input("clipboard-toast-trigger-btn", "n_clicks"),
    State("clipboard-toast", "is_open"),
    prevent_initial_call=True,
)
def show_clipboard_toast(n_clicks, is_open):
    if is_print:
        print("show_clipboard_toast")
    """Show toast notification when clipboard is clicked."""
    if n_clicks:
        return True
    return False


@app.callback(
    Output(component_id="download-dataframe-xlsx", component_property="data"),
    [
        Input(component_id="btn_xlsx", component_property="n_clicks"),
        Input(component_id="summary-memory", component_property="data"),
        Input(component_id="thermal-results-mass", component_property="data"),
        Input(component_id="thermal-results-time", component_property="data"),
        Input(component_id="econ-results", component_property="data"),
    ],
    [
        State(component_id="model-select", component_property="value"),
        State(component_id="sbt-params-store", component_property="data"),
    ],
    prevent_initial_call=True,
)
def generate_summary(n_clicks, df_tbl, df_mass_flow_rate, df_time, df_econ, model, sbt_params):
    if is_print:
        print("generate_summary")
    if "btn_xlsx" == ctx.triggered_id:
        write_excelsheet(
            df_summary=df_tbl,
            df_subsurf_res_mass=df_mass_flow_rate,
            df_subsurf_res_time=df_time,
            df_econ=df_econ,
            geoCLUSTER_results_pathname=geoCLUSTER_results_pathname,
            model=model,
            sbt_params=sbt_params or {},
        )
        return dcc.send_file(geoCLUSTER_results_pathname)
    else:
        raise PreventUpdate


@app.callback(
    Output(component_id="collapse", component_property="is_open"),
    [Input(component_id="collapse-button", component_property="n_clicks")],
    [State(component_id="collapse", component_property="is_open")],
)
def toggle_collapse(n, is_open):
    if is_print:
        print("toggle_collapse 1")
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
    if is_print:
        print("toggle_collapse 2")
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
    if is_print:
        print("toggle_collapse 3")
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
    if is_print:
        print("toggle_collapse 4")
    # ----------------------------------------------
    # Economic results graph guidance.
    # ----------------------------------------------

    if n:
        return not is_open
    return is_open


@app.callback(
    [
        Output(component_id="model-select", component_property="options"),
        Output(component_id="model-select", component_property="value", allow_duplicate=True),
    ],
    [
        Input(component_id="fluid-select", component_property="value"),
        Input(component_id="model-select", component_property="value"),
    ],
    prevent_initial_call=True,
)
def update_model_based_on_fluid_and_selection(fluid, model_selection):
    if is_print:
        print("update_model_based_on_fluid_and_selection")
    """
    Update model-select options and value based on fluid selection and model selection.
    - Database (HDF5): Always stays as HDF5 regardless of fluid
    - Simulator: H2O → SBT V1.0, sCO2/All → SBT V2.0
    """
    
    # Determine which input triggered the callback
    trigger_id = ctx.triggered_id if ctx.triggered_id else None
    
    # Get current fluid value
    current_fluid = fluid
    
    # Get updated options based on current fluid
    new_options = get_model_options(current_fluid)
    
    if model_selection == "HDF5":
        new_value = "HDF5"
    elif model_selection in ("SBT V1.0", "SBT V2.0"):
        if current_fluid == "H2O":
            new_value = "SBT V1.0"
        else:
            new_value = "SBT V2.0"
    else:
        new_value = model_selection
    
    return new_options, new_value


@app.callback(
    [
        Output(component_id="mesh-div", component_property="style", allow_duplicate=True),
        Output(component_id="pipe-roughness-container", component_property="style", allow_duplicate=True),
        Output(component_id="hyperparam5-container", component_property="style", allow_duplicate=True),
        Output(component_id="lat-flow-container", component_property="style", allow_duplicate=True),
        Output(component_id="lat-allo-container", component_property="style", allow_duplicate=True),
        Output(component_id="component-performance-div", component_property="style", allow_duplicate=True),
        Output(component_id="see-all-params-text", component_property="children"),
        Output(component_id="see-all-params-chevron", component_property="children"),
    ],
    [
        Input(component_id="see-all-params-button", component_property="n_clicks"),
        Input(component_id="model-select", component_property="value"),
        Input(component_id="case-select", component_property="value"),
        Input(component_id="tabs", component_property="value"),
    ],
    [
        State(component_id="see-all-params-text", component_property="children"),
        State(component_id="see-all-params-state", component_property="data"),
    ],
    prevent_initial_call='initial_duplicate',
)
def toggle_see_all_params(n_clicks, model, case, tab, current_button_text, params_state):
    if is_print:
        print("toggle_see_all_params")
    callback_ctx = ctx
    style_none = {"display": "none"}
    style_block = {"display": "block"}
    
    # Component Performance (Insulation Thermal Conductivity) only relevant for coaxial
    component_perf_style = style_none if case == "utube" else style_block
    
    # Initialize params_state if not present
    if params_state is None:
        params_state = {"expanded": False}
    
    # Get the trigger ID
    if not callback_ctx.triggered:
        # Initial state - check cache first, otherwise use defaults
        is_expanded = params_state.get("expanded", False)
        if is_expanded:
            if model == "HDF5":
                return style_none, style_none, style_none, style_none, style_none, style_none, "See less parameters", " ▲"
            elif model == "SBT V1.0":
                return style_block, style_block, style_block, style_block, style_none, component_perf_style, "See less parameters", " ▲"
            else:  # SBT V2.0
                return style_block, style_block, style_block, style_block, style_none, component_perf_style, "See less parameters", " ▲"
        else:
            # Use default collapsed state
            if model == "HDF5":
                return style_none, style_none, style_none, style_none, style_none, style_none, "See all parameters", " ▼"
            elif model == "SBT V1.0":
                return style_none, style_none, style_none, style_none, style_none, style_none, "See all parameters", " ▼"
            else:  # SBT V2.0
                return style_none, style_block, style_block, style_none, style_none, style_none, "See all parameters", " ▼"
    
    trigger_id = callback_ctx.triggered[0]["prop_id"].split(".")[0]
    
    # Handle tab change - restore cached state
    if trigger_id == "tabs":
        is_expanded = params_state.get("expanded", False)
        if is_expanded:
            if model == "HDF5":
                return style_none, style_none, style_none, style_none, style_none, style_none, "See less parameters", " ▲"
            elif model == "SBT V1.0":
                return style_block, style_block, style_block, style_block, style_none, component_perf_style, "See less parameters", " ▲"
            else:  # SBT V2.0
                return style_block, style_block, style_block, style_block, style_none, component_perf_style, "See less parameters", " ▲"
        else:
            # Use default collapsed state
            if model == "HDF5":
                return style_none, style_none, style_none, style_none, style_none, style_none, "See all parameters", " ▼"
            elif model == "SBT V1.0":
                return style_none, style_none, style_none, style_none, style_none, style_none, "See all parameters", " ▼"
            else:  # SBT V2.0
                return style_none, style_block, style_block, style_none, style_none, style_none, "See all parameters", " ▼"
    
    # Handle model or case change - restore cached state instead of resetting
    if trigger_id == "model-select" or trigger_id == "case-select":
        is_expanded = params_state.get("expanded", False)
        if is_expanded:
            if model == "HDF5":
                return style_none, style_none, style_none, style_none, style_none, style_none, "See less parameters", " ▲"
            elif model == "SBT V1.0":
                return style_block, style_block, style_block, style_block, style_none, component_perf_style, "See less parameters", " ▲"
            else:  # SBT V2.0
                return style_block, style_block, style_block, style_block, style_none, component_perf_style, "See less parameters", " ▲"
        else:
            # Use default collapsed state
            if model == "HDF5":
                return style_none, style_none, style_none, style_none, style_none, style_none, "See all parameters", " ▼"
            elif model == "SBT V1.0":
                return style_none, style_none, style_none, style_none, style_none, style_none, "See all parameters", " ▼"
            else:  # SBT V2.0
                return style_none, style_block, style_block, style_none, style_none, style_none, "See all parameters", " ▼"
    
    # Handle button click
    if trigger_id == "see-all-params-button" and n_clicks is not None and n_clicks > 0:
        # Check if currently showing "See less" (meaning params are visible)
        is_visible = current_button_text == "See less parameters" if current_button_text else False
        if is_visible:
            # Hide all parameters
            if model == "HDF5":
                return style_none, style_none, style_none, style_none, style_none, style_none, "See all parameters", " ▼"
            elif model == "SBT V1.0":
                # Hide all parameters for SBT V1.0
                return style_none, style_none, style_none, style_none, style_none, style_none, "See all parameters", " ▼"
            else:  # SBT V2.0
                # Hide: mesh, lat-flow (Lateral Flow Multiplier), lat-allo, component-performance
                # Keep pipe-roughness and hyperparam5 visible (they're always visible for SBT V2.0)
                return style_none, style_block, style_block, style_none, style_none, style_none, "See all parameters", " ▼"
        else:
            # Show all parameters - cache will be updated by separate callback
            if model == "HDF5":
                # Component Performance unavailable for HDF5 (database mode)
                return style_none, style_none, style_none, style_none, style_none, style_none, "See less parameters", " ▲"
            elif model == "SBT V1.0":
                # For SBT V1.0, show the SBT-specific parameters including component-performance (with Insulation Thermal Conductivity)
                # lat-allo (Lateral Flow Allocation) is unavailable for simulator
                # Component Performance only for coaxial
                return style_block, style_block, style_block, style_block, style_none, component_perf_style, "See less parameters", " ▲"
            else:  # SBT V2.0
                # For SBT V2.0, pipe-roughness and hyperparam5 are already visible, so keep them visible
                # Show: mesh, lat-flow (Lateral Flow Multiplier), component-performance (lat-allo is unavailable for simulator)
                # Component Performance only for coaxial
                return style_block, style_block, style_block, style_block, style_none, component_perf_style, "See less parameters", " ▲"
    
    # Fallback - check cache for state
    is_expanded = params_state.get("expanded", False)
    if is_expanded:
        if model == "HDF5":
            return style_none, style_none, style_none, style_none, style_none, style_none, "See less parameters", " ▲"
        elif model == "SBT V1.0":
            return style_block, style_block, style_block, style_block, style_none, component_perf_style, "See less parameters", " ▲"
        else:  # SBT V2.0
            return style_block, style_block, style_block, style_block, style_none, component_perf_style, "See less parameters", " ▲"
    else:
        if model == "HDF5":
            return style_none, style_none, style_none, style_none, style_none, style_none, "See all parameters", " ▼"
        elif model == "SBT V1.0":
            return style_none, style_none, style_none, style_none, style_none, style_none, "See all parameters", " ▼"
        else:  # SBT V2.0
            return style_none, style_block, style_block, style_none, style_none, style_none, "See all parameters", " ▼"


@app.callback(
    Output(component_id="see-all-params-state", component_property="data"),
    [
        Input(component_id="see-all-params-button", component_property="n_clicks"),
    ],
    [
        State(component_id="see-all-params-text", component_property="children"),
        State(component_id="see-all-params-state", component_property="data"),
    ],
    prevent_initial_call=True,
)
def update_see_all_params_cache(n_clicks, current_button_text, current_state):
    if is_print:
        print("update_see_all_params_cache")
    """Update cache when 'See all parameters' button is clicked"""
    if n_clicks is None or n_clicks == 0:
        raise PreventUpdate
    
    if current_state is None:
        current_state = {"expanded": False}
    
    is_expanded = current_button_text == "See less parameters" if current_button_text else False
    new_state = {"expanded": not is_expanded}
    
    return new_state






@app.callback(
    [
        Output(component_id="scenario1-div", component_property="style"),
        Output(component_id="scenario3-div", component_property="style"),
        Output(component_id="hr-break1", component_property="style"),
    ],
    [Input(component_id="model-select", component_property="value")],
    prevent_initial_call=True,
)
def update_scenario_buttons_visibility(selected_model):
    if is_print:
        print("update_scenario_buttons_visibility")
    if selected_model == "HDF5":
        return {"display": "block"}, {"display": "block"}, {"display": "block"}

    elif selected_model == "SBT V1.0" or selected_model == "SBT V2.0":
        # TODO: update tabs styline
        return {"display": "none"}, {"display": "none"}, {"display": "none"}

    else:
        raise PreventUpdate


@app.callback(
    Output(component_id="graphics-parent", component_property="children"),
    Input(component_id="model-select", component_property="value"),
    prevent_initial_call=True,
)
def update_graphics_container(selected_model):
    if is_print:
        print("update_graphics_container")
    if selected_model == "HDF5":
        return html.Div(
            id="graphics-container",
            children=[
                dcc.Graph(id="geothermal_time_plots", config=plotly_config),
                dbc.RadioItems(
                    options=[
                        {"label": "Auto Scale", "value": 1},
                        {"label": "Full Scale", "value": 2},
                    ],
                    value=1,
                    id="radio-graphic-control3",
                ),
            ],
        )

    elif selected_model == "SBT V1.0" or selected_model == "SBT V2.0":
        return html.Div(
            id="graphics-container",
            className="simulator-mode",
            children=[
                dcc.Graph(id="geothermal_time_plots", config=plotly_config),
                dbc.RadioItems(
                    options=[
                        {"label": "Auto Scale", "value": 1},
                        {"label": "Full Scale", "value": 2},
                    ],
                    value=1,
                    id="radio-graphic-control3",
                ),
            ],
        )
    else:
        raise PreventUpdate


@app.callback(
    Output(component_id="contour-tab", component_property="style"),
    [Input(component_id="model-select", component_property="value")],
    prevent_initial_call=True,
)
def update_contour_tab_visibility(selected_model):
    if is_print:
        print("update_contour_tab_visibility")
    if selected_model == "HDF5":
        return {"display": "block"}
    elif selected_model == "SBT V1.0" or selected_model == "SBT V2.0":
        return {"display": "none"}
    else:
        return {"display": "none"}


@app.callback(
    Output(component_id="tabs", component_property="value", allow_duplicate=True),
    [
        Input(component_id="model-select", component_property="value"),
    ],
    [
        State(component_id="tabs", component_property="value"),
    ],
    prevent_initial_call=True,
)
def switch_tab_on_model_change(model, current_tab):
    if is_print:
        print("switch_tab_on_model_change")
    """
    Automatically switch from subsurface contours tab (energy-tab) to subsurface results tab (energy-time-tab)
    when switching from Database to Simulator.
    """
    # Only trigger if switching from Database to Simulator
    if model in ("SBT V1.0", "SBT V2.0") and current_tab == "energy-tab":
        return "energy-time-tab"
    else:
        raise PreventUpdate


@app.callback(
    Output(component_id="end-use-select", component_property="value"),
    [
        Input(component_id="tabs", component_property="value"),
        Input(component_id="btn-nclicks-1", component_property="n_clicks"),
        Input(component_id="btn-nclicks-3", component_property="n_clicks"),
        Input(component_id="end-use-select", component_property="value"),
    ],
    prevent_initial_call=True,
)
def flip_to_tab(tab, btn1, btn3, end_use):
    if is_print:
        print("flip_to_tab")
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
    [
        Output(component_id="fluid-select", component_property="options"),
        Output(component_id="fluid-select", component_property="value"),
        Output(component_id="interpolation-select", component_property="options"),
        Output(component_id="interpolation-select", component_property="value"),
        Output(component_id="fluid-selection-store", component_property="data"),
    ],
    [
        Input(component_id="tabs", component_property="value"),
        Input(component_id="fluid-select", component_property="value"),
        Input(component_id="model-select", component_property="value"),
        Input(component_id="case-select", component_property="value"),
        Input(component_id="btn-nclicks-1", component_property="n_clicks"),
        Input(component_id="btn-nclicks-3", component_property="n_clicks"),
    ],
    [
        State(component_id="fluid-selection-store", component_property="data"),
    ],
    prevent_initial_call=True,
)
def change_dropdown(at, fluid, model, case, btn1, btn3, fluid_store):
    if is_print:
        print("change_dropdown")
    if model not in ("HDF5", "SBT V2.0", "SBT V1.0"):
        raise PreventUpdate

    trigger_id = ctx.triggered_id if ctx.triggered_id else None

    if fluid_store is None:
        fluid_store = {"preferred": "All", "last_specific": "H2O"}

    if trigger_id in ["btn-nclicks-1", "btn-nclicks-3"]:
        if fluid == "All" or fluid is None:
            preferred = fluid_store.get("preferred", "All")
            last_specific = fluid_store.get("last_specific", "H2O") or "H2O"
            new_store = {"preferred": preferred, "last_specific": last_specific}
            fluid_list = ["All", "H2O", "sCO2"]
            interpol_list = ["True"]
            return (
                [{"label": i, "value": i} for i in fluid_list],
                "H2O",
                [{"label": i, "value": i} for i in interpol_list],
                interpol_list[0],
                new_store,
            )

    if at == "energy-tab":
        last_specific = fluid_store.get("last_specific", "H2O") or "H2O"
        selected_fluid = fluid if fluid in ["H2O", "sCO2"] else last_specific
        if selected_fluid not in ["H2O", "sCO2"]:
            selected_fluid = "H2O"
        fluid_options = [{"label": "H2O", "value": "H2O"}, {"label": "sCO2", "value": "sCO2"}]
        interpol_options = [{"label": "True", "value": "True"}]
        new_store = {"preferred": fluid_store.get("preferred", "All"), "last_specific": last_specific}
        return (fluid_options, selected_fluid, interpol_options, "True", new_store)

    if trigger_id == "model-select" and at != "energy-tab":
        preferred = fluid_store.get("preferred", "All")
        last_specific = fluid_store.get("last_specific", "H2O") or "H2O"
        fluid_list = ["All", "H2O", "sCO2"]
        if preferred in fluid_list:
            selected_fluid = preferred
        elif last_specific in fluid_list:
            selected_fluid = last_specific
        else:
            selected_fluid = fluid_list[0] if fluid_list else "H2O"
        interpol_list = ["True"]
        new_store = {"preferred": preferred, "last_specific": last_specific}
        return (
            [{"label": i, "value": i} for i in fluid_list],
            selected_fluid,
            [{"label": i, "value": i} for i in interpol_list],
            interpol_list[0],
            new_store,
        )

    preferred = fluid_store.get("preferred", "All")
    last_specific = fluid_store.get("last_specific", "H2O") or "H2O"
    new_store = {"preferred": preferred, "last_specific": last_specific}

    if trigger_id == "fluid-select" and fluid is not None:
        new_store["preferred"] = fluid
        if fluid != "All":
            new_store["last_specific"] = fluid
        preferred = new_store["preferred"]
        last_specific = new_store["last_specific"]

    fluid_list = ["All", "H2O", "sCO2"]
    interpol_list = ["True"]
    selected_fluid = fluid if fluid is not None else preferred

    if preferred in fluid_list:
        selected_fluid = preferred

    fluid_options = [{"label": i, "value": i} for i in fluid_list]
    interpol_options = [{"label": i, "value": i} for i in interpol_list]

    if trigger_id == "tabs":
        if selected_fluid == fluid:
            current_fluid_options = [{"label": i, "value": i} for i in (["H2O", "sCO2"] if at == "energy-tab" else ["All", "H2O", "sCO2"])]
            if fluid_options == current_fluid_options:
                raise PreventUpdate

    return (
        fluid_options,
        selected_fluid,
        interpol_options,
        interpol_list[0],
        new_store,
    )


@app.callback(
    Output(component_id="tabs", component_property="value"),
    [
        Input(component_id="tabs", component_property="value"),
        Input(component_id="btn-nclicks-1", component_property="n_clicks"),
        # Input(component_id='btn-nclicks-2', component_property='n_clicks'),
        Input(component_id="btn-nclicks-3", component_property="n_clicks"),
    ],
)
def flip_to_tab(tab, btn1, btn3):
    if is_print:
        print("flip_to_tab")
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
    [
        Output(component_id="mdot-select", component_property="value"),
        Output(component_id="L2-select", component_property="value"),
        Output(component_id="L1-select", component_property="value"),
        Output(component_id="grad-select", component_property="value"),
        Output(component_id="diameter-select", component_property="value"),
        Output(component_id="Tinj-select", component_property="value"),
        Output(component_id="k-select", component_property="value"),
        Output(component_id="drillcost-select", component_property="value"),
        Output(component_id="discount-rate-select", component_property="value"),
        Output(component_id="lifetime-select", component_property="value"),
        Output(component_id="kwt-select", component_property="value"),
        Output(component_id="kwe-select", component_property="value"),
        Output(component_id="precool-select", component_property="value"),
        Output(component_id="turb-pout-select", component_property="value"),
        Output(component_id="Tsurf-select", component_property="value"),
        Output(component_id="c-select", component_property="value"),
        Output(component_id="rho-select", component_property="value"),
        Output(component_id="n-laterals-select", component_property="value"),
        Output(component_id="lateral-flow-select", component_property="value"),
        Output(component_id="lateral-multiplier-select", component_property="value"),
    ],
    [
        Input(component_id="btn-nclicks-1", component_property="n_clicks"),
        Input(component_id="btn-nclicks-3", component_property="n_clicks"),
        Input(component_id="tabs", component_property="value"),
        Input(component_id="case-select", component_property="value"),
        Input(component_id="fluid-select", component_property="value"),
        Input(component_id="end-use-select", component_property="value"),
        Input(component_id="model-select", component_property="value"),
    ],
    prevent_initial_call=True,
)
def update_slider_with_btn(btn1, btn3, at, case, fluid, end_use, model):
    if is_print:
        print("update_slider_with_btn")
    N_OUTPUTS = 20
    trigger = ctx.triggered_id

    if trigger not in ("btn-nclicks-1", "btn-nclicks-3"):
        raise PreventUpdate

    if model != "HDF5":
        raise PreventUpdate

    case = (case or "").strip()
    fluid = (fluid or "").strip()
    end_use = (end_use or "").strip()

    default_tail = (
        25,    # Tsurf
        790,   # c
        2750,  # rho
        1,     # n-laterals
        1,     # lateral-flow
        1,     # lateral-multiplier
    )

    if trigger == "btn-nclicks-1":
        head = (
            24,     # mdot
            10000,  # L2
            2500,   # L1
            0.050,  # grad
            0.30,   # diameter
            40,     # Tinj
            3.5,    # k
            1000,   # drillcost
            7.0,    # discount rate
            40,     # lifetime
            100,    # kwt
            3000,   # kwe
            13,     # precool
            80,     # turb-pout
        )
        return head + default_tail

    presets = {
        ("coaxial", "H2O", "Electricity"): (39.2, 20000, 5000, 0.070, 0.444, 60, 4.5, 1000, 7.0, 40, 100, 3000, 13, 80),
        ("coaxial", "H2O", "Heating"):     (73.4, 13000, 5000, 0.070, 0.444, 30, 4.5, 1000, 7.0, 40, 100, 3000, 13, 80),
        ("coaxial", "H2O", "All"):          (73.4, 13000, 5000, 0.070, 0.444, 30, 4.5, 1000, 7.0, 40, 100, 3000, 13, 80),
        ("utube",   "H2O", "Electricity"):  (43,   20000, 5000, 0.070, 0.44,  60, 4.5, 1000, 7.0, 40, 100, 3000, 13, 80),
        ("utube",   "H2O", "Heating"):     (100,  20000, 5000, 0.070, 0.44,  30, 4.5, 1000, 7.0, 40, 100, 3000, 13, 80),
        ("utube",   "H2O", "All"):         (100,  20000, 5000, 0.070, 0.44,  30, 4.5, 1000, 7.0, 40, 100, 3000, 13, 80),
        ("coaxial", "sCO2", "Electricity"): (69.6, 13000, 5000, 0.070, 0.44,  60, 4.5, 1000, 7.0, 40, 100, 3000, 13, 80),
        ("coaxial", "sCO2", "Heating"):     (69.6, 6000,  5000, 0.070, 0.44,  45, 4.5, 1000, 7.0, 40, 100, 3000, 13, 80),
        ("coaxial", "sCO2", "All"):         (69.6, 6000,  5000, 0.070, 0.44,  45, 4.5, 1000, 7.0, 40, 100, 3000, 13, 80),
        ("utube",   "sCO2", "Electricity"): (100,  20000, 5000, 0.070, 0.44,  60, 4.5, 1000, 7.0, 40, 100, 3000, 13, 80),
        ("utube",   "sCO2", "Heating"):     (100,  11000, 5000, 0.070, 0.44,  45, 4.5, 1000, 7.0, 40, 100, 3000, 13, 80),
        ("utube",   "sCO2", "All"):         (100,  11000, 5000, 0.070, 0.44,  45, 4.5, 1000, 7.0, 40, 100, 3000, 13, 80),
    }

    fluid_key = "H2O" if fluid in ("All", "H2O") else fluid

    head = presets.get((case, fluid_key, end_use))
    if head is None:
        return (no_update,) * N_OUTPUTS

    return head + default_tail


@app.callback(
    [
        Output(component_id="model-params-div", component_property="style"),
        Output(component_id="Tsurf-select-div", component_property="style"),
        Output(component_id="c-select-div", component_property="style"),
        Output(component_id="rho-select-div", component_property="style"),
        Output(component_id="diameter-vertical-select-div", component_property="style"),
        Output(component_id="diameter-lateral-select-div", component_property="style"),
        Output(
            component_id="diameter-select-div", component_property="style"
        ),  # show HDF5, hide SBT
        Output(
            component_id="fluid-mode-div", component_property="style"
        ),  # hide HDF5 and v1, show v2
        Output(component_id="num-lat-div", component_property="style"),
        Output(component_id="lat-allo-container", component_property="style"),
        # Output(component_id='lateral-flow-select-div', component_property='style'),
        Output(component_id="pipe-roughness-container", component_property="style"),
        Output(component_id="component-performance-div", component_property="style"),
        Output(component_id="see-all-params-button-container", component_property="style"),
        Output(component_id="coaxial-flow-type-wrapper", component_property="style"),
    ],
    [
        Input(component_id="model-select", component_property="value"),
        Input(component_id="case-select", component_property="value"),
    ],
    prevent_initial_call='initial_duplicate',
)
def show_model_params(model, case):
    if is_print:
        print("show_model_params")
    b = {"display": "block"}
    n = {"display": "none"}
    see_all_button_style = {"textAlign": "center", "marginTop": "20px", "marginBottom": "10px", "display": "block"}

    # Component Performance (Insulation Thermal Conductivity) only relevant for coaxial
    component_perf_style = n if case == "utube" else n  # Hidden by default, shown via "See all parameters" for coaxial
    
    coaxial_flow_type_style = b if (model in ["SBT V1.0", "SBT V2.0"] and case == "coaxial") else n
    
    if model == "HDF5":
        return n, n, n, n, n, n, b, n, n, n, n, n, n, n

    if model == "SBT V1.0":
        # lat-allo (Lateral Flow Allocation) unavailable for simulator
        return b, b, b, b, b, b, n, n, b, n, n, component_perf_style, see_all_button_style, coaxial_flow_type_style

    if model == "SBT V2.0":
        # lat-allo (Lateral Flow Allocation) unavailable for simulator
        return b, b, b, b, b, b, n, b, b, n, b, component_perf_style, see_all_button_style, coaxial_flow_type_style


@app.callback(
    [
        Output(
            component_id="economics-params-div",
            component_property="style",
            allow_duplicate=True,
        ),
        Output(component_id="kwt-div", component_property="style"),
        Output(component_id="kwe-div", component_property="style"),
    ],
    [
        Input(component_id="tabs", component_property="value"),
        # Input(component_id="model-select", component_property="value"),
        Input(component_id="fluid-select", component_property="value"),
        Input(component_id="end-use-select", component_property="value"),
    ],
    prevent_initial_call=True,
)
def econ_sliders_visibility(tab, fluid, end_use):
    if is_print:
        print("econ_sliders_visibility")
    b = {"display": "block"}
    n = {"display": "none"}

    econ_parms_div_style = {
        "display": "block",
        "border": "solid 3px #c4752f",
        "border-radius": "10px",
        "margin-bottom": "5px",
        "margin-right": "5px",
        "padding-bottom": "5px",
    }

    econ_parms_div_style_heating = {
        "display": "block",
        "border": "solid 3px #c4752f",
        "border-top": "solid 3px #c4752f",
        "border-left": "solid 3px #c4752f",
        "border-right": "solid 3px #c4752f",
        "border-bottom": "none",
        "border-radius": "10px 10px 0px 0px",
        "margin-bottom": "0px",
        "margin-right": "5px",
        "padding-bottom": "5px",
    }

    econ_parms_div_style_2 = {
        "display": "block",
        "border": "solid 3px #c4752f",
        "border-top": "solid 3px #c4752f",
        "border-left": "solid 3px #c4752f",
        "border-right": "solid 3px #c4752f",
        "border-bottom": "none",
        "border-radius": "10px 10px 0px 0px",
        "margin-bottom": "5px",
        "margin-right": "5px",
        "padding-bottom": "5px",
    }

    if tab == "energy-time-tab" or tab == "energy-tab":
        return n, n, n
    else:
        if fluid == "All" and end_use == "All":
            return econ_parms_div_style_2, b, b
        if fluid == "All" and end_use == "Heating":
            return econ_parms_div_style_heating, b, n
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
            return econ_parms_div_style_heating, b, n
        elif fluid == "sCO2" and end_use == "Electricity":
            return econ_parms_div_style_2, n, b
        else:
            return b, b, b


@app.callback(
    [
        Output(component_id="sCO2-card", component_property="style"),
        Output(component_id="sCO2-text", component_property="style"),
        Output(component_id="check-visual-card", component_property="style"),
    ],
    [
        Input(component_id="tabs", component_property="value"),
        Input(component_id="fluid-select", component_property="value"),
        Input(component_id="end-use-select", component_property="value"),
        Input(component_id="checklist", component_property="value"),
    ],
)
def show_hide_detailed_card(tab, fluid, end_use, checklist):
    if is_print:
        print("show_hide_detailed_card")
    if tab == "energy-time-tab" or tab == "energy-tab":
        return {"border": "solid 0px white"}, {"display": "none"}, {"display": "none"}

    if tab == "economics-time-tab" or tab == "about-tab" or tab == "summary-tab":
        if fluid == "H2O":
            return (
                {"border": "solid 0px white", "display": "none"},
                {"display": "none"},
                {"display": "none"},
            )
        elif end_use == "Heating":
            sCO2_card_style = {
                "border": "solid 3px #c4752f",
                "display": "block",
                "margin-top": "0px",
                "border-radius": "0px 0px 10px 10px",
                "border-top": "none",
                "border-right": "solid 3px #c4752f",
                "border-bottom": "solid 3px #c4752f",
                "border-left": "solid 3px #c4752f"
            }
            return (
                sCO2_card_style,
                {"display": "none"},
                {"display": "none"},
            )
        else:
            sCO2_card_style = {
                "border": "solid 3px #c4752f",
                "display": "block",
                "margin-top": "-8px",
                "border-radius": "0px 0px 10px 10px",
                "border-top": "none",
                "border-right": "solid 3px #c4752f",
                "border-bottom": "solid 3px #c4752f",
                "border-left": "solid 3px #c4752f"
            }
            check_visual_style = {
                "display": "inline-block",
                "border-top": "solid 2px #c4752f",
                "margin-top": "10px",
                "width": "auto",
                "height": "auto",
                "overflow": "visible",
                "clear": "both"
            }
            return (
                sCO2_card_style,
                {"display": "block"},
                check_visual_style,
            )
    
    return {"border": "solid 0px white", "display": "none"}, {"display": "none"}, {"display": "none"}


@app.callback(
    [
        Output(component_id="mdot-select-div", component_property="style"),
        Output(component_id="L2-select-div", component_property="style"),
        Output(component_id="L1-select-div", component_property="style"),
        Output(component_id="grad-select-div", component_property="style"),
        Output(
            component_id="diameter-select-div",
            component_property="style",
            allow_duplicate=True,
        ),
        Output(
            component_id="Tinj-select-div", component_property="style"
        ),  # AB DEBUG HERE
        Output(component_id="k-select-div", component_property="style"),
        Output(component_id="drillcost-div", component_property="style"),
        Output(component_id="discount-rate-div", component_property="style"),
        Output(component_id="lifetime-div", component_property="style"),
        Output(component_id="precool-div", component_property="style"),
        Output(component_id="turb-pout-div", component_property="style"),
        Output(component_id="wellbore-params-div", component_property="style"),
    ],
    [
        Input(component_id="param-select", component_property="value"),
        Input(component_id="tabs", component_property="value"),
        Input(component_id="fluid-select", component_property="value"),
        Input(component_id="end-use-select", component_property="value"),
        Input(component_id="model-select", component_property="value"),
    ],
    prevent_initial_call=True,
)
def show_hide_element(visibility_state, tab, fluid, end_use, model):
    if is_print:
        print("show_hide_element")
    # ----------------------------------------------------------------------------------------------
    # Reveals or hides sliders depending on which tab selected and which dropdowns.
    # ----------------------------------------------------------------------------------------------

    # print("show_hide_element: ", model, tab, fluid, end_use, visibility_state)
    # fluid getting changed makes this run "twice" but it's working as expected.

    b = {"display": "block"}
    n = {"display": "none"}

    if model == "HDF5":
        if tab == "about-tab":
            if fluid == "H2O" or end_use == "Heating":
                return b, b, b, b, b, b, b, b, b, b, n, n, b

            else:
                return b, b, b, b, b, b, b, b, b, b, b, b, b

        elif tab == "energy-time-tab":
            return b, b, b, b, b, b, b, n, n, n, n, n, b

        elif tab == "energy-tab":
            if visibility_state == param_list[0]:
                return n, n, b, b, b, b, b, n, n, n, n, n, b
            if visibility_state == param_list[1]:
                return n, b, n, b, b, b, b, n, n, n, n, n, b
            if visibility_state == param_list[2]:
                return n, b, b, n, b, b, b, n, n, n, n, n, b
            if visibility_state == param_list[3]:
                return n, b, b, b, n, b, b, n, n, n, n, n, b
            if visibility_state == param_list[4]:  # "Injection Temperature (˚C)"
                return n, b, b, b, b, n, b, n, n, n, n, n, n
            if visibility_state == param_list[5]:
                return n, b, b, b, b, b, n, n, n, n, n, n, b

        elif tab == "economics-time-tab":
            if fluid == "H2O":
                if end_use == "All":
                    return b, b, b, b, b, b, b, b, b, b, n, n, b
                if end_use == "Heating":
                    return b, b, b, b, b, b, b, b, b, b, n, n, b
                if end_use == "Electricity":
                    return b, b, b, b, b, b, b, b, b, b, n, n, b

            else:
                if end_use == "All":
                    return b, b, b, b, b, b, b, b, b, b, b, b, b
                if end_use == "Heating":
                    return b, b, b, b, b, b, b, b, b, b, n, n, b
                if end_use == "Electricity":
                    return b, b, b, b, b, b, b, b, b, b, b, b, b

        elif tab == "summary-tab":
            if fluid == "H2O":
                return b, b, b, b, b, b, b, b, b, b, n, n, b
            else:
                return b, b, b, b, b, b, b, b, b, b, b, b, b

    elif model == "SBT V1.0" or model == "SBT V2.0":
        if tab == "about-tab":
            if fluid == "H2O" or end_use == "Heating":
                return b, b, b, b, n, b, b, b, b, b, n, n, b
            else:
                return b, b, b, b, n, b, b, b, b, b, b, b, b

        elif tab == "energy-time-tab":
            return b, b, b, b, n, b, b, n, n, n, n, n, b

        elif tab == "energy-tab":
            if visibility_state == param_list[0]:
                return n, n, b, b, n, b, b, n, n, n, n, n, b
            if visibility_state == param_list[1]:
                return n, b, n, b, n, b, b, n, n, n, n, n, b
            if visibility_state == param_list[2]:
                return n, b, b, n, n, b, b, n, n, n, n, n, b
            if visibility_state == param_list[3]:
                return n, b, b, b, n, b, b, n, n, n, n, n, b
            if visibility_state == param_list[4]:
                return n, b, b, b, n, n, b, n, n, n, n, n, n
            if visibility_state == param_list[5]:
                return n, b, b, b, n, b, n, n, n, n, n, n, b

        elif tab == "economics-time-tab":
            if fluid == "H2O":
                if end_use == "All":
                    return b, b, b, b, n, b, b, b, b, b, n, n, b
                if end_use == "Heating":
                    return b, b, b, b, n, b, b, b, b, b, n, n, b
                if end_use == "Electricity":
                    return b, b, b, b, n, b, b, b, b, b, n, n, b
            else:
                if end_use == "All":
                    return b, b, b, b, n, b, b, b, b, b, b, b, b
                if end_use == "Heating":
                    return b, b, b, b, n, b, b, b, b, b, n, n, b
                if end_use == "Electricity":
                    return b, b, b, b, n, b, b, b, b, b, b, b, b

        elif tab == "summary-tab":
            if fluid == "H2O":
                return b, b, b, b, n, b, b, b, b, b, n, n, b
            else:
                return b, b, b, b, n, b, b, b, b, b, b, b, b
    else:
        raise PreventUpdate


@app.callback(
    [
        Output(
            component_id="grad-container",
            component_property="children",
            allow_duplicate=True,
        ),
        Output(
            component_id="k-container",
            component_property="children",
            allow_duplicate=True,
        ),
        Output(
            component_id="Tinj-container",
            component_property="children",
            allow_duplicate=True,
        ),
        Output(
            component_id="mdot-container",
            component_property="children",
            allow_duplicate=True,
        ),
        Output(component_id="diameter-container", component_property="children"),
        Output(component_id="L2-container", component_property="children"),
        Output(component_id="L1-container", component_property="children"),
    ],
    [
        Input(component_id="model-select", component_property="value"),
        Input(component_id="case-select", component_property="value"),
    ],
    [
        State(component_id="slider-values-store", component_property="data"),
    ],
    prevent_initial_call='initial_duplicate',
)
def update_slider_ranges(model, case, store_data):
    if is_print:
        print("update_slider_ranges")
    grad_dict = create_steps(
        arg_arr=u_sCO2.grad, str_round_place="{:.2f}", val_round_place=2
    )
    k_dict = create_steps(arg_arr=u_sCO2.k, str_round_place="{:.1f}", val_round_place=1)
    D_dict = create_steps(arg_arr=u_sCO2.D, str_round_place="{:.4f}", val_round_place=4)
    D_dict[0.444] = '0.444'

    if store_data is None:
        store_data = {}
    saved_values = store_data.get(model, {}).get(case, {})

    if model == "HDF5":  # hide the other params (happens in the next callback)
        Tinj_dict = {30: "30", 45: "", 59: "59"}
        mdot_dict = create_steps(
            arg_arr=u_sCO2.mdot, str_round_place="{:.1f}", val_round_place=1
        )
        L2_dict = create_steps(
            arg_arr=u_sCO2.L2, str_round_place="{:.0f}", val_round_place=0
        )
        L1_dict = create_steps(
            arg_arr=u_sCO2.L1, str_round_place="{:.0f}", val_round_place=0
        )

        grad_container = slider2(
            DivID="grad-select-div",
            ID="grad-select",
            ptitle="Geothermal Gradient (°C/m)",
            min_v=u_sCO2.grad[0],
            max_v=u_sCO2.grad[-1],
            mark_dict=grad_dict,
            start_v=saved_values.get("grad", start_vals_d["grad"]),
            div_style=div_block_style,
            parameter_name="Geothermal Gradient (°C/m)",
        )
        k_container = slider2(
            DivID="k-select-div",
            ID="k-select",
            ptitle="Rock Thermal Conductivity (W/m-°C)",
            min_v=u_sCO2.k[0],
            max_v=u_sCO2.k[-1],
            mark_dict=k_dict,
            start_v=saved_values.get("k", start_vals_d["k"]),
            div_style=div_block_style,
            parameter_name="Rock Thermal Conductivity (W/m-°C)",
            custom_title=True,
        )
        Tinj_container = slider2(
            DivID="Tinj-select-div",
            ID="Tinj-select",
            ptitle="Injection Temperature (˚C)",
            min_v=u_sCO2.Tinj[0] - 273.15,
            max_v=59.0,
            mark_dict=Tinj_dict,
            start_v=saved_values.get("Tinj", 55.0),
            div_style=div_block_style,
            parameter_name="Injection Temperature (˚C)",
        )
        mdot_container = slider2(
            DivID="mdot-select-div",
            ID="mdot-select",
            ptitle="Mass Flow Rate (kg/s)",
            min_v=u_sCO2.mdot[0],
            max_v=u_sCO2.mdot[-1],
            mark_dict=mdot_dict,
            start_v=saved_values.get("mdot", start_vals_d["mdot"]),
            div_style=div_block_style,
            parameter_name="Mass Flow Rate (kg/s)",
        )
        diameter_container = slider1(
            DivID="diameter-select-div",
            ID="diameter-select",
            ptitle="Borehole Diameter (m)",
            min_v=0.2159,
            max_v=0.444,
            mark_dict=D_dict,
            step_i=0.002,
            start_v=saved_values.get("D", start_vals_d["D"]),
            div_style=div_block_style,
            parameter_name="Borehole Diameter (m)",
        )
        L2_container = slider2(
            DivID="L2-select-div",
            ID="L2-select",
            ptitle="Horizontal Extent (m)",
            min_v=u_sCO2.L2[0],
            max_v=u_sCO2.L2[-1],
            mark_dict=L2_dict,
            start_v=saved_values.get("L2", start_vals_d["L2"]),
            div_style=div_block_style,
            parameter_name="Horizontal Extent (m)",
        )
        L1_container = slider2(
            DivID="L1-select-div",
            ID="L1-select",
            ptitle="Drilling Depth (m)",
            min_v=u_sCO2.L1[0],
            max_v=u_sCO2.L1[-1],
            mark_dict=L1_dict,
            start_v=saved_values.get("L1", start_vals_d["L1"]),
            div_style=div_block_style,
            parameter_name="Drilling Depth (m)",
        )

        return (
            grad_container,
            k_container,
            Tinj_container,
            mdot_container,
            diameter_container,
            L2_container,
            L1_container,
        )

    elif model == "SBT V1.0" or model == "SBT V2.0":
        # Tsurf_dict = {0: '0', 40: '40'}
        c_dict = {700: "700", 1200: "1200"}
        rho_dict = {400: "400", 4000: "4000"}

        Tinj_dict = {30: "30", 65: "", 100: "100"}
        grad_dict = {0.015: "0.015", 0.2: "0.200"}
        mdot_dict = {5: "5", 300: "300"}
        L2_dict = {1000: "1k", 50000: "50k"}
        L1_dict = {1000: "1k", 10000: "10k"}
        k_dict = {0.4: "0.4", 5: "5.0"}

        diameter_vertical_dict = {0.2159: "0.2159", 0.4445: "0.4445"}
        diameter_lateral_dict = {0.2159: "0.2159", 0.4445: "0.4445"}

        # Set defaults for coaxial SBT models to match database defaults
        if case == "coaxial":
            coaxial_default_mdot = 30.0
            coaxial_default_tinj = 55.0
            coaxial_default_grad = start_vals_d["grad"]
            coaxial_default_k = start_vals_d["k"]
        else:
            # Use normal defaults for U-tube
            coaxial_default_mdot = None
            coaxial_default_tinj = None
            coaxial_default_grad = None
            coaxial_default_k = None
        
        grad_container = slider2(
            DivID="grad-select-div",
            ID="grad-select",
            ptitle="Geothermal Gradient (°C/m)",  # min_v=0.01, max_v=0.1,
            min_v=0.015,
            max_v=0.200,
            mark_dict=grad_dict,
            start_v=coaxial_default_grad if (case == "coaxial" and coaxial_default_grad is not None) else (saved_values.get("grad", coaxial_default_grad if coaxial_default_grad is not None else start_vals_d["grad"])),
            div_style=div_block_style,
            parameter_name="Geothermal Gradient (°C/m)",
            step_i=0.001,
        )
        k_container = slider2(
            DivID="k-select-div",
            ID="k-select",
            ptitle="Rock Thermal Conductivity (W/m-°C)",  # min_v=0.1, max_v=7.0,
            min_v=0.4,
            max_v=5.0,
            mark_dict=k_dict,
            start_v=coaxial_default_k if (case == "coaxial" and coaxial_default_k is not None) else (saved_values.get("k", coaxial_default_k if coaxial_default_k is not None else start_vals_d["k"])),
            div_style=div_block_style,
            parameter_name="Rock Thermal Conductivity (W/m-°C)",
            custom_title=True,
        )
        Tinj_container = slider2(
            DivID="Tinj-select-div",
            ID="Tinj-select",
            ptitle="Injection Temperature (˚C)",
            min_v=30.0,
            max_v=100.0,
            # min_v=20.0, max_v=200.0,
            mark_dict=Tinj_dict,
            start_v=saved_values.get("Tinj", coaxial_default_tinj if (case == "coaxial" and coaxial_default_tinj is not None) else start_vals_d["Tinj"]),
            div_style=div_block_style,
            parameter_name="Injection Temperature (˚C)",
        )
        mdot_step = 1
        mdot_container = slider2(
            DivID="mdot-select-div",
            ID="mdot-select",
            ptitle="Mass Flow Rate (kg/s)",
            min_v=5,
            max_v=300,
            # min_v=u_sCO2.mdot[0], max_v=u_sCO2.mdot[-1],
            mark_dict=mdot_dict,
            start_v=saved_values.get("mdot", coaxial_default_mdot if (case == "coaxial" and coaxial_default_mdot is not None) else start_vals_d["mdot"]),
            div_style=div_block_style,
            parameter_name="Mass Flow Rate (kg/s)",
            step_i=mdot_step,
        )
        diameter_container = slider1(
            DivID="diameter-select-div",
            ID="diameter-select",
            ptitle="Borehole Diameter (m)",
            min_v=0.2159,
            max_v=0.444,
            mark_dict=D_dict,
            step_i=0.002,
            start_v=saved_values.get("D", start_vals_d["D"]),
            div_style=div_none_style,
        )
        L2_container = slider2(
            DivID="L2-select-div",
            ID="L2-select",
            ptitle="Horizontal Extent (m)",
            min_v=1000,
            max_v=50000,
            # min_v=u_sCO2.L2[0], max_v=u_sCO2.L2[-1],
            mark_dict=L2_dict,
            start_v=saved_values.get("L2", start_vals_d["L2"]),
            div_style=div_block_style,
            parameter_name="Horizontal Extent (m)",
        )
        L1_container = slider2(
            DivID="L1-select-div",
            ID="L1-select",
            ptitle="Drilling Depth (m)",
            min_v=1000,
            max_v=10000,
            # min_v=u_sCO2.L1[0], max_v=u_sCO2.L1[-1],
            mark_dict=L1_dict,
            start_v=saved_values.get("L1", start_vals_d["L1"]),
            div_style=div_block_style,
            parameter_name="Drilling Depth (m)",
        )

        return (
            grad_container,
            k_container,
            Tinj_container,
            mdot_container,
            diameter_container,
            L2_container,
            L1_container,
        )
    else:
        raise PreventUpdate


@app.callback(
    Output(component_id="slider-values-store", component_property="data"),
    [
        Input(component_id="mdot-select", component_property="value"),
        Input(component_id="L2-select", component_property="value"),
        Input(component_id="L1-select", component_property="value"),
        Input(component_id="grad-select", component_property="value"),
        Input(component_id="diameter-select", component_property="value"),
        Input(component_id="Tinj-select", component_property="value"),
        Input(component_id="k-select", component_property="value"),
        Input(component_id="Tsurf-select", component_property="value"),
        Input(component_id="c-select", component_property="value"),
        Input(component_id="rho-select", component_property="value"),
        Input(component_id="diameter-vertical-select", component_property="value"),
        Input(component_id="diameter-lateral-select", component_property="value"),
        Input(component_id="n-laterals-select", component_property="value"),
        Input(component_id="lateral-flow-select", component_property="value"),
        Input(component_id="insulation-thermal-conductivity-select", component_property="value"),
        Input(component_id="lateral-multiplier-select", component_property="value"),
        Input(component_id="mass-mode-select", component_property="data"),
        Input(component_id="inlet-pressure-select", component_property="value"),
        Input(component_id="model-select", component_property="value"),
        Input(component_id="case-select", component_property="value"),
    ],
    [
        State(component_id="slider-values-store", component_property="data"),
    ],
    prevent_initial_call=True,
)
def save_slider_values(mdot, L2, L1, grad, D, Tinj, k, Tsurf, c, rho, diameter_vertical, diameter_lateral, n_laterals, lateral_flow, insulation_k, lateral_multiplier, mass_mode, inlet_pressure, model, case, store_data):
    if is_print:
        print("save_slider_values")
    """Save slider values to store, keyed by (model, case) to prevent cross-contamination"""
    if model is None or case is None:
        raise PreventUpdate
    
    if store_data is None:
        store_data = {}
    
    model_bucket = store_data.get(model, {})
    case_bucket = model_bucket.get(case, {}).copy()
    
    case_bucket.update({
        "mdot": mdot,
        "L2": L2,
        "L1": L1,
        "grad": grad,
        "D": D,
        "Tinj": Tinj,
        "k": k,
        "Tsurf": Tsurf,
        "c": c,
        "rho": rho,
        "diameter-vertical": diameter_vertical,
        "diameter-lateral": diameter_lateral,
        "n-laterals": n_laterals,
        "lateral-multiplier": lateral_multiplier,
    })
    
    if model == "SBT V2.0":
        case_bucket["inlet-pressure"] = inlet_pressure
    else:
        case_bucket["inlet-pressure"] = mass_mode
    
    if case != "coaxial":
        case_bucket["lateral-flow"] = lateral_flow
    elif case == "coaxial" and model in ["SBT V1.0", "SBT V2.0"]:
        if insulation_k is not None and isinstance(insulation_k, (int, float)):
            case_bucket["coaxial-insulation-k"] = insulation_k
        if n_laterals is not None:
            case_bucket["thicknesscenterpipe"] = n_laterals
        if lateral_multiplier is not None and isinstance(lateral_multiplier, str):
            case_bucket["coaxial-flow-type"] = lateral_multiplier
    
    model_bucket[case] = case_bucket
    store_data[model] = model_bucket
    
    return store_data


@app.callback(
    [
        Output(component_id="Tsurf-select", component_property="value", allow_duplicate=True),
        Output(component_id="c-select", component_property="value", allow_duplicate=True),
        Output(component_id="rho-select", component_property="value", allow_duplicate=True),
    ],
    [
        Input(component_id="model-select", component_property="value"),
    ],
    [
        State(component_id="slider-values-store", component_property="data"),
        State(component_id="Tsurf-select", component_property="value"),
        State(component_id="c-select", component_property="value"),
        State(component_id="rho-select", component_property="value"),
    ],
    prevent_initial_call=True,
)
def restore_always_visible_sliders(model, store_data, current_tsurf, current_c, current_rho):
    if is_print:
        print("restore_always_visible_sliders")
    """Restore Tsurf, c, and rho when model changes - these sliders are always present"""
    if model is None:
        raise PreventUpdate
    
    if store_data is None:
        store_data = {}
    
    # Helper function to get value from current model or fall back to other models
    def get_value(key, current_value, default):
        # First check current model
        if model in store_data and key in store_data[model]:
            val = store_data[model][key]
            if val is not None:
                return val
        
        # If not found, check other models (prioritize SBT models since Tsurf is visible there)
        for other_model in ["SBT V2.0", "SBT V1.0", "HDF5"]:
            if other_model in store_data and key in store_data[other_model]:
                val = store_data[other_model][key]
                if val is not None:
                    return val
        
        # If no saved value found, preserve current value if it exists, otherwise use default
        if current_value is not None:
            return current_value
        return default
    
    # Restore these values, checking current model first, then other models
    # If no saved value found, preserve current slider value
    Tsurf = get_value("Tsurf", current_tsurf, start_vals_hdf5.get("Tsurf", 25))
    c = get_value("c", current_c, start_vals_hdf5.get("c", 790.0))
    rho = get_value("rho", current_rho, start_vals_hdf5.get("rho", 2800))
    
    # Only update if values have changed to prevent unnecessary cascading updates
    if current_tsurf == Tsurf and current_c == c and current_rho == rho:
        raise PreventUpdate
    
    return Tsurf, c, rho


# COMMENTED OUT: restore_slider_values was removed because it's redundant with update_slider_ranges
# which already sets saved values via start_v when creating sliders. This callback was causing
# React unmounting errors by trying to update slider values while update_slider_ranges was
# recreating them. The batching benefit is no longer needed since values are set during creation.
#
# @app.callback(
#     [
#         Output(component_id="mdot-select", component_property="value", allow_duplicate=True),
#         Output(component_id="L2-select", component_property="value", allow_duplicate=True),
#         Output(component_id="L1-select", component_property="value", allow_duplicate=True),
#         Output(component_id="grad-select", component_property="value", allow_duplicate=True),
#         Output(component_id="diameter-select", component_property="value", allow_duplicate=True),
#         Output(component_id="Tinj-select", component_property="value", allow_duplicate=True),
#         Output(component_id="k-select", component_property="value", allow_duplicate=True),
#     ],
#     [
#         Input(component_id="model-select", component_property="value"),  # Trigger only on model switches
#         Input(component_id="case-select", component_property="value"),
#     ],
#     [
#         State(component_id="slider-values-store", component_property="data"),
#         State(component_id="mdot-select", component_property="value"),
#         State(component_id="L2-select", component_property="value"),
#         State(component_id="L1-select", component_property="value"),
#         State(component_id="grad-select", component_property="value"),
#         State(component_id="diameter-select", component_property="value"),
#         State(component_id="Tinj-select", component_property="value"),
#         State(component_id="k-select", component_property="value"),
#     ],
#     prevent_initial_call=True,
# )
# def restore_slider_values(model, case, store_data, current_mdot, current_L2, current_L1, current_grad, current_D, current_Tinj, current_k):
#     """Restore slider values from store when model switches (not on initial load)"""
#     if is_print:
#         print("restore_slider_values")
#     if model is None:
#         raise PreventUpdate
#     
#     # Only restore when model-select actually changed (not on initial load)
#     if not ctx.triggered or ctx.triggered_id != "model-select":
#         raise PreventUpdate
#     
#     # Don't restore values for coaxial - let users set their own values
#     # This prevents the callback from overriding user-set slider values
#     if case == "coaxial" and model in ["SBT V1.0", "SBT V2.0"]:
#         raise PreventUpdate
#     
#     if store_data is None:
#         store_data = {}
#     
#     # Helper function to get value from current model or fall back to other models
#     # Returns (value, found) tuple where found indicates if value was found in store
#     def get_value(key):
#         # First check current model
#         if model in store_data and key in store_data[model]:
#             val = store_data[model][key]
#             if val is not None:
#                 return val, True
#         
#         # If not found, check other models
#         for other_model in ["HDF5", "SBT V2.0", "SBT V1.0"]:
#             if other_model in store_data and key in store_data[other_model]:
#                 val = store_data[other_model][key]
#                 if val is not None:
#                     return val, True
#         
#         # Return None if not found anywhere
#         return None, False
#     
#     # Try to restore values from store for all cases
#     # Don't force defaults - let users set their own values
#     mdot, mdot_found = get_value("mdot")
#     L2, L2_found = get_value("L2")
#     L1, L1_found = get_value("L1")
#     grad, grad_found = get_value("grad")
#     D, D_found = get_value("D")
#     Tinj, Tinj_found = get_value("Tinj")
#     k, k_found = get_value("k")
#     
#     if not any([mdot_found, L2_found, L1_found, grad_found, D_found, Tinj_found, k_found]):
#         raise PreventUpdate
#     
#     # Use defaults for any None values
#     mdot = mdot if mdot is not None else start_vals_d.get("mdot", 30.0)
#     L2 = L2 if L2 is not None else start_vals_d.get("L2", 10000)
#     L1 = L1 if L1 is not None else start_vals_d.get("L1", 3500)
#     grad = grad if grad is not None else start_vals_d.get("grad", 0.065)
#     D = D if D is not None else start_vals_d.get("D", 0.3500)
#     Tinj = Tinj if Tinj is not None else start_vals_d.get("Tinj", 60.0)
#     k = k if k is not None else start_vals_d.get("k", 3.0)
#     
#     # Only update if values have changed to prevent unnecessary cascading updates
#     if (current_mdot == mdot and current_L2 == L2 and current_L1 == L1 and 
#         current_grad == grad and current_D == D and current_Tinj == Tinj and current_k == k):
#         raise PreventUpdate
#     
#     return mdot, L2, L1, grad, D, Tinj, k


@app.callback(
    [
        Output(
            component_id="Diameter1-container", component_property="children"
        ),  # Diameter1
        Output(
            component_id="Diameter2-container", component_property="children"
        ),  # Diameter2
        Output(
            component_id="num-lat-container", component_property="children"
        ),  # PipeParam3
        Output(
            component_id="lat-allo-container", component_property="children"
        ),  # PipeParam4
        Output(
            component_id="lat-flow-container", component_property="children"
        ),  # PipeParam5
    ],
    [
        Input(component_id="model-select", component_property="value"),
        Input(component_id="case-select", component_property="value"),
    ],
    [
        State(component_id="slider-values-store", component_property="data"),
    ],
    prevent_initial_call=True,
)
def update_sliders_heat_exchanger(model, case, store_data):
    if is_print:
        print("update_sliders_heat_exchanger")
    if store_data is None:
        store_data = {}
    saved_values = store_data.get(model, {}).get(case, {})
    
    def get_saved_value(key, default_dict, default_key):
        if key in saved_values and saved_values[key] is not None:
            return saved_values[key]
        for other_model in ["SBT V2.0", "SBT V1.0"]:
            if other_model in store_data and case in store_data[other_model] and key in store_data[other_model][case]:
                val = store_data[other_model][case][key]
                if val is not None:
                    return val
        return default_dict.get(default_key)
    
    if model == "SBT V1.0" or model == "SBT V2.0":
        if case == "utube":
            diameter_vertical = slider1(
                DivID="diameter-vertical-select-div",
                ID="diameter-vertical-select",
                ptitle="Wellbore Diameter Vertical (m)",
                min_v=0.2159,
                max_v=0.4445,
                mark_dict=diameter_vertical_dict,
                step_i=0.002,
                start_v=get_saved_value("diameter-vertical", start_vals_sbt, "diameter-vertical"),
                div_style=div_block_style,
                parameter_name="Wellbore Diameter Vertical (m)",
            )
            diameter_lateral = slider1(
                DivID="diameter-lateral-select-div",
                ID="diameter-lateral-select",
                ptitle="Wellbore Diameter Lateral (m)",
                min_v=0.2159,
                max_v=0.4445,
                mark_dict=diameter_lateral_dict,
                step_i=0.002,
                start_v=get_saved_value("diameter-lateral", start_vals_sbt, "diameter-lateral"),
                div_style=div_block_style,
                parameter_name="Wellbore Diameter Lateral (m)",
            )
            n_laterals_val = get_saved_value("n-laterals", start_vals_hdf5, "n-laterals")
            n_laterals = input_box(
                DivID="num-lat-div",
                ID="n-laterals-select",
                ptitle="Number of Laterals",
                min_v=0,
                max_v=3,
                start_v=n_laterals_val if n_laterals_val is not None and n_laterals_val >= 0 else 0,
                step_i=1,
                div_style=div_block_style,
                parameter_name="Number of Laterals",
                horizontal=True,
                input_width="60px",
            )
            lateral_flow_value = get_saved_value("lateral-flow", start_vals_hdf5, "lateral-flow")
            # Lateral Flow Allocation is unavailable for simulator models, always return hidden placeholder
            lateral_flow = lateral_flow_placeholder(value=lateral_flow_value, style=div_none_style)
            # Calculate multiplier based on number of laterals: 1/n_laterals
            multiplier_val = 1.0 / n_laterals_val if n_laterals_val and n_laterals_val > 0 else 1.0
            lateral_multiplier = input_box(
                DivID="lat-flow-mul-div",
                ID="lateral-multiplier-select",
                ptitle="Lateral Flow Multiplier",
                min_v=0,
                max_v=1,
                start_v=multiplier_val,
                step_i=0.05,
                div_style=div_block_style,
                parameter_name="Lateral Flow Multiplier",
                horizontal=True,
                input_width="60px",
            )

            return (
                diameter_vertical,
                diameter_lateral,
                n_laterals,
                lateral_flow,
                lateral_multiplier,
            )

        elif case == "coaxial":
            insulation_thermal_k_dict = {0.025: "0.025", 0.50: "0.5"}
            radius_centerpipe_dict = {0.127: "0.127", 0.348: "0.348"}

            coaxial_defaults = {"diameter-vertical": 0.4445, "diameter-lateral": 0.2}
            # Get both values first - filtering will handle invalid values
            wellbore_diameter = get_saved_value("diameter-vertical", coaxial_defaults, "diameter-vertical")
            centerpipe_diameter = get_saved_value("diameter-lateral", coaxial_defaults, "diameter-lateral")
            
            # Only use defaults if we truly have no valid values (filtered out)
            # Otherwise use the slider values, even if they need adjustment
            if wellbore_diameter is None:
                wellbore_diameter = coaxial_defaults["diameter-vertical"]
            if centerpipe_diameter is None:
                centerpipe_diameter = coaxial_defaults["diameter-lateral"]

            radius = slider1(
                DivID="diameter-vertical-select-div",
                ID="diameter-vertical-select",
                ptitle="Wellbore Diameter (m)",
                min_v=0.2159,
                max_v=0.4445,
                mark_dict=diameter_vertical_dict,
                step_i=0.002,
                start_v=wellbore_diameter,
                div_style=div_block_style,
                parameter_name="Wellbore Diameter (m)",
            )

            radiuscenterpipe = slider1(
                DivID="diameter-lateral-select-div",
                ID="diameter-lateral-select",
                ptitle="Center Pipe Diameter (m)",
                min_v=0.127,
                max_v=0.348,
                mark_dict=radius_centerpipe_dict,
                step_i=0.002,
                start_v=centerpipe_diameter,
                div_style=div_block_style,
                parameter_name="Center Pipe Diameter (m)",
            )

            thicknesscenterpipe = slider1(
                DivID="num-lat-div",
                ID="n-laterals-select",
                ptitle="Center Pipe Thickness (m)",
                min_v=0.005,
                max_v=0.025,
                mark_dict=thickness_centerpipe_dict,
                step_i=0.001,
                start_v=get_saved_value("thicknesscenterpipe", start_vals_sbt, "thicknesscenterpipe"),
                div_style=div_block_style,
            )

            k_center_pipe = lateral_flow_placeholder(value=1.0, style=div_none_style)
            
            saved_flow_type = get_saved_value("coaxial-flow-type", {"coaxial-flow-type": "Inject in Annulus"}, "coaxial-flow-type")
            if saved_flow_type not in ["Inject in Annulus", "Inject in Center Pipe"]:
                saved_flow_type = "Inject in Annulus"
            
            # Create dropdown with explicit value
            flow_type_options = ["Inject in Annulus", "Inject in Center Pipe"]
            coaxialflowtype = html.Div(
                id="lat-flow-mul-div",
                className="name-input-container-dd",
                style=div_block_style,
                children=[
                    html.P("Coaxial Flow Type", className="input-title"),
                    dcc.Dropdown(
                        id="lateral-multiplier-select",
                        options=flow_type_options,
                        value=saved_flow_type,
                        clearable=False,
                        searchable=False,
                        disabled=False,
                        className="select-dropdown"
                    ),
                ])

            return (
                radius,
                radiuscenterpipe,
                thicknesscenterpipe,
                k_center_pipe,
                coaxialflowtype,
            )

    elif model == "HDF5":
        # HDF5 doesn't support laterals - but components must exist for callbacks
        diameter_vertical = slider1(
            DivID="diameter-vertical-select-div",
            ID="diameter-vertical-select",
            ptitle="Wellbore Diameter Vertical (m)",
            min_v=0.4,
            max_v=1.2,
            mark_dict=diameter_vertical_dict,
            step_i=0.002,
            start_v=start_vals_sbt["diameter-vertical"],
            div_style=div_none_style,
            parameter_name="Wellbore Diameter Vertical (m)",
        )
        diameter_lateral = slider1(
            DivID="diameter-lateral-select-div",
            ID="diameter-lateral-select",
            ptitle="Wellbore Diameter Lateral (m)",
            min_v=0.4,
            max_v=1.2,
            mark_dict=diameter_lateral_dict,
            step_i=0.002,
            start_v=start_vals_sbt["diameter-lateral"],
            div_style=div_none_style,
            parameter_name="Wellbore Diameter Lateral (m)",
        )
        # HDF5 doesn't use laterals - but create hidden components for callbacks
        n_laterals = input_box(
            DivID="num-lat-div",
            ID="n-laterals-select",
            ptitle="Number of Laterals",
            min_v=0,
            max_v=20,
            start_v=1,
            step_i=1,
            div_style={"position": "absolute", "left": "-9999px", "visibility": "hidden"},
            parameter_name="Number of Laterals",
        )
        lateral_flow = input_box(
            DivID="lat-allocation-div",
            ID="lateral-flow-select",
            ptitle="Lateral Flow Allocation",
            min_v=0,
            max_v=1,
            start_v=1.0,
            step_i=0.01,
            div_style={"position": "absolute", "left": "-9999px", "visibility": "hidden"},
            parameter_name="Lateral Flow Allocation",
        )
        lateral_multiplier = input_box(
            DivID="lat-flow-mul-div",
            ID="lateral-multiplier-select",
            ptitle="Lateral Flow Multiplier",
            min_v=0,
            max_v=1,
            start_v=1.0,
            step_i=0.05,
            div_style={"position": "absolute", "left": "-9999px", "visibility": "hidden"},
            parameter_name="Lateral Flow Multiplier",
            horizontal=True,
            input_width="60px",
        )

        return (
            diameter_vertical,
            diameter_lateral,
            n_laterals,
            lateral_flow,
            lateral_multiplier,
        )

    else:
        raise PreventUpdate


@app.callback(
    [
        Output(component_id="diameter-vertical-select", component_property="value", allow_duplicate=True),
        Output(component_id="diameter-lateral-select", component_property="value", allow_duplicate=True),
        Output(component_id="n-laterals-select", component_property="value", allow_duplicate=True),
        Output(component_id="lateral-flow-select", component_property="value", allow_duplicate=True),
        Output(component_id="lateral-multiplier-select", component_property="value", allow_duplicate=True),
    ],
    [
        Input(component_id="model-select", component_property="value"),
        Input(component_id="case-select", component_property="value"),
    ],
    [
        State(component_id="slider-values-store", component_property="data"),
        State(component_id="diameter-vertical-select", component_property="value"),
        State(component_id="diameter-lateral-select", component_property="value"),
        State(component_id="n-laterals-select", component_property="value"),
        State(component_id="lateral-flow-select", component_property="value"),
        State(component_id="lateral-multiplier-select", component_property="value"),
    ],
    prevent_initial_call=True,
)
def restore_sbt_sliders(model, case, store_data, current_diameter_vertical, current_diameter_lateral, current_n_laterals, current_lateral_flow, current_lateral_multiplier):
    if is_print:
        print("restore_sbt_sliders")
    """Restore SBT-specific slider values when model or case switches"""
    if model is None:
        raise PreventUpdate
    
    # Only restore for SBT models
    if model not in ["SBT V1.0", "SBT V2.0"]:
        raise PreventUpdate
    
    # Check what triggered this callback
    triggered_id = ctx.triggered[0]["prop_id"].split(".")[0] if ctx.triggered else None
    
    if case == "coaxial" and model in ["SBT V1.0", "SBT V2.0"]:
        raise PreventUpdate
    
    # If case-select triggered, restore values for the new case
    if triggered_id == "case-select":
        if store_data is None:
            store_data = {}
        
        saved_values = store_data.get(model, {}).get(case, {})
        
        def get_value(key, current_value, default_dict, default_key):
            if key in saved_values and saved_values[key] is not None:
                return saved_values[key]
            for other_model in ["SBT V2.0", "SBT V1.0"]:
                if other_model in store_data and case in store_data[other_model] and key in store_data[other_model][case]:
                    val = store_data[other_model][case][key]
                    if val is not None:
                        return val
            if current_value is not None:
                return current_value
            return default_dict.get(default_key)
        
        if case == "coaxial":
            saved_flow_type = get_value("coaxial-flow-type", None, {"coaxial-flow-type": "Inject in Annulus"}, "coaxial-flow-type")
            if saved_flow_type not in ["Inject in Annulus", "Inject in Center Pipe"]:
                saved_flow_type = "Inject in Annulus"
            lateral_multiplier = saved_flow_type
            
            if current_lateral_multiplier == lateral_multiplier:
                raise PreventUpdate
            
            return current_diameter_vertical, current_diameter_lateral, current_n_laterals, no_update, lateral_multiplier
        else:
            # For utube, restore lateral flow allocation and multiplier
            lateral_flow = get_value("lateral-flow", current_lateral_flow, start_vals_hdf5, "lateral-flow")
            n_laterals = get_value("n-laterals", current_n_laterals, start_vals_hdf5, "n-laterals")
            if n_laterals is not None and n_laterals > 0:
                lateral_multiplier = 1.0 / n_laterals
            else:
                lateral_multiplier = get_value("lateral-multiplier", current_lateral_multiplier, start_vals_hdf5, "lateral-multiplier")
            
            # Only update lateral-flow and lateral-multiplier, keep others unchanged
            if current_lateral_flow == lateral_flow and current_lateral_multiplier == lateral_multiplier:
                raise PreventUpdate
            
            return current_diameter_vertical, current_diameter_lateral, current_n_laterals, lateral_flow, lateral_multiplier
    
    # For model-select changes, restore all sliders (existing behavior)
    if triggered_id != "model-select":
        raise PreventUpdate
    
    if store_data is None:
        store_data = {}
    
    saved_values = store_data.get(model, {}).get(case, {})
    
    def get_value(key, current_value, default_dict, default_key):
        if key in saved_values and saved_values[key] is not None:
            return saved_values[key]
        for other_model in ["SBT V2.0", "SBT V1.0"]:
            if other_model in store_data and case in store_data[other_model] and key in store_data[other_model][case]:
                val = store_data[other_model][case][key]
                if val is not None:
                    return val
        if current_value is not None:
            return current_value
        return default_dict.get(default_key)
    
    if case == "coaxial":
        coaxial_defaults = {"diameter-vertical": 0.4445, "diameter-lateral": 0.2}
        diameter_vertical = get_value("diameter-vertical", current_diameter_vertical, coaxial_defaults, "diameter-vertical")
        diameter_lateral = get_value("diameter-lateral", current_diameter_lateral, coaxial_defaults, "diameter-lateral")
    else:
        diameter_vertical = get_value("diameter-vertical", current_diameter_vertical, start_vals_sbt, "diameter-vertical")
        diameter_lateral = get_value("diameter-lateral", current_diameter_lateral, start_vals_sbt, "diameter-lateral")
    
    n_laterals = get_value("n-laterals", current_n_laterals, start_vals_hdf5, "n-laterals")
    
    if case == "coaxial":
        
        saved_flow_type = get_value("coaxial-flow-type", None, {"coaxial-flow-type": "Inject in Annulus"}, "coaxial-flow-type")
        if saved_flow_type not in ["Inject in Annulus", "Inject in Center Pipe"]:
            saved_flow_type = "Inject in Annulus"
        lateral_multiplier = saved_flow_type
        
        # Only update if flow type changed
        if current_lateral_multiplier == lateral_multiplier:
            raise PreventUpdate
        
        return diameter_vertical, diameter_lateral, n_laterals, no_update, lateral_multiplier
    else:
        lateral_flow = get_value("lateral-flow", current_lateral_flow, start_vals_hdf5, "lateral-flow")
        if n_laterals is not None and n_laterals > 0:
            lateral_multiplier = 1.0 / n_laterals
        else:
            lateral_multiplier = get_value("lateral-multiplier", current_lateral_multiplier, start_vals_hdf5, "lateral-multiplier")
    
    # Only update if values have changed to prevent unnecessary cascading updates
    if (current_diameter_vertical == diameter_vertical and current_diameter_lateral == diameter_lateral and 
        current_n_laterals == n_laterals and current_lateral_flow == lateral_flow and 
        current_lateral_multiplier == lateral_multiplier):
        raise PreventUpdate
    
    return diameter_vertical, diameter_lateral, n_laterals, lateral_flow, lateral_multiplier


@app.callback(
    Output(component_id="lateral-multiplier-select", component_property="value", allow_duplicate=True),
    [
        Input(component_id="n-laterals-select", component_property="value"),
        Input(component_id="model-select", component_property="value"),
        Input(component_id="case-select", component_property="value"),
    ],
    prevent_initial_call=True,
)
def update_lateral_multiplier(n_laterals, model, case):
    if is_print:
        print("update_lateral_multiplier")
    """Automatically update lateral flow multiplier based on number of laterals.
    
    The multiplier should be 1/n_laterals:
    - 1 lateral: multiplier = 1
    - 2 laterals: multiplier = 1/2 = 0.5
    - 3 laterals: multiplier = 1/3 ≈ 0.333
    
    For coaxial, this component is used for the flow type dropdown, so don't update it.
    """
    if model is None or model not in ["SBT V1.0", "SBT V2.0"]:
        raise PreventUpdate
    
    if case == "coaxial":
        raise PreventUpdate
    
    if n_laterals is None or n_laterals <= 0:
        raise PreventUpdate
    
    multiplier = 1.0 / n_laterals
    
    return multiplier


@app.callback(
    Output(component_id="lateral-flow-select", component_property="value", allow_duplicate=True),
    [
        Input({"type": "lateral-flow-input", "index": ALL}, component_property="value"),
        Input(component_id="model-select", component_property="value"),
        Input(component_id="case-select", component_property="value"),
    ],
    prevent_initial_call=True,
)
def aggregate_lateral_flow_inputs(input_values, model, case):
    if is_print:
        print("aggregate_lateral_flow_inputs")
    """Aggregate multiple lateral flow allocation inputs into a single value for callbacks.
    
    For SBT models with utube case, we show multiple inputs (one per lateral).
    This callback aggregates them into a single value that other callbacks can use.
    The value is the average allocation per lateral (for backward compatibility).
    Handles both fraction strings (e.g., "1/3") and decimal values.
    """
    if model is None or model not in ["SBT V1.0", "SBT V2.0"]:
        raise PreventUpdate
    if case != "utube":
        raise PreventUpdate
    
    if not input_values or len(input_values) == 0:
        raise PreventUpdate
    
    # Convert fraction strings to decimals
    valid_decimal_values = []
    for v in input_values:
        if v is None:
            continue
        # Try to parse as fraction (e.g., "1/3") or decimal
        try:
            if isinstance(v, str) and '/' in v:
                # Parse fraction
                frac = Fraction(v)
                valid_decimal_values.append(float(frac))
            else:
                # Parse as decimal (handles both string and numeric)
                valid_decimal_values.append(float(v))
        except (ValueError, ZeroDivisionError):
            continue
    
    if len(valid_decimal_values) == 0:
        raise PreventUpdate
    
    # Return average allocation per lateral (for backward compatibility with single input)
    avg_allocation = sum(valid_decimal_values) / len(valid_decimal_values)
    return avg_allocation


@app.callback(
    Output(component_id="lat-allo-container", component_property="children", allow_duplicate=True),
    Input(component_id="n-laterals-select", component_property="value"),
    Input(component_id="model-select", component_property="value"),
    Input(component_id="case-select", component_property="value"),
    State(component_id="slider-values-store", component_property="data"),
    prevent_initial_call="initial_duplicate",
)
def update_lateral_flow_allocation_inputs(n_laterals, model, case, store_data):
    if is_print:
        print("update_lateral_flow_allocation_inputs")
    """Update the lateral flow allocation inputs when number of laterals changes.
    
    Creates multiple input fields (one per lateral) with equal distribution.
    Displays fractions (e.g., "1/3") in the inputs.
    
    Note: Lateral Flow Allocation is unavailable for simulator (SBT) models.
    """
    if case == "coaxial":
        raise PreventUpdate
    
    # Lateral Flow Allocation is unavailable for all models (always return placeholder)
    return lateral_flow_placeholder(style=div_none_style)


@app.callback(
    [
        Output(component_id="hyperparam1-container", component_property="children"),
        Output(component_id="hyperparam3-container", component_property="children"),
        Output(component_id="hyperparam5-container", component_property="children"),
        Output(component_id="pipe-roughness-container", component_property="children"),
        Output(component_id="inlet-pressure-container", component_property="children"),
    ],
    [
        Input(component_id="model-select", component_property="value"),
    ],
    [
        State(component_id="slider-values-store", component_property="data"),
    ],
    prevent_initial_call=False,
)
def update_sliders_hyperparms(model, store_data):
    if is_print:
        print("update_sliders_hyperparms")
    if store_data is None:
        store_data = {}
    
    if model == "SBT V1.0":
        hyperparam1 = html.Div(id="hyperparam1-container", children=[])
        hyperparam3 = html.Div(id="hyperparam3-container", children=[])
        hyperparam5 = dropdown_box(
            DivID="fluid-mode-div",
            ID="fluid-mode-select",
            ptitle="Fluid Properties Mode",
            options=["Constant", "Variable"],
            disabled=False,
            div_style=div_none_style,
        )
        pipe_roughness = slider1(
            DivID="pipe-roughness-div",
            ID="pipe-roughness-select",
            ptitle="Pipe Roughness (µm)",
            min_v=1,
            max_v=3,
            mark_dict=pipe_roughness_um_dict,
            step_i=0.1,
            start_v=start_vals_sbt["piperoughness"],
            div_style=div_none_style,  # Hidden for SBT V1.0
            parameter_name="Pipe Roughness (µm)",
        )
        inlet_pressure = slider1(
            DivID="inlet-pressure-div",
            ID="inlet-pressure-select",
            ptitle="Inlet Pressure (MPa)",
            min_v=5,
            max_v=20,
            mark_dict=inlet_pressure_dict,
            step_i=0.1,
            start_v=start_vals_sbt["inletpressure"],
            div_style=div_none_style,  # Hidden for SBT V1.0
            parameter_name="Inlet Pressure (MPa)",
        )
        return hyperparam1, hyperparam3, hyperparam5, pipe_roughness, inlet_pressure

    elif model == "SBT V2.0":
        # Get saved inlet pressure value from store, checking current model first, then other models
        saved_inlet_pressure = None
        if model in store_data and "inlet-pressure" in store_data[model]:
            saved_inlet_pressure = store_data[model]["inlet-pressure"]
        else:
            # Check other models for saved value
            for other_model in ["SBT V1.0", "HDF5"]:
                if other_model in store_data and "inlet-pressure" in store_data[other_model]:
                    val = store_data[other_model]["inlet-pressure"]
                    # Only use if it's a numeric value (for SBT V2.0, inlet pressure is a number)
                    if isinstance(val, (int, float)):
                        saved_inlet_pressure = val
                        break
        inlet_pressure_start = saved_inlet_pressure if (saved_inlet_pressure is not None and isinstance(saved_inlet_pressure, (int, float))) else start_vals_sbt["inletpressure"]
        
        inlet_pressure = slider1(
            DivID="inlet-pressure-div",
            ID="inlet-pressure-select",
            ptitle="Inlet Pressure (MPa)",
            min_v=5,
            max_v=20,
            mark_dict=inlet_pressure_dict,
            step_i=0.1,
            start_v=inlet_pressure_start,
            div_style=div_block_style,
            parameter_name="Inlet Pressure (MPa)",
        )
        hyperparam1 = html.Div()
        hyperparam3 = html.Div(id="hyperparam3-container", children=[])
        pipe_roughness = slider1(
            DivID="pipe-roughness-div",
            ID="pipe-roughness-select",
            ptitle="Pipe Roughness (µm)",
            min_v=1,
            max_v=3,
            mark_dict=pipe_roughness_um_dict,
            step_i=0.1,
            start_v=start_vals_sbt["piperoughness"],
            div_style=div_block_style,
            parameter_name="Pipe Roughness (µm)",
        )
        saved_fluid_mode = None
        if model in store_data and "fluid-mode" in store_data[model]:
            saved_fluid_mode = store_data[model]["fluid-mode"]
        
        fluid_mode_options = ["Constant", "Temperature–Pressure Dependent"]
        if saved_fluid_mode == "Variable":
            default_fluid_mode = "Temperature–Pressure Dependent"
        elif saved_fluid_mode in fluid_mode_options:
            default_fluid_mode = saved_fluid_mode
        else:
            default_fluid_mode = "Temperature–Pressure Dependent"
        
        hyperparam5 = html.Div(
            id="fluid-mode-div",
            className="name-input-container-dd",
            style=div_block_style,
            children=[
                html.P("Fluid Properties Mode", className="input-title"),
                dcc.Dropdown(
                    id="fluid-mode-select",
                    options=[{"label": opt, "value": opt} for opt in fluid_mode_options],
                    value=default_fluid_mode,
                    clearable=False,
                    searchable=False,
                    disabled=False,
                    className="select-dropdown"
                ),
            ])

        return hyperparam1, hyperparam3, hyperparam5, pipe_roughness, inlet_pressure

    elif model == "HDF5":
        hyperparam1 = html.Div(id="hyperparam1-container", children=[])
        hyperparam3 = html.Div(id="hyperparam3-container", children=[])
        hyperparam5 = dropdown_box(
            DivID="fluid-mode-div",
            ID="fluid-mode-select",
            ptitle="Fluid Properties Mode",
            options=["Constant", "Variable"],
            disabled=False,
            div_style=div_none_style,
        )
        pipe_roughness = slider1(
            DivID="pipe-roughness-div",
            ID="pipe-roughness-select",
            ptitle="Pipe Roughness (µm)",
            min_v=1,
            max_v=3,
            mark_dict=pipe_roughness_um_dict,
            step_i=0.1,
            start_v=start_vals_sbt["piperoughness"],
            div_style=div_none_style,  # Hidden for HDF5
            parameter_name="Pipe Roughness (µm)",
        )
        inlet_pressure = slider1(
            DivID="inlet-pressure-div",
            ID="inlet-pressure-select",
            ptitle="Inlet Pressure (MPa)",
            min_v=5,
            max_v=20,
            mark_dict=inlet_pressure_dict,
            step_i=0.1,
            start_v=start_vals_sbt["inletpressure"],
            div_style=div_none_style,  # Hidden for HDF5
            parameter_name="Inlet Pressure (MPa)",
        )
        return hyperparam1, hyperparam3, hyperparam5, pipe_roughness, inlet_pressure
    else:
        raise PreventUpdate


# -----------------------------------------------------------------------------
# Define dash app plotting callbacks.
# -----------------------------------------------------------------------------


@app.callback(
    Output(component_id="plot-params-store", component_property="data"),
    [
        Input(component_id="interpolation-select", component_property="value"),
        Input(component_id="fluid-select", component_property="value"),
        Input(component_id="case-select", component_property="value"),
        Input(component_id="mdot-select", component_property="value"),
        Input(component_id="L2-select", component_property="value"),
        Input(component_id="L1-select", component_property="value"),
        Input(component_id="grad-select", component_property="value"),
        Input(component_id="diameter-select", component_property="value"),
        Input(component_id="Tinj-select", component_property="value"),
        Input(component_id="k-select", component_property="value"),
        Input(component_id="radio-graphic-control3", component_property="value"),
        Input(component_id="model-select", component_property="value"),
        Input(component_id="Tsurf-select", component_property="value"),
        Input(component_id="c-select", component_property="value"),
        Input(component_id="rho-select", component_property="value"),
        Input(component_id="diameter-vertical-select", component_property="value"),
        Input(component_id="diameter-lateral-select", component_property="value"),
        Input(component_id="n-laterals-select", component_property="value"),
        Input(component_id="lateral-flow-select", component_property="value"),
        Input(component_id="insulation-thermal-conductivity-select", component_property="value"),
        Input(component_id="lateral-multiplier-select", component_property="value"),
        Input(component_id="mesh-select", component_property="value"),
        Input(component_id="accuracy-select", component_property="value"),
        Input(component_id="mass-mode-select", component_property="data"),
        Input(component_id="temp-mode-select", component_property="data"),
        Input(component_id="pipe-roughness-select", component_property="value"),
        Input(component_id="fluid-mode-select", component_property="value"),
        Input(component_id="inlet-pressure-select", component_property="value"),
        Input(component_id="slider-values-store", component_property="data"),
    ],
)
def build_plot_params(
    interp_time, fluid, case, mdot, L2, L1, grad, diameter, Tinj, k,
    radio3, model, Tsurf, c_m, rho_m, dia_vert, dia_lat, n_lat,
    lateral_flow, insulation_k, lateral_mult, mesh, accuracy, mass_mode, temp_mode,
    pipe_roughness, fluid_mode, inlet_pressure, slider_store
):
    
    if is_print:
        print("build_plot_params")
    triggered = ctx.triggered if ctx.triggered else []
    triggered_ids = [t["prop_id"].split(".")[0] for t in triggered] if triggered else []
    
    if case == "coaxial" and model in ["SBT V1.0", "SBT V2.0"]:
        param4_value = insulation_k if insulation_k is not None else lateral_flow
    else:
        param4_value = lateral_flow
    
    return {
        "vals": (
            interp_time, fluid, case, mdot, L2, L1, grad, diameter, Tinj, k,
            radio3, model, Tsurf, c_m, rho_m, dia_vert, dia_lat, n_lat,
            param4_value, lateral_mult, mesh, accuracy, mass_mode, temp_mode,
            pipe_roughness, fluid_mode, inlet_pressure,
        ),
        "slider_store": slider_store,
    }


@app.callback(
    Output(component_id="plot-run-trigger", component_property="data"),
    [
        Input(component_id="plot-params-store", component_property="data"),
    ],
    prevent_initial_call=True,
)
def trigger_plot_run(params):
    if is_print:
        print("trigger_plot_run")
    return time.time()


@app.callback(
    [
        Output(component_id="geothermal_time_plots", component_property="figure"),
        Output(component_id="thermal-memory", component_property="data"),
        Output(component_id="thermal-results-mass", component_property="data"),
        Output(component_id="thermal-results-time", component_property="data"),
        Output(component_id="thermal-results-errors", component_property="data"),
        Output(component_id="TandP-data", component_property="data"),
        Output(component_id="calculation-request-id", component_property="data"),
        Output(component_id="plot-inputs-cache", component_property="data"),
    ],
    [
        Input(component_id="plot-run-trigger", component_property="data"),
    ],
    [
        State(component_id="plot-params-store", component_property="data"),
        State(component_id="plot-inputs-cache", component_property="data"),
        State(component_id="calculation-request-id", component_property="data"),
        State(component_id="slider-values-store", component_property="data"),
    ],
    prevent_initial_call=True,
)
def update_subsurface_results_plots(
    _trigger,
    params_store,
    plot_inputs_cache,
    current_request_id,
    store_data,
):
    # -----------------------------------------------------------------------------
    # Creates and displays Plotly subplots of the subsurface results.
    # -----------------------------------------------------------------------------
    if is_print:
        print("update_subsurface_results_plots")
    if not params_store:
        return (no_update,) * 8

    vals = params_store["vals"]
    (
        interp_time, fluid, case, mdot, L2, L1, grad, D, Tinj, k_m,
        scale, model, Tsurf, c_m, rho_m, Diameter1, Diameter2, PipeParam3,
        PipeParam4, PipeParam5, mesh, accuracy, mass_mode, temp_mode,
        pipe_roughness, fluid_mode, inlet_pressure,
    ) = vals
    
    if model == "SBT V2.0":
        HyperParam1 = inlet_pressure
    else:
        HyperParam1 = mass_mode
    HyperParam3 = temp_mode
    HyperParam5 = fluid_mode

    # Convert pipe roughness from µm (UI) to meters (model)
    # pipe_roughness slider is in µm for the UI (1-3), model expects meters (1e-6 to 3e-6)
    if pipe_roughness is None:
        pipe_roughness = 1  # Default to 1 µm
    pipe_roughness_m = pipe_roughness * 1e-6
    
    # For SBT V2.0: HyperParam1 = Inlet Pressure, HyperParam3 = Pipe Roughness
    # For SBT V1.0: HyperParam1 = Mass Flow Rate Mode, HyperParam3 = Injection Temperature Mode
    if model == "SBT V2.0":
        # Ensure HyperParam1 (Inlet Pressure) is a float
        try:
            hyperparam1_value = float(HyperParam1) if HyperParam1 is not None else 20.0
        except (TypeError, ValueError):
            hyperparam1_value = 20.0
        # Use converted pipe roughness value in meters
        hyperparam3_value = pipe_roughness_m
    else:
        hyperparam1_value = HyperParam1
        hyperparam3_value = HyperParam3

    try:
        # Detect if multiple slider inputs changed at once (batch update from restore_slider_values)
        # If so, check cache first to avoid multiple runs
        triggered = ctx.triggered if ctx.triggered else []
        triggered_ids = [t["prop_id"].split(".")[0] for t in triggered] if triggered else []
        slider_inputs = ["mdot-select", "L2-select", "L1-select", "grad-select", "diameter-select", "Tinj-select", "k-select"]
        triggered_sliders = [t["prop_id"] for t in triggered if any(slider in t["prop_id"] for slider in slider_inputs)]
        
        # Canonicalization function to normalize numeric types
        def canon(v, nd=6):
            if v is None:
                return None
            if isinstance(v, bool):
                return v
            if isinstance(v, (int, float)):
                return round(float(v), nd)
            return v
        
        # Build current_inputs and canonicalize immediately
        current_inputs_raw = {
            "interp_time": interp_time,
            "fluid": fluid,
            "case": case,
            "mdot": mdot,
            "L2": L2,
            "L1": L1,
            "grad": grad,
            "D": D,
            "Tinj": Tinj,
            "k_m": k_m,
            "scale": scale,
            "model": model,
            "Tsurf": Tsurf,
            "c_m": c_m,
            "rho_m": rho_m,
            "Diameter1": Diameter1,
            "Diameter2": Diameter2,
            "PipeParam3": PipeParam3,
            "PipeParam4": PipeParam4,
            "PipeParam5": PipeParam5,
            "coaxial_flow_type": PipeParam5 if (case == "coaxial" and PipeParam5 is not None) else None,
            "mesh": mesh,
            "accuracy": accuracy,
            "HyperParam1": hyperparam1_value,
            "HyperParam3": hyperparam3_value,
            "pipe_roughness": pipe_roughness,
            "HyperParam5": HyperParam5,
        }
        current_inputs = {k: canon(v) for k, v in current_inputs_raw.items()}
        
        use_plot_cache = False
        
        if use_plot_cache and plot_inputs_cache and plot_inputs_cache.get("inputs"):
            cached_inputs_raw = plot_inputs_cache.get("inputs", {})
            cached_inputs = {k: canon(v) for k, v in cached_inputs_raw.items()}
            # Ensure model matches - don't return cache from different model
            if current_inputs == cached_inputs and cached_inputs.get("model") == model:
                cached_pack = plot_inputs_cache.get("outputs")
                if cached_pack and len(cached_pack) == 6:
                    cached_fig_dict, cached_mem, cached_mass, cached_time, cached_errs, cached_tandp = cached_pack
                    if cached_fig_dict and isinstance(cached_fig_dict, dict):
                        figure_data = cached_fig_dict.get("data", [])
                        if figure_data and len(figure_data) > 0:
                            fig = go.Figure(cached_fig_dict)
                            uirev = f"{L1}-{grad}-{Tinj}-{mdot}-{scale}-{model}-{case}-{fluid}-{current_request_id}"
                            fig.update_layout(
                                uirevision=uirev,
                                datarevision=current_request_id if current_request_id is not None else time.time()
                            )
                            new_fig_dict = fig.to_dict()
                            outputs6 = (new_fig_dict, cached_mem, cached_mass, cached_time, cached_errs, cached_tandp)
                            new_cache = {"inputs": current_inputs, "outputs": outputs6}
                            request_id = current_request_id if current_request_id is not None else 0
                            return (*outputs6, request_id, new_cache)
        
        (
            subplots,
            forty_yr_TPmeans_dict,
            df_mass_flow_rate,
            df_time,
            err_subres_dict,
            TandP_dict,
        ) = generate_subsurface_lineplots(
            interp_time,
            fluid,
            case,
            mdot,
            L2,
            L1,
            grad,
            D,
            Tinj,
            k_m,
            scale,
            model,
            Tsurf,
            c_m,
            rho_m,
            # diameter_vertical, diameter_lateral, n_laterals, lateral_flow, lateral_multiplier,
            Diameter1,
            Diameter2,
            PipeParam3,
            PipeParam4,
            PipeParam5,
            mesh,
            accuracy,
            hyperparam1_value,
            hyperparam3_value,
            HyperParam5,
        )

        # Build outputs6 explicitly: (figure, thermal_memory, mass, time, errors, TandP)
        outputs6 = (
            subplots,
            forty_yr_TPmeans_dict,
            df_mass_flow_rate,
            df_time,
            err_subres_dict,
            TandP_dict,
        )
        
        if subplots and isinstance(subplots, dict):
            figure_data = subplots.get("data", [])
            if not figure_data or len(figure_data) == 0:
                # Figure is empty, return empty figure to show blank state
                empty_fig = make_subplots(rows=2, cols=3)
                empty_outputs6 = (empty_fig.to_dict(), {}, {}, {}, {}, {})
                empty_cache = {"inputs": current_inputs, "outputs": empty_outputs6}
                request_id = current_request_id if current_request_id is not None else 0
                return (*empty_outputs6, request_id, empty_cache)
            
            # Force plotly.js to treat this as a new render every callback
            uirev = f"{L1}-{grad}-{Tinj}-{mdot}-{scale}-{model}-{case}-{fluid}-{current_request_id}"
            subplots.setdefault("layout", {})
            subplots["layout"]["uirevision"] = uirev
            subplots["layout"]["datarevision"] = current_request_id if current_request_id is not None else time.time()
            # Update outputs6 with modified subplots
            outputs6 = (subplots, forty_yr_TPmeans_dict, df_mass_flow_rate, df_time, err_subres_dict, TandP_dict)
        
        # Build cache with canonicalized inputs and outputs6
        new_cache = {
            "inputs": current_inputs,
            "outputs": outputs6
        }
        
        request_id = current_request_id if current_request_id is not None else 0
        return (*outputs6, request_id, new_cache)
    except PreventUpdate:
        # Re-raise PreventUpdate - it's not an error, it's a signal to skip updating
        raise
    except Exception as e:
        traceback.print_exc()
        # Return empty/default values on error
        empty_fig = make_subplots(rows=2, cols=3)
        empty_outputs6 = (empty_fig.to_dict(), {}, {}, {}, {}, {})
        error_inputs = current_inputs if 'current_inputs' in locals() else {}
        empty_cache = {"inputs": error_inputs, "outputs": empty_outputs6}
        request_id = current_request_id if current_request_id is not None else 0
        return (*empty_outputs6, request_id, empty_cache)






@app.callback(
    [
        Output(component_id="geothermal_plots", component_property="figure"),
        Output(component_id="thermal-contours-errors", component_property="data"),
    ],
    [
        Input(component_id="interpolation-select", component_property="value"),
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
    [
        State(component_id="plot-params-store", component_property="data"),
    ],
    prevent_initial_call=True,
)
def update_subsurface_contours_plots(
    interp_time, fluid, case, param, mdot, L2, L1, grad, D, Tinj, k_m, params_store
):
    if is_print:
        print("update_subsurface_contours_plots")
    # -----------------------------------------------------------------------------
    # Creates and displays Plotly subplots of the subsurface contours.
    # -----------------------------------------------------------------------------
    if not params_store or "vals" not in params_store:
        raise PreventUpdate
    
    vals = params_store["vals"]
    model = vals[11]
    
    # Contours are only available for HDF5 models, not SBT models
    if model is None or model in ["SBT V1.0", "SBT V2.0"]:
        raise PreventUpdate
    
    # Check for None values on initial load
    if fluid is None or case is None or param is None or mdot is None or L2 is None or L1 is None or Tinj is None or D is None or grad is None or k_m is None:
        raise PreventUpdate
    
    subplots, err_subcontour_dict = generate_subsurface_contours(
        interp_time, fluid, case, param, mdot, L2, L1, grad, D, Tinj, k_m
    )

    return subplots, err_subcontour_dict


@app.callback(
    [
        Output(component_id="econ_plots", component_property="figure"),
        Output(component_id="econ-memory", component_property="data"),
        Output(component_id="econ-results", component_property="data"),
        Output(component_id="econ-errors", component_property="data"),
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
        Input(component_id="radio-graphic-control4", component_property="value"),
        Input(component_id="checklist", component_property="value"),
        Input(component_id="model-select", component_property="value"),
        # SBT-specific inputs to trigger recalculation
        Input(component_id="Tsurf-select", component_property="value"),
        Input(component_id="c-select", component_property="value"),
        Input(component_id="rho-select", component_property="value"),
        Input(component_id="diameter-vertical-select", component_property="value"),
        Input(component_id="diameter-lateral-select", component_property="value"),
        Input(component_id="n-laterals-select", component_property="value"),
        Input(component_id="lateral-flow-select", component_property="value"),
        Input(component_id="lateral-multiplier-select", component_property="value"),
        Input(component_id="mesh-select", component_property="value"),
        Input(component_id="accuracy-select", component_property="value"),
        Input(component_id="mass-mode-select", component_property="data"),
        Input(component_id="temp-mode-select", component_property="data"),
        Input(component_id="pipe-roughness-select", component_property="value"),
        Input(component_id="fluid-mode-select", component_property="value"),
        Input(component_id="inlet-pressure-select", component_property="value"),
    ],
)
def update_econ_plots(
    TandP_dict,
    interp_time,
    fluid,
    case,
    end_use,
    mdot,
    L2,
    L1,
    grad,
    D,
    Tinj,
    k_m,
    Drilling_cost_per_m,
    Discount_rate,
    Lifetime,
    Direct_use_heat_cost_per_kWth,
    Power_plant_cost_per_kWe,
    Pre_Cooling_Delta_T,
    Turbine_outlet_pressure,
    scale,
    checklist,
    model,
    Tsurf,
    c_m,
    rho_m,
    Diameter1,
    Diameter2,
    PipeParam3,
    PipeParam4,
    PipeParam5,  # For U-tube: lateral multiplier, For coaxial: flow type (string from dropdown)
    mesh,
    accuracy,
    mass_mode,
    temp_mode,
    pipe_roughness,
    fluid_mode,
    inlet_pressure,
):
    if is_print:
        print("update_econ_plots")
    try:
        # Handle None values
        if TandP_dict is None:
            TandP_dict = {}
        if checklist is None:
            checklist = []
        
        # Convert pipe roughness from µm (UI) to meters (model)
        # pipe_roughness slider is in µm for the UI (1-3), model expects meters (1e-6 to 3e-6)
        if pipe_roughness is None:
            pipe_roughness = 1  # Default to 1 µm
        pipe_roughness_m = pipe_roughness * 1e-6
        
        # -----------------------------------------------------------------------------
        # Creates and displays Plotly subplots of the economic results.
        # -----------------------------------------------------------------------------
        
        if checklist == [" "]:
            is_plot_ts = True
        else:
            is_plot_ts = False
        # For SBT V2.0: HyperParam1 = Inlet Pressure (in MPa), HyperParam3 = Pipe Roughness
        # For SBT V1.0: HyperParam1 = Mass Flow Rate Mode, HyperParam3 = Injection Temperature Mode
        if model == "SBT V2.0":
            # Ensure HyperParam1 (Inlet Pressure) is a float (in MPa)
            try:
                hyperparam1_value = float(inlet_pressure) if inlet_pressure is not None else 10.0
            except (TypeError, ValueError):
                hyperparam1_value = 10.0
            # Use converted pipe roughness value in meters
            hyperparam3_value = pipe_roughness_m
        else:
            hyperparam1_value = None  # Not used for HDF5 or SBT V1.0
            hyperparam3_value = temp_mode
        
        economics_fig, econ_data_dict, econ_values_dict, err_econ_dict = (
            generate_econ_lineplots(
                TandP_dict,
                interp_time,
                case,
                end_use,
                fluid,
                mdot,
                L2,
                L1,
                grad,
                D,
                Tinj,
                k_m,
                Drilling_cost_per_m,
                Discount_rate,
                Lifetime,
                Direct_use_heat_cost_per_kWth,
                Power_plant_cost_per_kWe,
                Pre_Cooling_Delta_T,
                Turbine_outlet_pressure,
                scale,
                properties_H2O_pathname,
                properties_CO2v2_pathname,
                additional_properties_CO2v2_pathname,
                tmatrix_pathname,
                model,
                is_plot_ts_check=is_plot_ts,
                HyperParam1=hyperparam1_value,
            )
        )

        return economics_fig, econ_data_dict, econ_values_dict, err_econ_dict
    except Exception as e:
        traceback.print_exc()
        # Return empty figure and empty data on error
        empty_fig = go.Figure()
        return empty_fig, {}, {}, {}


@app.callback(
    Output(component_id="ts-text", component_property="style"),
    [
        Input(component_id="fluid-select", component_property="value"),
        Input(component_id="end-use-select", component_property="value"),
        Input(component_id="checklist", component_property="value"),
    ],
)
def update_plot_title(fluid, end_use, checklist):
    if is_print:
        print("update_plot_title")
    if checklist == [" "]:
        is_title_show = True
    else:
        is_title_show = False

    if not is_title_show:
        return {"display": "none"}

    if fluid == "H2O" or end_use == "Heating":
        return {"display": "none"}

    return {"display": "none"}


@app.callback(
    # Output(component_id="table", component_property="figure"),
    [
        Output(component_id="table", component_property="figure"),
        Output(component_id="summary-memory", component_property="data"),
        Output(component_id="sbt-params-store", component_property="data"),
    ],
    [
        Input(component_id="interpolation-select", component_property="value"),
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
        Input(component_id="econ-memory", component_property="data"),
        Input(component_id="thermal-memory", component_property="data"),
        Input(component_id="model-select", component_property="value"),
        Input(component_id="TandP-data", component_property="data"),
        # SBT-specific inputs
        Input(component_id="Tsurf-select", component_property="value"),
        Input(component_id="c-select", component_property="value"),
        Input(component_id="rho-select", component_property="value"),
        Input(component_id="diameter-vertical-select", component_property="value"),
        Input(component_id="diameter-lateral-select", component_property="value"),
        Input(component_id="n-laterals-select", component_property="value"),
        Input(component_id="lateral-flow-select", component_property="value"),
        Input(component_id="lateral-multiplier-select", component_property="value"),
        Input(component_id="mesh-select", component_property="value"),
        Input(component_id="accuracy-select", component_property="value"),
        Input(component_id="mass-mode-select", component_property="data"),
        Input(component_id="temp-mode-select", component_property="data"),
        Input(component_id="pipe-roughness-select", component_property="value"),
        Input(component_id="fluid-mode-select", component_property="value"),
        Input(component_id="inlet-pressure-select", component_property="value"),
    ],
)
def update_table(
    interp_time,
    fluid,
    case,
    mdot,
    L2,
    L1,
    grad,
    D,
    Tinj,
    k,
    Drilling_cost_per_m,
    Discount_rate,
    Lifetime,
    Direct_use_heat_cost_per_kWth,
    Power_plant_cost_per_kWe,
    Pre_Cooling_Delta_T,
    Turbine_outlet_pressure,
    econ_dict,
    thermal_dict,
    model,
    tandp_data,
    Tsurf,
    c_m,
    rho_m,
    Diameter1,
    Diameter2,
    PipeParam3,
    PipeParam4,
    PipeParam5,  # For U-tube: lateral multiplier, For coaxial: flow type (string from dropdown)
    mesh,
    accuracy,
    mass_mode,
    temp_mode,
    pipe_roughness,
    fluid_mode,
    inlet_pressure,
):
    if is_print:
        print("update_table")
    try:
        # Handle None values
        if econ_dict is None:
            econ_dict = {}
        if thermal_dict is None:
            thermal_dict = {}
        if tandp_data is None:
            tandp_data = None
        
        # Convert pipe roughness from µm (UI) to meters (model)
        # pipe_roughness slider is in µm for the UI (1-3), model expects meters (1e-6 to 3e-6)
        if pipe_roughness is None:
            pipe_roughness = 1  # Default to 1 µm
        pipe_roughness_m = pipe_roughness * 1e-6
        
        # For SBT V2.0: HyperParam1 = Inlet Pressure, HyperParam3 = Pipe Roughness
        # For SBT V1.0: HyperParam1 = Mass Flow Rate Mode, HyperParam3 = Injection Temperature Mode
        if model == "SBT V2.0":
            # Ensure HyperParam1 (Inlet Pressure) is a float
            try:
                hyperparam1_value = float(inlet_pressure) if inlet_pressure is not None else 20.0
            except (TypeError, ValueError):
                hyperparam1_value = 20.0
            # Use converted pipe roughness value in meters
            hyperparam3_value = pipe_roughness_m
        else:
            hyperparam1_value = mass_mode
            hyperparam3_value = temp_mode
        
        # Add TandP data to thermal_dict for SBT models
        if model != "HDF5" and tandp_data:
            thermal_dict["TandP-data"] = tandp_data

        sbt_params = {}
        if model in ["SBT V1.0", "SBT V2.0"]:
            lateral_flow_display = PipeParam4
            if isinstance(PipeParam4, list):
                lateral_flow_display = ", ".join([f"{x:.4f}" for x in PipeParam4])
            
            sbt_params = {
                "model": model,
                "Tsurf": Tsurf if Tsurf is not None else "-",
                "c_m": c_m if c_m is not None else "-",
                "rho_m": rho_m if rho_m is not None else "-",
                "Diameter1": Diameter1 if Diameter1 is not None else "-",
                "Diameter2": Diameter2 if Diameter2 is not None else "-",
                "n_laterals": PipeParam3 if (case == "utube" and PipeParam3 is not None) else "-",
                "lateral_flow_allocation": lateral_flow_display if case == "utube" else "-",
                "lateral_multiplier": PipeParam5 if (case == "utube" and PipeParam5 is not None) else "-",
                "center_pipe_thickness": PipeParam3 if (case == "coaxial" and PipeParam3 is not None) else "-",
                "insulation_thermal_conductivity": PipeParam4 if (case == "coaxial" and PipeParam4 is not None) else "-",
                "coaxial_flow_type": PipeParam5 if (case == "coaxial" and PipeParam5 is not None) else "-",
                "mesh": mesh if mesh is not None else "-",
                "accuracy": accuracy if accuracy is not None else "-",
                "pipe_roughness": pipe_roughness if pipe_roughness is not None else "-",
                "HyperParam1": hyperparam1_value if hyperparam1_value is not None else "-",
                "HyperParam3": hyperparam3_value if hyperparam3_value is not None else "-",
                "HyperParam5": fluid_mode if fluid_mode is not None else "-",
                "case": case if case is not None else "-",
                "fluid": fluid if fluid is not None else "-",
            }

        tbl, summary_dict = generate_summary_table(
            mdot,
            L2,
            L1,
            grad,
            D,
            Tinj,
            k,
            Drilling_cost_per_m,
            Discount_rate,
            Lifetime,
            Direct_use_heat_cost_per_kWth,
            Power_plant_cost_per_kWe,
            Pre_Cooling_Delta_T,
            Turbine_outlet_pressure,
            interp_time,
            case,
            fluid,
            model,
            thermal_dict,
            econ_dict,
        )

        return tbl, summary_dict, sbt_params
    except Exception as e:
        print(f"Error in update_table: {e}")
        # Return empty table and empty summary on error
        empty_table = go.Figure()
        return empty_table, {}, {}


@app.callback(
    [
        Output(component_id="error_block_div1", component_property="children"),
        Output(component_id="error_block_div2", component_property="children"),
        Output(component_id="error_block_div3", component_property="children"),
    ],
    [
        Input(component_id="thermal-results-errors", component_property="data"),
        Input(component_id="thermal-contours-errors", component_property="data"),
        Input(component_id="econ-errors", component_property="data"),
        Input(component_id="econ-results", component_property="data"),
    ],
)
def update_error_divs(err_sub_dict, err_contour_dict, err_econ_dict, econ_results_dict):
    if is_print:
        print("update_error_divs")
    # Handle None values
    if err_sub_dict is None:
        err_sub_dict = {}
    if err_contour_dict is None:
        err_contour_dict = {}
    if err_econ_dict is None:
        err_econ_dict = {}
    if econ_results_dict is None:
        econ_results_dict = {}

    err_div1 = html.Div(  # id="error_block_div1",
        style={"display": "none"}
    )

    err_div2 = html.Div(  # id="error_block_div2",
        style={"display": "none"}
    )

    err_div3 = html.Div(  # id="error_block_div3",
        style={"display": "none"}
    )

    error_style = {
        "display": "block",
        "width": "100%",
        "paddingTop": "8px",
        "paddingLeft": "20px",
        "paddingRight": "20px",
        "paddingBottom": "5px",
        "backgroundColor": lightbrown,
        "color": darkergrey,
    }
    if err_sub_dict != {}:
        try:
            error_message = next(iter(err_sub_dict.values()))
            if error_message is None:
                error_message = "Unknown error occurred"
        except (StopIteration, TypeError):
            error_message = "Unknown error occurred"

        if "No outputs" in error_message:
            error_message = "No outputs were able to be calculated because there are not enough data at these limits. Consider changing parameter value(s)."

        err_div1 = html.Div(  # id="error_block_div1",
            style=error_style,
            children=[
                html.Img(id="error-img1", src=app.get_asset_url("error.png")),
                dcc.Markdown(
                    "**Did not plot visual(s).**", style={"display": "inline-block"}
                ),
                html.P(error_message),
            ],
        )

    if err_contour_dict != {}:
        try:
            error_message = next(iter(err_contour_dict.values()))
            if error_message is None:
                error_message = "Unknown error occurred"
        except (StopIteration, TypeError):
            error_message = "Unknown error occurred"

        err_div2 = html.Div(  # id="error_block_div2",
            style=error_style,
            children=[
                html.Img(id="error-img2", src=app.get_asset_url("error.png")),
                dcc.Markdown(
                    "**Did not plot visual(s).**", style={"display": "inline-block"}
                ),
                html.P(error_message),
            ],
        )

    if err_econ_dict != {}:
        try:
            # Check if plots are actually being generated by looking at econ_results_dict
            plots_have_data = False
            if econ_results_dict and isinstance(econ_results_dict, dict):
                # Check if there are meaningful economic values (not just "-")
                lcoe_sco2 = econ_results_dict.get("LCOE sCO2", "")
                lcoe_h2o = econ_results_dict.get("LCOE H2O", "")
                lcoh_sco2 = econ_results_dict.get("LCOH sCO2", "")
                lcoh_h2o = econ_results_dict.get("LCOH H2O", "")

                # If any of these have meaningful values (not "-"), then plots are being generated
                plots_have_data = (
                    lcoe_sco2 != "-"
                    or lcoe_h2o != "-"
                    or lcoh_sco2 != "-"
                    or lcoh_h2o != "-"
                )

            # Only show warning if plots are NOT being generated
            if not plots_have_data:
                error_message = next(iter(err_econ_dict.values()))
                if error_message and error_message.strip():
                    if "object has no attribute" in error_message:
                        error_message = "No outputs could be calculated because there is not enough data at these limits. Consider changing parameter value(s)."
                    err_div3 = html.Div(  # id="error_block_div3",
                        style=error_style,
                        children=[
                            html.Img(
                                id="error-img3", src=app.get_asset_url("error.png")
                            ),
                            dcc.Markdown(
                                "**Did not plot visual(s).**",
                                style={"display": "inline-block"},
                            ),
                            html.P(
                                "No outputs could be calculated because there is not enough data at these limits. Consider changing parameter value(s)."
                            ),
                        ],
                    )
                else:
                    err_div3 = html.Div(style={"display": "none"})
            else:
                # Plots are being generated, don't show warning
                err_div3 = html.Div(style={"display": "none"})

        except Exception as e:
            err_div3 = html.Div(style={"display": "none"})

    return err_div1, err_div2, err_div3


@app.callback(
    Output(component_id="warning_block_div3", component_property="children"),
    [
        Input(component_id="econ-memory", component_property="data"),
        Input(component_id="econ-results", component_property="data"),
        Input(component_id="fluid-select", component_property="value"),
        Input(component_id="end-use-select", component_property="value"),
    ],
)
def update_warning_divs(levelized_cost_dict, econ_results_dict, fluid, end_use):
    if is_print:
        print("update_warning_divs")
    warning_div3 = html.Div(style={"display": "none"})

    # Handle None or missing data
    if levelized_cost_dict is None:
        return warning_div3
    
    if econ_results_dict is None:
        econ_results_dict = {}
        
    error_style = {
        "display": "block",
        "width": "100%",
        "paddingTop": "8px",
        "paddingLeft": "20px",
        "paddingRight": "20px",
        "paddingBottom": "5px",
        "backgroundColor": lightbrown,
        "color": darkergrey,
    }

    # Safely check for LCOE and LCOH values
    lcoe_sco2 = levelized_cost_dict.get("LCOE sCO2", "-")
    lcoe_h2o = levelized_cost_dict.get("LCOE H2O", "-")
    lcoh_sco2 = levelized_cost_dict.get("LCOH sCO2", "-")
    lcoh_h2o = levelized_cost_dict.get("LCOH H2O", "-")
    
    # Helper function to check if LCOE value is valid (numeric and reasonable)
    def is_valid_lcoe(value):
        if value is None or value == "-" or value == "Insufficient Inputs":
            return False
        if value == "9999.00" or value == "9999":
            return False
        try:
            # Check if it's a valid number
            num_val = float(value)
            # Valid if it's a positive number less than 9999
            return 0 < num_val < 9999
        except (ValueError, TypeError):
            return False
    
    # Helper function to check if LCOH value is valid (numeric and reasonable)
    def is_valid_lcoh(value):
        if value is None or value == "-" or value == "Insufficient Inputs":
            return False
        if value == "9999.00" or value == "9999":
            return False
        try:
            # Check if it's a valid number
            num_val = float(value)
            # Valid if it's a positive number less than 9999
            return 0 < num_val < 9999
        except (ValueError, TypeError):
            return False
    
    # Check if we have any valid LCOE values
    has_valid_lcoe = is_valid_lcoe(lcoe_sco2) or is_valid_lcoe(lcoe_h2o)
    
    # Check if we have any valid LCOH values
    has_valid_lcoh = is_valid_lcoh(lcoh_sco2) or is_valid_lcoh(lcoh_h2o)
    
    # Check if a calculation has been attempted (to avoid showing warning on initial load)
    # A calculation has been attempted if:
    # 1. econ_results_dict exists and has data, OR
    # 2. levelized_cost_dict has values that aren't just the default "-" (including "Insufficient Inputs")
    calculation_attempted = False
    if econ_results_dict and isinstance(econ_results_dict, dict) and len(econ_results_dict) > 0:
        calculation_attempted = True
    elif levelized_cost_dict and isinstance(levelized_cost_dict, dict):
        # Check if we have any non-default values (not just "-")
        # "Insufficient Inputs" means a calculation was attempted, so count it
        has_non_default = False
        for key, value in levelized_cost_dict.items():
            if value != "-" and value is not None:
                has_non_default = True
                break
        calculation_attempted = has_non_default
    
    # Check which LCOE values are invalid
    lcoe_sco2_is_insufficient = lcoe_sco2 == "Insufficient Inputs"
    lcoe_h2o_is_insufficient = lcoe_h2o == "Insufficient Inputs"
    lcoe_sco2_is_other_invalid = (
        lcoe_sco2 == "9999.00" or lcoe_sco2 == "-"
    )
    lcoe_h2o_is_other_invalid = (
        lcoe_h2o == "9999.00" or lcoe_h2o == "-"
    )
    
    lcoe_sco2_is_invalid = lcoe_sco2_is_insufficient or (lcoe_sco2_is_other_invalid and calculation_attempted)
    lcoe_h2o_is_invalid = lcoe_h2o_is_insufficient or (lcoe_h2o_is_other_invalid and calculation_attempted)
    
    # Determine which LCOE to check based on selected fluid
    if fluid == "H2O":
        should_check_lcoe = lcoe_h2o_is_invalid
    elif fluid == "sCO2":
        should_check_lcoe = lcoe_sco2_is_invalid
    else:  # fluid == "All"
        should_check_lcoe = lcoe_sco2_is_invalid or lcoe_h2o_is_invalid
    
    # Only show LCOE warning if end-use is "Electricity" or "All"
    # If end-use is "Heating", LCOE is not relevant, so don't show the warning
    if end_use == "Heating":
        return warning_div3  # Return empty warning div for heating end-use
    
    # Only show LCOH warning if end-use is "Heating" or "All"
    # If end-use is "Electricity", LCOH is not relevant, so don't show the warning
    # (Note: Currently this callback only shows LCOE warnings, but adding this check
    # for consistency and future-proofing in case LCOH warnings are added)
    if end_use == "Electricity":
        # Check if there are LCOH-related issues that should be ignored
        # For now, this callback only handles LCOE, but we return early to be consistent
        pass  # Continue to check LCOE warnings below
    
    # Show warning if the relevant LCOE is invalid and a calculation was attempted
    # This helps users understand why a specific fluid's LCOE isn't calculating
    if should_check_lcoe and calculation_attempted:
        # Different messages for different error types - provide actionable guidance
        # Check error codes to provide specific guidance, otherwise use generic "Insufficient Inputs" message
        error_codes = econ_results_dict.get("error_codes", [])
        
        # Check for specific error codes first
        has_error_6000 = 6000 in error_codes  # Zero electricity production
        has_error_7000 = 7000 in error_codes  # Negative LCOE
        has_error_1000 = 1000 in error_codes  # Production temp below injection temp
        has_error_2000 = 2000 in error_codes  # H2O temp outside ORC range
        has_error_3000 = 3000 in error_codes  # CO2 injection temp too low
        has_error_4000 = 4000 in error_codes  # CO2 turbine outlet temp too low
        
        if has_error_6000:
            # Zero electricity production - can't calculate LCOE
            warning_div3 = html.Div(
                style=error_style,
                children=[
                    html.Img(id="warning-img", src=app.get_asset_url("warning.png")),
                    dcc.Markdown(
                        "**Cannot calculate LCOE - zero electricity production.**", style={"display": "inline-block"}
                    ),
                    html.P(
                        "To get valid LCOE values, try increasing: Geothermal Gradient, Depth (L1), or Injection Temperature to raise the Outlet Temperature."
                    ),
                ],
            )
        elif has_error_7000:
            # Negative LCOE - system is losing energy
            warning_div3 = html.Div(
                style=error_style,
                children=[
                    html.Img(id="warning-img", src=app.get_asset_url("warning.png")),
                    dcc.Markdown(
                        "**Cannot calculate LCOE - negative net electricity production.**", style={"display": "inline-block"}
                    ),
                    html.P(
                        "The system is consuming more energy than it produces. Try increasing: Geothermal Gradient, Depth (L1), or Injection Temperature to increase Outlet Temperature and improve energy production. Alternatively, try reducing Mass Flow Rate to decrease pumping power requirements."
                    ),
                ],
            )
        elif has_error_1000:
                    # Production temperature dropped below injection temperature
                    warning_div3 = html.Div(
                        style=error_style,
                        children=[
                            html.Img(id="warning-img", src=app.get_asset_url("warning.png")),
                            dcc.Markdown(
                                "**Cannot calculate LCOE - production temperature too low.**", style={"display": "inline-block"}
                            ),
                            html.P(
                                "Production temperature dropped below injection temperature. Try increasing: Geothermal Gradient, Depth (L1), or reducing Mass Flow Rate to improve heat extraction."
                            ),
                        ],
                    )
        elif has_error_2000:
            # H2O temperature outside ORC range
            warning_div3 = html.Div(
                style=error_style,
                children=[
                    html.Img(id="warning-img", src=app.get_asset_url("warning.png")),
                    dcc.Markdown(
                        "**Cannot calculate LCOE - temperature outside ORC range.**", style={"display": "inline-block"}
                    ),
                    html.P(
                        "For H2O, injection temperature must be above 50°C and production temperature must be between 100-200°C. Try adjusting: Injection Temperature, Geothermal Gradient, or Depth (L1)."
                    ),
                ],
            )
        elif has_error_3000:
            # CO2 injection temperature too low
            warning_div3 = html.Div(
                style=error_style,
                children=[
                    html.Img(id="warning-img", src=app.get_asset_url("warning.png")),
                    dcc.Markdown(
                        "**Cannot calculate LCOE - CO2 injection temperature too low.**", style={"display": "inline-block"}
                    ),
                    html.P(
                        "CO2 injection temperature must be at least 32°C to remain supercritical. Increase Injection Temperature to at least 32°C."
                    ),
                ],
            )
        elif has_error_4000:
            # CO2 turbine outlet temperature too low
            warning_div3 = html.Div(
                style=error_style,
                children=[
                    html.Img(id="warning-img", src=app.get_asset_url("warning.png")),
                    dcc.Markdown(
                        "**Cannot calculate LCOE - CO2 turbine outlet temperature too low.**", style={"display": "inline-block"}
                    ),
                    html.P(
                        "CO2 turbine outlet temperature dropped below 37°C. Try increasing: Geothermal Gradient, Depth (L1), or Injection Temperature to raise the production temperature."
                    ),
                ],
            )
        else:
            # Generic message for cases without specific error codes
            # Use "Insufficient Inputs" message for consistency
            warning_div3 = html.Div(
                style=error_style,
                children=[
                    html.Img(id="warning-img", src=app.get_asset_url("warning.png")),
                    dcc.Markdown(
                        "**Insufficient inputs to calculate LCOE.**", style={"display": "inline-block"}
                    ),
                    html.P(
                        "To get valid LCOE values, try increasing: Outlet Temperature (increase Geothermal Gradient, Depth, or Injection Temperature), or reducing Mass Flow Rate."
                    ),
                ],
            )

    return warning_div3


# -----------------------------------------------------------------------------
# App runs here. Define configurations, proxies, etc.
# -----------------------------------------------------------------------------

server = app.server
# from app import server as application # in the wsgi.py file -- this targets the Flask server of Dash app

compress = Compress()
compress.init_app(app.server)  # gzip all static assets

app.server.config.update(
    SESSION_COOKIE_SECURE=True,
    SESSION_COOKIE_HTTPONLY=True,
    SESSION_COOKIE_SAMESITE="Lax",
)

"""" 
1. Dash/Flask code no longer issuing a Set-Cookie header.
2. No client-side script (such as the old Google-Analytics snippet) is writing a cookie anymore.
"""


@server.route("/.well-known/apple-app-site-association")
def aasa():
    return send_from_directory(
        "static",
        "apple-app-site-association",  # apple-app-site-association is purely a hint to iOS / iPadOS about when it may open a native app instead of Safari
        mimetype="application/json",
    )


if __name__ == "__main__":
    # test change.
    # app.run_server(port=8060, debug=True)
    app.run(
        # host="127.0.0.1",
        port=8060,
        debug=True,  # needs to be False in production
        # ssl_context="adhoc",
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
#  Output(component_id='diameter-vertical-select-div', component_property='style'),
#  Output(component_id='diameter-lateral-select-div', component_property='style'),
#  Output(component_id='n-laterals-select-div', component_property='style'),
#  Output(component_id='lateral-flow-select-div', component_property='style'),
#  Output(component_id='lateral-multiplier-select-div', component_property='style'),

# Clientside callback: resize econ plot when tab becomes active
app.clientside_callback(
    """
    function(tabValue, currentPing) {
      if (tabValue !== "economics-time-tab") {
        return currentPing || 0;
      }
      // Wait a tick so the tab is actually visible / laid out
      setTimeout(function() {
        try {
          const graphDiv = document.getElementById("econ_plots");
          if (!graphDiv) return;
          // The actual plotly graph is inside the dcc.Graph container
          const gd = graphDiv.querySelector(".js-plotly-plot");
          if (gd && typeof Plotly !== 'undefined' && Plotly && Plotly.Plots && Plotly.Plots.resize) {
            Plotly.Plots.resize(gd);
          }
        } catch(e) {
          console.log("econ resize error", e);
        }
      }, 80);
      return Date.now();
    }
    """,
    Output("econ-resize-ping", "data"),
    Input("tabs", "value"),
    State("econ-resize-ping", "data"),
)

# Clientside callback: resize econ plot when T-S diagram checkbox changes
app.clientside_callback(
    """
    function(val, currentPing) {
      setTimeout(function() {
        try {
          const graphDiv = document.getElementById("econ_plots");
          const gd = graphDiv && graphDiv.querySelector(".js-plotly-plot");
          if (gd && typeof Plotly !== 'undefined' && Plotly && Plotly.Plots && Plotly.Plots.resize) {
            Plotly.Plots.resize(gd);
          }
        } catch(e) {
          console.log("econ resize error (checkbox)", e);
        }
      }, 80);
      return currentPing || 0;
    }
    """,
    Output("econ-resize-ping", "data", allow_duplicate=True),
    Input("checklist", "value"),
    State("econ-resize-ping", "data"),
    prevent_initial_call=True,
)

