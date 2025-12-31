#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# standard libraries
import re

# data manipulation and plotting libraries
import numpy as np
import pandas as pd
import plotly.graph_objects as go

# --------------------------------------
# Colors for the contours.
# --------------------------------------

# colorscaleR = [[0, '#450500'], [0.5, 'red'], [1, 'lightsalmon']]
colorscaleR = [[0, '#2e2112'], [0.5, '#cf7304'], [1, 'lightsalmon']]
# colorscaleY = [[0, '#574200'], [0.5, 'gold'], [1, '#d0f000']]
colorscaleY = [[0, '#574200'], [0.5, 'gold'], [1, '#fff7d4']]
# colorscaleB = [[0, '#011f63'], [0.5, 'blue'], [1, '#02a2f2']]
colorscaleB = [[0, '#042140'], [0.5, '#0554ad'], [1, '#02a2f2']]
# colorscaleG = [[0, '#013801'], [0.5, 'green'], [1, '#07ed99']]
colorscaleG = [[0, '#043336'], [0.5, '#10a4ad'], [1, '#abebdc']]


slider_list = ["Mass Flow Rate (kg/s)", "Horizontal Extent (m)", "Drilling Depth (m)", "Geothermal Gradient (°C/m)",
                "Borehole Diameter (m)", "Injection Temperature (˚C)", "Rock Thermal Conductivity (W/m-°C)", 
                    "Drilling Cost ($/m)", "Discount Rate (%)", "Lifetime (years)", "Plant CAPEX ($/kWt)", 
                    "Plant CAPEX ($/kWe)", "Pre-cooling (˚C)", "Turbine Outlet Pressure (bar)"]



def blank_canvas(fig, row_n, col_n):

    # --------------------------------------------------------------------------------
    # Returns a blank figure with no x and y axes and unresponsive.
    # --------------------------------------------------------------------------------

    fig.update_xaxes(showticklabels=False, row=row_n, col=col_n)
    fig.update_yaxes(showticklabels=False, row=row_n, col=col_n)

    # fig.update_layout( xaxis =  { "visible": False }, yaxis = { "visible": False },
    #                     # annotations = [
    #                     #     {   
    #                     #         "text": "Please select v and m",
    #                     #         "xref": "paper",
    #                     #         "yref": "paper",
    #                     #         "showarrow": False,
    #                     #         "font": {
    #                     #             "size": 28
    #                     #         }
    #                     #     }],
    #                     row=row_n, col=col_n)


def blank_data():
    
    # ------------------------------
    # Returns blank pandas data.
    # ------------------------------

    sCO2_Tout = pd.Series(dtype=object)
    sCO2_Pout = pd.Series(dtype=object)
    H2O_Tout = pd.Series(dtype=object)
    H2O_Pout = pd.Series(dtype=object)
    sCO2_kWe = pd.Series(dtype=object)
    sCO2_kWt = pd.Series(dtype=object)
    H2O_kWe = pd.Series(dtype=object)
    H2O_kWt = pd.Series(dtype=object)

    return sCO2_Tout, sCO2_Pout, H2O_Tout, H2O_Pout, sCO2_kWe, sCO2_kWt, H2O_kWe,H2O_kWt

def blank_data_kW():

    # ------------------------------
    # Returns blank pandas data.
    # ------------------------------

    blank1 = pd.Series(dtype=object)
    blank2 = pd.Series(dtype=object)
    blank3 = pd.Series(dtype=object)
    blank4 = pd.Series(dtype=object)

    return blank1, blank2, blank3, blank4



def update_blank_econ(fig, nrow1, ncol1):

    fig.add_trace(go.Scatter(x=pd.Series(dtype=object), y=pd.Series(dtype=object),), row=nrow1, col=ncol1)
    blank_canvas(fig=fig, row_n=nrow1, col_n=ncol1)

    blank_val = "-"

    return fig, blank_val

def update_blank_econ2(fig, nrow1, ncol1, nrow2, ncol2):

    fig.add_trace(go.Scatter(x=pd.Series(dtype=object), y=pd.Series(dtype=object),), row=nrow1, col=ncol1)
    blank_canvas(fig=fig, row_n=nrow1, col_n=ncol1)

    fig.add_trace(go.Scatter(x=pd.Series(dtype=object), y=pd.Series(dtype=object),),row=nrow2, col=ncol2)
    blank_canvas(fig=fig, row_n=nrow2, col_n=ncol2)

    blank_val1 = "-"
    blank_val2 = "-"

    return fig, blank_val1 #, blank_val2


def parse_error_message(e, e_name, model=None):

    # print('\t', e)

    if e.__class__.__name__ == "ValueError":
        dim = re.findall(r'\d+', str(e))
        if dim != []:
            # For SBT models, show generic message instead of specific parameter
            if model and model in ["SBT V1.0", "SBT V2.0"]:
                error_message = "Cannot calculate results for these parameter values. Consider adjusting the parameter values."
            else:
                err_param = slider_list[int(dim[0])]
                error_message = f"{err_param} is out of bounds of possible values. Consider changing the value."
            # print('\t', e)
            # print(e_name, error_message)
        else:
            error_message = "Did not plot visual(s).\n\nNo outputs could be calculated because there is not enough data at these limits. Consider changing parameter value(s)."
            # print(e_name, e)
    else:
        error_message = "Did not plot visual(s).\n\nNo outputs could be calculated because there is not enough data at these limits. Consider changing parameter value(s)."
        # print(e_name, e)

    return error_message


def rename_for_contour(input1, input2, input3, input4):

    kWe_avg_flipped = input1
    kWt_avg_flipped = input2
    Tout_flipped = input3
    Pout_flipped = input4

    return kWe_avg_flipped, kWt_avg_flipped, Tout_flipped, Pout_flipped

def find_nearest(array, value):

    # --------------------------------------------------------------------------------
    # Returns the nearest value and nearest index from a given value and array.
    # --------------------------------------------------------------------------------

    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()

    return array[idx], idx



def update_layout_properties_subsurface_results(fig, m_dot, time, plot_scale, is_simulator=False):

    # --------------------------------------------------------------------------------------------------
    # Updates the x-axis, y-axis, and annotation properties of the subsurface results Plotly figure(s).
    # --------------------------------------------------------------------------------------------------

    fig.update_layout(hovermode="x unified", hoverlabel=dict(font_size=12), margin_pad=0)

    if hasattr(time, '__len__') and len(time) > 0:
        if hasattr(time, 'iloc'):  # pandas Series
            time_max = time.iloc[-1]
        else:  # numpy array or list
            time_max = time[-1]
    else:
        time_max = 40.0
    
    if is_simulator:
        fig.update_xaxes(visible=False, row=1, col=1)
        fig.update_xaxes(visible=False, row=1, col=2)
        fig.update_yaxes(visible=False, row=1, col=1)
        fig.update_yaxes(visible=False, row=1, col=2)
    else:
        if hasattr(m_dot, '__len__') and len(m_dot) > 0:
            if hasattr(m_dot, 'iloc'):  # pandas Series
                m_dot_max = m_dot.iloc[-1]
            else:  # numpy array or list
                m_dot_max = m_dot[-1]
        else:
            m_dot_max = 100.0
        
        fig.update_xaxes(title_text="Mass Flow Rate (kg/s)", showgrid=True, range=[0,m_dot_max], 
                            row=1, col=1, tickfont = dict(size=12), title_font=dict(size=14))
        fig.update_xaxes(title_text="Mass Flow Rate (kg/s)", showgrid=True, range=[0,m_dot_max], 
                            row=1, col=2, tickfont = dict(size=12), title_font=dict(size=14))
        
        if plot_scale == 2:
            fig.update_yaxes(title_text="40-Year Average of Exergy (kWe)<br>", range=[-1500,8000],
                                row=1, col=1, tickfont = dict(size=12), title_font=dict(size=14))
            fig.update_yaxes(title_text="40-Year Average of Available <br>Thermal Output (kWt)", range=[-7000,40000],
                                row=1, col=2, tickfont = dict(size=12), title_font=dict(size=14))
        else:
            fig.update_yaxes(title_text="40-Year Average of Exergy (kWe)<br>",
                                row=1, col=1, tickfont = dict(size=12), title_font=dict(size=14))
            fig.update_yaxes(title_text="40-Year Average of Available <br>Thermal Output (kWt)",
                                row=1, col=2, tickfont = dict(size=12), title_font=dict(size=14))
    
    fig.update_xaxes(title_text="Time (years)", showgrid=True, range=[0,time_max], 
                        row=2, col=1, tickfont = dict(size=12), title_font=dict(size=14))
    fig.update_xaxes(title_text="Time (years)", showgrid=True, range=[0,time_max], 
                        row=2, col=2, tickfont = dict(size=12), title_font=dict(size=14))
    fig.update_xaxes(title_text="Time (years)", showgrid=True, range=[0,time_max], 
                        row=2, col=3, tickfont = dict(size=12), title_font=dict(size=14))

    if plot_scale == 2:
        fig.update_yaxes(title_text="Exergy (kWe)", range=[0,13000],
                            row=2, col=1, tickfont = dict(size=12), title_font=dict(size=14))
        fig.update_yaxes(title_text="Outlet Temperature (˚C)", range=[0,600-273.15], 
                            row=2, col=2, tickfont = dict(size=12), title_font=dict(size=14))
        fig.update_yaxes(title_text="Outlet Pressure (MPa)", range=[0,45500000 /1000000],
                            row=2, col=3, tickfont = dict(size=12), title_font=dict(size=14))
    else:
        fig.update_yaxes(title_text="Exergy (kWe)",
                            row=2, col=1, tickfont = dict(size=12), title_font=dict(size=14))
        fig.update_yaxes(title_text="Outlet Temperature (˚C)",
                            row=2, col=2, tickfont = dict(size=12), title_font=dict(size=14))
        fig.update_yaxes(title_text="Outlet Pressure (MPa)",
                            row=2, col=3, tickfont = dict(size=12), title_font=dict(size=14))

    if is_simulator:
        fig.update_layout(
            legend_title_text='Working Fluid',
            legend=dict(y=1.05, yanchor='top'),
            template='none', margin=dict(l=70, r=70, t=50, b=70),
            title_text=f'<b>Output Over Time</b>', title_x=0.5, title_y=0.99,
            font=dict(size=10)
        )
        fig.update_yaxes(domain=[0.25, 0.75], row=2, col=1)
        fig.update_yaxes(domain=[0.25, 0.75], row=2, col=2)
        fig.update_yaxes(domain=[0.25, 0.75], row=2, col=3)
    else:
        fig.update_layout(legend_title_text='Working Fluid',
                          template='none', margin=dict(l=70, r=70, t=35, b=70),)
        fig.update_layout(title_text=f'<b>Output Over Mass Flow Rate</b>', title_x=0.35, title_y=0.99,
                            font=dict(size=10))
        fig.update_layout(annotations=[go.layout.Annotation(
                                        showarrow=False, text=f'<b>Output Over Time</b>',
                                        x=0.35,  y=0.44, xref='paper', yref='paper',
                                        font=dict(size=14))])

    return fig


def update_layout_properties_subsurface_contours(fig, param):

    fig.update_xaxes(title_text="Mass Flow Rate (kg/s)", row=1, col=1)
    fig.update_xaxes(title_text="Mass Flow Rate (kg/s)", row=1, col=2)
    fig.update_xaxes(title_text="Mass Flow Rate (kg/s)", row=2, col=1)
    fig.update_xaxes(title_text="Mass Flow Rate (kg/s)", row=2, col=2)


    # Update yaxis properties
    fig.update_yaxes(title_text=param, row=1, col=1)
    fig.update_yaxes(title_text=param, row=1, col=2)
    fig.update_yaxes(title_text=param, row=2, col=1)
    fig.update_yaxes(title_text=param, row=2, col=2)

    fig.update_layout(margin=dict(l=70, r=70, t=50, b=0)) 

    hover_replace0 = fig.data[0].__dict__['_parent'].__dict__['_data'][0]['hovertemplate'].replace("Y", param)
    fig.data[0].__dict__['_parent'].__dict__['_data'][0]['hovertemplate'] = hover_replace0

    hover_replace1 = fig.data[0].__dict__['_parent'].__dict__['_data'][1]['hovertemplate'].replace("Y", param)
    fig.data[0].__dict__['_parent'].__dict__['_data'][1]['hovertemplate'] = hover_replace1

    hover_replace2 = fig.data[0].__dict__['_parent'].__dict__['_data'][2]['hovertemplate'].replace("Y", param)
    fig.data[0].__dict__['_parent'].__dict__['_data'][2]['hovertemplate'] = hover_replace2

    hover_replace3 = fig.data[0].__dict__['_parent'].__dict__['_data'][3]['hovertemplate'].replace("Y", param)
    fig.data[0].__dict__['_parent'].__dict__['_data'][3]['hovertemplate'] = hover_replace3

    # print(fig.data[0].__dict__['_parent'].__dict__['_data_objs'])

    return fig



def update_layout_properties_econ_results(fig, end_use, plot_scale, is_plot_ts_check=False):

    if end_use == "Electricity":
        row_num = 1
    else:
        row_num = 2

    # Sync Axes
    fig.update_layout(hovermode="x unified", hoverlabel=dict(font_size=12))

    # X Axes
    fig.update_xaxes(title_text="Time (year)", showgrid=True, range=[0, 40], row=row_num, col=1,
                        tickfont = dict(size=12), title_font=dict(size=14))
    fig.update_xaxes(title_text="Time (year)", showgrid=True, range=[0, 40], row=row_num, col=3,
                        tickfont = dict(size=12), title_font=dict(size=14)) 
    fig.update_xaxes(title_text="Time (year)", showgrid=True, range=[0, 40], row=1, col=1, 
                        tickfont = dict(size=12), title_font=dict(size=14))
    fig.update_xaxes(title_text="Time (year)", showgrid=True, range=[0, 40], row=1, col=3,
                        tickfont = dict(size=12), title_font=dict(size=14))

    fig.update_yaxes(title_text="Heat Production (MWt)", # range=[0, 55], 
                        row=1, col=1, tickfont = dict(size=12), title_font=dict(size=14))
    fig.update_yaxes(title_text="Annual Heat Production (GWh)", # range=[0, 500], 
                    row=1, col=3, tickfont = dict(size=12), title_font=dict(size=14))
    fig.update_yaxes(title_text="Electricity Production (MWe)", # range=[-1.25, 7], 
                        row=row_num, col=1, tickfont = dict(size=12), title_font=dict(size=14)) # max(teaobj.Inst_Net_Electricity_production/1e3)+1
    fig.update_yaxes(title_text="Annual Electricity Production (GWe)", #range=[-10, 55], 
                        row=row_num, col=3, tickfont = dict(size=12), title_font=dict(size=14))

    if plot_scale == 2:

        fig.update_yaxes(title_text="Heat Production (MWt)", range=[0, 55], row=1, col=1,
                            tickfont = dict(size=12), title_font=dict(size=14))
        fig.update_yaxes(title_text="Annual Heat Production (GWh)", range=[0, 500], row=1, col=3,
                        tickfont = dict(size=12), title_font=dict(size=14))
        fig.update_yaxes(title_text="Electricity Production (MWe)", range=[-1.25, 7], row=row_num, col=1,
                            tickfont = dict(size=12), title_font=dict(size=14)) # max(teaobj.Inst_Net_Electricity_production/1e3)+1
        fig.update_yaxes(title_text="Annual Electricity Production (GWe)", range=[-10, 55], row=row_num, col=3,
                            tickfont = dict(size=12), title_font=dict(size=14))
    top_margin = 70 if is_plot_ts_check else 50
    figure_height = 1200 if is_plot_ts_check else 800
    layout_dict = {
        'legend_title_text': 'Working Fluid',
        'template': 'none',
        'margin': dict(l=70, r=70, t=top_margin, b=70),
        'legend': dict(
            y=0.98,
            yanchor='top'
        ),
        'height': figure_height
    }
    fig.update_layout(**layout_dict)

    # Subtitles
    if end_use == "All":
        fig.update_layout(title_text=f'<b>End-Use: Heating</b>', title_x=0.35, title_y=1,
                            font=dict(size=10)
                            )
        electricity_title_y = 0.673 if is_plot_ts_check else 0.50
        fig.update_layout(annotations=[go.layout.Annotation(
                                        showarrow=False, text=f'<b>End-Use: Electricity</b>',
                                        x=0.35, y=electricity_title_y, xref='paper', yref='paper',
                                        font=dict(size=14))])

    if end_use == "Heating":
        fig.update_layout(title_text=f'<b>End-Use: Heating</b>', title_x=0.35, title_y=1,
                            font=dict(size=10)
                            )

    if end_use == "Electricity":
        fig.update_layout(title_text=f'<b>End-Use: Electricity</b>', title_x=0.35, title_y=1,
                            font=dict(size=10)
                            )

    return fig


def update_lcoh_lcoe_table(fig, fluid, end_use, lcoh_sCO2, lcoh_H2O, lcoe_sCO2, lcoe_H2O):


    if end_use == "Electricity":
        row_num = 1
    else:
        row_num = 2
    
    if end_use == "Heating" or end_use == "All":

        if fluid == "All":
            fig.add_trace(go.Table(#header=dict(values=['Values']),
                     cells=dict(values=[[f'sCO2 LCOH: {lcoh_sCO2} $/MWh th', f'H2O LCOH: {lcoh_H2O} $/MWh th']])
                     ), row=1, col=5)

        if fluid == "sCO2":
            fig.add_trace(go.Table(
                     cells=dict(values=[[f'sCO2 LCOH: {lcoh_sCO2} $/MWh th']])
                     ), row=1, col=5)

        if fluid == "H2O":
            fig.add_trace(go.Table(
                     cells=dict(values=[[f'H2O LCOH: {lcoh_H2O} $/MWh th']])
                     ), row=1, col=5)


    if end_use == "Electricity" or end_use == "All":
        
        if fluid == "All":

            fig.add_trace(go.Table(
                     cells=dict(values=[[f'sCO2 LCOE: {lcoe_sCO2} $/MWh e', f'H2O LCOE: {lcoe_H2O} $/MWh e']])
                     ), row=row_num, col=5)

        if fluid == "sCO2":
            fig.add_trace(go.Table(
                     cells=dict(values=[[f'sCO2 LCOE: {lcoe_sCO2} $/MWh e']])
                     ), row=row_num, col=5)

        if fluid == "H2O":
            fig.add_trace(go.Table(
                     cells=dict(values=[[f'H2O LCOE: {lcoe_H2O} $/MWh e']])
                     ), row=row_num, col=5)

    return fig


