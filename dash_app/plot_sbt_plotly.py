#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import plotly.graph_objects as go
import numpy as np

def plot_borehole_geometry_plotly(clg_configuration, numberoflaterals, x, y, z, xinj, yinj, zinj, xprod, yprod, zprod, xlat, ylat, zlat):
    
    x = x.flatten()
    y = y.flatten()
    z = z.flatten()

    # Create a 3D figure for the borehole geometry
    fig = go.Figure()

    if clg_configuration == 1:  # co-axial geometry

        # Plot the borehole trajectory
        fig.add_trace(go.Scatter3d(x=x, y=y, z=z, 
                                mode='lines+markers', 
                                line=dict(color='black', width=4), 
                                marker=dict(symbol='circle', size=5)))
        
    elif clg_configuration == 2:  # U-loop geometry

        xinj = xinj.flatten()
        yinj = yinj.flatten()
        zinj = zinj.flatten()

        xprod = xprod.flatten()
        yprod = yprod.flatten()
        zprod = zprod.flatten()

        # xlat = xlat.flatten()
        # ylat = ylat.flatten()
        # zlat = zlat.flatten()

        # Plot injection well
        fig.add_trace(go.Scatter3d(x=xinj, y=yinj, z=zinj, 
                                    mode='lines+markers', 
                                    line=dict(color='blue', width=4), 
                                    marker=dict(symbol='circle', size=5),
                                    name='Injection Well'
                                    ))
        
        # Plot production well
        fig.add_trace(go.Scatter3d(x=xprod, y=yprod, z=zprod, 
                                    mode='lines+markers', 
                                    line=dict(color='red', width=4), 
                                    marker=dict(symbol='circle', size=5),
                                    name='Production Well'
                                    ))
        
        # Plot lateral wells
        for i in range(numberoflaterals):
            fig.add_trace(go.Scatter3d(x=xlat[:, i], y=ylat[:, i], z=zlat[:, i], 
                                       mode='lines+markers', 
                                       line=dict(color='black', width=4), 
                                       marker=dict(symbol='circle', size=5),
                                       name=f'Lateral(s)'
                                       ))
        
    # Set axes ranges
    fig.update_layout(
        scene=dict(
            xaxis=dict(range=[np.min(x) - 200, np.max(x) + 200], title='x (m)'),
            yaxis=dict(range=[np.min(y) - 200, np.max(y) + 200], title='y (m)'),
            zaxis=dict(range=[np.min(z) - 500, 0], title='Depth (m)'),
        ),
        margin=dict(pad=0), #(l=0, r=0, b=0, t=0),
        legend=dict(
            itemsizing='constant'
        )
    )
    # paper_bgcolor='rgba(255,255,255,0.10)', # or 0.40
    #                   plot_bgcolor=''

    RGB = "rgba(212, 163, 110, 1)"
    # "peru" or "sandybrown" can also be used
    fig.update_layout(#plot_bgcolor='rgb(12,163,135)',
                    #   paper_bgcolor='rgba(255,255,255,0.10)', # rgb(12,163,135)
                    #coloraxis={"colorbar": {"x": -0.2, "len": 0.5, "y": 0.8}}, #I think this is for contours
                    scene = dict(
                                xaxis = dict(
                                        backgroundcolor="peru",
                                        gridcolor="white",
                                        showbackground=True,
                                        zerolinecolor="white",),
                                yaxis = dict(
                                    backgroundcolor=RGB,
                                    gridcolor="white",
                                    showbackground=True,
                                    zerolinecolor="black"),
                                zaxis = dict(
                                    backgroundcolor=RGB,
                                    gridcolor="white",
                                    showbackground=True,
                                    zerolinecolor="white",),),
                    )
    # fig.update_layout()
    # Save the plot as an HTML file
    fig.write_html("borehole_geometry.html")
    
    # Show the plot
    # fig.show()