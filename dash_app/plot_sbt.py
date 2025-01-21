#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import numpy as np

def plot_borehole_geometry(clg_configuration, numberoflaterals, x, y, z, xinj, yinj, zinj, xprod, yprod, zprod, xlat, ylat, zlat):

    # Make 3D figure of borehole geometry to make sure it looks correct
    plt.close('all') #close all curent figures
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    if clg_configuration == 1: #co-axial geometry
        ax.plot(x, y, z, 'k-o', linewidth=2)
        ax.set_xlim([np.min(x) - 200, np.max(x) + 200])
        ax.set_ylim([np.min(y) - 200, np.max(y) + 200])
        ax.set_zlim([np.min(z) - 500, 0])
        ax.set_zlabel('Depth (m)')
        ax.set_xlabel('x (m)')
        ax.set_ylabel('y (m)')
    
    elif clg_configuration == 2: #U-loop geometry
        ax.plot(xinj, yinj, zinj, 'b-o', linewidth=2)
        #ax.axis('equal') # Uncomment this next line to set the plotted geometry to correct scale with equal axis unit spacing
        ax.plot(xprod, yprod, zprod, 'r-o', linewidth=2)
        for i in range(numberoflaterals):
            ax.plot(xlat[:, i], ylat[:, i], zlat[:, i], 'k-o', linewidth=2)
        ax.set_xlim([np.min(x) - 200, np.max(x) + 200])
        ax.set_ylim([np.min(y) - 200, np.max(y) + 200])
        ax.set_zlim([np.min(z) - 500, 0])
        ax.set_zlabel('Depth (m)')
        ax.set_xlabel('x (m)')
        ax.set_ylabel('y (m)')
        ax.legend(['Injection Well', 'Production Well', 'Lateral(s)'])

    plt.show()



def plot_final_fluid_temp_profile_v1(sbt_version, clg_configuration, 
                                    Tw_up_previous, Tw_down_previous, Tfluiddownnodes, 
                                    Deltaz, TwMatrix, 
                                    numberoflaterals, coaxialflowtype,
                                    interconnections,
                                    lateralflowallocation,
                                    xinj, xlat, xprod
                                    ):

    # Plot final fluid temperature profile
    
    if clg_configuration == 1: #co-axial geometry 
        
        plt.figure()
        if sbt_version == 1:
            plt.plot(Tw_down_previous,-np.cumsum(Deltaz))
            plt.plot(Tw_up_previous,-np.cumsum(Deltaz))
        elif sbt_version == 2:
            plt.plot(Tfluiddownnodes,-np.cumsum(np.concatenate(([0], Deltaz))))
            plt.plot(Tfluidupnodes,-np.cumsum(np.concatenate(([0], Deltaz))))        
        plt.grid(True)
        plt.xlabel('Fluid Temperature [°C]', fontsize=12)
        plt.ylabel('Measured Depth [m]', fontsize=12)
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        plt.title('Final Fluid Temperature')
        if coaxialflowtype == 1:
            plt.legend(['Annulus (Injection)', 'Center Pipe (Production)'], loc='upper left')
        elif coaxialflowtype == 2:
            plt.legend(['Center Pipe (Injection)', 'Annulus (Production)'], loc='upper left')
        plt.show()

    elif clg_configuration == 2: #U-loop geometry 

        Tw_final_injector = TwMatrix[-1, 0:(interconnections[0] - 1)]  # Final fluid temperature profile in injection well [°C]
        Tw_final_producer = TwMatrix[-1, interconnections[0]:interconnections[1] - 1]  # Final fluid temperature profile in production well [°C]
        Tw_final_lateral=np.empty((numberoflaterals,TwMatrix[-1, interconnections[1] - 1 : interconnections[2] - 2].shape[0]))
        plt.figure()
        plt.plot(range(1, interconnections[0]), Tw_final_injector, 'b-', linewidth=2)
        plt.grid(True)
        plt.plot([-2, -1], [-2, -1], 'k-', linewidth=2)  # Dummy plot for legend
        plt.plot([-2, -1], [-2, -1], 'r-', linewidth=2)  # Dummy plot for legend
        
        for kk in range(numberoflaterals):
            if kk < numberoflaterals-1:
                Tw_final_lateral[kk,:] = TwMatrix[-1, interconnections[kk+1] - kk - 1 : interconnections[kk+2] - kk - 2]
                
            else:
                Tw_final_lateral[kk,:] = TwMatrix[-1, interconnections[kk+1] - numberoflaterals :]
            if(kk==0):
                plt.plot(np.arange(len(xinj)-2, len(xinj) - 2 + len(xlat)), np.append(np.array([Tw_final_injector[-1]]), Tw_final_lateral[kk,:].reshape(1,len(Tw_final_lateral[kk,:]))), 'k-', linewidth=2)
            if (kk==1):
                plt.plot(np.arange(len(xinj)-2, len(xinj) - 2 + len(xlat)), np.append(np.array([Tw_final_injector[-1]]), Tw_final_lateral[kk,:].reshape(1,len(Tw_final_lateral[kk,:]))), 'm-', linewidth=2)
            if (kk==2):
                plt.plot(np.arange(len(xinj)-2, len(xinj) - 2 + len(xlat)), np.append(np.array([Tw_final_injector[-1]]), Tw_final_lateral[kk,:].reshape(1,len(Tw_final_lateral[kk,:]))), 'c-', linewidth=2)
        
        plt.plot((np.arange(interconnections[0] - 1, interconnections[1] - 1) + len(xlat) - 1), 
                np.append(np.sum(lateralflowallocation * Tw_final_lateral[:,-1]) , np.array([Tw_final_producer])), 'r-', linewidth=2)
        
        plt.ylabel('Fluid Temperature [°C]', fontsize=12)
        plt.xlabel('Position along flow path [-]', fontsize=12)
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        plt.title('Final Fluid Temperature')
        plt.legend(['Injection Well', 'Lateral(s)', 'Production Well'], loc='upper left')
        plt.axis([0, len(xinj) + len(xprod) + len(xlat) - 2, min(TwMatrix[-1, :]) - 1, max(TwMatrix[-1, :]) + 1])
        plt.show()  

def plot_final_fluid_temp_profile_v2(sbt_version, coaxialflowtype, Pfluiddownnodes, Pfluidupnodes, 
                                            Deltaz
                                            ):

    #Plot final fluid temperature profile (SBT v2 only)\
    if sbt_version == 2:
        plt.figure()
        plt.plot(Pfluiddownnodes,-np.cumsum(np.concatenate(([0], Deltaz))))
        plt.plot(Pfluidupnodes,-np.cumsum(np.concatenate(([0], Deltaz))))        
        plt.grid(True)
        plt.xlabel('Fluid Pressure [bar]', fontsize=12)
        plt.ylabel('Measured Depth [m]', fontsize=12)
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        plt.title('Final Fluid Pressure')
        if coaxialflowtype == 1:
            plt.legend(['Annulus (Injection)', 'Center Pipe (Production)'], loc='upper left')
        elif coaxialflowtype == 2:
            plt.legend(['Center Pipe (Injection)', 'Annulus (Production)'], loc='upper left')
        plt.show()    

def plot_heat_production(HeatProduction, times):

    # Plot heat production
    plt.figure()
    plt.plot(times[1:]/3600/24/365, HeatProduction[1:], linewidth=2, color='green')
    plt.axis([0, times[-1]/3600/24/365, 0, max(HeatProduction)])
    plt.xlabel('Time [years]', fontsize=12)
    plt.ylabel('Heat Production [MWt]', fontsize=12)
    plt.gca().set_facecolor((1, 1, 1))
    plt.grid(True)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.show()

def plot_production_temperature_linear(Toutput, Tinstore, times):

    # Plot production temperature with linear time scale
    plt.figure()
    plt.plot(times[1:]/365/24/3600, Toutput[1:], linewidth=2, color='blue', label='Production Temperature')
    plt.axis([0, times[-1]/3600/24/365, min(Tinstore)-10, max(Toutput)])
    plt.plot(times[1:]/365/24/3600, Tinstore[1:], linewidth=2, color='black', label='Injection Temperature')
    plt.xlabel('Time [years]', fontsize=12)
    plt.ylabel('Temperature [°C]', fontsize=12)
    plt.gca().set_facecolor((1, 1, 1))
    plt.grid(True)
    plt.legend()
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.show()

def plot_production_tempterature_log(Toutput, Tinstore, times):

    # Plot production temperature with logarithmic time scale
    plt.figure()
    plt.semilogx(times[1:], Toutput[1:], linewidth=2, color='blue')
    plt.axis([10**2, times[-1], 0, max(Toutput)+10])
    plt.xlabel('Time [s]', fontsize=12)
    plt.ylabel('Production Temperature [°C]', fontsize=12)
    plt.gca().set_facecolor((1, 1, 1))
    plt.grid(True)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.show()