#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import numpy as np

import plotly.graph_objects as go
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

    # elif clg_configuration == 2: #U-loop geometry 

    #     Tw_final_injector = TwMatrix[-1, 0:(interconnections[0] - 1)]  # Final fluid temperature profile in injection well [°C]
    #     Tw_final_producer = TwMatrix[-1, interconnections[0]:interconnections[1] - 1]  # Final fluid temperature profile in production well [°C]
    #     print(interconnections) # when more than one lateral [ 35  71  91 111]
    #     # one lateral: [35 71]
    #     if numberoflaterals == 1:
    #         Tw_final_lateral=np.empty((numberoflaterals,TwMatrix[-1, interconnections[1] - 1 :].shape[0]))
    #     else:
    #         Tw_final_lateral=np.empty((numberoflaterals,TwMatrix[-1, interconnections[1] - 1 : interconnections[2] - 2].shape[0]))
    #     plt.figure()
    #     plt.plot(range(1, interconnections[0]), Tw_final_injector, 'b-', linewidth=2)
    #     plt.grid(True)
    #     plt.plot([-2, -1], [-2, -1], 'k-', linewidth=2)  # Dummy plot for legend
    #     plt.plot([-2, -1], [-2, -1], 'r-', linewidth=2)  # Dummy plot for legend
        
    #     for kk in range(numberoflaterals):
    #         if kk < numberoflaterals-1:
    #             Tw_final_lateral[kk,:] = TwMatrix[-1, interconnections[kk+1] - kk - 1 : interconnections[kk+2] - kk - 2]
                
    #         else:
    #             Tw_final_lateral[kk,:] = TwMatrix[-1, interconnections[kk+1] - numberoflaterals :]
    #         if(kk==0):
    #             plt.plot(np.arange(len(xinj)-2, len(xinj) - 2 + len(xlat)), np.append(np.array([Tw_final_injector[-1]]), Tw_final_lateral[kk,:].reshape(1,len(Tw_final_lateral[kk,:]))), 'k-', linewidth=2)
    #         if (kk==1):
    #             plt.plot(np.arange(len(xinj)-2, len(xinj) - 2 + len(xlat)), np.append(np.array([Tw_final_injector[-1]]), Tw_final_lateral[kk,:].reshape(1,len(Tw_final_lateral[kk,:]))), 'm-', linewidth=2)
    #         if (kk==2):
    #             plt.plot(np.arange(len(xinj)-2, len(xinj) - 2 + len(xlat)), np.append(np.array([Tw_final_injector[-1]]), Tw_final_lateral[kk,:].reshape(1,len(Tw_final_lateral[kk,:]))), 'c-', linewidth=2)
        
    #     plt.plot((np.arange(interconnections[0] - 1, interconnections[1] - 1) + len(xlat) - 1), 
    #             np.append(np.sum(lateralflowallocation * Tw_final_lateral[:,-1]) , np.array([Tw_final_producer])), 'r-', linewidth=2)
        
    #     plt.ylabel('Fluid Temperature [°C]', fontsize=12)
    #     plt.xlabel('Position along flow path [-]', fontsize=12)
    #     plt.xticks(fontsize=12)
    #     plt.yticks(fontsize=12)
    #     plt.title('Final Fluid Temperature')
    #     plt.legend(['Injection Well', 'Lateral(s)', 'Production Well'], loc='upper left')
    #     plt.axis([0, len(xinj) + len(xprod) + len(xlat) - 2, min(TwMatrix[-1, :]) - 1, max(TwMatrix[-1, :]) + 1])
    #     plt.show()  

    elif clg_configuration == 2: #U-loop geometry 
        if sbt_version == 1:
            Tw_final_injector = TwMatrix[-1, 0:(interconnections[0] - 1)]  # Final fluid temperature profile in injection well [°C]
            Tw_final_producer = TwMatrix[-1, interconnections[0]:interconnections[1] - 1]  # Final fluid temperature profile in production well [°C]
            if numberoflaterals == 1:
                Tw_final_lateral=np.empty((numberoflaterals,TwMatrix[-1, interconnections[1] - 1 :].shape[0]))
            else:
                Tw_final_lateral=np.empty((numberoflaterals,TwMatrix[-1, interconnections[1] - 1 : interconnections[2] - 2].shape[0]))
            # Tw_final_lateral=np.empty((numberoflaterals,TwMatrix[-1, interconnections[1] - 1 : interconnections[2] - 2].shape[0]))
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
        elif sbt_version == 2:
            Tw_final_injector = Tfluidnodes[:interconnections[0]+1]  # Final fluid temperature profile in injection well [°C]
            Tw_final_producer = Tfluidnodes[interconnections[0]+1:interconnections[1]+1]  # Final fluid temperature profile in production well [°C]
            plt.figure()
            plt.plot(range(1, interconnections[0] + 2), Tw_final_injector, 'b-', linewidth=2, label='Injection Well')
            plt.grid(True)
            #Dummy plot for legend purposes
            plt.plot([-2, -1], [-2, -1], 'k-', linewidth=2, label='Lateral(s)')
            plt.plot([-2, -1], [-2, -1], 'r-', linewidth=2, label='Production Well')
            for dd in range(numberoflaterals):
                Tw_final_lateral = np.concatenate([
                    [Tfluidnodes[interconnections[0]]],
                    Tfluidnodes[lateralnodalstartpoints[dd]:lateralnodalendpoints[dd] + 1],
                    [Tfluidlateralexitstore[dd, -1]]
                ])  # Final fluid temperature profile in lateral dd [°C]
            
                lateral_x = range(len(xinj), len(xinj) + len(xlat))
                plt.plot(lateral_x, Tw_final_lateral, 'k-', linewidth=2)
            producer_x = range(interconnections[0] + len(xlat), interconnections[1] + len(xlat))
            plt.plot(producer_x, Tw_final_producer, 'r-', linewidth=2)
            plt.ylabel('Fluid Temperature [°C]', fontsize=12)
            plt.xlabel('Position along flow path [-]', fontsize=12)
            plt.title('Final Fluid Temperature', fontsize=14)
            plt.gca().tick_params(labelsize=12)
            plt.legend(loc='upper left')
            x_axis_limit = len(xinj) + len(xprod) + len(xlat) - 2
            y_min = min(Tfluidnodes) - 1
            y_max = max(Tfluidnodes) + 1
            plt.axis([0, x_axis_limit, y_min, y_max])

    #Plot final fluid pressure profile (SBT v2 only)
    if sbt_version == 2:
        if clg_configuration == 1: #co-axial geometry 
            plt.figure()
            plt.plot(Pfluiddownnodes/1e5,-np.cumsum(np.concatenate(([0], Deltaz))),'b-')
            plt.plot(Pfluidupnodes/1e5,-np.cumsum(np.concatenate(([0], Deltaz))),'r-')        
            plt.grid(True)
            plt.xlabel('Fluid Pressure [bar]', fontsize=12)
            plt.ylabel('Measured Depth [m]', fontsize=12)
            plt.xticks(fontsize=12)
            plt.yticks(fontsize=12)
            plt.title('Final Fluid Pressure')
            if coaxialflowtype == 1:
                plt.legend(['Annulus (Injection)', 'Center Pipe (Production)'], loc='best')
            elif coaxialflowtype == 2:
                plt.legend(['Center Pipe (Injection)', 'Annulus (Production)'], loc='best')
            plt.show()    
            
        elif clg_configuration == 2: #U-loop geometry 
            Pw_final_injector = np.array(Pfluidnodes[:interconnections[0]+1]) / 1e5  # Final fluid pressure profile in injection well [bar]
            Pw_final_producer = np.array(Pfluidnodes[interconnections[0]+1:interconnections[1]+1]) / 1e5  # Final fluid pressure profile in production well [bar]
            plt.figure()
            plt.plot(range(1, interconnections[0] + 2), Pw_final_injector, 'b-', linewidth=2, label='Injector')
            plt.grid(True)
            plt.plot([-2, -1], [-2, -1], 'k-', linewidth=2, label='Lateral')  # Dummy plot for legend
            plt.plot([-2, -1], [-2, -1], 'r-', linewidth=2, label='Producer')  # Dummy plot for legend
            for dd in range(numberoflaterals):
                Pw_final_lateral = np.concatenate([
                    [Pfluidnodes[interconnections[0]]],
                    Pfluidnodes[lateralnodalstartpoints[dd]:lateralnodalendpoints[dd] + 1],
                    [Pfluidlateralexit[dd]]
                ]) / 1e5  # Final fluid pressure profile in lateral dd [bar]
            
                lateral_x = range(len(xinj), len(xinj) + len(xlat))
                plt.plot(lateral_x, Pw_final_lateral, 'k-', linewidth=2)
            
            producer_x = range(interconnections[0] + len(xlat), interconnections[1] + len(xlat))
            plt.plot(producer_x, Pw_final_producer, 'r-', linewidth=2, label='Producer')
            plt.ylabel('Fluid Pressure [bar]', fontsize=12)
            plt.xlabel('Position along flow path [-]', fontsize=12)
            plt.title('Final Fluid Pressure', fontsize=14)
            plt.gca().tick_params(labelsize=12)
            plt.legend(['Injection Well', 'Lateral(s)', 'Production Well'], loc='best')
            x_axis_limit = len(xinj) + len(xprod) + len(xlat) - 2
            y_min = min(Pfluidnodes) / 1e5 - 10
            y_max = max(Pfluidnodes) / 1e5 + 10
            plt.axis([0, x_axis_limit, y_min, y_max])
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

def plot_heat_production(clg_configuration, HeatProduction, times):

    # Plot heat production
    # plt.figure()
    # plt.plot(times[1:]/3600/24/365, HeatProduction[1:], linewidth=2, color='green')
    # plt.axis([0, times[-1]/3600/24/365, 0, max(HeatProduction)])
    # plt.xlabel('Time [years]', fontsize=12)
    # plt.ylabel('Heat Production [MWt]', fontsize=12)
    # plt.gca().set_facecolor((1, 1, 1))
    # plt.grid(True)
    # plt.xticks(fontsize=12)
    # plt.yticks(fontsize=12)
    # plt.show()

    plt.figure()
    plt.plot(times[1:]/3600/24/365, HeatProduction[1:], linewidth=2, color='green')
    if clg_configuration == 1: #co-axial geometry
        plt.axis([0, times[-1]/3600/24/365, 0, 5*AverageHeatProduction])
    elif clg_configuration == 2: #U-loop geometry
        plt.axis([0, times[-1]/3600/24/365, 0, max(HeatProduction)])
    plt.xlabel('Time [years]', fontsize=12)
    plt.ylabel('Heat Production [MWt]', fontsize=12)
    plt.gca().set_facecolor((1, 1, 1))
    plt.grid(True)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.show()

def plot_production_temperature(sbt_version, Poutput, Pin, times):

    # Plot production pressure
    if sbt_version == 2:
        time_years = np.array(times[1:]) / (365 * 24 * 3600)
        plt.figure()
        plt.plot(time_years, Poutput[1:], linewidth=2, color='red', label='Production Pressure')
        plt.axis([0, times[-1] / (365 * 24 * 3600), Pin - 10, max(Poutput) + 10])
        plt.plot([0, times[-1] / (365 * 24 * 3600)], [Pin, Pin], linewidth=2, color='blue', label='Injection Pressure')
        plt.xlabel('Time [years]', fontsize=12)
        plt.ylabel('Pressure [bar]', fontsize=12)
        plt.grid(True)
        plt.gca().tick_params(labelsize=12)
        plt.legend(fontsize=12)
        plt.gcf().set_facecolor('white')
        plt.show()

def plot_production_temperature_linear(clg_configuration, Toutput, Tinstore, times):

    # Plot production temperature with linear time scale
    # plt.figure()
    # plt.plot(times[1:]/365/24/3600, Toutput[1:], linewidth=2, color='blue', label='Production Temperature')
    # plt.axis([0, times[-1]/3600/24/365, min(Tinstore)-10, max(Toutput)])
    # plt.plot(times[1:]/365/24/3600, Tinstore[1:], linewidth=2, color='black', label='Injection Temperature')
    # plt.xlabel('Time [years]', fontsize=12)
    # plt.ylabel('Temperature [°C]', fontsize=12)
    # plt.gca().set_facecolor((1, 1, 1))
    # plt.grid(True)
    # plt.legend()
    # plt.xticks(fontsize=12)
    # plt.yticks(fontsize=12)
    # plt.show()

    plt.figure()
    plt.plot(times[1:]/365/24/3600, Toutput[1:], linewidth=2, color='red', label='Production Temperature')
    if clg_configuration == 1: #co-axial geometry
        plt.axis([0, times[-1]/3600/24/365, min(Tinstore)-10, AverageProductionTemperature+5*(AverageProductionTemperature-min(Tinstore))])
    elif clg_configuration == 2: #U-loop geometry
        plt.axis([0, times[-1]/3600/24/365, min(Tinstore)-10, max(Toutput)])
    plt.plot(times[1:]/365/24/3600, Tinstore[1:], linewidth=2, color='blue', label='Injection Temperature')
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