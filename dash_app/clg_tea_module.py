#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# data manipulation libraries
import numpy as np
import scipy.io
from scipy.interpolate import interpn, interp1d
# import pdb

# sourced scripts
# import clgs as clgs_v2



def get_error_message(error_code: int) -> str:
    error_code_lookup = {
        1000: "production temperature drops below injection temperature",
        2000: "H2O Injection temperature too low. H20 Injection Temperature must be above 50C to fall within ORC correlations",
        2001: "H2O production temperature must be within 100C to 200C too fall within ORC correlations",
        3000: "Too low CO2 reinjection temperature. Injection temperature must be greater or equal to 32 degrees for it to be supercritical",
        4000: "Too low CO2 turbine outlet temperature. Turbine outlet temperature must be above 37 degrees",
    }

    error_message = error_code_lookup.get(error_code, None)
    if error_message:
        error_message = "• " + error_message
    return error_message

class TEA:
    def __init__(self, u_sCO2, u_H2O, c_sCO2, c_H2O,
                    Fluid, End_use, Configuration, Flow_user, Hor_length_user, Depth_user, Gradient_user, 
                    Diameter_user, Tin_user, krock_user, Drilling_cost_per_m, O_and_M_cost_plant, Discount_rate,
                    Pump_efficiency, Lifetime, Direct_use_heat_cost_per_kWth, Electricity_rate, 
                    Power_plant_cost_per_kWe, T0, P0, Turbine_isentropic_efficiency, Generator_efficiency,
                    Compressor_isentropic_efficiency, Pre_Cooling_Delta_T, Turbine_outlet_pressure,
                    properties_H2O_pathname, 
                    properties_CO2v2_pathname, 
                    additional_properties_CO2v2_pathname):
        
        self.Fluid = Fluid
        self.End_use = End_use
        self.Configuration = Configuration
        self.Flow_user = Flow_user
        self.Hor_length_user = Hor_length_user
        self.Depth_user = Depth_user
        self.Gradient_user = Gradient_user
        self.Diameter_user = Diameter_user
        self.Tin_user = Tin_user
        self.krock_user = krock_user
        self.Drilling_cost_per_m = Drilling_cost_per_m
        self.O_and_M_cost_plant = O_and_M_cost_plant
        self.Discount_rate = Discount_rate
        self.Pump_efficiency = Pump_efficiency
        self.Lifetime = Lifetime
        self.Direct_use_heat_cost_per_kWth = Direct_use_heat_cost_per_kWth
        self.Electricity_rate = Electricity_rate
        self.Power_plant_cost_per_kWe = Power_plant_cost_per_kWe
        self.T0 = T0
        self.P0 = P0
        self.Turbine_isentropic_efficiency = Turbine_isentropic_efficiency
        self.Generator_efficiency = Generator_efficiency
        self.Compressor_isentropic_efficiency = Compressor_isentropic_efficiency
        self.Pre_Cooling_Delta_T = Pre_Cooling_Delta_T
        self.Turbine_outlet_pressure = Turbine_outlet_pressure

        # self.filename = "../../geo-data/clgs_results_final.h5"               #Filename of h5 database with simulation results [-]
        self.Number_of_points_per_year = 4               #Number of time steps per year in database [-] (must be 4)
        
        self.point = (Flow_user, Hor_length_user, Depth_user, Gradient_user, Diameter_user, Tin_user, krock_user)
        
        self.P_in = 2e7         #Constant Injection pressure [Pa]
        self.T_in = Tin_user-273.15   #Injection temperature [deg.C]
        
    def verify(self): #Verify inputs are within allowable bounds
        self.error = 0 
        if self.Fluid != 1 and self.Fluid !=2:
            print("Error: Fluid must be 1 (H2O) or 2 (CO2). Simulation terminated.")
            self.error = 1
        if self.End_use != 1 and self.End_use !=2:
            print("Error: End_use must be 1 (Direct-Use) or 2 (Electricity). Simulation terminated.")
            self.error = 1
        if self.Flow_user < 5 or self.Flow_user > 100:
            print("Error: Flow rate must be between 5 and 100 kg/s. Simulation terminated.")
            self.error = 1
        if self.Hor_length_user < 1000 or self.Hor_length_user > 20000:
            print("Error: Horizontal length must be between 1,000 and 20,000 m. Simulation terminated.")
            self.error = 1
        if self.Depth_user < 1000 or self.Depth_user > 5000:
            print("Error: Vertical depth must be between 1,000 and 5,000 m. Simulation terminated.")
            self.error = 1
        if self.Gradient_user < 0.03 or self.Gradient_user > 0.07:
            print("Error: Geothermal gradient must be between 0.03 and 0.07 degrees C per m. Simulation terminated.")
            self.error = 1
        if self.Diameter_user < 0.2159 or self.Diameter_user > 0.4445:
            print("Error: Wellbore diameter must be between 0.2159 and 0.4445 m. Simulation terminated.")
            self.error = 1
        if self.Tin_user < 303.15 or self.Tin_user > 333.15:
            print("Error: Injection temperature must be between 303.15 and 333.15 K. Simulation terminated.")
            self.error = 1
        if self.krock_user < 1.5 or self.krock_user > 4.5:
            print("Error: Rock thermal conductivity must be between 1.5 and 4.5 W/m/K. Simulation terminated.")
            self.error = 1
        if self.Drilling_cost_per_m < 0 or self.Drilling_cost_per_m > 10000:
            print("Error: Drilling costs per m of measured depth must be between 0 and 10,000 $/m. Simulation terminated.")
            self.error = 1
        if self.O_and_M_cost_plant < 0 or self.O_and_M_cost_plant > 0.2:
            print("Error: Operation & maintance cost of surface plant (expressed as fraction of total surface plant capital cost) must be between 0 and 0.2. Simulation terminated.")
            self.error = 1
        if self.Discount_rate < 0 or self.Discount_rate > 0.2:
            print("Error: Discount rate must be between 0 and 0.2. Simulation terminated.")
            self.error = 1
        if self.Pump_efficiency < 0.5 or self.Pump_efficiency > 1:
            print("Error: Pump efficiency must be between 0.5 and 1. Simulation terminated.")
            self.error = 1
        if self.Lifetime < 5 or self.Lifetime > 40:
            print("Error: System lifetime must be between 5 and 40 years. Simulation terminated.")
            self.error = 1
        if isinstance(self.Lifetime, int) == False:
            print("Error: System lifetime must be integer. Simulation terminated.")
            self.error = 1
        if self.End_use == 1:    
            if self.Direct_use_heat_cost_per_kWth < 0 or self.Direct_use_heat_cost_per_kWth > 10000:
                print("Error: Capital cost for direct-use surface plant must be between 0 and 10,000 $/kWth. Simulation terminated.")
                self.error = 1
            if self.Electricity_rate < 0 or self.Electricity_rate > 0.5:
                print("Error: Electricity rate in direct-use for pumping power must be between 0 and 0.5 $/kWh. Simulation terminated.")
                self.error = 1
        if self.End_use == 2:    
            if self.Power_plant_cost_per_kWe < 0 or self.Power_plant_cost_per_kWe > 10000:
                print("Error: Power plant capital cost must be between 0 and 10,000 $/kWe. Simulation terminated.")
                self.error = 1
            if self.T0 < 278.15 or self.T0 > 303.15:
                print("Error: Dead-state temperature must be between 278.15 and 303.15 K. Simulation terminated.")
                self.error = 1
            if self.P0 < 0.8e5 or self.P0 > 1.1e5:
                print("Error: Dead state pressure must be between 0.8e5 and 1.1e5 Pa. Simulation terminated.")
                self.error = 1
        if self.Fluid == 2 and self.End_use == 2:
            if self.Turbine_isentropic_efficiency < 0.8 or self.Turbine_isentropic_efficiency > 1:
                print("Error: Turbine isentropic efficiency must be between 0.8 and 1. Simulation terminated.")
                self.error = 1
            if self.Generator_efficiency < 0.8 or self.Generator_efficiency > 1:
                print("Error: Generator efficiency must be between 0.8 and 1. Simulation terminated.")
                self.error = 1
            if self.Compressor_isentropic_efficiency < 0.8 or self.Compressor_isentropic_efficiency > 1:
                print("Error: Compressor isentropic efficiency must be between 0.8 and 1. Simulation terminated.")
                self.error = 1
            if self.Pre_Cooling_Delta_T < 0 or self.Pre_Cooling_Delta_T > 15:
                print("Error: CO2 temperature decline after turbine and before compressor must be between 0 and 15 degrees C. Simulation terminated.")
                self.error = 1
            if self.Turbine_outlet_pressure < 75 or self.Turbine_outlet_pressure > 200:
                print("Error: CO2 turbine outlet pressure must be between 75 and 200 bar. Simulation terminated.")
                self.error = 1
        return self.error
    
    def initialize(self, TandP_dict, u_sCO2, u_H2O, c_sCO2, c_H2O, properties_H2O_pathname, properties_CO2v2_pathname, additional_properties_CO2v2_pathname):
        
        if self.Fluid == 1:
            if self.Configuration == 1:
                self.u_H2O = u_H2O # clgs_v2.data(self.filename, "utube", "H2O")
            elif self.Configuration == 2:
                self.u_H2O = c_H2O # clgs_v2.data(self.filename, "coaxial", "H2O")
            # self.timearray = self.u_H2O.time
            self.timearray = self.u_H2O.time
            # self.timearray = np.array(TandP_dict["time"])
            self.FlowRateVector = self.u_H2O.mdot #length of 26
            self.HorizontalLengthVector = self.u_H2O.L2 #length of 20
            self.DepthVector = self.u_H2O.L1 #length of 9
            self.GradientVector = self.u_H2O.grad #length of 5
            self.DiameterVector = self.u_H2O.D #length of 3
            self.TinVector = self.u_H2O.Tinj #length of 3
            self.KrockVector = self.u_H2O.k #length of 3
            self.Fluid_name = 'Water'
        elif self.Fluid == 2:
            if self.Configuration == 1:
                self.u_sCO2 = u_sCO2 # clgs_v2.data(self.filename, "utube", "sCO2")
            elif self.Configuration == 2:
                self.u_sCO2 = c_sCO2 # clgs_v2.data(self.filename, "coaxial", "sCO2")
            self.timearray = self.u_sCO2.time
            # self.timearray = np.array(TandP_dict["time"])
            self.FlowRateVector = self.u_sCO2.mdot #length of 26
            self.HorizontalLengthVector = self.u_sCO2.L2 #length of 20
            self.DepthVector = self.u_sCO2.L1 #length of 9
            self.GradientVector = self.u_sCO2.grad #length of 5
            self.DiameterVector = self.u_sCO2.D #length of 3
            self.TinVector = self.u_sCO2.Tinj #length of 3
            self.KrockVector = self.u_sCO2.k #length of 3   
            self.Fluid_name = 'CarbonDioxide'     

   
            
        self.numberofcases = len(self.FlowRateVector)*len(self.HorizontalLengthVector)*len(self.DepthVector)*len(self.GradientVector)*len(self.DiameterVector)*len(self.TinVector)*len(self.KrockVector)
        
        
        self.Time_array = np.linspace(0,self.Lifetime*365*24*3600,1+self.Lifetime*self.Number_of_points_per_year) #[s]
        self.Linear_time_distribution = self.Time_array/365/24/3600
        self.TNOP = (self.Lifetime*self.Number_of_points_per_year+1)      #Total number of points for selected lifetime
        #Find closests lifetime
        closestlifetime = self.timearray.flat[np.abs(self.timearray - self.Lifetime).argmin()]    
        self.indexclosestlifetime = np.where(self.timearray == closestlifetime)[0][0]

        #load property data
        if self.Fluid == 1:
            mat = scipy.io.loadmat(properties_H2O_pathname) 
        else:
            mat = scipy.io.loadmat(properties_CO2v2_pathname) 
            additional_mat = scipy.io.loadmat(additional_properties_CO2v2_pathname)
        self.Pvector = mat['Pvector'][0]
        self.Tvector = mat['Tvector'][0]
        self.density = mat['density']
        self.enthalpy = mat['enthalpy']
        self.entropy = mat['entropy']
        if self.Fluid == 2:
            self.Pvector_ap = additional_mat['Pvector_ap'][0]
            self.hvector_ap = additional_mat['hvector_ap'][0]
            self.svector_ap = additional_mat['svector_ap'][0]
            self.TPh = additional_mat['TPh']
            self.hPs = additional_mat['hPs']
    
        #Define ORC power plant conversion efficiencies
        self.Utilization_efficiency_correlation_temperatures = np.array([100, 200]) #Linear correlation assumed here based on GEOPHIRES ORC correlation between 100 and 200 deg C [deg.C]
        self.Utilization_efficiency_correlation_conversion = np.array([0.2, 0.45])  #Efficiency of ORC conversion from production exergy to electricity based on GEOPHIRES correlation [-]
        self.Heat_to_power_efficiency_correlation_temperatures = np.array([100, 200]) #Linear correlation based on Chad Augustine's thesis [deg.C]
        self.Heat_to_power_efficiency_correlation_conversion = np.array([0.05, 0.14]) #Conversion from enthalpy to electricity [-]

        #Calculate dead-state enthalpy and entropy in case of electricity production
        if self.End_use == 2:   
            self.h_0 = interpn((self.Pvector,self.Tvector),self.enthalpy,np.array([self.P0,self.T0]))[0] #dead-state enthalpy [J/kg]
            self.s_0 = interpn((self.Pvector,self.Tvector),self.entropy,np.array([self.P0,self.T0]))[0] #dead-state entropy [J/kg/K]


        #Pre-populate specific heat capacity of air in case of electricity production
        if self.End_use == 2:
            self.Tair_for_cp_array = np.linspace(10,100,num=10)
            self.cp_air_array = np.array([1005.65818063, 1005.87727966, 1006.19281999, 1006.60616167, 1007.11890862, 1007.73265999, 1008.44882744, 1009.26850304, 1010.19236691, 1011.2206266])
              
        #Initialize heat/electricity arrays
        self.Instantaneous_production_enthalpy = np.zeros(len(self.Time_array))
        self.Instantaneous_temperature_after_isenthalpic_throttling = np.zeros(len(self.Time_array))
        self.Instantaneous_heat_production = np.zeros(len(self.Time_array))
        self.Annual_heat_production = np.zeros(self.Lifetime)
        self.Annual_pumping_power = np.zeros(self.Lifetime)
        self.Average_fluid_density = np.zeros(len(self.Time_array))
        if self.End_use == 2: #electricity generation
            self.Instantaneous_exergy_production = np.zeros(len(self.Time_array))  #Produced exergy only (independent from injection conditions)
            self.Instantaneous_exergy_extraction = np.zeros(len(self.Time_array))  #Difference between produced exergy and injected exergy
            self.Instantaneous_electricity_production_method_1 = np.zeros(len(self.Time_array)) #based on exergy produced (only for water)
            self.Instantaneous_electricity_production_method_2 = np.zeros(len(self.Time_array)) #based on exergy extracted
            self.Instantaneous_electricity_production_method_3 = np.zeros(len(self.Time_array)) #based on thermal efficiency
            self.Instantaneous_electricity_production_method_4 = np.zeros(len(self.Time_array)) #based on direct turbine expansion (for CO2)
            self.Instantaneous_utilization_efficiency_method_1 = np.zeros(len(self.Time_array)) #conversion from produced exergy to electricity
            self.Instantaneous_utilization_efficiency_method_2 = np.zeros(len(self.Time_array)) #conversion from extracted exergy to electricity
            self.Instantaneous_themal_efficiency = np.zeros(len(self.Time_array)) #conversion from enthalpy to electricity
            self.Annual_electricity_production = np.zeros(self.Lifetime)
        if self.Fluid == 2:
            self.Instantaneous_turbine_power = np.zeros(len(self.Time_array)) #Direct turbine expansion considered for systems using sCO2

        #Initialize error code
        self.error_codes = np.zeros(0)  #if error occurs, code will be assigned to this tag



    def getTandP(self, u_sCO2, u_H2O, c_sCO2, c_H2O, sbt_version, TandP_dict):
        
        hdf5_times = [0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.25, 
                        2.5, 2.75, 3, 3.25, 3.5, 3.75, 4, 4.25, 4.5, 4.75,
                        5, 5.25, 5.5, 5.75, 6, 6.25, 6.5, 6.75, 7, 7.25, 
                        7.5, 7.75, 8, 8.25, 8.5, 8.75, 9, 9.25, 9.5, 9.75, 
                        10, 10.25, 10.5, 10.75, 11, 11.25, 11.5, 11.75,
                        12, 12.25, 12.5, 12.75, 13, 13.25, 13.5, 13.75, 
                        14, 14.25, 14.5, 14.75, 15, 15.25, 15.5, 15.75, 
                        16, 16.25, 16.5, 16.75, 17, 17.25, 17.5, 17.75, 18, 
                        18.25, 18.5, 18.75, 19, 19.25, 19.5, 19.75, 20, 20.25, 
                        20.5, 20.75, 21, 21.25, 21.5, 21.75, 22, 22.25, 22.5, 
                        22.75, 23, 23.25, 23.5, 23.75, 24, 24.25, 24.5, 24.75, 
                        25, 25.25, 25.5, 25.75, 26, 26.25, 26.5, 26.75, 27, 
                        27.25, 27.5, 27.75, 28, 28.25, 28.5, 28.75, 29, 29.25, 29.5, 
                        29.75, 30, 30.25, 30.5, 30.75, 31, 31.25, 31.5, 31.75, 32, 32.25,
                        32.5, 32.75, 33, 33.25, 33.5, 33.75, 34, 34.25, 34.5, 34.75, 35,
                        35.25, 35.5, 35.75, 36, 36.25, 36.5, 36.75, 37, 37.25, 37.5,
                        37.75, 38, 38.25, 38.5, 38.75, 39, 39.25, 39.5, 39.75, 40]
  
        if self.Fluid == 1:
            # Check if H2O data exists and is not None
            if TandP_dict.get("H2O_Tout") is not None and TandP_dict.get("H2O_Pout") is not None and TandP_dict.get("time") is not None:
                time_data = np.array(TandP_dict["time"])
                temp_data = np.array(TandP_dict["H2O_Tout"])
                pressure_data = np.array(TandP_dict["H2O_Pout"])
                
                # Check array lengths and handle mismatches
                if len(time_data) != len(temp_data) or len(temp_data) == 0:
                    if len(temp_data) == 0:
                        # If temperature data is empty, use default values
                        self.Tout = np.ones(len(hdf5_times)) * (self.T_in + 273.15)
                        self.Pout = np.ones(len(hdf5_times)) * self.P_in
                        self.error_codes = np.append(self.error_codes, 8000)  # Missing thermal data
                        return
                    # Use the shorter length
                    min_len = min(len(time_data), len(temp_data))
                    time_data = time_data[:min_len]
                    temp_data = temp_data[:min_len]
                    pressure_data = pressure_data[:min_len] if len(pressure_data) >= min_len else pressure_data
                
                try:
                    f = interp1d(time_data, temp_data, fill_value="extrapolate")
                    self.Tout = f(np.array(hdf5_times))
                    self.Pout = pressure_data
                except Exception as e:
                    # Fallback to default values
                    self.Tout = np.ones(len(hdf5_times)) * (self.T_in + 273.15)
                    self.Pout = np.ones(len(hdf5_times)) * self.P_in
            else:
                # Use default values when data is None
                self.Tout = np.ones(len(hdf5_times)) * (self.T_in + 273.15)
                self.Pout = np.ones(len(hdf5_times)) * self.P_in
                # Mark this as invalid thermal data
                self.error_codes = np.append(self.error_codes, 8000)  # Missing thermal data

        elif self.Fluid == 2:
            # Check if sCO2 data exists and is not None
            if TandP_dict.get("sCO2_Tout") is not None and TandP_dict.get("sCO2_Pout") is not None and TandP_dict.get("time") is not None:
                time_data = np.array(TandP_dict["time"])
                temp_data = np.array(TandP_dict["sCO2_Tout"])
                pressure_data = np.array(TandP_dict["sCO2_Pout"])
                
                # Check array lengths and handle mismatches
                if len(time_data) != len(temp_data) or len(temp_data) == 0:
                    if len(temp_data) == 0:
                        # If temperature data is empty, use default values
                        self.Tout = np.ones(len(hdf5_times)) * (self.T_in + 273.15)
                        self.Pout = np.ones(len(hdf5_times)) * self.P_in
                        self.error_codes = np.append(self.error_codes, 8000)  # Missing thermal data
                        return
                    # Use the shorter length
                    min_len = min(len(time_data), len(temp_data))
                    time_data = time_data[:min_len]
                    temp_data = temp_data[:min_len]
                    pressure_data = pressure_data[:min_len] if len(pressure_data) >= min_len else pressure_data
                
                try:
                    f = interp1d(time_data, temp_data, fill_value="extrapolate")
                    self.Tout = f(np.array(hdf5_times))
                    self.Pout = pressure_data
                except Exception as e:
                    # Fallback to default values
                    self.Tout = np.ones(len(hdf5_times)) * (self.T_in + 273.15)
                    self.Pout = np.ones(len(hdf5_times)) * self.P_in
            else:
                # Use default values when data is None
                self.Tout = np.ones(len(hdf5_times)) * (self.T_in + 273.15)
                self.Pout = np.ones(len(hdf5_times)) * self.P_in
                # Mark this as invalid thermal data
                self.error_codes = np.append(self.error_codes, 8000)  # Missing thermal data

        #Initial time correction (Correct production temperature and pressure at time 0 (the value at time 0 [=initial condition] is not a good representation for the first few months)
        self.Tout[0] = self.Tout[1]
        self.Pout[0] = self.Pout[1]
        
        #Extract Tout and Pout over lifetime
        self.InterpolatedTemperatureArray = self.Tout[0:self.indexclosestlifetime+1]-273.15
        self.InterpolatedPressureArray = self.Pout[0:self.indexclosestlifetime+1]
        
    def calculateLC(self):
        self.Linear_production_temperature = self.InterpolatedTemperatureArray
        self.Linear_production_pressure = self.InterpolatedPressureArray
        self.AveProductionTemperature = np.average(self.Linear_production_temperature)

        # print(self.Linear_production_temperature[1:10]) # AB: way too low
        
        # if self.Fluid_name == "Water":
        #     # print(self.Linear_time_distribution) # ok
        #     # print(self.timearray) # ok
        #         # print(self.FlowRateVector) # same
        #         # print(self.HorizontalLengthVector, self.DepthVector, self.GradientVector, self.DiameterVector, self.TinVector, self.KrockVector) #same
        #     if self.End_use == 2:
        #         print("TEA START")
        #         print(self.timearray)
        #         # print(times)
        #         print("----------")
        #         print(self.Fluid_name, self.End_use)
        #         print("avg prod temp: ", self.AveProductionTemperature)
        #         print("Pout: ", self.Pout[1:10])
        #         print("Tout: ", self.Tout[1:10])
        #         print(self.Tout.shape)
        #         print(self.timearray.shape)
        #         print(self.indexclosestlifetime) # same
        #         print(self.timearray.flat[np.abs(self.timearray - self.Lifetime).argmin()]) # closest lifetime # same
        #         # print(self.Lifetime) # same 
        #         # self.Tout[0:self.indexclosestlifetime+1]-273.15 # same 
        #         print("Interp Temp Array: ", self.InterpolatedTemperatureArray[1:10])        

        self.AveProductionPressure = np.average(self.Linear_production_pressure)/1e5  #[bar]
        self.Flow_rate = self.Flow_user #Total flow rate [kg/s]
        self.calculatedrillinglength()
        if min(self.Linear_production_temperature) > self.T_in:
            self.calculateheatproduction()
            if self.End_use == 2:
                self.calculateelectricityproduction()
            self.calculatecapex()
            self.calculatopex()
            
            
            Discount_vector = 1./np.power(1+self.Discount_rate,np.linspace(0,self.Lifetime-1,self.Lifetime))
            if self.End_use == 1:   #direct-use heating
                self.LCOH = (self.TotalCAPEX + np.sum(self.OPEX_Plant*Discount_vector))*1e6/np.sum(self.Annual_heat_production/1e3*Discount_vector) #$/MWh
                # print(self.LCOH)
                # print("\n\n")

                if isinstance(self.LCOH, (int, float)) and self.LCOH<0:
                    self.LCOH = "Negative LCOH - System may be losing heat"
                    self.error_codes = np.append(self.error_codes,5000)
            elif self.End_use == 2: #electricity production
                if self.Average_electricity_production == 0:
                    self.LCOE = "No electricity production - Check system parameters"
                    self.error_codes = np.append(self.error_codes,6000)
                else:
                    self.LCOE = (self.TotalCAPEX + np.sum(self.OPEX_Plant*Discount_vector))*1e6/np.sum((self.Annual_electricity_production-self.Annual_pumping_power)/1e3*Discount_vector) #$/MWh
                if isinstance(self.LCOE, (int, float)) and self.LCOE<0:
                    self.LCOE = "Negative LCOE - System may be inefficient"
                    self.error_codes = np.append(self.error_codes,7000)
            
        else:  #Production temperature went below injection temperature # AB HERE *****
            self.error_codes = np.append(self.error_codes,1000)
            # Set meaningful error messages for when thermal data is invalid
            if self.End_use == 1:
                self.LCOH = "Invalid thermal data - Production temp below injection temp"
            elif self.End_use == 2:
                self.LCOE = "Invalid thermal data - Production temp below injection temp"
        
        # Check for missing thermal data (error code 8000)
        if 8000 in self.error_codes:
            if self.End_use == 1:
                self.LCOH = "Missing thermal data - Check system parameters"
            elif self.End_use == 2:
                self.LCOE = "Missing thermal data - Check system parameters"
    
 
    def calculatedrillinglength(self):
        if self.Configuration == 1:
            self.Drilling_length = self.Hor_length_user + 2*self.Depth_user  #Total drilling depth of both wells and lateral in U-loop [m]
        elif self.Configuration == 2:
            self.Drilling_length = self.Hor_length_user + self.Depth_user  #Total drilling depth of well and lateral in co-axial case [m]        
    
    def calculateheatproduction(self):
        #Calculate instantaneous heat production
        self.Average_fluid_density = interpn((self.Pvector,self.Tvector),self.density,np.dstack((0.5*self.P_in + 0.5*self.Linear_production_pressure,0.5*self.T_in + 0.5*self.Linear_production_temperature+273.15))[0])
        self.hprod = interpn((self.Pvector,self.Tvector),self.enthalpy,np.dstack((self.Linear_production_pressure,self.Linear_production_temperature+273.15))[0])
        self.hinj = interpn((self.Pvector,self.Tvector),self.enthalpy,np.array([self.P_in,self.T_in+273.15]))
        self.Instantaneous_heat_production = self.Flow_rate*(self.hprod - self.hinj)/1000 #Heat production based on produced minus injected enthalpy [kW]
        
        # print(self.Instantaneous_heat_production.shape)

        #Calculate annual heat production (kWh)
        self.Annual_heat_production = 8760/5*(self.Instantaneous_heat_production[0::4][0:-1]+self.Instantaneous_heat_production[1::4]+self.Instantaneous_heat_production[2::4]+self.Instantaneous_heat_production[3::4]+self.Instantaneous_heat_production[4::4])
      
        #Calculate average heat production
        self.AveAnnualHeatProduction = np.average(self.Annual_heat_production) #kWh
        self.AveInstHeatProduction = np.average(self.Instantaneous_heat_production) #kWth
        
        #Calculate average heat production and first year heat production
        self.Average_heat_production = np.average(self.Instantaneous_heat_production) #[kW]
        #Average_production_temperature = np.average(Linear_production_temperature) #[deg.C]
        self.FirstYearHeatProduction = self.Annual_heat_production[0] #kWh
        
        self.calculatepumpingpower()
        
        
    def calculateelectricityproduction(self):
        
        #Calculate instantaneous exergy production, exergy extraction, and electricity generation (MW) and annual electricity generation [kWh]
        self.h_prod = self.hprod #produced enthalpy [J/kg]
        self.h_inj = self.hinj #injected enthalpy [J/kg]
        self.s_prod = interpn((self.Pvector,self.Tvector),self.entropy,np.dstack((self.Linear_production_pressure,self.Linear_production_temperature+273.15))[0]) #produced entropy [J/kg/K]
        self.s_inj = interpn((self.Pvector,self.Tvector),self.entropy,np.array([self.P_in,self.T_in+273.15])) #injected entropy [J/kg/K]
            
        self.Instantaneous_exergy_production = (self.Flow_rate*(self.h_prod-self.h_0 - self.T0*(self.s_prod-self.s_0)))/1000 #[kW]
        self.Instantaneous_exergy_extraction = (self.Flow_rate*(self.h_prod-self.h_inj - self.T0*(self.s_prod-self.s_inj)))/1000 #[kW]     
            
        self.AverageInstNetExergyProduction = np.average(self.Instantaneous_exergy_production) #[kW]
        self.AverageInstNetExergyExtraction = np.average(self.Instantaneous_exergy_extraction) #[kW]
            
        if self.Fluid == 1:
            
            if self.T_in >= 50 and min(self.Linear_production_temperature) >= 100 and max(self.Linear_production_temperature) <= 385:
                self.Instantaneous_utilization_efficiency_method_1 = np.interp(self.Linear_production_temperature,self.Utilization_efficiency_correlation_temperatures,self.Utilization_efficiency_correlation_conversion,left = 0) #Utilization efficiency based on conversion of produced exergy to electricity
                self.Instantaneous_electricity_production_method_1 = self.Instantaneous_exergy_production*self.Instantaneous_utilization_efficiency_method_1 #[kW]
                self.Instantaneous_themal_efficiency = np.interp(self.Linear_production_temperature,self.Heat_to_power_efficiency_correlation_temperatures,self.Heat_to_power_efficiency_correlation_conversion,left = 0) #Utilization efficiency based on conversion of produced exergy to electricity
                self.Instantaneous_electricity_production_method_3 = self.Instantaneous_heat_production*self.Instantaneous_themal_efficiency #[kW]
                
            else: #Water injection temperature and/or production tempeature fall outside the range used in the correlations
                if self.T_in < 50: 
                    self.error_codes = np.append(self.error_codes,2000)
                if min(self.Linear_production_temperature) < 100 or max(self.Linear_production_temperature) > 385:
                    self.error_codes = np.append(self.error_codes,2001)
                self.Instantaneous_utilization_efficiency_method_1 = np.zeros(len(self.Time_array))
                self.Instantaneous_electricity_production_method_1 = np.zeros(len(self.Time_array))
                self.Instantaneous_themal_efficiency = np.zeros(len(self.Time_array))
                self.Instantaneous_electricity_production_method_3 = np.zeros(len(self.Time_array))
                    
            #based on method 1 for now (could be 50-50)
            self.Annual_electricity_production = 8760/5*(self.Instantaneous_electricity_production_method_1[0::4][0:-1]+self.Instantaneous_electricity_production_method_1[1::4]+self.Instantaneous_electricity_production_method_1[2::4]+self.Instantaneous_electricity_production_method_1[3::4]+self.Instantaneous_electricity_production_method_1[4::4])
            self.Inst_electricity_production = self.Instantaneous_electricity_production_method_1 #[kW]
            self.AveInstElectricityProduction = np.average(self.Instantaneous_electricity_production_method_1) #[kW]
            
    
        elif self.Fluid == 2:
            T_prod = self.Linear_production_temperature #Production temperature [deg.C]
            P_prod = self.Linear_production_pressure    #Production pressure [Pa]
    
            h_turbine_out_ideal = interpn((self.Pvector_ap,self.svector_ap),self.hPs,np.dstack((np.ones(self.TNOP)*self.Turbine_outlet_pressure*1e5,self.s_prod))[0]) 
            self.Instantaneous_turbine_power = self.Flow_rate*(self.h_prod-h_turbine_out_ideal)*self.Turbine_isentropic_efficiency/1000 #Turbine output [kW]
            h_turbine_out_actual = self.h_prod-self.Instantaneous_turbine_power/self.Flow_rate*1000 #Actual fluid enthalpy at turbine outlet [J/kg]
            self.T_turbine_out_actual = interpn((self.Pvector_ap,self.hvector_ap),self.TPh,np.dstack((np.ones(self.TNOP)*self.Turbine_outlet_pressure*1e5,h_turbine_out_actual))[0])-273.15 
                
            if min(self.T_turbine_out_actual) > 37 and self.T_in > 32:
                self.Pre_cooling_temperature = min(self.T_turbine_out_actual) - self.Pre_Cooling_Delta_T
                
                Pre_compressor_h = interpn((self.Pvector,self.Tvector),self.enthalpy,np.array([self.Turbine_outlet_pressure*1e5,self.Pre_cooling_temperature+273.15])) 
                    
                Pre_cooling = self.Flow_rate*(h_turbine_out_actual - Pre_compressor_h)/1e3 #Pre-compressor cooling [kWth]
                Pre_compressor_s = interpn((self.Pvector,self.Tvector),self.entropy,np.array([self.Turbine_outlet_pressure*1e5,self.Pre_cooling_temperature+273.15])) 
                    
                    
                Post_compressor_h_ideal = interpn((self.Pvector_ap,self.svector_ap),self.hPs,np.array([self.P_in,Pre_compressor_s[0]])) 
                Post_compressor_h_actual = Pre_compressor_h + (Post_compressor_h_ideal-Pre_compressor_h)/self.Compressor_isentropic_efficiency #Actual fluid enthalpy at compressor outlet [J/kg]
                self.Post_compressor_T_actual = interpn((self.Pvector_ap,self.hvector_ap),self.TPh,np.array([self.P_in,Post_compressor_h_actual[0]])) - 273.15
                Compressor_Work = self.Flow_rate*(Post_compressor_h_actual - Pre_compressor_h)/1e3 #[kWe]
                Post_cooling = self.Flow_rate*(Post_compressor_h_actual - self.h_inj)/1e3 #Fluid cooling after compression [kWth]
                    
                if Post_cooling<0:
                    ResistiveHeating = -Post_cooling
                    Post_cooling = 0
                else:
                    ResistiveHeating = 0
                   
                Total_cooling = Pre_cooling + Post_cooling #Total CO2 cooling requirements [kWth]
                
                T_air_in_pre_cooler = self.T0-273.15
                T_air_out_pre_cooler = (self.T_turbine_out_actual+self.Pre_cooling_temperature)/2 #Air outlet temperature in pre-cooler [deg.C]
                cp_air = np.interp(0.5*T_air_in_pre_cooler+0.5*T_air_out_pre_cooler,self.Tair_for_cp_array,self.cp_air_array) #Air specific heat capacity in pre-cooler [J/kg/K]
                m_air_pre_cooler = Pre_cooling*1000/(cp_air*(T_air_out_pre_cooler - T_air_in_pre_cooler)) #Air flow rate in pre-cooler [kg/s]
               
                T_air_in_post_cooler = self.T0-273.15
                T_air_out_post_cooler = (self.Post_compressor_T_actual+self.T_in)/2 #Air outlet temperature in post-cooler [deg.C]
                cp_air = np.interp(0.5*T_air_in_post_cooler+0.5*T_air_out_post_cooler,self.Tair_for_cp_array,self.cp_air_array) #Air specific heat capacity in post-cooler [J/kg/K]
                    
                m_air_post_cooler = Post_cooling*1000/(cp_air*(T_air_out_post_cooler - T_air_in_post_cooler)) #Air flow rate in post-cooler [kg/s]
                
                Air_cooling_power = (m_air_pre_cooler+m_air_post_cooler)*0.25  #Electricity for air-cooling, assuming 0.25 kWe per kg/s [kWe] 
                
                self.Instantaneous_electricity_production_method_4 = self.Instantaneous_turbine_power*self.Generator_efficiency-Compressor_Work-Air_cooling_power-ResistiveHeating # Net electricity using CO2 direct turbine expansion cycle [kWe]
                self.Inst_electricity_production = self.Instantaneous_electricity_production_method_4 #[kW]
                self.Annual_electricity_production = 8760/5*(self.Instantaneous_electricity_production_method_4[0::4][0:-1]+self.Instantaneous_electricity_production_method_4[1::4]+self.Instantaneous_electricity_production_method_4[2::4]+self.Instantaneous_electricity_production_method_4[3::4]+self.Instantaneous_electricity_production_method_4[4::4])
                self.AveInstElectricityProduction = np.average(self.Instantaneous_electricity_production_method_4) #[kW]
                #check if negative
                if min(self.Instantaneous_electricity_production_method_4)<0:
                    self.error_codes = np.append(self.error_codes,5500) #Calculated electricity generation is negative
                    
                    self.Annual_electricity_production = np.zeros(self.Lifetime)
                    self.Inst_electricity_production = np.zeros(self.TNOP)
                    self.AveInstElectricityProduction = 0
            else: #turbine outlet or reinjection temperature too low
                if (self.T_in <= 32):
                    self.error_codes = np.append(self.error_codes,3000)
                    
                if (min(self.T_turbine_out_actual)<=37):
                    self.error_codes = np.append(self.error_codes,4000)
                    
                self.Annual_electricity_production = np.zeros(self.Lifetime)
                self.Inst_electricity_production = np.zeros(self.TNOP)
                self.AveInstElectricityProduction = 0
        
        self.calculatepumpingpower()
        
        self.Average_electricity_production = np.average(self.Annual_electricity_production)/8760 #[kW]
        self.AveAnnualElectricityProduction = np.average(self.Annual_electricity_production) #[kWh]
        self.AveInstNetElectricityProduction = self.AveInstElectricityProduction - np.average(self.PumpingPower) #[kW]
        self.AveAnnualNetElectricityProduction = self.AveAnnualElectricityProduction - np.average(self.Annual_pumping_power) #kWh
        self.FirstYearElectricityProduction = self.Annual_electricity_production[0] #kWh
        self.Inst_Net_Electricity_production = self.Inst_electricity_production-self.PumpingPower #[kW]

    def calculatepumpingpower(self):
        #Calculate pumping power
        self.PumpingPower = (self.P_in-self.Linear_production_pressure)*self.Flow_rate/self.Average_fluid_density/self.Pump_efficiency/1e3 #Pumping power [kW]
        self.PumpingPower[self.PumpingPower<0] = 0 #Set negative values to zero (if the production pressure is above the injection pressure, we throttle the fluid)
        self.Annual_pumping_power = 8760/5*(self.PumpingPower[0::4][0:-1]+self.PumpingPower[1::4]+self.PumpingPower[2::4]+self.PumpingPower[3::4]+self.PumpingPower[4::4]) #kWh
            
    def calculatecapex(self):
        self.CAPEX_Drilling = self.Drilling_length*self.Drilling_cost_per_m/1e6 #Drilling capital cost [M$]
        if self.End_use == 1:   #direct-use heating
            self.CAPEX_Surface_Plant = np.max(self.Instantaneous_heat_production)*self.Direct_use_heat_cost_per_kWth/1e6 #[M$]
        elif self.End_use == 2: #electricity production
            if self.Fluid == 1:
                self.CAPEX_Surface_Plant = np.max(self.Instantaneous_electricity_production_method_1)*self.Power_plant_cost_per_kWe/1e6 #[M$]
            elif self.Fluid == 2:
                self.CAPEX_Surface_Plant = np.max(self.Instantaneous_electricity_production_method_4)*self.Power_plant_cost_per_kWe/1e6 #[M$]
        
        self.TotalCAPEX = self.CAPEX_Drilling + self.CAPEX_Surface_Plant           #Total system capital cost (only includes drilling and surface plant cost) [M$]
        
    def calculatopex(self):
        #Calculate OPEX
        if self.End_use == 1: #direct-use heating
            self.OPEX_Plant = self.O_and_M_cost_plant*self.CAPEX_Surface_Plant + self.Annual_pumping_power*self.Electricity_rate/1e6  #Annual plant O&M cost [M$/year]
        elif self.End_use == 2: #electricity production
            self.OPEX_Plant = self.O_and_M_cost_plant*self.CAPEX_Surface_Plant  #Annual plant O&M cost [M$/year]
        self.AverageOPEX_Plant = np.average(self.OPEX_Plant)
    
    def printresults(self):
        ### Print results to screen
        print('##################')
        print('Simulation Results')
        print('##################')
        print('Number of cases in database = ' + str(self.numberofcases))
        print(" ")
        print('### Configuration ###')
        if self.End_use == 1:
            print('End-Use = Direct-Use')
        elif self.End_use == 2:    
            print('End-Use = Electricity')
        if self.Fluid == 1:
            print('Fluid = water')
        elif self.Fluid == 2:
            print('Fluid = sCO2')     
        if self.Configuration == 1:
            print('Design = U-loop')
        elif self.Configuration == 2:
            print('Design = Co-axial')
            
        #Print conditions
        print("Flow rate = " + "{0:.1f}".format((self.Flow_user)) +" kg/s")
        print("Lateral Length = " + str(round(self.Hor_length_user)) +" m")
        print("Vertical Depth = " + str(round(self.Depth_user)) +" m")
        print("Geothermal Gradient = " + "{0:.1f}".format((self.Gradient_user*1000)) +" deg.C/km")
        print("Wellbore Diameter = " + "{0:.4f}".format((self.Diameter_user)) +" m")
        print("Injection Temperature = " + "{0:.1f}".format((self.Tin_user-273.15)) +" deg.C")
        print("Thermal Conductivity = " + "{0:.2f}".format((self.krock_user)) +" W/m/K")
        
        #check if error occured
        if len(self.error_codes)>0:
            print(" ")
            if np.in1d(1000,self.error_codes): #plot the temperature and pressure for these
                # TODO: ERROR HAPPENS HERE (AB)
                print("Error: production temperature drops below injection temperature. Simulation terminated.\n")
                
            if np.in1d(2000,self.error_codes): #plot the temperature and pressure for these
                print("Error: Water injection temperature and/or production tempeature fall outside the range of the ORC correlations. These correlations require injection temperature larger than 50 deg.C and production temperature in the range of 100 to 200 deg.C. Electricity production set to 0. \n")
                
            if np.in1d(3000,self.error_codes): #CO2 injection temperature cannot be below 32 degrees C (CO2 must be supercritical)
                print("Error: Too low CO2 reinjection temperature. Injection temperature must be greater or equal to 32 degrees for it to be supercritical.\n")

            if np.in1d(4000,self.error_codes): #Turbine outlet CO2 temperature dropped below 37 degrees C
                print("Error: Too low CO2 turbine outlet temperature. Turbine outlet temperature must be above 37 degrees C.\n")    
                
            if np.in1d(5000,self.error_codes): #Calculated LCOH was negative, set to $9999/MWh
                print("Error: Calculated LCOH was negative and has been reset to $9,999/MWh.\n")     
                
            if np.in1d(5500,self.error_codes): #Negative electricity calculated with CO2 cycle
                print("Error: Calculated net electricity generation is negative and reset to 0.\n")
                
            if np.in1d(6000,self.error_codes): #Zero electricity production. LCOE set to $9999/MWh
                print("Error: Calculated net electricity production was 0. LCOE reset to $9,999/MWh.\n")   
        
            if np.in1d(7000,self.error_codes): #Calculated LCOE was negative, set to $9999/MWh
                print("Error: Calculated LCOE was negative and has been reset to $9,999/MWh.\n")  
                
        
        #Print results for heating
        if self.End_use==1:
            print(" ")
            print('### Reservoir Simulation Results ###')
            print("Average Production Temperature = " + "{0:.1f}".format((self.AveProductionTemperature)) + " deg.C")
            print("Average Production Pressure = " + "{0:.1f}".format((self.AveProductionPressure)) + " bar")
            if np.in1d(1000,self.error_codes) == False:
                print("Average Heat Production = " + "{0:.1f}".format((self.AveInstHeatProduction)) +" kWth" )    
                print("First Year Heat Production = " + "{0:.1f}".format((self.FirstYearHeatProduction/1e3)) + " MWh")    
                print(" ")
                print('### Cost Results ###')
                print("Total CAPEX = " + "{0:.1f}".format((self.TotalCAPEX)) + " M$")
                print("Drilling Cost = " + "{0:.1f}".format((self.CAPEX_Drilling)) + " M$")
                print("Surface Plant Cost = " + "{0:.1f}".format((self.CAPEX_Surface_Plant)) + " M$")
                print("OPEX = " + "{0:.1f}".format((self.AverageOPEX_Plant*1000)) + " k$/year")
                print("LCOH = " + "{0:.1f}".format((self.LCOH)) +" $/MWh")
            
        if self.End_use == 2:
            print(" ")
            print('### Reservoir Simulation Results ###')
            print("Average Production Temperature = " + "{0:.1f}".format((self.AveProductionTemperature)) + " deg.C")
            print("Average Production Pressure = " + "{0:.1f}".format((self.AveProductionPressure)) + " bar")
            if np.in1d(1000,self.error_codes) == False:
                print("Average Heat Production = " + "{0:.1f}".format((self.AveInstHeatProduction)) +" kWth" )    
                print("Average Net Electricity Production = " + "{0:.1f}".format((self.AveInstNetElectricityProduction)) +" kWe" )    
                print("First Year Heat Production = " + "{0:.1f}".format((self.FirstYearHeatProduction/1e3)) + " MWh")    
                print("First Year Electricity Production = " + "{0:.1f}".format((self.FirstYearElectricityProduction/1e3)) + " MWh")
                print(" ")
                print('### Cost Results ###')
                print("Total CAPEX = " + "{0:.1f}".format((self.TotalCAPEX)) + " M$")
                print("Drilling Cost = " + "{0:.1f}".format((self.CAPEX_Drilling)) + " M$")
                print("Surface Plant Cost = " + "{0:.1f}".format((self.CAPEX_Surface_Plant)) + " M$")
                print("OPEX = " + "{0:.1f}".format((self.AverageOPEX_Plant*1000)) + " k$/year")
                print("LCOE = " + "{0:.1f}".format((self.LCOE)) +" $/MWh")
                