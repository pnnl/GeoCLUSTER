#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# data manipulation libraries
import numpy as np
import h5py
from scipy.interpolate import interpn, interp1d
import CoolProp.CoolProp as CP
import itertools as iter
from sbt_v27 import run_sbt as run_sbt_final
import time

class data:
    def __init__(self, fname, case, fluid):

        self.fluid = fluid
        self.case = case

        file = h5py.File(fname, "r") 
        fixed_loc = "/" + case + "/fixed_params/"
        input_loc = "/" + case + "/" + fluid + "/input/"
        output_loc = "/" + case + "/" + fluid + "/output/"

        # independent vars
        self.mdot = file[input_loc + "mdot"][:]  # i0
        self.L2 = file[input_loc + "L2"][:]  # i1
        self.L1 = file[input_loc + "L1"][:]  # i2
        self.grad = file[input_loc + "grad"][:]  # i3
        self.D = file[input_loc + "D"][:]  # i4
        self.Tinj = file[input_loc + "T_i"][:]  # i5
        self.k = file[input_loc + "k_rock"][:]  # i6
        self.time = file[input_loc + "time"][:]  # i7
        self.ivars = (
            self.mdot,
            self.L2,
            self.L1,
            self.grad,
            self.D,
            self.Tinj,
            self.k,
            self.time,
        )

        # fixed vars
        self.Pinj = file[fixed_loc + "Pinj"][()]
        self.Tamb = file[fixed_loc + "Tamb"][()]

        # dim = Mdot x L2 x L1 x grad x D x Tinj x k
        Wt = file[output_loc + "Wt"][:]  # int mdot * dh dt
        We = file[output_loc + "We"][:]  # int mdot * (dh - Too * ds) dt

        self.GWhr = 1e6 * 3600000.0

        self.kWe_avg = (
            We * self.GWhr / (1000.0 * self.time[-1] * 86400.0 * 365.0)
        )
        self.kWt_avg = (
            Wt * self.GWhr / (1000.0 * self.time[-1] * 86400.0 * 365.0)
        )

        # dim = Mdot x L2 x L1 x grad x D x Tinj x k x time
        self.shape = (
            len(self.mdot),
            len(self.L2),
            len(self.L1),
            len(self.grad),
            len(self.D),
            len(self.Tinj),
            len(self.k),
            len(self.time),
        )

        # if you get an error that these files don't exist, run the make_zarr.py file in the data directory to build these files!
        self.Tout = file[f"/{case}/{fluid}/output/Tout_chunked"]
        self.Pout = file[f"/{case}/{fluid}/output/Pout_chunked"]

        self.CP_fluid = "CO2"
        if fluid == "H2O":
            self.CP_fluid = "H2O"

    def get_parameter_indices(self, array, target):
        """
        returns the slice that would give you the points in the given array around target. Assumes that array is sorted.
        I.E if array is [1, 2, 3] and you give it 1.5, the points around 1.5 are [1, 2], so it will return the indices of [1,2], which are [0, 1].
        However, since we need to use this indices to slice into the 8D matrix, we return a slice object that will slice indices.

        if the array contains the target exactly, it'll just return a slice for getting specifically that point.

        If we pass in "all" as the target, it will return a slice that will slice the entirety of that array
        """
        
        # Handle empty array
        if len(array) == 0:
            print(f"Warning: Empty array passed to get_parameter_indices")
            return slice(0, 1)

        if target == "all":
            return slice(None)  # slice all of the points
        
        # Handle None values - return a slice for the middle of the array
        if target is None:
            print(f"Warning: None value passed to get_parameter_indices, using middle of array")
            middle_idx = len(array) // 2
            return slice(middle_idx, middle_idx + 1)
        
            
        # NOTE: PROBLEM CLAUSE FOR NOT ALLOWING GEOGRAD TO BE MORE THAN 0.7
        # Check if the target value is within the array bounds
        # For metric units, use strict validation
        # For imperial units, we need to check if this is a converted value
        if target < array[0] or target > array[-1]:
            # Check if this might be an imperial unit conversion
            # We can detect imperial by checking if the values are in typical imperial ranges
            is_likely_imperial = False
            
            # Check for common imperial patterns
            if target > 1000 and array[-1] < 100:  # Length: meters vs feet
                is_likely_imperial = True
            elif target < 1 and array[-1] > 10:  # Geothermal gradient: K/m vs F/ft
                is_likely_imperial = True
            elif target > 100 and array[-1] < 10:  # Thermal conductivity: W/m-K vs Btu/ft-h-F
                is_likely_imperial = True
            elif target < 1 and array[-1] > 100:  # Heat capacity: J/kg-K vs Btu/lb-F
                is_likely_imperial = True
            elif target > 50 and array[-1] < 10:  # Density: kg/m³ vs lb/ft³
                is_likely_imperial = True
            elif target < 1 and array[-1] > 100:  # Pressure: Pa vs psi
                is_likely_imperial = True
            
            if is_likely_imperial:
                # Don't print warnings for imperial units - this is expected
            else:
                # This is likely a real metric validation error
                lineprint = f"Warning: expected given value {target} to be between min and max of given array ({array[0], array[-1]})"
                print(lineprint)
                # Don't raise exception, but log the warning
            
            # For values outside bounds, clamp to the nearest valid range
            if target < array[0]:
                return slice(0, 1)  # Use first value
            else:  # target > array[-1]
                return slice(-1, None)  # Use last value
        
        # Normal case: find the appropriate slice for interpolation
        for i, value in enumerate(array):
            if value == target:
                return slice(i, i + 1)
            if value > target:
                return slice(i - 1, i + 1)
        
        # Fallback: if we get here, return the last slice
        return slice(-1, None)

    def read_values_around_point_for_interpolation(
        self, zarr_array, point, parameter_values
    ):
        """
        reads the values around the given point(s) for interpolation from the passed in zarr_array.
        If "all" is passed in as any of the coordinates of the point,
        then it'll slice across the entirety of that coordinate's dimension
        """
        if any(val is None for val in point):
        
        indices = [
            self.get_parameter_indices(params, value)
            for value, params in zip(point, parameter_values)
        ]
        
        if any(idx is None for idx in indices):
            # Replace None values with slice(None) as a fallback
            indices = [idx if idx is not None else slice(None) for idx in indices]
        
        data = zarr_array[
            tuple(indices)
        ]  # need to convert into into a tuple or zarr thinks we're doing fancy indexing
        return indices, data

    def interpolate_points(self, zarr_array, point_to_read_around, points):
        """
        Reads in only the data directly around points_to_read_around from zarr_array and interpolates the passed in points on that data.
        Essentially a wrapper around interpn where we only read the data we need.
        Alot of the data from the original zarr_array is not needed for interpolation since the default linear interpolation method of interpn only needs the corners of the hypercube around a given point to interpolate it's value
        Instead of the entire 8D parameter space.
        """
        indices, values_around_point = self.read_values_around_point_for_interpolation(
            zarr_array, point_to_read_around, self.ivars
        )
        grid = [
            params[these_indices] for these_indices, params in zip(indices, self.ivars)
        ]  # the grid is the values of the parameters at the points we're interpolating between
        interpolated_points = interpn(grid, values_around_point, points)
        return interpolated_points

    def reshape_output(self, tout):
        NEW_SIZE = 161
        original_indices = np.arange(tout.shape[0])
        new_indices = np.arange(NEW_SIZE + 1,)
        interpolator = interp1d(original_indices, tout, kind='quadratic')
        new_data = interpolator(new_indices)
        return new_data


    def interp_outlet_states(self, point, sbt_version, 
                    Tsurf, c_m, rho_m, 
                    # radius_vertical, radius_lateral, n_laterals, lateral_flow, lateral_multiplier,
                    Diameter1, Diameter2, PipeParam3, PipeParam4, PipeParam5,
                    mesh, accuracy, 
                    HyperParam1, HyperParam3, HyperParam5
                    # mass_mode, temp_mode
                    ): # needs to be a callback option
        """
        :param sbt_version: 0 if not using SBT, 1 if using SBT v1, 2 if using SBT v2 
        """
        if sbt_version == 0:
            points = list(
                iter.product(
                    (point[0],),
                    (point[1],),
                    (point[2],),
                    (point[3],),
                    (point[4],),
                    (point[5],),
                    (point[6],),
                    self.time,
                )
            )

            point_to_read_around = (
                *point,
                "all",
            )  # unpacking point and adding "all" to the end of it. This tells interpolate_points to read in all of the time dimension
            Tout = self.interpolate_points(self.Tout, point_to_read_around, points)
            Pout = self.interpolate_points(self.Pout, point_to_read_around, points)
            times = self.time
            # print("TEMP OUT: -----****")
            # print(Tout)

        else:
            mdot, L2, L1, grad, D , Tinj, k = point

            # AB UNIT CONVERSIONS AND RENAMING
            # DIAMETER NEEDS TO GO FROM M TO KM so divide by 1000
            L2_original = L2
            L1_original = L1
            Tinj_original = Tinj
            
            L2 = L2/1000
            L1 = L1/1000
            Tinj = Tinj-273.15

            if self.CP_fluid == "H2O":
                fluid = 1
            elif self.CP_fluid == "CO2":
                fluid = 2
            else:
                fluid = 1 # water
            
            if sbt_version == 1:

                mass_mode = HyperParam1
                temp_mode = HyperParam3

                # Convert string values to integers for SBT model
                if mass_mode == "Constant" or mass_mode == 0:
                    mass_mode_b = 0
                elif mass_mode == "Variable" or mass_mode == 1:
                    mass_mode_b = 1
                else:
                    mass_mode_b = 0  # Default to constant

                if temp_mode == "Constant" or temp_mode == 0:
                    temp_mode_b = 0
                elif temp_mode == "Variable" or temp_mode == 1:
                    temp_mode_b = 1
                else:
                    temp_mode_b = 0  # Default to constant
                
                hyperparam1 = mass_mode_b
                hyperparam2 = "MassFlowRate.xlsx"
                hyperparam3 = temp_mode_b
                hyperparam4 = "InjectionTemperatures.xlsx"
                hyperparam5 = None
            
            if sbt_version == 2:

                reltolerance =  1e-5 
                maxnumberofiterations =  15
                fluid_mode = HyperParam5
                

                if fluid_mode == "Constant" or fluid_mode == 0:
                    fluid_mode_b = 0
                elif fluid_mode == "Variable" or fluid_mode == 1:
                    fluid_mode_b = 1
                else:
                    fluid_mode_b = 1  # Default to variable
                
                hyperparam1 = HyperParam1*10 # Pin (convert MPa to bar)
                hyperparam2 = HyperParam3 # pipe roughness
                hyperparam3 = fluid_mode_b # fluid mode
                hyperparam4 = reltolerance
                hyperparam5 = maxnumberofiterations


            # TODO: !!! ***
            if self.case == "coaxial":
                case = 1
                if PipeParam5 == "Inject in Annulus":
                    PipeParam5 = 1
                if PipeParam5 == "Inject in Center Pipe":
                    PipeParam5 = 2
            #     Diameter1 = radius_vertical #radius # Diameter1/2
            #     Diameter2 = radius_lateral #radiuscenterpipe # Diameter2/2
            #     PipeParam3 = n_laterals #thicknesscenterpipe
            #     PipeParam4 = lateral_flow # k_center_pipe
            #     PipeParam5 = lateral_multiplier # coaxialflowtype
            if self.case == "utube":
                case = 2
                # Diameter1 = radius_vertical
                # Diameter2 = radius_lateral
                # PipeParam3 = n_laterals
                PipeParam4 = [PipeParam4]
                # PipeParam5 = lateral_multiplier

            start = time.time()
            
            # print(hyperparam1, hyperparam2, hyperparam3, hyperparam4, hyperparam5)

            try:
                times, Tout, Pout = run_sbt_final(
                        ## Model Specifications 
                        sbt_version=sbt_version, mesh_fineness=mesh, HYPERPARAM1=hyperparam1, HYPERPARAM2=hyperparam2, 
                        HYPERPARAM3=hyperparam3, HYPERPARAM4=hyperparam4, HYPERPARAM5=hyperparam5, 
                        accuracy=accuracy,

                        ## Operations
                        clg_configuration=case, mdot=mdot, Tinj=Tinj, fluid=fluid, ## Operations
                        DrillingDepth_L1=L1, HorizontalExtent_L2=L2, #BoreholeDiameter=D, ## Wellbore Geometry
                        Diameter1=Diameter1, Diameter2=Diameter2, 
                        PipeParam3=PipeParam3, PipeParam4=PipeParam4, 
                        # PipeParam3=3, PipeParam4=[1/3,1/3,1/3], 
                        PipeParam5=PipeParam5, ## Tube Geometry

                        ## Geologic Properties
                        Tsurf=Tsurf, GeoGradient=grad, k_m=k, c_m=c_m, rho_m=rho_m, 
                        # Tsurf=20, GeoGradient=grad, k_m=k, c_m=825, rho_m=2875, 
                )
                
            except Exception as e:
                times, Tout, Pout = None, None, None
                
            if Pout is None:
                constant_pressure = 2e7 # 200 Bar in pascal || 2.09e7 
                # constant_pressure = 22228604.37405011
                Pout = constant_pressure * np.ones_like(Tout)
                
            end = time.time()

            if times is not None and Tout is not None and Pout is not None:
                times = times[14:]
                Tout = Tout[14:]
                Pout = Pout[14:]
            else:
                # Handle case where outputs are None
                pass

        return Tout, Pout, times

    def interp_outlet_states_contour(self, param, point, time_index=None):
        """
        Build a 2D (varied_param × mdot) slice of Tout, Pout at a single time.
        point = (mdot, L2, L1, grad, D, Tinj, k) in SI units.
        time_index: optional index into self.time (default: last snapshot).
        """
        from utils_labels import canonicalize_param_label
        
        # Canonicalize the parameter label to handle both metric and imperial units
        param_key = canonicalize_param_label(param) if param else ""
        
        # Map canonical keys to axis indices in self.ivars (0..7). time is 7.
        axis_by_key = {
            "L2": 1,    # Horizontal Extent
            "L1": 2,    # Vertical Extent  
            "grad": 3,  # Geothermal Gradient
            "D": 4,     # Borehole Diameter
            "Tinj": 5,  # Injection Temperature
            "k": 6,     # Rock Thermal Conductivity
        }
        
        if param_key not in axis_by_key:
            raise ValueError(f"Unknown parameter label: {param} (canonicalized: {param_key})")

        var_axis = axis_by_key[param_key]
        # Choose a time index for the contours (use end-of-life by default)
        if time_index is None:
            time_index = len(self.time) - 1
        t_fixed = self.time[time_index]

        # 2) Build the interpolation grid "points": (mdot_grid, varied_grid, fixed scalars..., time_fixed)
        varied_grid = self.ivars[var_axis]
        mdot_grid = self.mdot

        # unpack point in the same order as self.ivars[:-1]
        mdot, L2, L1, grad, D, Tinj, k = point

        def fixed_val(axis):
            return {
                1: L2, 2: L1, 3: grad, 4: D, 5: Tinj, 6: k
            }[axis]

        points = list(iter.product(
            mdot_grid,                                 # axis 0 (mdot) -> vary (x)
            (varied_grid if var_axis == 1 else (L2,)), # axis 1 (L2)
            (varied_grid if var_axis == 2 else (L1,)), # axis 2 (L1)
            (varied_grid if var_axis == 3 else (grad,)),# axis 3 (grad)
            (varied_grid if var_axis == 4 else (D,)),  # axis 4 (D)
            (varied_grid if var_axis == 5 else (Tinj,)),# axis 5 (Tinj in K)
            (varied_grid if var_axis == 6 else (k,)),  # axis 6 (k)
            (t_fixed,)                                  # axis 7 (time) -> fixed
        ))

        # 3) Tell interpolate_points to read "all" for mdot and the varied axis, and fix the rest:
        indices = []
        for axis, arr in enumerate(self.ivars):
            if axis == 0 or axis == var_axis:
                indices.append(slice(None))            # read all mdot & varied
            elif axis == 7:
                # time fixed to nearest index
                # use a tiny helper from get_parameter_indices to clamp/locate time
                idx_slice = self.get_parameter_indices(self.time, t_fixed)
                indices.append(idx_slice)
            else:
                # fixed values from 'point'
                val = [mdot, L2, L1, grad, D, Tinj, k][axis]
                idx_slice = self.get_parameter_indices(self.ivars[axis], val)
                indices.append(idx_slice)

        # 4) Read the corner hypercube and interpn on that reduced grid
        def interpn_on(zarr):
            data = zarr[tuple(indices)]
            grid = [arr[idx] for idx, arr in zip(indices, self.ivars)]
            return interpn(grid, data, points)

        try:
            Tout = interpn_on(self.Tout)
            Pout = interpn_on(self.Pout)
            # reshape to (len(varied), len(mdot)) then transpose to (mdot, varied)
            Tout = np.reshape(Tout, (len(varied_grid), len(mdot_grid))).T
            Pout = np.reshape(Pout, (len(varied_grid), len(mdot_grid))).T
        except Exception as e:
            print(f"Error in interp_outlet_states_contour: {e}")
            # Fallback dummy data to avoid callback crash (and make the error visible)
            Tout = np.full((len(mdot_grid), len(varied_grid)), 365.0)
            Pout = np.full((len(mdot_grid), len(varied_grid)), 2.2228604e7)

        return Tout, Pout

    def interp_kW(
        self,
        point,
        Tout,
        Pout,
        Tamb=300.0, # AB: is this surface temperature? What is Tamb?
    ):

        mdot = point[0]
        Tinj = point[5]

        enthalpy_out = CP.PropsSI("H", "T", Tout, "P", Pout, self.CP_fluid)
        enthalpy_in = CP.PropsSI("H", "T", Tinj, "P", self.Pinj, self.CP_fluid)
        entropy_out = CP.PropsSI("S", "T", Tout, "P", Pout, self.CP_fluid)
        entropy_in = CP.PropsSI("S", "T", Tinj, "P", self.Pinj, self.CP_fluid)

        kWe = (
            mdot
            * (enthalpy_out - enthalpy_in - Tamb * (entropy_out - entropy_in))
            / 1000.0
        )
        kWt = mdot * (enthalpy_out - enthalpy_in) / 1000.0

        return kWe, kWt

    def interp_kW_contour(self, param, point, Tout, Pout, Tamb=300.0):

        mdot = self.mdot
        Tinj = point[4]

        accept_one_dim_only = False

        try:  # PropsSI accepts multidim for Python 3.10 only

            enthalpy_out = CP.PropsSI("H", "T", Tout, "P", Pout, self.CP_fluid)
            enthalpy_in = CP.PropsSI("H", "T", Tinj, "P", self.Pinj, self.CP_fluid)
            entropy_out = CP.PropsSI("S", "T", Tout, "P", Pout, self.CP_fluid)
            entropy_in = CP.PropsSI("S", "T", Tinj, "P", self.Pinj, self.CP_fluid)

        except Exception:  # PropsSI accepts only one dim for Python 3.8, 3.9. 3.11

            accept_one_dim_only = True
            dim_one = Tout.shape[0]
            dim_two = Tout.shape[1]

            Tout = Tout.flatten(order="C")
            Pout = Pout.flatten(order="C")
            enthalpy_out = CP.PropsSI("H", "T", Tout, "P", Pout, self.CP_fluid)
            enthalpy_in = CP.PropsSI("H", "T", Tinj, "P", self.Pinj, self.CP_fluid)
            entropy_out = CP.PropsSI("S", "T", Tout, "P", Pout, self.CP_fluid)
            entropy_in = CP.PropsSI("S", "T", Tinj, "P", self.Pinj, self.CP_fluid)

        if not accept_one_dim_only:

            kWe = (
                mdot
                * (enthalpy_out - enthalpy_in - Tamb * (entropy_out - entropy_in))
                / 1000.0
            )
            kWt = mdot * (enthalpy_out - enthalpy_in) / 1000.0

        else:

            calc = (
                enthalpy_out - enthalpy_in - Tamb * (entropy_out - entropy_in)
            ) / 1000.0
            calc = np.reshape(calc, (dim_one, dim_two), order="C")
            kWe = mdot * calc
            calc = (enthalpy_out - enthalpy_in) / 1000.0
            calc = np.reshape(calc, (dim_one, dim_two), order="C")
            kWt = mdot * calc

        return kWe, kWt

    def interp_kWe_avg(self, point):
        ivars = self.ivars[:-1]
        return (
            self.GWhr
            * interpn(ivars, self.We, point)
            / (1000.0 * self.time[-1] * 86400.0 * 365.0)
        )

    def interp_kWt_avg(self, point):
        ivars = self.ivars[:-1]
        return (
            self.GWhr
            * interpn(ivars, self.Wt, point)
            / (1000.0 * self.time[-1] * 86400.0 * 365.0)
        )

    def compute_kW(self, index, Tamb=300.0):
        mdot = self.mdot[index[0]]
        Tinj = self.Tinj[index[5]]
        point = (
            index[0],
            index[1],
            index[2],
            index[3],
            index[4],
            index[5],
            index[6],
            ...,
        )
        Tout = self.Tout[point]
        Pout = self.Pout[point]
        enthalpy_out = CP.PropsSI("H", "T", Tout, "P", Pout, self.CP_fluid)
        enthalpy_in = CP.PropsSI("H", "T", Tinj, "P", self.Pinj, self.CP_fluid)
        entropy_out = CP.PropsSI("S", "T", Tout, "P", Pout, self.CP_fluid)
        entropy_in = CP.PropsSI("S", "T", Tinj, "P", self.Pinj, self.CP_fluid)
        kWe = (
            mdot
            * (enthalpy_out - enthalpy_in - Tamb * (entropy_out - entropy_in))
            / 1000.0
        )
        kWt = mdot * (enthalpy_out - enthalpy_in) / 1000.0

        return kWe, kWt