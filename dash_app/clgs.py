#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# data manipulation libraries
import numpy as np
import h5py
from scipy.interpolate import interpn, interp1d
import CoolProp.CoolProp as CP
import itertools
from sbt_v27 import run_sbt as run_sbt_final
from sbt_coaxial_support import check_coaxial_diameters
import time
from functools import lru_cache
import hashlib
import pickle
# import uuid

# Global cache for SBT calculations
_sbt_cache = {}
_cache_max_size = 50  # Limit cache size to prevent memory issues

def _make_cache_key(*args, **kwargs):
    """Create a hash key from function arguments for caching"""
    # Convert all arguments to a tuple of hashable values
    key_data = []
    for arg in args:
        if isinstance(arg, (int, float, str, bool, type(None))):
            key_data.append(arg)
        elif isinstance(arg, np.ndarray):
            key_data.append(tuple(arg.flatten()[:10]))  # Use first 10 elements for arrays
        else:
            key_data.append(str(arg))
    for k, v in sorted(kwargs.items()):
        if isinstance(v, (int, float, str, bool, type(None))):
            key_data.append((k, v))
        elif isinstance(v, np.ndarray):
            key_data.append((k, tuple(v.flatten()[:10])))
        else:
            key_data.append((k, str(v)))
    return hashlib.md5(pickle.dumps(tuple(key_data))).hexdigest()

def _get_cached_result(cache_key):
    """Get cached result if available"""
    return _sbt_cache.get(cache_key)

def _set_cached_result(cache_key, result):
    """Store result in cache, with size limit"""
    # prif"[_set_cached_result] Called with cache_key: {cache_key}, type: {type(cache_key)}", flush=True)
    cache_size_before = len(_sbt_cache)
    cache_keys_before = list(_sbt_cache.keys())[:5]  # First 5 keys
    if len(_sbt_cache) >= _cache_max_size:
        # Remove oldest entry (simple FIFO)
        oldest_key = next(iter(_sbt_cache))
        # prif"[CACHE EVICT] Cache full ({cache_size_before} entries), removing oldest key: {oldest_key}", flush=True)
        del _sbt_cache[oldest_key]
    _sbt_cache[cache_key] = result
    cache_size_after = len(_sbt_cache)
    cache_keys_after = list(_sbt_cache.keys())[:5]  # First 5 keys
    # prif"[_set_cached_result] Stored cache key: {cache_key}, cache size: {cache_size_before} -> {cache_size_after}", flush=True)
    # prif"[_set_cached_result] Cache keys before: {cache_keys_before}, after: {cache_keys_after}", flush=True)
    # Verify the key is actually in the cache
    # if cache_key in _sbt_cache:
    #     print(f"[_set_cached_result] Verified: cache key {cache_key} is now in cache", flush=True)
    # else:
    #     print(f"[_set_cached_result] ERROR: cache key {cache_key} was NOT stored in cache!", flush=True)

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

        self.kWe_avg = (We * self.GWhr / (1000.0 * self.time[-1] * 86400.0 * 365.0))
        self.kWt_avg = (Wt * self.GWhr / (1000.0 * self.time[-1] * 86400.0 * 365.0))

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

        if target == "all":
            return slice(None)  # slice all of the points
        # NOTE: PROBLEM CLAUSE FOR NOT ALLOWING GEOGRAD TO BE MORE THAN 0.7
        # Suppress warning for SBT models which use different parameter ranges
        # Only warn if value is slightly outside range (interpolation case), not way outside (SBT case)
        if target < array[0] or target > array[-1]:
            # Check if value is way outside range (likely SBT model) - suppress warning
            array_range = array[-1] - array[0]
            if target < array[0] - array_range * 0.1 or target > array[-1] + array_range * 0.1:
                # Value is significantly outside range, likely SBT model - don't warn
                pass
            else:
                # Value is slightly outside range, might be interpolation - warn
                lineprint = f"Warning: expected given value {target} to be between min and max of given array ({array[0], array[-1]})"
            # raise Exception(
            #     f"expected given value {target} to be between min and max of given array ({array[0], array[-1]})"
            # )
        for i, value in enumerate(array):
            if value == target:
                return slice(i, i + 1)
            if value > target:
                return slice(i - 1, i + 1)

    def read_values_around_point_for_interpolation(
        self, zarr_array, point, parameter_values
    ):
        """
        reads the values around the given point(s) for interpolation from the passed in zarr_array.
        If "all" is passed in as any of the coordinates of the point,
        then it'll slice across the entirety of that coordinate's dimension
        """
        indices = [
            self.get_parameter_indices(params, value)
            for value, params in zip(point, parameter_values)
        ]
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
        
        # call_id = str(uuid.uuid4())[:8]  # Unique ID for this call
        # print(f"[CALL START] {call_id} - interp_outlet_states called: fluid={self.fluid}, case={self.case}, Tinj={point[5] if len(point) > 5 else 'N/A'}K", flush=True)
        """
        :param sbt_version: 0 if not using SBT, 1 if using SBT v1, 2 if using SBT v2 
        """
        if sbt_version == 0:
            points = list(
                itertools.product(
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
            try:
                Tout = self.interpolate_points(self.Tout, point_to_read_around, points)
                Pout = self.interpolate_points(self.Pout, point_to_read_around, points)
            except ValueError as e:
                if "out of bounds" in str(e):
                    raise ValueError(f"Parameter values out of bounds for HDF5 database interpolation: {e}")
                raise
            times = self.time

            mdot, L2, L1, grad, D , Tinj, k = point
            print(mdot, L2, L1, grad, D , Tinj, k)

        else:
            mdot, L2, L1, grad, D , Tinj, k = point

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

                if mass_mode == "Constant":
                    mass_mode_b = 0
                elif mass_mode == "Variable":
                    mass_mode_b = 1

                if temp_mode == "Constant":
                    temp_mode_b = 0
                elif temp_mode == "Variable":
                    temp_mode_b = 1
                
                hyperparam1 = mass_mode_b
                hyperparam2 = "MassFlowRate.xlsx"
                hyperparam3 = temp_mode_b
                hyperparam4 = "InjectionTemperatures.xlsx"
                hyperparam5 = None
            
            if sbt_version == 2:

                reltolerance =  1e-5 
                maxnumberofiterations =  15
                fluid_mode = HyperParam5
                
                if fluid_mode == "Constant":
                    fluid_mode_b = 0
                elif fluid_mode == "Variable" or fluid_mode == "Temperature–Pressure Dependent":
                    fluid_mode_b = 1
                
                # Ensure HyperParam1 and HyperParam3 are floats 
                hyperparam1 = float(HyperParam1)*10 # Pin (convert MPa to bar) | (Inlet Pressure in MPa)
                hyperparam2 = float(HyperParam3) # pipe roughness, default is 1e-6
                hyperparam3 = fluid_mode_b # fluid mode
                hyperparam4 = reltolerance
                hyperparam5 = maxnumberofiterations

            if self.case == "coaxial":
                case = 1
                if PipeParam5 == "Inject in Annulus":
                    PipeParam5 = 1
                if PipeParam5 == "Inject in Center Pipe":
                    PipeParam5 = 2
                
                # TODO: THIS FUNCTION is NON-FUNCTIONING: check_coaxial_diameters() -- and should likely be called elsewhere

            if self.case == "utube":
                case = 2
                num_laterals = int(PipeParam3) if PipeParam3 is not None else 1
                allocation_per_lateral = 1 / num_laterals
                PipeParam4 = [allocation_per_lateral] * num_laterals # PipeParam4 is lateralflowallocation in this case | e.g., [1/3,1/3,1/3]

            # Create cache key from all parameters
            cache_key = _make_cache_key(
                sbt_version=sbt_version, mesh=mesh, hyperparam1=hyperparam1, hyperparam2=hyperparam2,
                hyperparam3=hyperparam3, hyperparam4=hyperparam4, hyperparam5=hyperparam5,
                accuracy=accuracy, case=case, mdot=mdot, Tinj=Tinj, fluid=fluid,
                L1=L1, L2=L2, Diameter1=Diameter1, Diameter2=Diameter2,
                PipeParam3=PipeParam3, PipeParam4=PipeParam4, PipeParam5=PipeParam5,
                Tsurf=Tsurf, grad=grad, k=k, c_m=c_m, rho_m=rho_m
            )
            cache_key_saved = str(cache_key) 
            # if isinstance(cache_key, tuple):
            #     cache_key_str = str(cache_key)[:200]  # First 200 chars
            #     print(f"[DEBUG] Cache key (first 200 chars): {cache_key_str}...", flush=True)
            # else:
            #     print(f"[DEBUG] Cache key type: {type(cache_key)}, value: {cache_key}, saved as: {cache_key_saved}", flush=True)

            cached_result = _get_cached_result(cache_key_saved)
            if cached_result is not None:
                tout_max_cached = np.max(cached_result[1]) - 273.15
                times, Tout, Pout = cached_result
            else:
                start = time.time()

                try:
                    print("SBT VERSION: ", sbt_version)
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

                    # if case == 1:
                    #     if fluid == 1:
                            # print(" ************* ")
                            # print(f"[run-sbt-params] {call_id} - All parameters for SBT:", flush=True)
                            # print(f"[run-sbt-params] {call_id} -   sbt_version={sbt_version}, mesh_fineness={mesh}, accuracy={accuracy}", flush=True)
                            # print(f"[run-sbt-params] {call_id} -   variableflowrate={hyperparam1}, flowratefilename={hyperparam2}, variableinjectiontemperature={hyperparam3}, injectiontemperaturefilename={hyperparam4}, HYPERPARAM5={hyperparam5}", flush=True)
                            # print(f"[run-sbt-params] {call_id} -   clg_configuration={case}, mdot={mdot}, Tinj={Tinj}, fluid={fluid}", flush=True)
                            # print(f"[run-sbt-params] {call_id} -   DrillingDepth_L1={L1}, HorizontalExtent_L2={L2}", flush=True)
                            # print(f"[run-sbt-params] {call_id} -   Diameter1(Annulus)={Diameter1}, Diameter2 (Center Pipe Radius)={Diameter2}, PipeParam3 (thicknesscenterpipe)={PipeParam3}, PipeParam4 (k_center_pipe)={PipeParam4}, PipeParam5(coaxialflowtype)={PipeParam5}", flush=True)
                            # print(f"[run-sbt-params] {call_id} -   Tsurf={Tsurf}, GeoGradient={grad}, k_m={k}, c_m={c_m}, rho_m={rho_m}", flush=True)

                    # print(f"[DEBUG] Calling SBT with: sbt_version={sbt_version}, mesh={mesh}, accuracy={accuracy}", flush=True)
                    # print(f"[DEBUG]   HYPERPARAM1={hyperparam1}, HYPERPARAM2={hyperparam2}, HYPERPARAM3={hyperparam3}, HYPERPARAM4={hyperparam4}, HYPERPARAM5={hyperparam5}", flush=True)
                    # print(f"[DEBUG]   mdot={mdot}, Tinj={Tinj:.1f}°C (already converted from K to °C), fluid={fluid} (1=H2O, 2=sCO2), case={case} (utube/coaxial)", flush=True)
                    # print(f"[DEBUG]   L1={L1} km, L2={L2} km, grad={grad}°C/m", flush=True)
                    # print(f"[DEBUG]   Diameter1={Diameter1}, Diameter2={Diameter2}, PipeParam3={PipeParam3}, PipeParam4={PipeParam4}, PipeParam5={PipeParam5}", flush=True)
                    # print(f"[DEBUG]   Tsurf={Tsurf}°C, k_m={k}, c_m={c_m}, rho_m={rho_m}", flush=True)
                    print(len(times))
                    print(Tout[-1])
                    print(Pout)

                    # Cache the result if valid (before validation)
                    if Tout is not None and len(Tout) > 0:
                        tout_max_before_cache = np.max(Tout) - 273.15
                        final_cache_key = str(cache_key_saved)
                        _set_cached_result(final_cache_key, (times, Tout, Pout))
                        cached_verify = _get_cached_result(final_cache_key)
                        end = time.time()
                except Exception as e:
                    # Re-raise the exception so the calling code knows it failed
                    raise
            
            # Validate Tout immediately after simulation (for both cached and new results)
            if Tout is not None and len(Tout) > 0:
                Tout_arr = np.array(Tout)
                if (np.any(np.isnan(Tout_arr)) or np.any(np.isinf(Tout_arr)) or 
                    np.any(Tout_arr < 200) or np.any(Tout_arr > 1000)):
                    # Remove invalid result from cache
                    if cache_key in _sbt_cache:
                        del _sbt_cache[cache_key]
                    # Return empty arrays to signal failure
                    return np.array([]), np.array([]), np.array([])
            
            if Pout is None:
                constant_pressure = 2e7 # 200 Bar in pascal || 2.09e7 
                # constant_pressure = 22228604.37405011
                Pout = constant_pressure * np.ones_like(Tout)
            
            end = time.time()
            
            times = times[14:]
            Tout = Tout[14:]
            Pout = Pout[14:]
            
            # Note: SBT model already returns times in years (converted from seconds)
            # Debug: Check time range for SBT2
            if sbt_version == 2 and len(times) > 0:
                pass

        return Tout, Pout, times

    def interp_outlet_states_contour(self, param, point):

        try:
            var_index = None
            if param == "Horizontal Extent (m)":
                points = list(
                    itertools.product(
                        self.mdot,
                        self.L2,
                        (point[1],),
                        (point[2],),
                        (point[3],),
                        (point[4],),
                        (point[5],),
                        (point[6],),
                    )
                )
                var_index = 1
            if param == "Vertical Extent (m)":
                points = list(
                    itertools.product(
                        self.mdot,
                        (point[0],),
                        self.L1,
                        (point[2],),
                        (point[3],),
                        (point[4],),
                        (point[5],),
                        (point[6],),
                    )
                )
                var_index = 2
            if param == "Geothermal Gradient (°C/m)":
                points = list(
                    itertools.product(
                        self.mdot,
                        (point[0],),
                        (point[1],),
                        self.grad,
                        (point[3],),
                        (point[4],),
                        (point[5],),
                        (point[6],),
                    )
                )
                var_index = 3
            if param == "Borehole Diameter (m)":
                points = list(
                    itertools.product(
                        self.mdot,
                        (point[0],),
                        (point[1],),
                        (point[2],),
                        self.D,
                        (point[4],),
                        (point[5],),
                        (point[6],),
                    )
                )
                var_index = 4
            if param == "Injection Temperature (˚C)":
                points = list(
                    itertools.product(
                        self.mdot,
                        (point[0],),
                        (point[1],),
                        (point[2],),
                        (point[3],),
                        self.Tinj,
                        (point[5],),
                        (point[6],),
                    )
                )
                var_index = 5
            if param == "Rock Thermal Conductivity (W/m-°C)":
                points = list(
                    itertools.product(
                        self.mdot,
                        (point[0],),
                        (point[1],),
                        (point[2],),
                        (point[3],),
                        (point[4],),
                        self.k,
                        (point[6],),
                    )
                )
                var_index = 6

            N_DIMENSIONS = 8
            points_to_read_around = [None] * N_DIMENSIONS
            point_index = 0
            # the passed in point contains the coordinate of the value we're slicing over in conjuction with mass flow rate.
            # we don't want to include that point the points to fetch, so we're removing it from the list of parameters to read at a point
            # and replacing that with "all", so the function will slice over that entire dimension
            parameters_to_read_from_point = list(point)
            del parameters_to_read_from_point[var_index - 1]
            for i in range(len(points_to_read_around)):
                if i == 0 or i == var_index:
                    points_to_read_around[i] = "all"
                else:
                    points_to_read_around[i] = parameters_to_read_from_point[point_index]
                    point_index += 1
            points_to_read_around = tuple(
                points_to_read_around
            )  # interpolate_points expects this to be a tuple so converting it

            Tout = self.interpolate_points(self.Tout, points_to_read_around, points)
            Pout = self.interpolate_points(self.Pout, points_to_read_around, points)

            Tout = np.transpose(
                np.reshape(Tout, (len(self.mdot), len(self.ivars[var_index])))
            )  
            Pout = np.transpose(
                np.reshape(Pout, (len(self.mdot), len(self.ivars[var_index])))
            )
        except Exception:
            print("Flag: Check if SBT model selected")
            # AB: Dummy data for the contours:
            Tout = np.full((20, 26), 365)
            Pout = np.full((20, 26), 22228604)

        return Tout, Pout

    def interp_kW(
        self,
        point,
        Tout,
        Pout,
        Tamb=300.0, # AB: is this surface temperature? What is Tamb?
        Pinj=None,  # Optional inlet pressure override (for SBT models)
    ):

        mdot = point[0]
        Tinj = point[5]
        
        # Use provided Pinj if given (for SBT models), otherwise use database value
        inlet_pressure = Pinj if Pinj is not None else self.Pinj

        enthalpy_out = CP.PropsSI("H", "T", Tout, "P", Pout, self.CP_fluid)
        enthalpy_in = CP.PropsSI("H", "T", Tinj, "P", inlet_pressure, self.CP_fluid)
        entropy_out = CP.PropsSI("S", "T", Tout, "P", Pout, self.CP_fluid)
        entropy_in = CP.PropsSI("S", "T", Tinj, "P", inlet_pressure, self.CP_fluid)

        kWe = (
            mdot
            * (enthalpy_out - enthalpy_in - Tamb * (entropy_out - entropy_in))
            / 1000.0
        )
        kWt = mdot * (enthalpy_out - enthalpy_in) / 1000.0

        return kWe, kWt

    def interp_kW_contour(self, param, point, Tout, Pout, Tamb=300.0, Pinj=None):

        mdot = self.mdot
        Tinj = point[4]
        
        # TODO AB: CHECK ! ***********
        # Use provided Pinj if given (for SBT models), otherwise use database value
        inlet_pressure = Pinj if Pinj is not None else self.Pinj

        accept_one_dim_only = False

        try:  # PropsSI accepts multidim for Python 3.10 only

            enthalpy_out = CP.PropsSI("H", "T", Tout, "P", Pout, self.CP_fluid)
            enthalpy_in = CP.PropsSI("H", "T", Tinj, "P", inlet_pressure, self.CP_fluid)
            entropy_out = CP.PropsSI("S", "T", Tout, "P", Pout, self.CP_fluid)
            entropy_in = CP.PropsSI("S", "T", Tinj, "P", inlet_pressure, self.CP_fluid)

        except Exception:  # PropsSI accepts only one dim for Python 3.8, 3.9. 3.11

            accept_one_dim_only = True
            dim_one = Tout.shape[0]
            dim_two = Tout.shape[1]

            Tout = Tout.flatten(order="C")
            Pout = Pout.flatten(order="C")
            enthalpy_out = CP.PropsSI("H", "T", Tout, "P", Pout, self.CP_fluid)
            enthalpy_in = CP.PropsSI("H", "T", Tinj, "P", inlet_pressure, self.CP_fluid)
            entropy_out = CP.PropsSI("S", "T", Tout, "P", Pout, self.CP_fluid)
            entropy_in = CP.PropsSI("S", "T", Tinj, "P", inlet_pressure, self.CP_fluid)

        if not accept_one_dim_only:

            kWe = (
                mdot
                * (enthalpy_out - enthalpy_in - Tamb * (entropy_out - entropy_in))
                / 1000.0
            )
            kWt = mdot * (enthalpy_out - enthalpy_in) / 1000.0

        else:

            # Suppress warnings for NaN values in calculations (will be handled by contour NaN filling)
            with np.errstate(invalid='ignore'):
                calc = (
                    enthalpy_out - enthalpy_in - Tamb * (entropy_out - entropy_in)
                ) / 1000.0
            calc = np.reshape(calc, (dim_one, dim_two), order="C")
            kWe = mdot * calc
            
            # TODO: AB CHECK by with with np.errstate is added?
            with np.errstate(invalid='ignore'):
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
