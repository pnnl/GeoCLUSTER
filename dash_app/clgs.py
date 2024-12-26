#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# data manipulation libraries
import numpy as np
import h5py
from scipy.interpolate import interpn
import CoolProp.CoolProp as CP
import itertools as iter
import zarr
from paths import absolute_path


class data:
    def __init__(self, fname, case, fluid):

        self.fluid = fluid
        self.case = case

        with h5py.File(fname, "r") as file:
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
            self.Wt = file[output_loc + "Wt"][:]  # int mdot * dh dt
            self.We = file[output_loc + "We"][:]  # int mdot * (dh - Too * ds) dt

            self.GWhr = 1e6 * 3600000.0

            self.kWe_avg = (
                self.We * self.GWhr / (1000.0 * self.time[-1] * 86400.0 * 365.0)
            )
            self.kWt_avg = (
                self.Wt * self.GWhr / (1000.0 * self.time[-1] * 86400.0 * 365.0)
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

            # if you get an error that these files don't exist, run the zarr.ipynb notebook in the data directory to build these files!
            self.Tout = zarr.open(f"{absolute_path}/data/{case}_{fluid}_Tout.zarr", mode="r")
            self.Pout = zarr.open(f"{absolute_path}/data/{case}_{fluid}_Pout.zarr", mode="r")

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
        if target < array[0] or target > array[-1]:
            raise Exception(
                f"expected given value {target} to be between min and max of given array ({array[0], array[-1]})"
            )
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
        print(point)
        print(indices)
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

    def interp_outlet_states(self, point):
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

        return Tout, Pout

    def interp_outlet_states_contour(self, param, point):
        var_index = None
        if param == "Horizontal Extent (m)":
            points = list(
                iter.product(
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
                iter.product(
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
        if param == "Geothermal Gradient (K/m)":
            points = list(
                iter.product(
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
                iter.product(
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
                iter.product(
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
        if param == "Rock Thermal Conductivity (W/m-K)":
            points = list(
                iter.product(
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

        if param == "Horizontal Extent (m)":
            Tout = np.transpose(
                np.reshape(Tout, (self.Tout.shape[0], self.Tout.shape[1]))
            )  # 20,26
            Pout = np.transpose(
                np.reshape(Pout, (self.Tout.shape[0], self.Tout.shape[1]))
            )

        if param == "Vertical Extent (m)":
            Tout = np.transpose(
                np.reshape(Tout, (self.Tout.shape[0], self.Tout.shape[2]))
            )  # 20,9
            Pout = np.transpose(
                np.reshape(Pout, (self.Tout.shape[0], self.Tout.shape[2]))
            )

        if param == "Geothermal Gradient (K/m)":
            Tout = np.transpose(
                np.reshape(Tout, (self.Tout.shape[0], self.Tout.shape[3]))
            )  # 20,5
            Pout = np.transpose(
                np.reshape(Pout, (self.Tout.shape[0], self.Tout.shape[3]))
            )

        if param == "Borehole Diameter (m)":
            Tout = np.transpose(
                np.reshape(Tout, (self.Tout.shape[0], self.Tout.shape[4]))
            )  # 20,3
            Pout = np.transpose(
                np.reshape(Pout, (self.Tout.shape[0], self.Tout.shape[4]))
            )

        if param == "Injection Temperature (˚C)":
            Tout = np.transpose(
                np.reshape(Tout, (self.Tout.shape[0], self.Tout.shape[5]))
            )  # 20,3
            Pout = np.transpose(
                np.reshape(Pout, (self.Tout.shape[0], self.Tout.shape[5]))
            )

        if param == "Rock Thermal Conductivity (W/m-K)":
            Tout = np.transpose(
                np.reshape(Tout, (self.Tout.shape[0], self.Tout.shape[6]))
            )  # 20,3
            Pout = np.transpose(
                np.reshape(Pout, (self.Tout.shape[0], self.Tout.shape[6]))
            )

        return Tout, Pout

    def interp_kW(
        self,
        point,
        Tout,
        Pout,
        Tamb=300.0,
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
