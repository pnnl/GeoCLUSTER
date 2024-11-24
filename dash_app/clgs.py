#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# data manipulation libraries
import numpy as np
import h5py
from scipy.interpolate import interpn
import CoolProp.CoolProp as CP
import itertools as iter
from pympler import asizeof
import zarr
import time

class data:
  def __init__(self, fname, case, fluid):

    self.fluid = fluid
    self.case = case

    with h5py.File(fname, 'r') as file:
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
      self.ivars = (self.mdot, self.L2, self.L1, self.grad, self.D, self.Tinj, self.k, self.time)
      
      # fixed vars
      self.Pinj = file[fixed_loc + "Pinj"][()]
      self.Tamb = file[fixed_loc + "Tamb"][()]

      # dim = Mdot x L2 x L1 x grad x D x Tinj x k
      self.Wt = file[output_loc + "Wt"][:]  # int mdot * dh dt
      self.We = file[output_loc + "We"][:]  # int mdot * (dh - Too * ds) dt

      self.GWhr = 1e6 * 3600000.0

      self.kWe_avg = self.We * self.GWhr / (1000. * self.time[-1] * 86400. * 365.)
      self.kWt_avg = self.Wt * self.GWhr / (1000. * self.time[-1] * 86400. * 365.)

      # dim = Mdot x L2 x L1 x grad x D x Tinj x k x time
      self.shape = (
          len(self.mdot),
          len(self.L2),
          len(self.L1),
          len(self.grad),
          len(self.D),
          len(self.Tinj),
          len(self.k),
          len(self.time))
      print(case, fluid)
      self.Tout = zarr.open(f"data/{case}_{fluid}_Tout.zarr", mode="r")
      self.Pout = zarr.open(f"data/{case}_{fluid}_Pout.zarr", mode="r")

      print("sizes", self.kWe_avg.shape, self.Tout.shape)
      print("dtypes", self.kWe_avg.dtype, self.Tout.dtype)

    self.CP_fluid = "CO2"
    if (fluid == "H2O"):
      self.CP_fluid = "H2O"
      
  def get_parameter_indices(self, array, target):
    """
    returns a slice of the parameter indices around a given point.
    I.E array = [1, 2, 3] and target = 1.5, returns slice(0, 2)

    """
    
    if target == "all":
        return slice(None)
    if target < array[0] or target > array[-1]:
        raise Exception(f"expected given value {target} to be between min and max of given array ({array[0], array[-1]})")
    for i, value in enumerate(array):
        if value == target:
            return slice(i, i+1)
        if value > target:
            return slice(i-1, i+1)

  def read_values_around_point_for_interpolation(self, zarr_array, point, parameter_values):
    indices = [self.get_parameter_indices(params, value) for value, params in zip(point, parameter_values)]
    data = zarr_array[tuple(indices)]
    return indices, data


  def __uncompress(self, file, output_loc, state):
    U = file[output_loc + state + "/" + "U"][:]
    sigma = file[output_loc + state + "/" + "sigma"][:]
    Vt = file[output_loc + state + "/" + "Vt"][:]

    M_k = np.dot(U, np.dot(np.diag(sigma), Vt))

    shape = self.shape
    valid_runs = np.argwhere(np.isfinite(self.We.flatten()))[:, 0]
    M_k_full = np.full((shape[-1], np.prod(shape[:-1])), np.nan)
    M_k_full[:, valid_runs] = M_k
    ans = np.reshape(M_k_full.T, shape)
    #cleaning up
    del U, sigma, Vt, M_k
    return ans

  def interpolate_points(self, zarr_array, point_to_read_around, points):
    def read_values_around_points_for_interpolation(zarr_array, point_to_read_around, parameter_values):
        indices = [self.get_parameter_indices(params, value) for value, params in zip(point_to_read_around, parameter_values)]
        read_in_data = zarr_array[tuple(indices)]
        return indices, read_in_data

    indices, values_around_point = read_values_around_points_for_interpolation(zarr_array, point_to_read_around, self.ivars)
    grid = [params[these_indices] for these_indices, params in zip(indices, self.ivars)]
    interpolated_points = interpn(grid, values_around_point, points)
    return interpolated_points

  def interp_outlet_states(self, point):
    start_time = time.time()

    points = list(iter.product(
            (point[0],),
            (point[1],),
            (point[2],),
            (point[3],),
            (point[4],),
            (point[5],),
            (point[6],),
            self.time))
    

  #   print("--------")
  #   print(self.ivars)
  #   print(point)
  #   Tout = interpn(self.ivars, self.Tout, points)
  #   Pout = interpn(self.ivars, self.Pout, points)
  #   print("------")
  #   print("Tout", Tout[-1])
  #   print(len(self.time), len(Tout))

    point_to_read_around =(*point, "all") # unpacking point
    Tout = self.interpolate_points(self.Tout, point_to_read_around, points)
    Pout = self.interpolate_points(self.Pout, point_to_read_around, points)

    print("--- %s seconds ---" % (time.time() - start_time))
    return Tout, Pout

  def interp_outlet_states_contour(self, param, point):
    print("param", param, point)

    var_index = None
    if param == "Horizontal Extent (m)":
      points = list(iter.product(
            self.mdot,
            self.L2,
            (point[1],),
            (point[2],),
            (point[3],),
            (point[4],),
            (point[5],),
            (point[6],)
            )) 
      var_index = 1
    if param == "Vertical Extent (m)":
      points = list(iter.product(
            self.mdot,
            (point[0],),
            self.L1,
            (point[2],),
            (point[3],),
            (point[4],),
            (point[5],),
            (point[6],)
            )) 
      var_index = 2
    if param == "Geothermal Gradient (K/m)":
      points = list(iter.product(
            self.mdot,
            (point[0],),
            (point[1],),
            self.grad,
            (point[3],),
            (point[4],),
            (point[5],),
            (point[6],)
            )) 
      var_index = 3
    if param == "Borehole Diameter (m)":
      points = list(iter.product(
            self.mdot,
            (point[0],),
            (point[1],),
            (point[2],),
            self.D,
            (point[4],),
            (point[5],),
            (point[6],)
            )) 
      var_index = 4
    if param == "Injection Temperature (˚C)":
      points = list(iter.product(
            self.mdot,
            (point[0],),
            (point[1],),
            (point[2],),
            (point[3],),
            self.Tinj,
            (point[5],),
            (point[6],)
            )) 
      var_index = 5
    if param == "Rock Thermal Conductivity (W/m-K)":
      points = list(iter.product(
            self.mdot,
            (point[0],),
            (point[1],),
            (point[2],),
            (point[3],),
            (point[4],),
            self.k,
            (point[6],)
            )) 
      var_index = 6

    # converting to a list so we can mutate it
    N_DIMENSIONS = 8
    points_to_read_around = [None] * N_DIMENSIONS
    point_index = 0
    # the passed in point contains the coordinate of the value we're slicing over in conjuction with mass flow rate.
    # we don't want to include that point the points to fetch, so we're removing it from the list of parameters to read at a point
    parameters_to_read_from_point = list(point)
    del parameters_to_read_from_point[var_index - 1]
    for i in range(len(points_to_read_around)):
      if i == 0 or i == var_index:
        points_to_read_around[i] = "all"
      else:
        points_to_read_around[i] = parameters_to_read_from_point[point_index]
        point_index += 1
    points_to_read_around = tuple(points_to_read_around)

    print("points_to_read_around", points_to_read_around)
    print("og point", point)
    
    Tout = self.interpolate_points(self.Tout, points_to_read_around, points)
    Pout = self.interpolate_points(self.Pout, points_to_read_around, points)


    if param == "Horizontal Extent (m)":
      Tout = np.transpose(np.reshape(Tout, (self.Tout.shape[0], self.Tout.shape[1]))) # 20,26
      Pout = np.transpose(np.reshape(Pout, (self.Tout.shape[0], self.Tout.shape[1])))

    if param == "Vertical Extent (m)":
      Tout = np.transpose(np.reshape(Tout, (self.Tout.shape[0], self.Tout.shape[2]))) # 20,9
      Pout = np.transpose(np.reshape(Pout, (self.Tout.shape[0], self.Tout.shape[2])))

    if param == "Geothermal Gradient (K/m)":
      Tout = np.transpose(np.reshape(Tout, (self.Tout.shape[0], self.Tout.shape[3]))) # 20,5
      Pout = np.transpose(np.reshape(Pout, (self.Tout.shape[0], self.Tout.shape[3])))
    
    if param == "Borehole Diameter (m)":
      Tout = np.transpose(np.reshape(Tout, (self.Tout.shape[0], self.Tout.shape[4]))) # 20,3
      Pout = np.transpose(np.reshape(Pout, (self.Tout.shape[0], self.Tout.shape[4])))

    if param == "Injection Temperature (˚C)":
      Tout = np.transpose(np.reshape(Tout, (self.Tout.shape[0], self.Tout.shape[5]))) # 20,3
      Pout = np.transpose(np.reshape(Pout, (self.Tout.shape[0], self.Tout.shape[5])))

    if param == "Rock Thermal Conductivity (W/m-K)":
      Tout = np.transpose(np.reshape(Tout, (self.Tout.shape[0], self.Tout.shape[6]))) # 20,3
      Pout = np.transpose(np.reshape(Pout, (self.Tout.shape[0], self.Tout.shape[6])))


    print("return tout from intepr outlet states shape", Tout.shape)
    return Tout, Pout


  def interp_kW(self, point,  Tout, Pout, Tamb=300.0,):
    
    mdot = point[0]
    Tinj = point[5]
    enthalpy_out = CP.PropsSI('H', 'T', Tout, 'P', Pout, self.CP_fluid)
    enthalpy_in = CP.PropsSI('H', 'T', Tinj, 'P', self.Pinj, self.CP_fluid)
    entropy_out = CP.PropsSI('S', 'T', Tout, 'P', Pout, self.CP_fluid)
    entropy_in = CP.PropsSI('S', 'T', Tinj, 'P', self.Pinj, self.CP_fluid)
    kWe = mdot * (enthalpy_out - enthalpy_in - Tamb * (entropy_out - entropy_in)) / 1000.
    kWt = mdot * (enthalpy_out - enthalpy_in) / 1000.
    
    return kWe, kWt

  def interp_kW_contour(self, param, point, Tout, Pout, Tamb=300.0):
    
    mdot = self.mdot
    Tinj = point[4]

    accept_one_dim_only = False
    
    try: # PropsSI accepts multidim for Python 3.10 only

      enthalpy_out = CP.PropsSI('H', 'T', Tout, 'P', Pout, self.CP_fluid)
      enthalpy_in = CP.PropsSI('H', 'T', Tinj, 'P', self.Pinj, self.CP_fluid)
      entropy_out = CP.PropsSI('S', 'T', Tout, 'P', Pout, self.CP_fluid)
      entropy_in = CP.PropsSI('S', 'T', Tinj, 'P', self.Pinj, self.CP_fluid)
    
    except Exception: # PropsSI accepts only one dim for Python 3.8, 3.9. 3.11

      accept_one_dim_only = True
      dim_one = Tout.shape[0]
      dim_two = Tout.shape[1]

      Tout = Tout.flatten(order='C')
      Pout = Pout.flatten(order='C')
      enthalpy_out = CP.PropsSI('H', 'T', Tout, 'P', Pout, self.CP_fluid)
      enthalpy_in = CP.PropsSI('H', 'T', Tinj, 'P', self.Pinj, self.CP_fluid)
      entropy_out = CP.PropsSI('S', 'T', Tout, 'P', Pout, self.CP_fluid)
      entropy_in = CP.PropsSI('S', 'T', Tinj, 'P', self.Pinj, self.CP_fluid)

    if not accept_one_dim_only:

      kWe = mdot * (enthalpy_out - enthalpy_in - Tamb * (entropy_out - entropy_in)) / 1000.
      kWt = mdot * (enthalpy_out - enthalpy_in) / 1000.
    
    else:

      calc = (enthalpy_out - enthalpy_in - Tamb * (entropy_out - entropy_in)) / 1000.
      calc = np.reshape(calc, (dim_one, dim_two), order='C') 
      kWe = mdot * calc
      calc = (enthalpy_out - enthalpy_in) / 1000.
      calc = np.reshape(calc, (dim_one, dim_two), order='C') 
      kWt = mdot * calc

    return kWe, kWt

  def interp_kWe_avg(self, point):
    ivars = self.ivars[:-1]
    return self.GWhr * interpn(ivars, self.We, point) / (1000. * self.time[-1] * 86400. * 365.)

  def interp_kWt_avg(self, point):
    ivars = self.ivars[:-1]
    return self.GWhr * interpn(ivars, self.Wt, point) / (1000. * self.time[-1] * 86400. * 365.)

  def compute_kW(self, index, Tamb = 300.0):
    mdot = self.mdot[index[0]]
    Tinj = self.Tinj[index[5]]
    point = (index[0], index[1], index[2], index[3], index[4], index[5], index[6], ...)
    Tout = self.Tout[point]
    Pout = self.Pout[point]
    enthalpy_out = CP.PropsSI('H', 'T', Tout, 'P', Pout, self.CP_fluid)
    enthalpy_in = CP.PropsSI('H', 'T', Tinj, 'P', self.Pinj, self.CP_fluid)
    entropy_out = CP.PropsSI('S', 'T', Tout, 'P', Pout, self.CP_fluid)
    entropy_in = CP.PropsSI('S', 'T', Tinj, 'P', self.Pinj, self.CP_fluid)
    kWe = mdot * (enthalpy_out - enthalpy_in - Tamb * (entropy_out - entropy_in)) / 1000.
    kWt = mdot * (enthalpy_out - enthalpy_in) / 1000.

    return kWe, kWt