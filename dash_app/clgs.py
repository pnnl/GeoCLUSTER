#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# data manipulation libraries
import numpy as np
import h5py
from scipy.interpolate import interpn
import CoolProp.CoolProp as CP
import itertools as iter

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
      self.Tout = self.__uncompress(file, output_loc, "Tout")
      self.Pout = self.__uncompress(file, output_loc, "Pout")

    self.CP_fluid = "CO2"
    if (fluid == "H2O"):
      self.CP_fluid = "H2O"

  def __uncompress(self, file, output_loc, state):
    U = file[output_loc + state + "/" + "U"][:]
    sigma = file[output_loc + state + "/" + "sigma"][:]
    Vt = file[output_loc + state + "/" + "Vt"][:]
    M_k = np.dot(U, np.dot(np.diag(sigma), Vt))

    shape = self.shape
    valid_runs = np.argwhere(np.isfinite(self.We.flatten()))[:, 0]
    M_k_full = np.full((shape[-1], np.prod(shape[:-1])), np.nan)
    M_k_full[:, valid_runs] = M_k
    return np.reshape(M_k_full.T, shape)

  def interp_outlet_states(self, point):

    points = list(iter.product(
            (point[0],),
            (point[1],),
            (point[2],),
            (point[3],),
            (point[4],),
            (point[5],),
            (point[6],),
            self.time))
    Tout = interpn(self.ivars, self.Tout, points)
    Pout = interpn(self.ivars, self.Pout, points)

    return Tout, Pout

  def interp_outlet_states_contour(self, param, point):

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

    # print(self.ivars) # tuple
    # print(self.Tout.shape) # e.g., (26, 20, 9, 5, 3, 3, 3, 161)
    # print(points) # list

    # points=ivars (tuple)values=Tout, points=the coordinates 
    Tout = interpn(self.ivars, self.Tout, points) # mass flow vs. param
    Pout = interpn(self.ivars, self.Pout, points)

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


    return Tout, Pout


  def interp_kW(self, point, Tamb=300.0):
    
    mdot = point[0]
    Tinj = point[5]
    Tout, Pout = self.interp_outlet_states(point) # one dimensional output
    enthalpy_out = CP.PropsSI('H', 'T', Tout, 'P', Pout, self.CP_fluid)
    enthalpy_in = CP.PropsSI('H', 'T', Tinj, 'P', self.Pinj, self.CP_fluid)
    entropy_out = CP.PropsSI('S', 'T', Tout, 'P', Pout, self.CP_fluid)
    entropy_in = CP.PropsSI('S', 'T', Tinj, 'P', self.Pinj, self.CP_fluid)
    kWe = mdot * (enthalpy_out - enthalpy_in - Tamb * (entropy_out - entropy_in)) / 1000.
    kWt = mdot * (enthalpy_out - enthalpy_in) / 1000.
    
    return kWe, kWt

  def interp_kW_contour(self, param, point, Tamb=300.0):
    
    mdot = self.mdot
    Tinj = point[4]
    Tout, Pout = self.interp_outlet_states_contour(param, point) # multidimensional output

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