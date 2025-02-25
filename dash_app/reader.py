#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# sourced scripts
import clgs as clgs
from paths import inpath_dict
from data.decompress_hdf5 import decompress__hdf5
import os

def initialize_data():

	filename = inpath_dict["decompressed_h5_filepath"]
	need_to_make_decompressed_hdf5 = not (os.path.exists(filename))
	if need_to_make_decompressed_hdf5:
		#if running locally, decompress it
		if os.getenv("deployment_type", "aws") == "local":
			decompress__hdf5()

	u_sCO2 = clgs.data(filename, "utube", "sCO2")
	u_H2O = clgs.data(filename, "utube", "H2O")
	c_sCO2 = clgs.data(filename, "coaxial", "sCO2")
	c_H2O = clgs.data(filename, "coaxial", "H2O")

	# print('Read in data.')

	return u_sCO2, u_H2O, c_sCO2, c_H2O


def data_dict(u_sCO2, u_H2O, c_sCO2, c_H2O):

	param_dict = {("utube", "sCO2", "mdot"): u_sCO2.mdot,
	              ("utube", "sCO2", "Horizontal Extent (m)"): u_sCO2.L2,
	              ("utube", "sCO2", "Vertical Extent (m)"): u_sCO2.L1,
	              ("utube", "sCO2", "Geothermal Gradient (K/m)"): u_sCO2.grad,
	              ("utube", "sCO2", "Borehole Diameter (m)"): u_sCO2.D,
	              ("utube", "sCO2", "Injection Temperature (˚C)"): u_sCO2.Tinj,
	              ("utube", "sCO2", "Rock Thermal Conductivity (W/m-K)"): u_sCO2.k,

	              ("utube", "H2O", "mdot"): u_H2O.mdot,
	              ("utube", "H2O", "Horizontal Extent (m)"): u_H2O.L2,
	              ("utube", "H2O", "Vertical Extent (m)"): u_H2O.L1,
	              ("utube", "H2O", "Geothermal Gradient (K/m)"): u_H2O.grad,
	              ("utube", "H2O", "Borehole Diameter (m)"): u_H2O.D,
	              ("utube", "H2O", "Injection Temperature (˚C)"): u_H2O.Tinj,
	              ("utube", "H2O", "Rock Thermal Conductivity (W/m-K)"): u_H2O.k,

	              ("coaxial", "sCO2", "mdot"): c_sCO2.mdot,
	              ("coaxial", "sCO2", "Horizontal Extent (m)"): c_sCO2.L2,
	              ("coaxial", "sCO2", "Vertical Extent (m)"): c_sCO2.L1,
	              ("coaxial", "sCO2", "Geothermal Gradient (K/m)"): c_sCO2.grad,
	              ("coaxial", "sCO2", "Borehole Diameter (m)"): c_sCO2.D,
	              ("coaxial", "sCO2", "Injection Temperature (˚C)"): c_sCO2.Tinj,
	              ("coaxial", "sCO2", "Rock Thermal Conductivity (W/m-K)"): c_sCO2.k,

	              ("coaxial", "H2O", "mdot"): c_H2O.mdot,
	              ("coaxial", "H2O", "Horizontal Extent (m)"): c_H2O.L2,
	              ("coaxial", "H2O", "Vertical Extent (m)"): c_H2O.L1,
	              ("coaxial", "H2O", "Geothermal Gradient (K/m)"): c_H2O.grad,
	              ("coaxial", "H2O", "Borehole Diameter (m)"): c_H2O.D,
	              ("coaxial", "H2O", "Injection Temperature (˚C)"): c_H2O.Tinj,
	              ("coaxial", "H2O", "Rock Thermal Conductivity (W/m-K)"): c_H2O.k,          
	            }
	
	return param_dict




