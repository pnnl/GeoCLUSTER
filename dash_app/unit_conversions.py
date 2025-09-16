#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Unit Conversion Module for GeoCLUSTER
Handles conversion between different unit systems for all parameters
"""

TEMPERATURE_CONVERSIONS = {
    'C': 1.0,
    'F': lambda x: (x * 9/5) + 32
}

TEMPERATURE_CONVERSIONS_INVERSE = {
    'C': 1.0,
    'F': lambda x: (x - 32) * 5/9
}

LENGTH_CONVERSIONS = {
    'm': 1.0,
    'ft': 3.28084,
    'km': 0.001,
    'mi': 0.000621371,
    'cm': 100,
    'mm': 1000,
    'in': 39.3701,
    'yd': 1.09361
}

LENGTH_CONVERSIONS_INVERSE = {
    'm': 1.0,
    'ft': 0.3048,
    'km': 1000,
    'mi': 1609.34,
    'cm': 0.01,
    'mm': 0.001,
    'in': 0.0254,
    'yd': 0.9144
}

MASS_FLOW_CONVERSIONS = {
    'kg/s': 1.0,
    'lb/s': 2.20462,
    'kg/h': 3600,
    'lb/h': 7936.64,
    'g/s': 1000,
    'ton/h': 3.6
}

MASS_FLOW_CONVERSIONS_INVERSE = {
    'kg/s': 1.0,
    'lb/s': 0.453592,
    'kg/h': 0.000277778,
    'lb/h': 0.000125998,
    'g/s': 0.001,
    'ton/h': 0.277778
}

THERMAL_CONDUCTIVITY_CONVERSIONS = {
    'W/m-K': 1.0,
    'W/m-C': 1.0,
    'Btu/ft-h-F': 0.577789,
    'Btu/yd-h-F': 1.733367,
    'cal/cm-s-C': 0.002388,
    'kcal/m-h-C': 0.859845
}

THERMAL_CONDUCTIVITY_CONVERSIONS_INVERSE = {
    'W/m-K': 1.0,
    'W/m-C': 1.0,
    'Btu/ft-h-F': 1.730735,
    'Btu/yd-h-F': 0.576912,
    'cal/cm-s-C': 418.68,
    'kcal/m-h-C': 1.163
}

HEAT_CAPACITY_CONVERSIONS = {
    'J/kg-K': 1.0,
    'J/kg-C': 1.0,
    'Btu/lb-F': 0.000238846,
    'cal/g-C': 0.000238846,
    'kJ/kg-K': 0.001,
    'kcal/kg-C': 0.000238846
}

HEAT_CAPACITY_CONVERSIONS_INVERSE = {
    'J/kg-K': 1.0,
    'J/kg-C': 1.0,
    'Btu/lb-F': 4186.8,
    'cal/g-C': 4186.8,
    'kJ/kg-K': 1000,
    'kcal/kg-C': 4186.8
}

DENSITY_CONVERSIONS = {
    'kg/m3': 1.0,
    'g/cm3': 0.001,
    'lb/ft3': 0.062428,
    'lb/yd3': 0.002309,
    'lb/gal': 0.008345,
    'g/L': 0.001,
    't/m3': 0.001
}

DENSITY_CONVERSIONS_INVERSE = {
    'kg/m3': 1.0,
    'g/cm3': 1000,
    'lb/ft3': 16.0185,
    'lb/yd3': 433.0,
    'lb/gal': 119.826,
    'g/L': 1000,
    't/m3': 1000
}

PRESSURE_CONVERSIONS = {
    'Pa': 1.0,
    'bar': 1e-5,
    'MPa': 1e-6,
    'psi': 0.000145038,
    'atm': 9.86923e-6,
    'kPa': 1e-3,
    'torr': 7.50062e-3
}

PRESSURE_CONVERSIONS_INVERSE = {
    'Pa': 1.0,
    'bar': 1e5,
    'MPa': 1e6,
    'psi': 6894.76,
    'atm': 101325,
    'kPa': 1e3,
    'torr': 133.322
}

GEOTHERMAL_GRADIENT_CONVERSIONS = {
    'K/m': 1.0,
    'C/m': 1.0,
    'C/km': 1000,
    'F/ft': 0.54864,
    'F/yd': 1.64592,
    'F/100ft': 54.864,
    'K/km': 1000
}

GEOTHERMAL_GRADIENT_CONVERSIONS_INVERSE = {
    'K/m': 1.0,
    'C/m': 1.0,
    'C/km': 0.001,
    'F/ft': 1.82269,
    'F/yd': 0.60756,
    'F/100ft': 0.018227,
    'K/km': 0.001
}

# Default unit preferences
DEFAULT_UNITS = {
    'temperature': 'C',
    'length': 'm',
    'mass_flow': 'kg/s',
    'thermal_conductivity': 'W/m-K',
    'heat_capacity': 'J/kg-K',
    'density': 'kg/m3',
    'pressure': 'Pa',
    'geothermal_gradient': 'K/m'
}

class UnitConverter:
    """Main unit conversion class for GeoCLUSTER"""
    
    def __init__(self, user_preferences=None):
        """
        Initialize unit converter with user preferences
        
        Args:
            user_preferences (dict): User's preferred units for each parameter type
        """
        self.user_preferences = user_preferences or DEFAULT_UNITS.copy()
    
    def convert_temperature(self, value, from_unit, to_unit):
        """Convert temperature between different units"""
        if from_unit == to_unit:
            return value
        
        if from_unit != 'C':
            if from_unit in TEMPERATURE_CONVERSIONS_INVERSE:
                if callable(TEMPERATURE_CONVERSIONS_INVERSE[from_unit]):
                    value = TEMPERATURE_CONVERSIONS_INVERSE[from_unit](value)
                else:
                    value = value * TEMPERATURE_CONVERSIONS_INVERSE[from_unit]
        
        if to_unit != 'C':
            if to_unit in TEMPERATURE_CONVERSIONS:
                if callable(TEMPERATURE_CONVERSIONS[to_unit]):
                    value = TEMPERATURE_CONVERSIONS[to_unit](value)
                else:
                    value = value * TEMPERATURE_CONVERSIONS[to_unit]
        
        return value
    
    def convert_length(self, value, from_unit, to_unit):
        """Convert length between different units"""
        if from_unit == to_unit:
            return value
        
        if from_unit != 'm':
            value = value * LENGTH_CONVERSIONS_INVERSE[from_unit]
        
        if to_unit != 'm':
            value = value * LENGTH_CONVERSIONS[to_unit]
        
        return value
    
    def convert_mass_flow(self, value, from_unit, to_unit):
        """Convert mass flow rate between different units"""
        if from_unit == to_unit:
            return value
        
        if from_unit != 'kg/s':
            value = value * MASS_FLOW_CONVERSIONS_INVERSE[from_unit]
        
        if to_unit != 'kg/s':
            value = value * MASS_FLOW_CONVERSIONS[to_unit]
        
        return value
    
    def convert_thermal_conductivity(self, value, from_unit, to_unit):
        """Convert thermal conductivity between different units"""
        if from_unit == to_unit:
            return value
        
        if from_unit != 'W/m-K':
            value = value * THERMAL_CONDUCTIVITY_CONVERSIONS_INVERSE[from_unit]
        
        if to_unit != 'W/m-K':
            value = value * THERMAL_CONDUCTIVITY_CONVERSIONS[to_unit]
        
        return value
    
    def convert_heat_capacity(self, value, from_unit, to_unit):
        """Convert heat capacity between different units"""
        if from_unit == to_unit:
            return value
        
        if from_unit != 'J/kg-K':
            value = value * HEAT_CAPACITY_CONVERSIONS_INVERSE[from_unit]
        
        if to_unit != 'J/kg-K':
            value = value * HEAT_CAPACITY_CONVERSIONS[to_unit]
        
        return value
    
    def convert_density(self, value, from_unit, to_unit):
        """Convert density between different units"""
        if from_unit == to_unit:
            return value
        
        if from_unit != 'kg/m3':
            value = value * DENSITY_CONVERSIONS_INVERSE[from_unit]
        
        if to_unit != 'kg/m3':
            value = value * DENSITY_CONVERSIONS[to_unit]
        
        return value
    
    def convert_pressure(self, value, from_unit, to_unit):
        """Convert pressure between different units"""
        if from_unit == to_unit:
            return value
        
        if from_unit != 'Pa':
            value = value * PRESSURE_CONVERSIONS_INVERSE[from_unit]
        
        if to_unit != 'Pa':
            value = value * PRESSURE_CONVERSIONS[to_unit]
        
        return value
    
    def convert_geothermal_gradient(self, value, from_unit, to_unit):
        """Convert geothermal gradient between different units"""
        if from_unit == to_unit:
            return value
        
        if from_unit != 'K/m':
            value = value * GEOTHERMAL_GRADIENT_CONVERSIONS_INVERSE[from_unit]
        
        if to_unit != 'K/m':
            value = value * GEOTHERMAL_GRADIENT_CONVERSIONS[to_unit]
        
        return value
    
    def convert_parameter(self, value, parameter_type, from_unit=None, to_unit=None):
        """
        Convert any parameter based on its type
        
        Args:
            value: The value to convert
            parameter_type: Type of parameter (temperature, length, etc.)
            from_unit: Source unit (if None, uses base unit)
            to_unit: Target unit (if None, uses user preference)
        
        Returns:
            Converted value
        """
        if to_unit is None:
            to_unit = self.user_preferences.get(parameter_type, DEFAULT_UNITS[parameter_type])
        
        if from_unit is None:
            from_unit = DEFAULT_UNITS[parameter_type]
        
        conversion_methods = {
            'temperature': self.convert_temperature,
            'length': self.convert_length,
            'mass_flow': self.convert_mass_flow,
            'thermal_conductivity': self.convert_thermal_conductivity,
            'heat_capacity': self.convert_heat_capacity,
            'density': self.convert_density,
            'pressure': self.convert_pressure,
            'geothermal_gradient': self.convert_geothermal_gradient
        }
        
        if parameter_type in conversion_methods:
            return conversion_methods[parameter_type](value, from_unit, to_unit)
        else:
            return value  # No conversion available
    
    def update_preferences(self, new_preferences):
        """Update user unit preferences"""
        self.user_preferences.update(new_preferences)
    
    def get_available_units(self, parameter_type):
        """Get available units for a parameter type"""
        unit_maps = {
            'temperature': TEMPERATURE_CONVERSIONS.keys(),
            'length': LENGTH_CONVERSIONS.keys(),
            'mass_flow': MASS_FLOW_CONVERSIONS.keys(),
            'thermal_conductivity': THERMAL_CONDUCTIVITY_CONVERSIONS.keys(),
            'heat_capacity': HEAT_CAPACITY_CONVERSIONS.keys(),
            'density': DENSITY_CONVERSIONS.keys(),
            'pressure': PRESSURE_CONVERSIONS.keys(),
            'geothermal_gradient': GEOTHERMAL_GRADIENT_CONVERSIONS.keys()
        }
        return list(unit_maps.get(parameter_type, []))
    
    def get_unit_symbol(self, unit):
        """Get the proper symbol for a unit"""
        unit_symbols = {
            'C': '°C',
            'F': '°F',
            'm': 'm',
            'ft': 'ft',
            'km': 'km',
            'mi': 'mi',
            'cm': 'cm',
            'mm': 'mm',
            'in': 'in',
            'yd': 'yd',
            'kg/s': 'kg/s',
            'lb/s': 'lb/s',
            'kg/h': 'kg/h',
            'lb/h': 'lb/h',
            'g/s': 'g/s',
            'ton/h': 'ton/h',
            'W/m-K': 'W/m·K',
            'W/m-C': 'W/m·°C',
            'Btu/ft-h-F': 'Btu/ft·h·°F',
            'Btu/yd-h-F': 'Btu/yd·h·°F',
            'cal/cm-s-C': 'cal/cm·s·°C',
            'kcal/m-h-C': 'kcal/m·h·°C',
            'J/kg-K': 'J/kg·K',
            'J/kg-C': 'J/kg·°C',
            'Btu/lb-F': 'Btu/lb·°F',
            'cal/g-C': 'cal/g·°C',
            'kJ/kg-K': 'kJ/kg·K',
            'kcal/kg-C': 'kcal/kg·°C',
            'kg/m3': 'kg/m³',
            'g/cm3': 'g/cm³',
            'lb/ft3': 'lb/ft³',
            'lb/yd3': 'lb/yd³',
            'lb/gal': 'lb/gal',
            'g/L': 'g/L',
            't/m3': 't/m³',
            'Pa': 'Pa',
            'bar': 'bar',
            'MPa': 'MPa',
            'psi': 'psi',
            'atm': 'atm',
            'kPa': 'kPa',
            'torr': 'torr',
            'K/m': 'K/m',
            'C/m': '°C/m',
            'C/km': '°C/km',
            'F/ft': '°F/ft',
            'F/yd': '°F/yd',
            'F/100ft': '°F/100ft',
            'K/km': 'K/km'
        }
        return unit_symbols.get(unit, unit)

# Global unit converter instance
unit_converter = UnitConverter()

def convert_value(value, parameter_type, from_unit=None, to_unit=None):
    """Global function for easy unit conversion"""
    return unit_converter.convert_parameter(value, parameter_type, from_unit, to_unit)

def get_available_units(parameter_type):
    """Global function to get available units for a parameter type"""
    return unit_converter.get_available_units(parameter_type)

def get_unit_symbol(unit):
    """Global function to get unit symbol"""
    return unit_converter.get_unit_symbol(unit)
