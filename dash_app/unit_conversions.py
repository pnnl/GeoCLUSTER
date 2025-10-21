#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Unit Conversion Module for GeoCLUSTER
Handles conversion between different unit systems for all parameters
"""

# Unit conversion constants
TEMPERATURE_CONVERSIONS = {
    'C': 1.0,           # Celsius (base unit)
    'F': lambda x: (x * 9/5) + 32   # Fahrenheit
}

TEMPERATURE_CONVERSIONS_INVERSE = {
    'C': 1.0,           # Celsius (base unit)
    'F': lambda x: (x - 32) * 5/9   # Fahrenheit to Celsius
}

LENGTH_CONVERSIONS = {
    'm': 1.0,           # Meters (base unit)
    'ft': 3.28084,      # Feet
    'km': 0.001,        # Kilometers
    'mi': 0.000621371,  # Miles
    'cm': 100,          # Centimeters
    'mm': 1000,         # Millimeters
    'in': 39.3701,      # Inches
    'yd': 1.09361       # Yards
}

LENGTH_CONVERSIONS_INVERSE = {
    'm': 1.0,           # Meters (base unit)
    'ft': 0.3048,       # Feet to meters
    'km': 1000,         # Kilometers to meters
    'mi': 1609.34,      # Miles to meters
    'cm': 0.01,         # Centimeters to meters
    'mm': 0.001,        # Millimeters to meters
    'in': 0.0254,       # Inches to meters
    'yd': 0.9144        # Yards to meters
}

MASS_FLOW_CONVERSIONS = {
    'kg/s': 1.0,        # kg/s (base unit)
    'lb/s': 2.20462,    # Pounds per second
    'kg/h': 3600,       # Kilograms per hour
    'lb/h': 7936.64,    # Pounds per hour
    'g/s': 1000,        # Grams per second
    'ton/h': 3.6        # Metric tons per hour
}

MASS_FLOW_CONVERSIONS_INVERSE = {
    'kg/s': 1.0,        # kg/s (base unit)
    'lb/s': 0.453592,   # Pounds per second to kg/s
    'kg/h': 0.000277778, # Kilograms per hour to kg/s
    'lb/h': 0.000125998, # Pounds per hour to kg/s
    'g/s': 0.001,       # Grams per second to kg/s
    'ton/h': 0.277778    # Metric tons per hour to kg/s
}

THERMAL_CONDUCTIVITY_CONVERSIONS = {
    'W/m-K': 1.0,       # W/m-K (base unit)
    'W/m-C': 1.0,       # W/m-C (same as W/m-K)
    'Btu/ft-h-F': 0.577789,  # BTU per foot-hour-Fahrenheit
    'Btu/yd-h-F': 1.733367,  # BTU per yard-hour-Fahrenheit
    'cal/cm-s-C': 0.002388,  # Calories per cm-second-Celsius
    'kcal/m-h-C': 0.859845   # Kilocalories per meter-hour-Celsius
}

THERMAL_CONDUCTIVITY_CONVERSIONS_INVERSE = {
    'W/m-K': 1.0,       # W/m-K (base unit)
    'W/m-C': 1.0,       # W/m-C (same as W/m-K)
    'Btu/ft-h-F': 1.730735,  # BTU per foot-hour-Fahrenheit to W/m-K
    'Btu/yd-h-F': 0.576912,  # BTU per yard-hour-Fahrenheit to W/m-K
    'cal/cm-s-C': 418.68,    # Calories per cm-second-Celsius to W/m-K
    'kcal/m-h-C': 1.163      # Kilocalories per meter-hour-Celsius to W/m-K
}

HEAT_CAPACITY_CONVERSIONS = {
    'J/kg-K': 1.0,      # J/kg-K (base unit)
    'J/kg-C': 1.0,      # J/kg-C (same as J/kg-K)
    'Btu/lb-F': 0.000238846, # BTU per pound-Fahrenheit
    'cal/g-C': 0.000238846,  # Calories per gram-Celsius
    'kJ/kg-K': 0.001,   # Kilojoules per kilogram-Kelvin
    'kcal/kg-C': 0.000238846 # Kilocalories per kilogram-Celsius
}

HEAT_CAPACITY_CONVERSIONS_INVERSE = {
    'J/kg-K': 1.0,      # J/kg-K (base unit)
    'J/kg-C': 1.0,      # J/kg-C (same as J/kg-K)
    'Btu/lb-F': 4186.8, # BTU per pound-Fahrenheit to J/kg-K
    'cal/g-C': 4186.8,  # Calories per gram-Celsius to J/kg-K
    'kJ/kg-K': 1000,    # Kilojoules per kilogram-Kelvin to J/kg-K
    'kcal/kg-C': 4186.8 # Kilocalories per kilogram-Celsius to J/kg-K
}

DENSITY_CONVERSIONS = {
    'kg/m3': 1.0,       # kg/m³ (base unit)
    'g/cm3': 0.001,     # g/cm³
    'lb/ft3': 0.062428, # lb/ft³
    'lb/yd3': 0.002309, # lb/yd³
    'lb/gal': 0.008345, # lb/gal (US)
    'g/L': 0.001,       # g/L
    't/m3': 0.001       # metric tons/m³
}

DENSITY_CONVERSIONS_INVERSE = {
    'kg/m3': 1.0,       # kg/m³ (base unit)
    'g/cm3': 1000,      # g/cm³ to kg/m³
    'lb/ft3': 16.0185,  # lb/ft³ to kg/m³
    'lb/yd3': 433.0,    # lb/yd³ to kg/m³
    'lb/gal': 119.826,  # lb/gal (US) to kg/m³
    'g/L': 1000,        # g/L to kg/m³
    't/m3': 1000        # metric tons/m³ to kg/m³
}

PRESSURE_CONVERSIONS = {
    'Pa': 1.0,          # Pascal (base unit)
    'bar': 1e-5,        # Bar
    'MPa': 1e-6,        # Megapascal
    'psi': 0.000145038, # Pounds per square inch
    'atm': 9.86923e-6,  # Atmospheres
    'kPa': 1e-3,        # Kilopascal
    'torr': 7.50062e-3  # Torr
}

PRESSURE_CONVERSIONS_INVERSE = {
    'Pa': 1.0,          # Pascal (base unit)
    'bar': 1e5,         # Bar to Pa
    'MPa': 1e6,         # Megapascal to Pa
    'psi': 6894.76,     # Pounds per square inch to Pa
    'atm': 101325,      # Atmospheres to Pa
    'kPa': 1e3,         # Kilopascal to Pa
    'torr': 133.322     # Torr to Pa
}

GEOTHERMAL_GRADIENT_CONVERSIONS = {
    'K/m': 1.0,         # K/m (base unit)
    'C/m': 1.0,         # C/m (same as K/m)
    'C/km': 1000,       # C/km
    'F/ft': 0.54864,    # F/ft
    'F/yd': 1.64592,    # F/yd
    'F/100ft': 54.864,  # F/100ft
    'K/km': 1000        # K/km
}

GEOTHERMAL_GRADIENT_CONVERSIONS_INVERSE = {
    'K/m': 1.0,         # K/m (base unit)
    'C/m': 1.0,         # C/m (same as K/m)
    'C/km': 0.001,      # C/km to K/m
    'F/ft': 1.82269,    # F/ft to K/m
    'F/yd': 0.60756,    # F/yd to K/m
    'F/100ft': 0.018227, # F/100ft to K/m
    'K/km': 0.001       # K/km to K/m
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
        
        # Convert to Celsius first if needed
        if from_unit != 'C':
            if from_unit in TEMPERATURE_CONVERSIONS_INVERSE:
                if callable(TEMPERATURE_CONVERSIONS_INVERSE[from_unit]):
                    value = TEMPERATURE_CONVERSIONS_INVERSE[from_unit](value)
                else:
                    value = value * TEMPERATURE_CONVERSIONS_INVERSE[from_unit]
        
        # Convert from Celsius to target unit
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
        
        # Convert to meters first
        if from_unit != 'm':
            value = value * LENGTH_CONVERSIONS_INVERSE[from_unit]
        
        # Convert from meters to target unit
        if to_unit != 'm':
            value = value * LENGTH_CONVERSIONS[to_unit]
        
        return value
    
    def convert_mass_flow(self, value, from_unit, to_unit):
        """Convert mass flow rate between different units"""
        if from_unit == to_unit:
            return value
        
        # Convert to kg/s first
        if from_unit != 'kg/s':
            value = value * MASS_FLOW_CONVERSIONS_INVERSE[from_unit]
        
        # Convert from kg/s to target unit
        if to_unit != 'kg/s':
            value = value * MASS_FLOW_CONVERSIONS[to_unit]
        
        return value
    
    def convert_thermal_conductivity(self, value, from_unit, to_unit):
        """Convert thermal conductivity between different units"""
        if from_unit == to_unit:
            return value
        
        # Convert to W/m-K first
        if from_unit != 'W/m-K':
            value = value * THERMAL_CONDUCTIVITY_CONVERSIONS_INVERSE[from_unit]
        
        # Convert from W/m-K to target unit
        if to_unit != 'W/m-K':
            value = value * THERMAL_CONDUCTIVITY_CONVERSIONS[to_unit]
        
        return value
    
    def convert_heat_capacity(self, value, from_unit, to_unit):
        """Convert heat capacity between different units"""
        if from_unit == to_unit:
            return value
        
        # Convert to J/kg-K first
        if from_unit != 'J/kg-K':
            value = value * HEAT_CAPACITY_CONVERSIONS_INVERSE[from_unit]
        
        # Convert from J/kg-K to target unit
        if to_unit != 'J/kg-K':
            value = value * HEAT_CAPACITY_CONVERSIONS[to_unit]
        
        return value
    
    def convert_density(self, value, from_unit, to_unit):
        """Convert density between different units"""
        if from_unit == to_unit:
            return value
        
        # Convert to kg/m³ first
        if from_unit != 'kg/m3':
            value = value * DENSITY_CONVERSIONS_INVERSE[from_unit]
        
        # Convert from kg/m³ to target unit
        if to_unit != 'kg/m3':
            value = value * DENSITY_CONVERSIONS[to_unit]
        
        return value
    
    def convert_pressure(self, value, from_unit, to_unit):
        """Convert pressure between different units"""
        if from_unit == to_unit:
            return value
        
        # Convert to Pa first
        if from_unit != 'Pa':
            value = value * PRESSURE_CONVERSIONS_INVERSE[from_unit]
        
        # Convert from Pa to target unit
        if to_unit != 'Pa':
            value = value * PRESSURE_CONVERSIONS[to_unit]
        
        return value
    
    def convert_geothermal_gradient(self, value, from_unit, to_unit):
        """Convert geothermal gradient between different units"""
        if from_unit == to_unit:
            return value
        
        # Convert to K/m first
        if from_unit != 'K/m':
            value = value * GEOTHERMAL_GRADIENT_CONVERSIONS_INVERSE[from_unit]
        
        # Convert from K/m to target unit
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
