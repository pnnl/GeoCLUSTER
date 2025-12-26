# GeoCLUSTER API

A FastAPI-based REST API for simulating geothermal energy systems with closed-loop configurations. This API provides subsurface thermal modeling and optional techno-economic analysis (TEA) for geothermal energy production.

## Table of Contents

- [Features](#features)
- [Getting Started](#getting-started)
- [API Endpoints](#api-endpoints)
- [Request Parameters](#request-parameters)
- [Example Usage](#example-usage)
- [Response Format](#response-format)
- [Error Handling](#error-handling)

## Features

- **Subsurface Thermal Modeling**: Simulate temperature and pressure profiles for closed-loop geothermal systems
- **Multiple Configurations**: 
  - Tube shapes: U-tube and Coaxial
  - Working fluids: Water (H2O) and Supercritical CO2 (sCO2)
- **Techno-Economic Analysis**: Optional economic modeling for heating and electricity generation
- **Interactive Documentation**: Auto-generated API docs at `/docs`
- **Health Monitoring**: Built-in health check endpoint

## Getting Started

### Prerequisites

- Python 3.12

### Installation

1. Create and activate a virtual environment (recommended):
```bash
# Create virtual environment
python -m venv .venv

# Activate virtual environment
# On Windows PowerShell:
.\.venv\Scripts\Activate.ps1

# On Windows Command Prompt:
.\.venv\Scripts\activate.bat

# On macOS/Linux:
source .venv/bin/activate
```

2. Install dependencies using `uv` (recommended) or `pip`:
```bash
# Using uv
uv sync

# Or using pip
pip install .

```

3. Run the server:
```bash
# On Windows Command Prompt:
fastapi dev .\src\geocluster_api\main.py
```

```bash
# On macOS/Linux:
fastapi dev src/geocluster_api/main.py
```

The API will be available at `http://localhost:8000`

### Quick Health Check

```bash
curl http://localhost:8000/health
```

## Editing Base Path
To deploy the API on a different domain (i.e. "/some_other_domain/sub_domain"), set the `BASE_PATH` variable 
in the .env file like so:
```
BASE_PATH=/example_path/sub_domain
```

## API Endpoints

### `GET /`
Home endpoint with basic API information and status.

**Response:**
```json
{
  "message": "Welcome to the GeoCLUSTER API! Please view the /docs endpoint for more information",
  "status": "online",
  "endpoints": {
    "docs": "/docs",
    "HDF5 simulation": "/HDF5",
    "SBT v1 simulation": "/SBT/v1",
    "SBT v2 simulation": "/SBT/v2",
  }
}
```

### `GET /health`
Health check endpoint showing API status and available configurations.

**Response:**
```json
{
  "status": "healthy",
  "data_loaded": true,
  "available_configurations": {
    "tube_shapes": ["utube", "coaxial"],
    "fluids": ["H2O", "sCO2"]
  }
}
```

### `POST /HDF5`
Run a geothermal simulation with optional economic analysis.

**Content-Type:** `application/json`

### `POST /SBT/V1`
Run a simulation using the SBT V1 model. Accepts standard simulation parameters plus specific SBT configuration.

**Content-Type:** `application/json`

### `POST /SBT/V2`
Run a simulation using the SBT V2 model. Accepts standard simulation parameters plus specific SBT configuration including pressure and roughness.

**Content-Type:** `application/json`

## Request Parameters

### Subsurface Parameters (Required)

```json
{
  "subsurface_params": {
    "mass_flow_rate": 24.0,           // kg/s - Mass flow rate of working fluid
    "horizontal_extent": 10000.0,     // m - Horizontal length of wellbore
    "drilling_depth": 3500.0,         // m - Vertical depth of well
    "geothermal_gradient": 0.05,      // K/m - Temperature increase per meter depth
    "borehole_diameter": 0.35,        // m - Diameter of borehole
    "injection_temp_kelvin": 303.15,  // K - Injection temperature (30°C)
    "rock_thermal_conductivity": 3.0  // W/m-K - Thermal conductivity of rock
  }
}
```

### Configuration Parameters (Required)

```json
{
  "fluid": "H2O",          // "H2O" or "sCO2" - Working fluid
  "tube_shape": "utube"    // "utube" or "coaxial" - Well configuration
}
```

### Economic Parameters (Optional)

#### For Heating Applications:

```json
{
  "economic_params": {
    "drilling_cost_per_m": 1000.0,    // $/m - Cost per meter of drilling
    "discount_rate": 7.0,             // % - Discount rate for NPV calculations
    "lifetime": 40,                   // years - Project lifetime
    "end_use": "heating",
    "heating_params": {
      "capex": 100.0                  // $/kWt - Capital cost per kW thermal
    }
  }
}
```

#### For Electricity Generation:

```json
{
  "economic_params": {
    "drilling_cost_per_m": 1000.0,
    "discount_rate": 7.0,
    "lifetime": 40,
    "end_use": "electric",
    "electrical_params": {
      "capex": 3000.0,                // $/kWe - Capital cost per kW electric
      "pre_cooling": 13.0,            // °C - Pre-cooling temperature (for sCO2)
      "outlet_pressure": 80.0         // bar - Turbine outlet pressure (for sCO2)
    }
  }
}
```

### SBT Parameters (Required for /SBTV1 and /SBTV2)

#### For SBTV1:

```json
{
  "sbt_params": {
    "surface_temperature_c": 25.0,        // Surface temperature (°C)
    "rock_specific_heat_capacity": 790.0, // J/kg-K
    "rock_density": 2750.0,               // kg/m^3
    "well_bore_radius_vertical": 0.35,    // m
    "well_bore_radius_horizontal": 0.175, // m
    "mesh_fineness": 0,                   // Mesh refinement level
    "accuracy": 1                         // Accuracy level
  }
}
```

#### For SBTV2 (Adds to SBTV1):

```json
{
  "sbt_params": {
    // ... All SBTV1 parameters plus:
    "inlet_pressure": 10,                 // MPa - Inlet pressure
    "pipe_roughess": 0.0000001            // m - Pipe roughness
  }
}
```

## Example Usage

### Example 1: Basic Subsurface Simulation (No Economics)

**Python:**
```python
import requests
import json

url = "http://localhost:8000/HDF5"

payload = {
    "subsurface_params": {
        "mass_flow_rate": 24.0,
        "horizontal_extent": 10000.0,
        "drilling_depth": 4000.0,
        "geothermal_gradient": 0.05,
        "borehole_diameter": 0.25,
        "injection_temp_kelvin": 323.15,
        "rock_thermal_conductivity": 3.0
    },
    "fluid": "sCO2",
    "tube_shape": "utube"
}

response = requests.post(url, json=payload)
results = response.json()

print(f"Status: {response.status_code}")
print(json.dumps(results, indent=2))
```

**cURL:**
```bash
curl -X POST "http://localhost:8000/HDF5" \
  -H "Content-Type: application/json" \
  -d '{
    "subsurface_params": {
      "mass_flow_rate": 25.0,
      "horizontal_extent": 2000.0,
      "drilling_depth": 4000.0,
      "geothermal_gradient": 0.05,
      "borehole_diameter": 0.25,
      "injection_temp_kelvin": 323.15,
      "rock_thermal_conductivity": 3.0
    },
    "fluid": "sCO2",
    "tube_shape": "utube"
  }'
```

### Example 2: Simulation with Heating Economics

**Python:**
```python
import requests

url = "http://localhost:8000/HDF5"

payload = {
    "subsurface_params": {
        "mass_flow_rate": 30.0,
        "horizontal_extent": 5000.0,
        "drilling_depth": 3000.0,
        "geothermal_gradient": 0.06,
        "borehole_diameter": 0.3,
        "injection_temp_kelvin": 303.15,
        "rock_thermal_conductivity": 2.5
    },
    "fluid": "H2O",
    "tube_shape": "coaxial",
    "economic_params": {
        "drilling_cost_per_m": 1200.0,
        "discount_rate": 5.0,
        "lifetime": 30,
        "end_use": "heating",
        "heating_params": {
            "capex": 150.0
        }
    }
}

response = requests.post(url, json=payload)
results = response.json()

# Access economic results
if "economic_results" in results["data"]:
    econ = results["data"]["economic_results"]
    print(f"LCOH: {econ['LCOH']}")
    print(f"Total CAPEX: {econ['Total CAPEX']}")
    print(f"Average Heat Production: {econ['average_heat_production']}")
```

### Example 3: Electricity Generation with sCO2

**Python:**
```python
import requests

url = "http://localhost:8000/HDF5"

payload = {
    "subsurface_params": {
        "mass_flow_rate": 50.0,
        "horizontal_extent": 8000.0,
        "drilling_depth": 5000.0,
        "geothermal_gradient": 0.07,
        "borehole_diameter": 0.4,
        "injection_temp_kelvin": 323.15,
        "rock_thermal_conductivity": 3.5
    },
    "fluid": "sCO2",
    "tube_shape": "utube",
    "economic_params": {
        "drilling_cost_per_m": 1500.0,
        "discount_rate": 8.0,
        "lifetime": 40,
        "end_use": "electric",
        "electrical_params": {
            "capex": 3500.0,
            "pre_cooling": 15.0,
            "outlet_pressure": 85.0
        }
    }
}

response = requests.post(url, json=payload)
results = response.json()

# Access results
if response.status_code == 200:
    econ = results["data"]["economic_results"]
    print(f"LCOE: {econ['LCOE']}")
    print(f"Average Net Electricity: {econ['average_net_electricity_production']}")
else:
    print(f"Error: {response.status_code} - {results['detail']}")
```

### Example 4: JavaScript/TypeScript (Node.js)

```javascript
const axios = require('axios');

const payload = {
  subsurface_params: {
    mass_flow_rate: 25.0,
    horizontal_extent: 2000.0,
    drilling_depth: 4000.0,
    geothermal_gradient: 0.05,
    borehole_diameter: 0.35,
    injection_temp_kelvin: 323.15,
    rock_thermal_conductivity: 3.0
  },
  fluid: "sCO2",
  tube_shape: "utube"
};

axios.post('http://localhost:8000/HDF5', payload)
  .then(response => {
    console.log('Simulation Results:', response.data);
    const subsurface = response.data.data.subsurface_results;
    console.log('Output Temperatures:', subsurface.Tout);
    console.log('Electric Power:', subsurface.kWe);
  })
  .catch(error => {
    console.error('Error:', error.response?.data || error.message);
  });
```

### Example 5: SBT V1 Simulation

**Python:**
```python
import requests
import json

url = "http://localhost:8000/SBT/V1"

payload = {
    "subsurface_params": {
        "mass_flow_rate": 24.0,
        "horizontal_extent": 10000.0,
        "drilling_depth": 3500.0,
        "geothermal_gradient": 0.05,
        "borehole_diameter": 0.35,
        "injection_temp_kelvin": 303.15,
        "rock_thermal_conductivity": 3.0
    },
    "fluid": "H2O",
    "tube_shape": "utube",
    "sbt_params": {
        "surface_temperature_c": 25.0,
        "rock_specific_heat_capacity": 790.0,
        "rock_density": 2750.0,
        "well_bore_radius_vertical": 0.35,
        "well_bore_radius_horizontal": 0.175,
        "mesh_fineness": 0,
        "accuracy": 1
    }
}

response = requests.post(url, json=payload)
results = response.json()

print(f"Status: {response.status_code}")
if response.status_code == 200:
    print(json.dumps(results, indent=2))
else:
    print(f"Error: {results}")
```

## Response Format

### Successful Response (200 OK)

```json
{
  "data": {
    "subsurface_results": {
      "Tout": [/* Array of outlet temperatures (K) over time */],
      "Pout": [/* Array of outlet pressures (Pa) over time */],
      "times": [/* Array of time points (years) */],
      "kWe": [/* Array of electric power output (kW) over time */],
      "kWt": [/* Array of thermal power output (kW) over time */]
    },
    "economic_results": {  // Only present if economic_params provided
      "average_production_temporate": "150.5 C",
      "average_production_pressure": "45.2 bar",
      "average_heat_production": "1250.3 kWth",
      "first_year_heat_production": "10500.2 MWh",
      "Total CAPEX": "12.5 M$",
      "Drilling Cost": "5.2 M$",
      "Surface Plant Cost": "7.3 M$",
      "OPEX": "150.0 k$/year",
      "LCOH": "45.5 $/MWh"  // For heating
      // OR
      "average_net_electricity_production": "450.2 kWe",
      "first_year_electricity_production": "3500.8 MWh",
      "LCOE": "85.2 $/MWh"  // For electricity
    }
  }
}
```

## Error Handling

The API uses standard HTTP status codes and returns detailed error messages:

### 400 Bad Request
Invalid input parameters or validation errors.

```json
{
  "detail": "Invalid tube_shape: invalid_shape. Must be 'utube' or 'coaxial'."
}
```

### 500 Internal Server Error
Unexpected server errors during simulation.

```json
{
  "detail": "An error occurred during simulation: [error description]"
}
```

### 503 Service Unavailable
API is running but data files are not loaded.

```json
{
  "detail": "Server data not properly initialized. Please contact the administrator."
}
```

## Parameter Validation

The API validates:
- ✅ Required fields are present
- ✅ Tube shape is either "utube" or "coaxial"
- ✅ Fluid is either "H2O" or "sCO2"
- ✅ Numeric values are valid (not NaN or Inf)
- ✅ Parameters are within physically reasonable bounds

### Typical Valid Ranges:
- **mass_flow_rate**: 5 - 100 kg/s
- **horizontal_extent**: 1,000 - 20,000 m
- **drilling_depth**: 1,000 - 5,000 m
- **geothermal_gradient**: 0.03 - 0.7 K/m
- **borehole_diameter**: 0.2159 - 0.4445 m
- **injection_temp_kelvin**: 303.15 - 333.15 K (30°C - 60°C)
- **rock_thermal_conductivity**: 1.5 - 4.5 W/m-K

## Interactive Documentation

FastAPI provides automatic interactive API documentation:

- **Swagger UI**: [http://localhost:8000/docs](http://localhost:8000/docs)
- **ReDoc**: [http://localhost:8000/redoc](http://localhost:8000/redoc)

These interfaces allow you to:
- View all endpoints and their parameters
- Test API calls directly from the browser
- See request/response schemas
- Download OpenAPI specification

## Logging

The API logs important events including:
- Successful data file loading
- Simulation completions
- Errors and exceptions

Logs are output to the console with INFO level by default.

## Development

### Running in Development Mode

```bash
uvicorn geocluster_api.main:app --reload --host 0.0.0.0 --port 8000
```

The `--reload` flag enables auto-reloading when code changes are detected.

## Acknowledgments
This API is built using:
- [FastAPI](https://fastapi.tiangolo.com/) - Modern web framework
- [NumPy](https://numpy.org/) - Numerical computing
- [SciPy](https://scipy.org/) - Scientific computing
- [CoolProp](http://www.coolprop.org/) - Thermophysical property calculations
- [h5py](https://www.h5py.org/) - HDF5 data storage
