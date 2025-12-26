import logging
from os import getenv
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Literal, Union
from dotenv import load_dotenv

import CoolProp.CoolProp as CP
import numpy as np
from fastapi import Depends, FastAPI, HTTPException, Request, status, APIRouter
from fastapi.responses import JSONResponse
from fastapi.security import APIKeyHeader
from numpy.typing import NDArray
from pydantic import BaseModel, Field

from .geocluster_code.clgs import Data, interp_kW_standalone
from .geocluster_code.data.decompress_hdf5 import decompress__hdf5
from .geocluster_code.econ.econ_wrapper import Econ
from .geocluster_code.sbt.sbt_v27 import run_sbt
from .types import Payload, Subsurface

# Configure logging
load_dotenv()
default_path = ""
base_path = getenv("BASE_PATH", default_path)
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

app = FastAPI(
    title="GeoCLUSTER API", description="API for geothermal energy simulations",
    docs_url=f"{base_path}/docs",
    redoc_url=f"{base_path}/redoc",
)
prefix_router = APIRouter(prefix=base_path)

parent_path = Path(__file__).resolve().parent
# region declaring global constants
data_path = (
    parent_path / "geocluster_code/data/decompressed_clgs_results_final_float32.h5"
)
h2o_props_path = parent_path / "geocluster_code/data/properties_H2O.mat"
co2_props_path = parent_path / "geocluster_code/data/properties_CO2v2.mat"
co2_additional_props_path = (
    parent_path / "geocluster_code/data/additional_properties_CO2v2.mat"
)
props_paths = [h2o_props_path, co2_props_path, co2_additional_props_path]

if not Path(data_path).exists():
    logger.info("decompressing source HDF5 data")
    compressed_file_name = "clgs_results_final_float32.h5"
    compressed_hdf5_path = data_path.parent / compressed_file_name
    decompress__hdf5(compressed_hdf5_path, data_path)
logger.info("source HDF5 exists, continuing...")


# Validate that data files exist
try:
    if not data_path.exists():
        raise FileNotFoundError(f"Data file not found: {data_path}")
    for path in props_paths:
        if not path.exists():
            raise FileNotFoundError(f"Properties file not found: {path}")

    clgs_dict = {
        "utube": {
            "H2O": Data(data_path, "utube", "H2O"),
            "sCO2": Data(data_path, "utube", "sCO2"),
        },
        "coaxial": {
            "H2O": Data(data_path, "coaxial", "H2O"),
            "sCO2": Data(data_path, "coaxial", "sCO2"),
        },
    }
    logger.info("Successfully loaded all data files")
except Exception as e:
    logger.error(f"Failed to initialize data: {e}")
    clgs_dict = None
# endregion

SECRET_HEADER = "HTTP_X_API_UMBRELLA_REQUEST_ID"
api_key_header = APIKeyHeader(name=SECRET_HEADER)

@prefix_router.post("/HDF5")
def run_hdf5_simulation(body: Payload, request: Request):
    """
    Run a geothermal simulation with optional economic analysis.

    Raises:
        HTTPException: If data is not loaded or simulation fails
    """

    try:
        do_economic = body.economic_params is not None

        subsurface_params = extract_subsurface_params(body.subsurface_params)

        target_clgs = clgs_dict[body.tube_shape][body.fluid]

        # Get the results
        subsurface_results = get_outlet_states(target_clgs, subsurface_params)
        subsurface_results = prepare_outlet_states_to_return(subsurface_results)

        ans = {"subsurface_results": subsurface_results}

        if do_economic:
            TandP = get_T_and_P(body.tube_shape, clgs_dict, subsurface_params)
            clgs_s = get_clgs_s(clgs_dict)
            sbt_version = "0"
            end_use = get_end_use(body)
            econ_params = get_economic_params(body)
            tea = Econ(
                TandP,
                clgs_s,
                [body.tube_shape, end_use, body.fluid, sbt_version],
                subsurface_params,
                econ_params,
                props_paths,
                is_heating=end_use == "heating",
            ).tea
            ans["economic_results"] = prepare_economic_results(tea, end_use)

        logger.info(
            f"Successfully completed simulation for {body.tube_shape}/{body.fluid}"
        )
        return {"data": ans}

    except HTTPException:
        # Re-raise HTTP exceptions
        raise
    except ValueError as e:
        logger.error(f"Validation error in simulation: {e}")
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail=f"Invalid input parameters: {str(e)}",
        )
    except Exception as e:
        logger.error(f"Unexpected error in simulation: {e}", exc_info=True)
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=f"An error occurred during simulation: {str(e)}",
        )


def extract_subsurface_params(params: Subsurface):
    """
    Extract and convert subsurface parameters to numpy float32 array.

    Args:
        params: Subsurface parameters object

    Returns:
        List of numpy float32 values

    Raises:
        ValueError: If any parameter is invalid or out of expected range
    """
    try:
        # relies on tuple being ordered in the same order as SubSurface class
        params = tuple(params.model_dump().values())
        params_as_np_float = [np.float32(val) for val in params]

        # Basic validation
        if any(np.isnan(p) or np.isinf(p) for p in params_as_np_float):
            raise ValueError("Parameters contain NaN or Inf values")

        return params_as_np_float
    except Exception as e:
        logger.error(f"Error extracting subsurface params: {e}")
        raise ValueError(f"Failed to process subsurface parameters: {str(e)}")


def get_end_use(body: Payload):
    return body.economic_params.end_use


def get_economic_params(body: Payload):
    params = body.economic_params
    # the code that this gets passed into expects us to return None for unused params
    return [
        params.drilling_cost_per_m,
        params.discount_rate,
        params.lifetime,
        # TODO: clean this up lol
        getattr(getattr(params, "heating_params", None), "capex", None),
        getattr(getattr(params, "electrical_params", None), "capex", None),
        getattr(getattr(params, "electrical_params", None), "pre_cooling", None),
        getattr(getattr(params, "electrical_params", None), "outlet_pressure", None),
    ]


def prepare_economic_results(tea, end_use):
    ans = {
        "average_production_temporate": "{0:.1f}".format(tea.AveProductionTemperature)
        + " C",
        "average_production_pressure": "{0:.1f}".format(tea.AveProductionPressure)
        + " bar",
    }
    # TODO: better error handling
    if np.any(np.isin(tea.error_codes, [1000])):
        ans["error"] = (
            "production temperature drops below injection temperature. Simulation terminated"
        )
        return ans

    if end_use == "heating":
        ans["average_heat_production"] = (
            "{0:.1f}".format((tea.AveInstHeatProduction)) + " kWth"
        )
        ans["first_year_heat_production"] = (
            "{0:.1f}".format((tea.FirstYearHeatProduction / 1e3)) + " MWh"
        )
        ans["OPEX"] = "{0:.1f}".format((tea.AverageOPEX_Plant * 1000)) + " k$/year"
        ans["LCOH"] = "{0:.1f}".format((tea.LCOH)) + " $/MWh"
    # it's electricity
    else:
        ans["average_heat_production"] = (
            "{0:.1f}".format((tea.AveInstHeatProduction)) + " kWth"
        )
        ans["average_net_electricity_production"] = (
            "{0:.1f}".format((tea.AveInstNetElectricityProduction)) + " kWe"
        )
        ans["first_year_heat_production"] = (
            "{0:.1f}".format((tea.FirstYearHeatProduction / 1e3)) + " MWh"
        )
        ans["first_year_electricity_production"] = (
            "{0:.1f}".format((tea.FirstYearElectricityProduction / 1e3)) + " MWh"
        )
        ans["OPEX"] = "{0:.1f}".format((tea.AverageOPEX_Plant * 1000)) + " k$/year"
        ans["LCOE"] = "{0:.1f}".format((tea.LCOE * 1000)) + " $/MWh"

    ans["Total CAPEX"] = "{0:.1f}".format(tea.TotalCAPEX) + " M$"
    ans["Drilling Cost"] = "{0:.1f}".format((tea.CAPEX_Drilling)) + " M$"
    ans["Surface Plant Cost"] = "{0:.1f}".format((tea.CAPEX_Surface_Plant)) + " M$"
    return ans


@dataclass
class OutletStates:
    Tout: NDArray
    Pout: NDArray
    times: NDArray
    kWe: NDArray
    kWt: NDArray


def get_outlet_states(clgs: Data, point, Tamb=300) -> OutletStates:
    """
    Calculate outlet states from CLGS data.

    Args:
        clgs: Data object containing simulation data
        point: Input parameters for interpolation
        Tamb: Ambient temperature (default: 300K)

    Returns:
        OutletStates object with temperature, pressure, time, and power data

    Raises:
        ValueError: If interpolation fails or point is out of bounds
    """
    try:
        Tout, Pout, time = clgs.interp_outlet_states(point)

        if Tout is None or Pout is None or time is None:
            raise ValueError("Interpolation returned None values")

        kWe, kWt = clgs.interp_kW(point, Tout, Pout, Tamb)

        return OutletStates(Tout, Pout, time, kWe, kWt)
    except ValueError as e:
        if "out of bounds" in str(e):
            raise ValueError(f"Input parameters are outside valid range: {str(e)}")
        raise
    except Exception as e:
        logger.error(f"Error calculating outlet states: {e}")
        raise ValueError(f"Failed to calculate outlet states: {str(e)}")


def prepare_outlet_states_to_return(outlet: OutletStates):
    """
    Convert OutletStates numpy arrays to JSON-serializable lists.

    Args:
        outlet: OutletStates object with numpy arrays

    Returns:
        Dictionary with lists instead of numpy arrays

    Raises:
        ValueError: If conversion fails
    """
    try:
        # we can't jsonify np arrays so we gotta make them list
        ans = asdict(outlet)
        ans = {key: value.tolist() for key, value in ans.items()}
        return ans
    except Exception as e:
        logger.error(f"Error preparing outlet states: {e}")
        raise ValueError(f"Failed to prepare output data: {str(e)}")


def get_m_dot(point):
    return point[0]


def get_T_and_P(tube_shape: str, clgs_dict, point, Tamb=300):
    clgs_s_of_tube = clgs_dict[tube_shape]
    h2o = get_outlet_states(clgs_s_of_tube["H2O"], point, Tamb)
    sco2 = get_outlet_states(clgs_s_of_tube["sCO2"], point, Tamb)
    mdot = get_m_dot(point)
    TandP_dict = {
        "sCO2_Tout": sco2.Tout,
        "sCO2_Pout": sco2.Pout,
        "H2O_Tout": h2o.Tout,
        "H2O_Pout": h2o.Pout,
        "sCO2_kWe": sco2.kWe,
        "sCO2_kWt": sco2.kWt,
        "H2O_kWe": h2o.kWe,
        "H2O_kWt": h2o.kWt,
        "time": sco2.times,
        "mdot": mdot,
    }
    return TandP_dict


def get_clgs_s(clgs_dict):
    return [
        clgs_dict["utube"]["H2O"],
        clgs_dict["utube"]["sCO2"],
        clgs_dict["coaxial"]["H2O"],
        clgs_dict["coaxial"]["sCO2"],
    ]


@prefix_router.get("/")
def home():
    """Health check endpoint."""
    return {
        "message": "Welcome to the GeoCLUSTER API! Please view the /docs endpoint for more information",
        "status": "online" if clgs_dict is not None else "data_not_loaded",
        "endpoints": {
            "docs": "/docs",
            "hdf5 simulation": "/HDF5",
            "sbt v1 simulation": "/SBT/v1",
            "sbt v2 simulation": "/SBT/v2",
        },
    }


@prefix_router.get("/health")
def health_check():
    """
    Check if the API is running and data is properly loaded.

    Returns:
        Health status with data availability information
    """
    data_loaded = clgs_dict is not None

    if data_loaded:
        return {
            "status": "healthy",
            "data_loaded": True,
            "available_configurations": {
                "tube_shapes": list(clgs_dict.keys()),
                "fluids": ["H2O", "sCO2"],
            },
        }
    else:
        return JSONResponse(
            status_code=status.HTTP_503_SERVICE_UNAVAILABLE,
            content={
                "status": "unhealthy",
                "data_loaded": False,
                "error": "Data files not properly initialized",
            },
        )


class SBTV1Params(BaseModel):
    surface_temperature_c: float = 25.0
    rock_specific_heat_capacity: float = 790.0
    rock_density: float = 2750.0
    well_bore_radius_vertical: float = 0.35
    well_bore_radius_horizontal: float = 0.175
    mesh_fineness: int = 0
    accuracy: int = 1


class SBTV2Params(SBTV1Params):
    # MPA
    inlet_pressure: int = 10
    # m
    pipe_roughess: float = 0.0000001


class SBTV1Payload(Payload):
    sbt_params: SBTV1Params = SBTV1Params()


class SBTV2Payload(Payload):
    sbt_params: SBTV2Params = SBTV2Params()


class SBTPayload(Payload):
    sbt_params: Union[SBTV1Params, SBTV2Params]


@prefix_router.post("/SBT/V1")
def run_sbt_v1_endpoint(body: SBTV1Payload, request: Request):
    """
    Run simulation using SBT V1 model.
    """
    return sbt_wrapper(body, 1)


@prefix_router.post("/SBT/V2")
def run_sbt_v2_endpoint(body: SBTV2Payload, request: Request):
    """
    Run simulation using SBT V2 model.
    """
    return sbt_wrapper(body, 2)


def calculate_power_sbt(mdot, Tinj, Pinj, Tout, Pout, fluid_name, Tamb=300.0):
    """
    Calculate power output for SBT simulation.

    Args:
        mdot: Mass flow rate [kg/s]
        Tinj: Injection temperature [K]
        Pinj: Injection pressure [Pa]
        Tout: Outlet temperature array [K]
        Pout: Outlet pressure array [Pa]
        fluid_name: 'H2O' or 'CO2'
        Tamb: Ambient temperature [K]

    Returns:
        Tuple of (kWe, kWt) arrays
    """
    try:
        enthalpy_out = CP.PropsSI("H", "T", Tout, "P", Pout, fluid_name)
        enthalpy_in = CP.PropsSI("H", "T", Tinj, "P", Pinj, fluid_name)
        entropy_out = CP.PropsSI("S", "T", Tout, "P", Pout, fluid_name)
        entropy_in = CP.PropsSI("S", "T", Tinj, "P", Pinj, fluid_name)

        kWe = (
            mdot
            * (enthalpy_out - enthalpy_in - Tamb * (entropy_out - entropy_in))
            / 1000.0
        )
        kWt = mdot * (enthalpy_out - enthalpy_in) / 1000.0

        return kWe, kWt
    except Exception as e:
        logger.error(f"Error calculating power in SBT: {e}")
        # Return zeros in case of failure
        return np.zeros_like(Tout), np.zeros_like(Tout)


def run_sbt_simulation(body: SBTPayload, sbt_version=Literal[1, 2]):
    """
    Run SBT simulation and return outlet states.
    """
    sub = body.subsurface_params
    sbt_params = body.sbt_params

    # Fluid mapping
    fluid_map = {"H2O": 1, "sCO2": 2}
    fluid_id = fluid_map.get(body.fluid, 1)
    cp_fluid = "H2O" if body.fluid == "H2O" else "CO2"

    # Configuration mapping
    config_map = {"coaxial": 1, "utube": 2}
    clg_config = config_map.get(body.tube_shape, 2)

    Pinj = 10 if "inlet_pressure" not in sbt_params else sbt_params.inlet_pressure
    pipe_roughness = (
        1e-6 if "pipe_roughness" not in sbt_params else sbt_params.pipe_roughness
    )

    # Run SBT
    # Note: Converting units to match SBT expectations
    # L1, L2: m -> km (because is_app=True in sbt_v27.py multiplies by 1000)
    # Tinj: K -> C
    try:
        times, Tout, Pout = run_sbt(
            sbt_version=1,
            clg_configuration=clg_config,
            mdot=sub.mass_flow_rate,
            Tinj=sub.injection_temp_kelvin - 273.15,
            fluid=fluid_id,
            DrillingDepth_L1=sub.drilling_depth / 1000.0,
            HorizontalExtent_L2=sub.horizontal_extent / 1000.0,
            Diameter1=sbt_params.well_bore_radius_vertical,
            Diameter2=sbt_params.well_bore_radius_horizontal,
            Tsurf=sbt_params.surface_temperature_c,
            GeoGradient=sub.geothermal_gradient,
            k_m=sub.rock_thermal_conductivity,
            c_m=sbt_params.rock_specific_heat_capacity,
            rho_m=sbt_params.rock_density,
            mesh_fineness=sbt_params.mesh_fineness,
            accuracy=sbt_params.accuracy,
            HYPERPARAM1=Pinj,
            HYPERPARAM2=pipe_roughness,
        )
    except Exception as e:
        logger.error(f"SBT execution failed: {e}")
        raise HTTPException(
            status_code=500, detail=f"SBT model execution failed: {str(e)}"
        )

    # In SBT V1, Pout is always None, only SBT V2 Computes POUT
    if Pout is None:
        constant_pressure = 2e7  # 200 Bar in pascal || 2.09e7
        # constant_pressure = 22228604.37405011
        Pout = constant_pressure * np.ones_like(Tout)

    times = times[14:]
    Tout = Tout[14:]
    Pout = Pout[14:]

    if times is None or Tout is None:
        raise HTTPException(status_code=500, detail="SBT model returned no data")

    # # Temperature: C -> K (for CoolProp and consistency)
    # Tout_K = Tout + 273.15

    # Calculate Power
    kWe, kWt = interp_kW_standalone(
        sub.mass_flow_rate, sub.injection_temp_kelvin, Tout, Pout, Pinj, cp_fluid
    )

    return OutletStates(Tout, Pout, times, kWe, kWt), Pinj


def sbt_wrapper(body, sbt_version: Literal[1, 2]):
    try:
        outlet_states, Pinj = run_sbt_simulation(body, sbt_version)

        subsurface_results = prepare_outlet_states_to_return(outlet_states)
        ans = {"subsurface_results": subsurface_results}

        if body.economic_params is not None:
            # Prepare TandP dict for Econ
            # Zero-fill inactive fluid
            zeros = np.zeros_like(outlet_states.times)

            if body.fluid == "H2O":
                TandP = {
                    "H2O_Tout": outlet_states.Tout,
                    "H2O_Pout": outlet_states.Pout,
                    "H2O_kWe": outlet_states.kWe,
                    "H2O_kWt": outlet_states.kWt,
                    "sCO2_Tout": zeros,
                    "sCO2_Pout": zeros,
                    "sCO2_kWe": zeros,
                    "sCO2_kWt": zeros,
                }
            else:
                TandP = {
                    "sCO2_Tout": outlet_states.Tout,
                    "sCO2_Pout": outlet_states.Pout,
                    "sCO2_kWe": outlet_states.kWe,
                    "sCO2_kWt": outlet_states.kWt,
                    "H2O_Tout": zeros,
                    "H2O_Pout": zeros,
                    "H2O_kWe": zeros,
                    "H2O_kWt": zeros,
                }

            TandP["time"] = outlet_states.times
            TandP["mdot"] = body.subsurface_params.mass_flow_rate

            clgs_s = get_clgs_s(clgs_dict) if clgs_dict else []

            sbt_version = "2"
            end_use = get_end_use(body)
            econ_params = get_economic_params(body)

            tea = Econ(
                TandP,
                clgs_s,
                [body.tube_shape, end_use, body.fluid, sbt_version],
                extract_subsurface_params(body.subsurface_params),
                econ_params,
                props_paths,
                is_heating=end_use == "heating",
            ).tea

            ans["economic_results"] = prepare_economic_results(tea, end_use)

        return {"data": ans}

    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Error in SBT simulation: {e}", exc_info=True)
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=f"An error occurred during SBT simulation: {str(e)}",
        )
app.include_router(prefix_router)
