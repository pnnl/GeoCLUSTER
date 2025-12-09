from typing import Literal

from numpy.typing import NDArray
from pydantic import BaseModel, ConfigDict


class Subsurface(BaseModel):
    # kg/s
    mass_flow_rate: float = 24
    # m
    horizontal_extent: float = 10_000
    # m
    drilling_depth: float = 3500
    # K/m
    geothermal_gradient: float = 0.03
    # m
    borehole_diameter: float = 0.35

    injection_temp_kelvin: float = 323.15
    # W/m-k
    rock_thermal_conductivity: float = 3


class Heating(BaseModel):
    # $/kWt
    capex: float = 100


class Electrical(BaseModel):
    # $/kWe
    capex: float = 3_000


class SCO2_Turbine(BaseModel):
    # $ / kWe
    capex: float = 3_000
    # C
    pre_cooling: float = 13
    # bar
    outlet_pressure: float = 80


class Economics(BaseModel):
    drilling_cost_per_m: float = 1_000
    discount_rate: float = 7
    # years
    lifetime: int = 40
    heating_params: Heating | None = Heating()
    electrical_params: Electrical | SCO2_Turbine | None = None
    end_use: Literal["heating", "electric"] = "heating"


class Payload(BaseModel):
    subsurface_params: Subsurface = Subsurface()
    fluid: Literal["H2O", "sCO2"] = "H2O"
    tube_shape: Literal["utube", "coaxial"] = "utube"
    economic_params: Economics | None = Economics()


class SubSurfaceResponse(BaseModel):
    Tout: NDArray
    Pout: NDArray
    kWe: NDArray
    kWt: NDArray
    times: NDArray

    model_config = ConfigDict(arbitrary_types_allowed=True)


class Response(BaseModel):
    sub_surface: SubSurfaceResponse
