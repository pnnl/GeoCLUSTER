import sys
import json
from os import path
from pathlib import Path
import pytest
from pprint import pprint

import numpy as np
import numpy.testing as npt
from fastapi.testclient import TestClient

import geocluster_api.main as api
from geocluster_api.types import Payload
from geocluster_api.geocluster_code.econ.econ_wrapper import Econ

client = TestClient(api.app)

payload = Payload()
print("payload:")
print(payload.model_dump_json())
payload_dict = payload.model_dump()
data_path = Path(path.join(__file__, "..", "data"))
clgs_dict = api.clgs_dict

def test_hdf5_endpoint_unauthorized():
    response = client.post("HDF5", json=payload_dict)
    assert response.status_code == 401


def test_subsurface_paraes_are_same():
    assert payload.subsurface_params.mass_flow_rate == 24
    assert payload.subsurface_params.horizontal_extent == 10_000
    assert payload.subsurface_params.drilling_depth == 3500
    assert payload.subsurface_params.geothermal_gradient == 0.03
    assert payload.subsurface_params.borehole_diameter == 0.35
    assert payload.subsurface_params.injection_temp_kelvin == 323.15
    assert payload.subsurface_params.rock_thermal_conductivity == 3


def test_hdf5_outlet_states():
    subsurface_params = api.extract_subsurface_params(payload.subsurface_params)
    target_clgs = api.clgs_dict[payload.tube_shape][payload.fluid]
    subsurface_results = api.get_outlet_states(target_clgs, subsurface_params)

    expected_pout = np.load(data_path / "Pout.npy")
    expected_tout = np.load(data_path / "Tout.npy")
    expected_kWe = np.load(data_path / "kWe.npy")
    expected_kWt = np.load(data_path / "kWt.npy")
    npt.assert_almost_equal(subsurface_results.Pout, expected_pout, decimal=2)
    npt.assert_almost_equal(subsurface_results.Tout, expected_tout, decimal=2)
    npt.assert_almost_equal(subsurface_results.kWe, expected_kWe, decimal=2)
    npt.assert_almost_equal(subsurface_results.kWt, expected_kWt, decimal=2)

def test_hdf5_economic_result():
    with open("./data/tea_results_1.json") as f:
        expected_result_vars = json.load(f)

    subsurface_params = api.extract_subsurface_params(payload.subsurface_params)

    target_clgs = api.clgs_dict[payload.tube_shape][payload.fluid]

    # Get the results
    subsurface_results = api.get_outlet_states(target_clgs, subsurface_params)
    subsurface_results = api.prepare_outlet_states_to_return(subsurface_results)

    ans = {"subsurface_results": subsurface_results}

    TandP = api.get_T_and_P(payload.tube_shape, clgs_dict, subsurface_params)
    clgs_s = api.get_clgs_s(clgs_dict)
    sbt_version = "0"
    end_use = api.get_end_use(payload)
    econ_params = api.get_economic_params(payload)
    tea = Econ(
        TandP,
        clgs_s,
        [payload.tube_shape, end_use, payload.fluid, sbt_version],
        subsurface_params,
        econ_params,
        api.props_paths,
        is_heating=end_use == "heating",
    ).tea
    economic_results = tea.extract_result_vars()

    for key, actual in economic_results.items():
        expected = expected_result_vars[key]
        assert actual == pytest.approx(expected)