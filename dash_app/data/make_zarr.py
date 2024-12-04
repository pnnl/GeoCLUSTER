import h5py
import numpy as np
import os
import zarr
import h5py
import numpy as np

filename = "./clgs_results_final_float32.h5"


def get_parameters():
    with h5py.File(filename, "r") as file:
        case = "coaxial"
        fluid = "sCO2"
        input_loc = "/" + case + "/" + fluid + "/input/"
        mdot = file[input_loc + "mdot"][:]  # i0
        L2 = file[input_loc + "L2"][:]  # i1
        L1 = file[input_loc + "L1"][:]  # i2
        grad = file[input_loc + "grad"][:]  # i3
        D = file[input_loc + "D"][:]  # i4
        Tinj = file[input_loc + "T_i"][:]  # i5
        k = file[input_loc + "k_rock"][:]  # i6
        time = file[input_loc + "time"][:]  # i7
        return [mdot, L2, L1, grad, D, Tinj, k, time]


parameter_values = get_parameters()

filename = "./clgs_results_final_float32.h5"
shape = (26, 20, 9, 5, 3, 3, 3, 161)
shapes = ["utube", "coaxial"]
fluids = ["H2O", "sCO2"]
types = ["Tout", "Pout"]
for tube_shape in shapes:
    for fluid in fluids:
        for type in types:
            with h5py.File(filename, "r") as file:
                We = file[f"/{tube_shape}/{fluid}/output/We"][:]
                U = file[f"/{tube_shape}/{fluid}/output/{type}" + "/" + "U"][:]
                sigma = file[f"/{tube_shape}/{fluid}/output/{type}" + "/" + "sigma"][:]
                Vt = file[f"/{tube_shape}/{fluid}/output/{type}" + "/" + "Vt"][:]

                M_k = np.dot(U, np.dot(np.diag(sigma), Vt))

                valid_runs = np.argwhere(np.isfinite(We.flatten()))[:, 0]
                M_k_full = np.full((shape[-1], np.prod(shape[:-1])), np.nan)
                M_k_full[:, valid_runs] = M_k
                ans = np.reshape(M_k_full.T, shape)
                ans = ans.astype(np.float32)
                print(ans.shape)

                chunk_size = (26, 1, 1, 1, 1, 1, 1, 161)
                z = zarr.creation.array(
                    ans,
                    chunks=chunk_size,
                    store=f"{tube_shape}_{fluid}_{type}.zarr",
                    overwrite=True,
                )
