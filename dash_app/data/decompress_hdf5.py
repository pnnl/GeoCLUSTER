import h5py
from paths import inpath_dict
import numpy as np
import os

def copy_hdf5(source_path, destination_path):
    """Copies an HDF5 file from source to destination."""
    with h5py.File(source_path, 'r') as source_file, h5py.File(destination_path, 'w') as destination_file:
        for name in source_file:
            source_file.copy(name, destination_file)

def decompress__hdf5():
    #uncompressing hdf5 files and saving them as chunked files if they dont' exist
    print("decompressing hdf5...")

    in_filename = inpath_dict["h5_filepath"]
    out_filename = inpath_dict["decompressed_h5_filepath"]
    
    if not os.path.exists(in_filename):
        alt_paths = [
            "dash_app/data/clgs_results_final_float32.h5",
            "data/clgs_results_final_float32.h5",
            os.path.join(os.path.dirname(os.path.dirname(__file__)), "data", "clgs_results_final_float32.h5")
        ]
        for alt_path in alt_paths:
            if os.path.exists(alt_path):
                in_filename = alt_path
                print(f"Found HDF5 file at: {in_filename}")
                break
        else:
            raise FileNotFoundError(f"Could not find HDF5 file. Tried: {inpath_dict['h5_filepath']} and {alt_paths}")
    
    copy_hdf5(source_path=in_filename, destination_path=out_filename)

    shape = (26, 20, 9, 5, 3, 3, 3, 161)
    shapes = ["utube", "coaxial"]
    fluids = ["H2O", "sCO2"]
    types = ["Tout", "Pout"]
    with h5py.File(out_filename, "r+") as file:
        for tube_shape in shapes:
            for fluid in fluids:
                for output_type in types:
                    output_location = f"/{tube_shape}/{fluid}/output/{output_type}_chunked"

                    We = file[f"/{tube_shape}/{fluid}/output/We"][:]
                    U = file[f"/{tube_shape}/{fluid}/output/{output_type}" + "/" + "U"][:]
                    sigma = file[f"/{tube_shape}/{fluid}/output/{output_type}" + "/" + "sigma"][:]
                    Vt = file[f"/{tube_shape}/{fluid}/output/{output_type}" + "/" + "Vt"][:]

                    M_k = np.dot(U, np.dot(np.diag(sigma), Vt))

                    valid_runs = np.argwhere(np.isfinite(We.flatten()))[:, 0]
                    M_k_full = np.full((shape[-1], np.prod(shape[:-1])), np.nan)
                    M_k_full[:, valid_runs] = M_k
                    ans = np.reshape(M_k_full.T, shape)
                    ans = ans.astype(np.float32)

                    chunk_size = (26, 1, 1, 1, 1, 1, 1, 161)
                    file.create_dataset(output_location, ans.shape, chunks=chunk_size, data=ans)
        file.close()
    print("finished decompressing hdf5!")

if __name__ == "__main__":
    decompress__hdf5()