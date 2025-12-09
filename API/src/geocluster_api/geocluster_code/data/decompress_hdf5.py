import h5py
import numpy as np


def copy_hdf5(source_path, destination_path):
    """Copies an HDF5 file from source to destination."""
    with (
        h5py.File(source_path, "r") as source_file,
        h5py.File(destination_path, "w") as destination_file,
    ):
        for name in source_file:
            source_file.copy(name, destination_file)


def decompress__hdf5(source_path, out_path):
    """
    uncompressing hdf5 files and saving them as chunked files
    """

    # copy hdf5
    copy_hdf5(source_path=source_path, destination_path=out_path)

    shape = (26, 20, 9, 5, 3, 3, 3, 161)
    shapes = ["utube", "coaxial"]
    fluids = ["H2O", "sCO2"]
    types = ["Tout", "Pout"]
    with h5py.File(out_path, "r+") as file:
        for tube_shape in shapes:
            for fluid in fluids:
                for output_type in types:
                    output_location = (
                        f"/{tube_shape}/{fluid}/output/{output_type}_chunked"
                    )

                    We = file[f"/{tube_shape}/{fluid}/output/We"][:]
                    U = file[f"/{tube_shape}/{fluid}/output/{output_type}" + "/" + "U"][
                        :
                    ]
                    sigma = file[
                        f"/{tube_shape}/{fluid}/output/{output_type}" + "/" + "sigma"
                    ][:]
                    Vt = file[
                        f"/{tube_shape}/{fluid}/output/{output_type}" + "/" + "Vt"
                    ][:]

                    M_k = np.dot(U, np.dot(np.diag(sigma), Vt))

                    valid_runs = np.argwhere(np.isfinite(We.flatten()))[:, 0]
                    M_k_full = np.full((shape[-1], np.prod(shape[:-1])), np.nan)
                    M_k_full[:, valid_runs] = M_k
                    ans = np.reshape(M_k_full.T, shape)
                    ans = ans.astype(np.float32)

                    chunk_size = (26, 1, 1, 1, 1, 1, 1, 161)
                    file.create_dataset(
                        output_location, ans.shape, chunks=chunk_size, data=ans
                    )
        file.close()
    print("finished decompressing hdf5!")
