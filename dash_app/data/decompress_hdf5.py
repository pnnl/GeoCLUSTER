def main():
    #uncompressing hdf5 files and saving them as chunked files if they dont' exist
    import numpy as np
    import h5py
    filename = "./clgs_results_final_float32.h5"
    shape = (26, 20, 9, 5, 3, 3, 3, 161)
    shapes = ["utube", "coaxial"]
    fluids = ["H2O", "sCO2"]
    types = ["Tout", "Pout"]
    with h5py.File(filename, "r+") as file:
        for tube_shape in shapes:
            for fluid in fluids:
                for type in types:
                    output_location = f"/{tube_shape}/{fluid}/output/{type}_chunked"
                    if output_location in file:
                        print("hdf5 file already decompressed")
                        exit()

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

                    chunk_size = (26, 1, 1, 1, 1, 1, 1, 161)
                    file.create_dataset(output_location, ans.shape, chunks=chunk_size, data=ans)
        file.close()

if __name__ == "__main__":
    main()