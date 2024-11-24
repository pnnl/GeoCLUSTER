GEOCLUSTER uses a precomputed simulation to predict the output heat and exergy of a Closed Loop Geothermal System (CLGS) modeled by the heat-exchanger type, working fluid, and 7 other independent parameters (the ones in the left-hand column). This output is combined with some modifiable economic parameters and is fed into an economic model which computes the levelized cost of heat (LCOF) and the levelized cost of energy (LCOE) of the CLGS in that economic environment.

## Below-Ground Model

The simulation output of every combination of parameters is saved to a compressed HDF5 file, clgs_results_final_float32.h5 which is ~90MB. On the initial load, this HDF5 file is decompressed, read into memory, and loaded into an numpy ndarray. When parameters are selected, it uses those parameters to index into that array and return the resulting simulation output, which is the output heat and exergy of the modeled CLGS over time.

## Above Ground Economic Model

The economic model takes the output of the below-ground model and combines it with 4 modifiable economic parameters (the ones in the right-hand column) to compute the LCOH and the LCOE. This computation happens at runtime, and isn't pre-computed like the below ground model is.

## Computational Architecture

Currently, the webserver and HDF5 files are on the same EC2 instance. The webserver functions as both the view and the controller, with the hdf5 file as the model. GeoCluster can also be ran locally by downloading this repository and following the readme.
