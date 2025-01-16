class GeoCLUSTER {
    constructor() {
      this.initPromise = this.init();
    }
  
    init = async () => {
      this.pyodide = await loadPyodide();
  
      await this.pyodide.loadPackage("h5py");
      await this.pyodide.loadPackage("zarr");
  
      // Use Pyodide's API to define and run Python code
    //   await this.pyodide.runPythonAsync(`
    //         from pyodide.http import pyfetch
    //         import h5py
    //         import zarr
    //         import io
        
    //         # Fetch the HDF5 file using pyodide's pyfetch
    //         async def load_zarr_file():
    //             # response = await pyfetch("/assets/clgs_results_final_float32.h5", method="GET", mode="cors")
    //             response = await pyfetch("/zarr/utube_sCO2_Pout.zarr", method="GET", mode="cors")
    //             zarr_results_bytes = await response.bytes()
    //             print("GOT ZARR_RESULTS_BYTES")
        
    //             # Create a byte buffer from the downloaded bytes
    //             byte_buffer = io.BytesIO(zarr_results_bytes)
    //             print("CREATED BYTE BUFFER")
        
    //             # Open the HDF5 file from the byte buffer
    //             # return h5py.File(byte_buffer, 'r')
    //             store = zarr.storage.KVStore(byte_buffer)
    //             print(f'KVStore KEYS: {list(store.keys())}')
    //             return zarr.open(store, mode='r')

    //         # Call the async function in your browser environment
    //         clgs_load_zarr = await load_zarr_file()
    //         print("ZARR FILES LOADED")
    //     `);
    // };
      console.log("IM DOING NOTHING IN GEOCLUSTER.JS FOR NOW");
      // console.log("ABOUT TO RUN load_zarr_directory()");
      // await this.load_zarr_directory("utube_sCO2_Pout.zarr");
};

    load_zarr_directory = async (dirName) => {
      // await this.initPromise;
      const url = `/zarr_list/${dirName}`;
      const response = await fetch(url);
      const file_list = await response.json();

      const store = new Map();

      console.log("ABOUT TO LOOP AND FETCH ALL INDIVIDUAL ZARR FILES");
      for (const fileName of file_list) {
        const fileResponse = await fetch(`/zarr/${fileName}`);
        const fileBuffer = await fileResponse.arrayBuffer();
        const relativePath = fileName.replace(`${dirName}/`, "");
        store.set(relativePath, new Uint8Array(fileBuffer));
      }

      console.log("FINISHED MAKING STORE");

      // Use the bytes in store to create the Zarr group
      await this.pyodide.globals.set("file_store", store);
      await this.pyodide.runPythonAsync(`
          import zarr

          # Create a Zarr store from the file bytes
          store = zarr.storage.KVStore({a for a in file_store})
          print(f'KVStore Keys: {list(store.keys())}')
          # group = zarr.open(store, mode='r')
      `);

      return "Zarr directory loaded";
    };
    
    // printGroupsToConsole = async () => {
    //   await this.initPromise;
    //   await this.pyodide.runPythonAsync(`
    //         def list_groups(name, node):
    //             if isinstance(node, h5py.Group):
    //                 print(f"Group: {name}")
        
    //         clgs_hdf5.visititems(list_groups)
    //     `);
    // };
  
    // getGroups = async () => {
    //   await this.initPromise;
    //   return this.pyodide.runPython(
    //     `
    //         group_names = []
    //         def visit_func(name, obj):
    //             if isinstance(obj, h5py.Group):
    //                 group_names.append(name)
            
    //         # Traverse the file and collect group names
    //         clgs_hdf5.visititems(visit_func)
    //         group_names
    //         `
    //   );
    // };
  }
  
  window.geoCluster = new GeoCLUSTER();