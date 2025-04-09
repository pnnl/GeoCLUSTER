import requests
import os


def download_file(url, local_filename):
    """
    Downloads a file from a URL to a local file.

    Args:
        url (str): The URL of the file to download.
        local_filename (str): The local filename to save the downloaded file to.
    """
    try:
        with requests.get(url, stream=True) as r:
            r.raise_for_status()  # Raise an exception for HTTP errors
            with open(local_filename, 'wb') as f:
                for chunk in r.iter_content(chunk_size=8192):
                    f.write(chunk)
        print(f"Downloaded {url} to {local_filename}")
    except requests.exceptions.RequestException as e:
        print(f"An error occurred: {e}")

def download_hdf5():
    h5_filepath = "./decompressed_clgs_results_final_float32.h5"
    if os.path.exists(h5_filepath):
        print("Decompressed h5 file already exists, not downloading")
        return 
    
    print("downloading h5 file")
    download_file("https://geocluster-data.s3.us-west-2.amazonaws.com/development/chunked_clgs_results_final_float32.h5", h5_filepath)

if __name__ == "__main__":
    download_hdf5()