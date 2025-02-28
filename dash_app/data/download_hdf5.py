# Don't run this locally! This is only for the CICD pipeline to download the hdf5 file. 
# If you want to get the updated hdf5 file, just run decompress_hdf5.py script
import boto3
import botocore
import os

def download_s3_file(bucket_name, object_key, file_name):
    """
    Download a file from S3 and save it locally
    :param bucket_name: S3 bucket name
    :param object_key: S3 object key (full path)
    :param file_name: Local file name to save as
    """
    s3 = boto3.client('s3')
    
    try:
        s3.download_file(bucket_name, object_key, file_name)
        print(f"File downloaded successfully to {file_name}")
    except botocore.exceptions.ClientError as e:
        if e.response['Error']['Code'] == "404":
            print("The object does not exist.")
        else:
            raise
    except Exception as e:
        print(f"Error downloading file: {str(e)}")
        raise

def download_hdf5():
    h5_filepath = "./decompressed_clgs_results_final_float32.h5"
    if os.path.exists(h5_filepath):
        print("Decompressed h5 file already exists, not downloading")
        return 

    download_s3_file(
        bucket_name='geocluster-data',
        object_key='development/chunked_clgs_results_final_float32.h5',
        file_name=h5_filepath
    )

if __name__ == "__main__":
    download_hdf5()