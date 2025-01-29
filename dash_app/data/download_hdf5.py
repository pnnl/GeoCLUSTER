import boto3
import botocore

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

if __name__ == "__main__":
    # Usage example
    download_s3_file(
        bucket_name='geocluster-data',
        object_key='development/chunked_clgs_results_final_float32.h5',
        file_name='./clgs_results_final_float32.h5'
    )