import os
import sys
from tqdm import tqdm
import requests
from beaupy.spinners import Spinner, DOTS
from beaupy import confirm, select, select_multiple
import argparse


def setup_args():
    args = argparse.ArgumentParser()

    args.add_argument('-f', '--file', required=True, action='store', help='Filename to download')
    args.add_argument('-t', '--token', required=True, action='store', help='Token to token')

    parser = args.parse_args()

    return parser


def download_file(file_id: str, token: str):
    with open(file_id, "wb") as f:
        link = f"https://api.gdc.cancer.gov/data/{file_id}"
        link = "https://api.gdc.cancer.gov/data/fd89bfa5-b3a7-4079-bf90-709580c006e5"
        print(f"Downloading {file_id}")
        temp_headers = {
            "X-Auth-Token": token
        }
        response = requests.get(link, stream=True, headers=temp_headers)
        total_length = response.headers.get('content-length')

        if total_length is None: # no content length header
            f.write(response.content)
        else:
            dl = 0
            total_length = int(total_length)
            for data in tqdm(response.iter_content(chunk_size=4096), "Downloading...", total=total_length, unit='bytes'):
                dl += len(data)
                f.write(data)
                done = int(50 * dl / total_length)


def main():
    args = setup_args()
    file_id = args.file
    token = args.token
    
    spinner = Spinner(DOTS, "Testing download...")
    #spinner.start()
    
    download_file(file_id, token)
    #spinner.stop()
    print(f"Download completed with code {r.status_code} - {r.text}")


if __name__ in ['__main__', 'builtin', '__buitlins__']:
    main()