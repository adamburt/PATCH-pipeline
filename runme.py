import os
import sys
import requests
from beaupy.spinners import Spinner, DOTS
from beaupy import confirm, select, select_multiple
import argparse


def setup_args():
    args = argparse.ArgumentParser()

    args.add_argument('-u', '--url', action='store', help='URL to download')
    args.add_argument('-t', '--token', action='store', help='Token to token')

    parser = args.parse_args()

    return parser


def main():
    args = setup_args()
    url = args.url
    token = args.token
    
    spinner = Spinner(DOTS, "Testing download...")
    spinner.start()
    url = ""
    temp_headers = {
        "X-Auth-Token": token
    }
    r = requests.get(url, headers=temp_headers)
    spinner.stop()
    print(f"Download completed with code {r.status_code}")


if __name__ in ['__main__', 'builtin', '__buitlins__']:
    main()