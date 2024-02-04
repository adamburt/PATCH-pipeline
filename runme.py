import os
import sys
import requests
from beaupy.spinners import Spinner, DOTS
from beaupy import confirm, select, select_multiple


def main():
    args = sys.argv
    try:
        token = args[1]
    except:
        print('You must provide a token')
        sys.exit(-1)
    
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