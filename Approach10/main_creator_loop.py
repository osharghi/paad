import collections
import json
import main_creator
import merge_transform
import os
import pandas as pd
import requests
import sys
import time


def start():
    percent = 15.5
    for i in range (0,1):
        user_input = [0, 'none', percent, percent, 1]
        main_creator.main(user_input)
        percent += 3
#        print("Sleeping for 20 mins to cool down")
##        time.sleep(1200)
    print(percent)


if __name__ == "__main__":
    start()
