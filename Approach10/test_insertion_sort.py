import requests
import csv
import json
import pandas as pd
import os
from sklearn.preprocessing import LabelEncoder


def insertion_sort(list, color):
    if(len(list) < 1):
        list.append(color)
        return list
    else:
        list.append(color)
        i = len(list)-1
        if(list[i]<list[i-1]):
            key = color
            j = i-1
            while(j>=0 and key<list[j]):
                list[j+1] = list[j]
                j = j-1
            list[j+1] = key





def main():
    list = ["green"]
    insertion_sort(list, "blue")
    insertion_sort(list, "black")
    insertion_sort(list, "cyan")
    insertion_sort(list, "gray")
    insertion_sort(list, "red")
    print(list)




if __name__ == "__main__":
    main()
