import requests
import csv
import json
from itertools import chain
from itertools import combinations
import numpy as np
import pandas as pd
import os
from sklearn.preprocessing import LabelEncoder


def build_df():
#    data = [["TP53", "T1"], ["KRAS", "K1"], ["TP53", "T2"], ["KRAS", "K2"]]
#    set = list(powerset(data))
#    for subset in set:
#        print("subset")
#        for item in subset:
#            print(item)

    genes = ["TP53", "SMAD4", "KRAS"]
    mutations = [1,2,3]
    sort(genes, mutations)
    print(genes)
    print(mutations)

def sort(genes, mutations):
    if(len(genes)>1):
        for i in range(1, len(genes)):
            if(genes[i]<genes[i-1]):
                g_key = genes[i]
                m_key = mutations[i]
                j = i-1
                while(j>=0 and g_key<genes[j]):
                    genes[j+1] = genes[j]
                    mutations[j+1] = mutations[j]
                    j = j-1
                genes[j+1] = g_key
                mutations[j+1] = m_key


def powerset(iterable):
    print('here')
    "powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(len(s)+1))

def main():
    df = build_df()




if __name__ == "__main__":
    main()
