import requests
import csv
import json
import pandas as pd
import os

#Use this module to get a list of genes for pancreatic cancer and output results to xlsx

url = "https://api.gdc.cancer.gov/genes"

filt = {
    "op": "and",
    "content":[ {"op": "=","content":{"field": "case.project.project_id", "value": "TCGA-PAAD"}}]
}

#fields = ["symbol", "case.case_id", "case.project.project_id"]
fields = ["symbol", "gene_id"]


fields = ",".join(fields)

params = {"filters":json.dumps(filt),"fields": fields,"format": "JSON", "size": "32000"}

response = requests.get(url, params = params)
result = json.loads(response.content.decode("utf-8"))

#print(response.content.decode("utf-8"))

columns = ["gene", "gene_id", "cases"]
gene_df  = pd.DataFrame(columns = columns)

for item in result.get("data").get("hits"):
    gene_df.loc[len(gene_df)] = [item["symbol"], item["id"], 0]

gene_df.to_excel("init_genes_list.xlsx")
#print(gene_df.to_string())



