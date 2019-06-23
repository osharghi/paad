import requests
import csv
import json
import pandas as pd
import os

#This module loads the list of genes xlsx and queries the cases associated with each gene
#and keeps a count of the number of cases for each gene that is associated with PAAD and
#saves result to genes_cases_output.xlsx


gene_df = pd.read_excel('init_genes_list.xlsx', index_col=0)

for (index, row) in gene_df.iterrows():
    symbol = row["gene"]

    url = "https://api.gdc.cancer.gov/ssms"

    filt = {
        "op": "and",
            "content":[
                       {
                       "op": "=",
                       "content":{
                       "field": "occurrence.case.project.project_id",
                       "value": "TCGA-PAAD"
                       }
                       },
                       {
                       "op": "=",
                       "content":{
                       "field": "consequence.transcript.gene.symbol",
                       "value": symbol
                       }
                       }
                       ]
    }

    fields = ["occurrence.case.case_id", "occurrence.case.project.project_id"]

    fields = ",".join(fields)

    params = {"filters":json.dumps(filt),"fields": fields,"format": "JSON", "size": "32000"}

    response = requests.get(url, params = params)
    result = json.loads(response.content.decode("utf-8"))

    #print(response.content.decode("utf-8"))

    columns = ["cases"]
    case_df  = pd.DataFrame(columns = columns)

    for item in result.get("data").get("hits"):
        occurrence = item["occurrence"]
        for sub_item in occurrence:
            case = sub_item["case"]
            if case["project"]["project_id"] == "TCGA-PAAD":
                case_df.loc[len(case_df)] = [case["case_id"]]


#    case_df.drop_duplicates(keep=False, inplace=True)
    case_df = case_df.drop_duplicates()
    i = int(index)
    gene_df.ix[i, 'cases'] = case_df.shape[0]
    row = gene_df['cases'].values[i]

gene_df.to_excel("init_genes_case_count.xlsx")


