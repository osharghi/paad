import requests
import csv
import json
import pandas as pd
import os

#Test file

def main():

    gene = "KRAS"
    dna_change = "chr12:g.25245350C>T"
    
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
                       "field": "genomic_dna_change",
                       "value": dna_change
                       }
                       }
                       ]
    }

    fields = ["occurrence.case.case_id", "occurrence.case.project.project_id"]

    fields = ",".join(fields)

    params = {"filters":json.dumps(filt),"fields": fields,"format": "JSON", "size": "32000"}

    response = requests.get(url, params = params)
    result = json.loads(response.content.decode("utf-8"))

    print(response.content.decode("utf-8"))

    columns = ["gene", "mutation", "project_id", "case_id"]
    cases_df  = pd.DataFrame(columns = columns)

    for item in result.get("data").get("hits"):
        occurrence = item["occurrence"]
        for cases in occurrence:
            case = cases["case"]
            project_id = case["project"]["project_id"]
            case_id = case["case_id"]
            cases_df.loc[len(cases_df)] = [gene, dna_change, project_id, case_id]

    return cases_df



if __name__ == "__main__":
    main()


