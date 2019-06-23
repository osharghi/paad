import collections
import json
import merge_transform
import os
import pandas as pd
import requests
import sys
import time

def get_user_input(lines):
    input=[]
    
    #upload file
    step_1 = lines[1].replace('value:', '')
    step_1 = step_1.replace(' ', '')
    step_1 = int(step_1)
    input.append(step_1)
    
    #file path
    step_2 = lines[4].replace('path:', '')
    step_2 = step_2.replace(' ', '')
    input.append(step_2)
    
    #gene percentage
    step_3_genes = lines[7].replace('genes:', '')
    step_3_genes = step_3_genes.replace(' ', '')
    step_3_genes = float(step_3_genes)
    input.append(step_3_genes)
    
    #mutation percentage
    step_3_mutations = lines[8].replace('mutations:', '')
    step_3_mutations = step_3_mutations.replace(' ', '')
    step_3_mutations = float(step_3_mutations)
    input.append(step_3_mutations)

    #Save query results
    step_4 = lines[11].replace('value:', '')
    step_4 = step_4.replace(' ', '')
    step_4 = int(step_4)
    input.append(step_4)

    return input

def get_title(gene_flt, mut_flt):
    gene_str = str(gene_flt)
    mut_str = str(mut_flt)
    
#    gene_str.replace('.', ',')
#    mut_str.replace('.', ',')

    title = "G" + gene_str + "M" + mut_str

    return title

def save_time(path, title, time):
    time_str = str(time)
    text_file = open(path + title + "_" + "time.txt", "w")
    text_file.write("Total time: %s" % time)
    text_file.close()

def trim_df(df, percentage):
    row_count = df.shape[0]
    rows_needed = row_count*(percentage/100)
    if(rows_needed<1):
        return df
    else:
        cut_off = int(-1*(row_count-rows_needed))
        df = df[:cut_off]
        return df

def join_dfs(df_list):
    
    if(len(df_list)>1):
        result = df_list[0]
        for i in range(1, len(df_list)):
            result = result.append(df_list[i])
        return result
    elif len(df_list) == 1:
        return df_list[0]
    else:
        return None

def get_top_mutations(gene_df, mutations_percent):
    mutations_df_list = []
    for (index, row) in gene_df.iterrows():
        mutation_df = get_mutations_for_gene(row['gene'])
        mutation_df = trim_df(mutation_df, mutations_percent)
        mutations_df_list.append(mutation_df)
    return mutations_df_list

def get_cases_for_mutation(df_mutations):
    cases_df_list = []
    for (index, row) in df_mutations.iterrows():
        case_df = query_cases_for_mutation(row['gene'], row['mutation'])
        cases_df_list.append(case_df)

    return cases_df_list


def get_mutations_for_gene(gene):
    cases_endpt = "https://api.gdc.cancer.gov/ssms"
    
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
                   "value": gene
                   }
                   }
                   ]
    }

    fields = ["genomic_dna_change", "occurrence.case.case_id", "occurrence.case.project.project_id"]

    fields = ",".join(fields)

    params = {
    "filters":json.dumps(filt),
    "fields": fields,
        "format": "JSON",
        "size": "32985"
    }

    response = requests.get(cases_endpt, params = params)
    result = json.loads(response.content.decode("utf-8"))

#    print(response.content.decode("utf-8"))
#    print(json.dumps(response.json(), indent=2))

    mutation_dict = {}
    columns = ["gene", "mutation", "cases"]

    mutation_df = pd.DataFrame(columns = columns)
    
    for item in result.get("data").get("hits"):
        dna_change = item.get("genomic_dna_change")
        occurence_list = item.get("occurrence")
        for sub_item in occurence_list:
            project_id = sub_item['case']['project']['project_id']
            if(project_id == "TCGA-PAAD"):
                if dna_change in mutation_dict.keys():
                    count = mutation_dict[dna_change]
                    count = count + 1
                    mutation_dict[dna_change] = count
                else:
                    mutation_dict[dna_change] = 1

    for key, value in mutation_dict.items():
        mutation_df.loc[len(mutation_df)] = [gene, key, value]

    mutation_df = mutation_df.sort_values('cases', ascending=False)
    return mutation_df

def query_cases_for_mutation(gene, dna_change):
    
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

    #    print(response.content.decode("utf-8"))

#    columns = ["gene", "mutation", "project_id", "case_id"]
    columns = ["gene", gene, "project_id", "case_id"]

    cases_df  = pd.DataFrame(columns = columns)

    for item in result.get("data").get("hits"):
        occurrence = item["occurrence"]
        for cases in occurrence:
            case = cases["case"]
            project_id = case["project"]["project_id"]
            case_id = case["case_id"]
            cases_df.loc[len(cases_df)] = [gene, dna_change, project_id, case_id]
        
    return cases_df


def main():
    
    lines = []
    user_input = []
    
    try:
        lines = [line.rstrip('\n') for line in open('menu.txt')]
    except:
        print('error')
        sys.exit()
        
    user_input = get_user_input(lines)
    if(user_input[0] == 0):
        pd.set_option('display.max_columns', None)
        gene_percent = user_input[2]
        mutations_percent = user_input[3]
        title = get_title(gene_percent, mutations_percent)
        if not os.path.exists(title):
            os.makedirs(title)
        path = os.path.dirname(os.path.abspath(__file__))
        path = path + "/" + title + "/"
        start = time.time()
        gene_df = pd.read_excel('init_genes_case_count.xlsx', index_col=0)
        #Sort dataframe by case count in descending order.
        gene_df = gene_df.sort_values('cases', ascending=False)
        gene_df = trim_df(gene_df, gene_percent)
        print("Querying and calculating mutations for genes...")
        mutations_df_list = get_top_mutations(gene_df, mutations_percent)
        print("Merging mutation dataframes...")
        mutations_df = join_dfs(mutations_df_list)
        if(user_input[4] == 1):
            print("Saving mutations dataframe to file...")
            print(mutations_df.to_string())
            file_name_mutations = path + title + "_" + "output_top_mutations.xlsx"
            mutations_df.to_excel(file_name_mutations)
            print("Top mutations saved to: " + file_name_mutations)
        print("Getting cases...")
        cases_df_list = get_cases_for_mutation(mutations_df)
        print("Merging cases dataframes...")
        cases_df = join_dfs(cases_df_list)
        cases_df = cases_df.fillna(0)
        if(user_input[4] == 1):
            print("Saving cases dataframe to file...")
            file_name_cases = path +title + "_" + "output_cases_for_mutations.xlsx"
            cases_df.to_excel(file_name_cases)
            print("Cases for mutations saved to: " + file_name_cases)
        final = merge_transform.merge_and_transform(cases_df)
        file_name_transformed = path + title + "_" + "output_transformed.xlsx"
        final.to_excel(file_name_transformed)
        print("Transformation results saved to: " + file_name_transformed)
        end = time.time()
        diff = end-start
        save_time(path, title, diff)
    elif(user_input[0] == 1):
        if(user_input[1] == 'none'):
            print("You must provide a path to the cases file")
        else:
            cases_df = pd.read_excel('output_cases_for_mutations.xlsx', index_col=0)
            final = merge_transform.merge_and_transform(cases_df)
            file_name_transformed = "output_transformed.xlsx"
            final.to_excel(file_name_transformed)
            print("Transformation results saved to: " + file_name_transformed)


if __name__ == "__main__":
    main()
