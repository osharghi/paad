import requests
import csv
from itertools import chain
from itertools import combinations
import json
import pandas as pd
import os
from sklearn.preprocessing import LabelEncoder


def get_dups(df):
    #Get DF of just data of duplicate cases
    dups = df[df.duplicated(['case_id'], keep=False)]
    return dups

def get_dups_list(dups):
    #Get list of duplicate case ID's
    unique_cases = dups["case_id"].unique()
    return unique_cases


def merger_new(dup_case_list, dup_df):
    #iterate through list of duplicate case_id's and
    #create a df for each unique case_id from the df of duplicate cases
    df_list = []
    
    for case in dup_case_list:
        one_case_df = dup_df[dup_df['case_id'] == case] #all data of a single duplicate case
        df_list.append(one_case_df)

    #create empty df
    columns = list(dup_df.columns.values)
    final_df = pd.DataFrame(columns = columns)
    final_df['mutations'] = ""
    final_df['gene_mutations'] = ""

    #for each df of duplicate case_id's, extract case_id and project_id
    #and call function that returns df of powerset
    for i in range(0, len(df_list)):
        df = df_list[i]
        #Setting a cap on the number of possible mutation combinations so we don't have 2^81 combos again
        if len(df)<11:
            case_list = df['case_id'].unique()
            project_list = df['project_id'].unique()
            case = case_list[0]
            project = project_list[0]
            merged_rows_df  = merge_df_rows_new(case, project, df, i, len(df_list))
            final_df = final_df.append(merged_rows_df, sort=True)
        
    print("FINISHED POWER SETS AND MERGING ALL DUP CASES")
#    print(final_df.to_string())
    return final_df

def merge_df_rows_new(case, project, df, case_index, num_of_cases):

    #will be an array of gene, mutation pairs [[gene, mutation]]
    row_gene_mutation_list = []
    
    for (index, row) in df.iterrows():
        #for each row, create a list of [gene,  mutation]
        #add list to row_gene_mutation_list
        row_list = []
        gene = row['gene']
        mutation = row[gene]
        row_list.append(gene)
        row_list.append(mutation)
        row_gene_mutation_list.append(row_list)

    #Create new df to store results of power set
    columns = list(df.columns.values)
    powerset_df = pd.DataFrame(columns = columns)
    powerset_df['mutations'] = ""
    powerset_df['gene_mutations'] = ""

    #PowerSet
    
    set = powerset(row_gene_mutation_list)
#    set = list(powerset(row_gene_mutation_list))

    subset_index = 0
    total_subsets = pow(2, len(row_gene_mutation_list))

    for subset in set:
        #get subset in each powerset
        print("Case_ID: " + case + " Case_Index: " + str(case_index) + " : " + str(num_of_cases-1) + " Subset: " + str(subset_index) + " : " + str(total_subsets-1))
        gene_mut_dict = {}
        for item in subset:
            #get each gene mutation pair in subset and add to dictionary
            #dict structure = gene:[mutations]
            gene = item[0]
            mut = str(item[1])
            if gene in gene_mut_dict:
#                gene_mut_dict[gene].append(mut)
                insertion_sort(gene_mut_dict[gene], mut)
            else:
                mut_list = []
                mut_list.append(mut)
                gene_mut_dict[gene]=mut_list
    
        if len(gene_mut_dict) != 0:
            #initalize df row with case_id and project_id
            powerset_df.loc[len(powerset_df), 'case_id'] = case
            powerset_df.loc[len(powerset_df)-1, 'project_id'] = project
            
            #initalize empty list to store genes and a dictionary of sorted mutations list for each gene
            gene_list = []
            gene_mut_dict_sorted = {}
            
            for key in gene_mut_dict:
                #get mutation list for each gene, sort, and then add to new sorted dictionary
#                gene_list.append(key) Commented out because now using insertion sort
                insertion_sort(gene_list, key)
                mut_list_for_gene = gene_mut_dict[key]
#                mut_list_for_gene.sort() Commented out because now using insertion sort
                gene_mut_dict_sorted[key] = mut_list_for_gene
                
                #build string mut1-mut2-mut3 for each gene and add to df row initalized above
                mut_total_str = ""
                for i in range(len(mut_list_for_gene)):
                    if(i==0):
                        mut_total_str = mut_list_for_gene[i]
                    else:
                        mut_total_str = mut_total_str + "-" + mut_list_for_gene[i]
                
                gene_key = key
                powerset_df.loc[len(powerset_df)-1, gene_key] = mut_total_str

            #sort list of genes
#            gene_list.sort() Commented out because now using insertion sort

            #initalize string for gen1-mut1-gene2-mut2 format
            genes_muts_combo_str = ""
            
            #intalize empty list to hold all mutations for all genes in subset to build mut1-mut2-mut3 format
            mut_list_total = []
            
            for j in range(len(gene_list)):
                #build gene1-mut1-gene2-mut2 string
                #add sorted mutation list for each gene to master list to build mut1-mut2-mut3 string further down
                gene = gene_list[j]
                mut_list_for_gene = gene_mut_dict_sorted[gene]
                mut_list_total.extend(mut_list_for_gene)
                
                for i in range(0, len(mut_list_for_gene)):
                    if j==0 and i==0:
                        genes_muts_combo_str = gene + "-" + mut_list_for_gene[i]
                    else:
                        genes_muts_combo_str = genes_muts_combo_str + "-" + gene + "-" + mut_list_for_gene[i]

            powerset_df.loc[len(powerset_df)-1, 'gene_mutations'] = genes_muts_combo_str
            genes_str = "-"
            mutations_str = "-"
            
            #build gene1-gene2-gene3 string
            genes_str = genes_str.join(gene_list)
            mut_list_total.sort()
            
            #build mut1-mut2-mut3 string for all mutations in subset
            mutations_str = mutations_str.join(mut_list_total)
            
            #update df row
            powerset_df.loc[len(powerset_df)-1, 'mutations'] = mutations_str
            powerset_df.loc[len(powerset_df)-1, 'gene'] = genes_str
            subset_index = subset_index+1
                
    print("POWERSET TOTAL: ")
    return powerset_df

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
    "powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(len(s)+1))


def delete_dups_from_original_df(df_orig_encoded, dup_list):
    #delete duplicate case id's from original df to prepare for joining unique case dfs
    df = df_orig_encoded.copy()
    for case in dup_list:
        df = df[df["case_id"] != case]
    return df

def insert_mutations_gene_column_original_df(df):
    df['mutations'] = ""
    df['gene_mutations'] = ""
    for (index, row) in df.iterrows():
        gene = row['gene']
        mutation = row[gene]
        row['mutations'] = mutation
        row['gene_mutations'] = gene + "-" + str(mutation)
    return df

def add_pancreatic_cancer_column(df):
    df['pancreatic_cancer'] = ""
    for (index, row) in df.iterrows():
        if(row['project_id'] == "TCGA-PAAD"):
            row['pancreatic_cancer'] = 1
        else:
            row['pancreatic_cancer'] = 0
    df = df.drop(columns=['project_id'])
    return df

def insertion_sort(list, str):
    if(len(list) < 1):
        list.append(str)
        return list
    else:
        list.append(str)
        i = len(list)-1
        if(list[i]<list[i-1]):
            key = str
            j = i-1
            while(j>=0 and key<list[j]):
                list[j+1] = list[j]
                j = j-1
            list[j+1] = key

def merge_and_transform(df):
    dups_df = get_dups(df)
    dup_case_list = get_dups_list(dups_df)
    merged_dups_df = merger_new(dup_case_list, dups_df)
    unique_case_df = delete_dups_from_original_df(df, dup_case_list)
    unique_case_df = insert_mutations_gene_column_original_df(unique_case_df)
    unique_case_df = unique_case_df.append(merged_dups_df, sort=True)
    unique_case_df = unique_case_df.reset_index(drop=True)
    unique_case_df = unique_case_df.fillna(0)
    final = add_pancreatic_cancer_column(unique_case_df)
#    final.to_excel("merged_output.xlsx")
    return final


def main():
    df = pd.read_excel('run_output.xlsx', index_col=0)
    merge_and_transform(df)


if __name__ == "__main__":
    main()
