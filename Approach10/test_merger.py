import requests
import csv
import json
import pandas as pd
import os
from sklearn.preprocessing import LabelEncoder

##NOT USED ANYMORE

def build_df():
    columns = ["case", "gene", "KRAS", "TP53", "SMAD4"]
    df  = pd.DataFrame(columns = columns)

    df.loc[len(df)] = ["Omid", "KRAS", 1, 0, 0]
    df.loc[len(df)] = ["John", "TP53", 0, 1, 0]
    df.loc[len(df)] = ["Omid", "TP53", 0, 2, 0]
    df.loc[len(df)] = ["Jake", "SMAD4", 0, 0, 5]
    df.loc[len(df)] = ["Chris", "KRAS", 3, 0, 0]
    df.loc[len(df)] = ["David", "KRAS", 1, 0, 0]
    df.loc[len(df)] = ["Elton", "TP53", 0, 1, 0]
    df.loc[len(df)] = ["John", "KRAS", 3, 0, 0]
    df.loc[len(df)] = ["Jake", "KRAS", 1, 0, 0]
    df.loc[len(df)] = ["Roger", "SMAD4", 0, 0, 2]

    print(df.to_string())
    return df

def encode_df(df):
    
    #Encode all columns except for gene
    encoder_list = []
    columns = list(df.columns.values)

    
    for i in range(len(columns)):
        encoder_list.append(LabelEncoder())
    
    df[columns[0]] = encoder_list[0].fit_transform(df[columns[0]])
    for i in range(1, len(columns)):
        if(columns[i] != "gene"):
            df[columns[i]] = encoder_list[i].fit_transform(df[columns[i]])

    print("ENCODE:")
    print(df.to_string())
    return df

def get_dups(df):
    #Get DF of just duplicate cases
    dups = df[df.duplicated(['case'], keep=False)]
    print("Dups:")
    print(dups.to_string())
    return dups

def get_dups_list(dups):
    #Get list of duplicate case ID's
    unique_cases = dups["case"].unique()
    return unique_cases


def merger(dup_case_list, dup_df):
    
    #iterate through list of duplicate case_id's and
    #create a df for each unique case_id from the df of duplicate cases
    df_list = []
    for case in dup_case_list:
        one_case_df = dup_df[dup_df['case'] == case]
        df_list.append(one_case_df)

    #create empty df
    columns = list(dup_df.columns.values)
    final_df = pd.DataFrame(columns = columns)

    #iterate through list of unique case dfs and merge rows to get a single row df
    #merge dfs to get final df with only unqieu cases
    for df in df_list:
        case = df.loc[df.index[0], 'case']
        merged_rows_df  = merge_df_rows(case, df)
        final_df = final_df.append(merged_rows_df, sort=True)

    print("FINAL MERGED")
    print(final_df.to_string())
    return final_df

def merge_df_rows(case, df):
    final = df.copy()
    new_df = pd.DataFrame()
    new_df['gene'] = final.groupby('case')['gene'].apply(lambda x: sorted(x.unique())).transform(lambda x: '->'.join(x))
    new_df['case'] = case
    cols_to_agg = [col for col in final if col!="gene" and col!="case"]
    print(cols_to_agg)
    new_df[cols_to_agg] = final.groupby('case')[cols_to_agg].agg(max)
    print("MERGED:")
    print(new_df.to_string())
    return new_df

def delete_dups_from_original_df(df_orig_encoded, dup_list):
    #delete duplicate case id's from original df to prepare for joining unique case dfs
    df = df_orig_encoded.copy()
    for case in dup_list:
        df = df[df["case"] != case]
    return df


def main():
    df = build_df()
    df_encode = encode_df(df)
    dups_df = get_dups(df_encode)
    dup_case_list = get_dups_list(dups_df)
    merged_dups_df = merger(dup_case_list, dups_df)
    unique_case_df = delete_dups_from_original_df(df_encode, dup_case_list)
    unique_case_df = unique_case_df.append(merged_dups_df, sort=True)
    print('END')
    print(unique_case_df.to_string())



if __name__ == "__main__":
    main()
