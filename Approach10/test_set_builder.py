
import argparse
import collections
import json
import numpy as np
import os
import pandas as pd
from sklearn.preprocessing import LabelEncoder
import sys


def main():

    #directory with Gene data
    directory = '/Users/omidsharghi/Desktop/CIS_Classes/CS297/data/Approach3/PAAD Mutated Gene Details/'

    #iterate through each gene folder in the directory and create dataframes for each sequence and add data frame to list.
    data_frames = []
    for dir1 in os.listdir(directory):
        if os.path.isdir(directory+dir1):
            for dir2 in os.listdir(directory+dir1):
                if os.path.isdir(directory+dir1+'/'+dir2):
                    for file in os.listdir(directory+dir1+'/'+dir2+'/'):
                        if file == 'clinical.xls':
                            df = pd.read_excel(directory+dir1+'/'+dir2+'/'+file)
                            data_frames.append(df)


    #merge data frames
    result = pd.concat(data_frames, axis=0, join='outer', join_axes=None, ignore_index=True, sort=False)

    #drop  columns
    drop_columns = ['gender', 'year_of_birth', 'race', 'ethnicity', 'submitter_id',  'site_of_resection_or_biopsy', 'project_id', 'days_to_birth', 'vital_status', 'days_to_death',
                    'year_of_death', 'classification_of_tumor', 'last_known_disease_status', 'age_at_diagnosis',
                    'tumor_stage', 'hpv_positive_type', 'days_to_last_follow_up', 'tumor_grade', 'progression_or_recurrence', 'prior_malignancy', 'ann_arbor_b_symptoms',
                    'year_of_diagnosis', 'method_of_diagnosis', 'laterality', 'residual_disease', 'ann_arbor_clinical_stage', 'morphology',
                    'new_event_anatomic_site', 'days_to_last_known_disease_status', 'perineural_invasion_present', 'ajcc_clinical_n', 'prior_treatment',
                    'ajcc_clinical_m', 'colon_polyps_history', 'ajcc_pathologic_m', 'ajcc_pathologic_n', 'ann_arbor_pathologic_stage', 'days_to_recurrence',
                    'figo_stage', 'cause_of_death', 'lymphatic_invasion_present', 'ajcc_clinical_t', 'hpv_status', 'vascular_invasion_present', 'new_event_type',
                    'ajcc_pathologic_stage', 'burkitt_lymphoma_clinical_variant', 'circumferential_resection_margin', 'ldh_normal_range_upper',
                    'ann_arbor_extranodal_involvement', 'lymph_nodes_positive', 'ajcc_pathologic_t', 'days_to_hiv_diagnosis', 'ajcc_clinical_stage',
                    'days_to_new_event', 'hiv_positive']
    result = result.drop(drop_columns, axis=1)

    result = result.fillna(0)

    #For patients with multiple mutations, we want every combination. Therefore we need to perform the following:
    joined = pd.merge(result, result, on='case_id', suffixes=('_left', '_right')) #merge table with itself on case ID's
    filtered = joined[joined['Gene_left'] != joined['Gene_right']] #filter out rows where patient is being merged with the same data.
    transformed = filtered.apply(agg_genes, axis=1) #Remaining rows are patients with multiple mutations. Perform concatentation to format 'Gene1->Gene2'
    transformed['KRAS'] = np.where(transformed['KRAS_left'] != 0, transformed['KRAS_left'], transformed['KRAS_right']) #We know either _left or _right is zero, take the one that isn't zero
    transformed['TP53'] = np.where(transformed['TP53_left'] != 0, transformed['TP53_left'], transformed['TP53_right']) #We know either _left or _right is zero, take the one that isn't zero
    transformed['primary_diagnosis'] = transformed['primary_diagnosis_left'] #We want to preserve this column to know what kind of cancer. Left and right data will be the same.
    transformed['tissue_or_organ_of_origin'] = transformed['tissue_or_organ_of_origin_left'] #We want to preserve this column to know which organ. Left and right data will be the same.

    transformed = transformed.drop(['Gene_left', 'Gene_right', 'KRAS_left', 'KRAS_right',
                                    'TP53_left', 'TP53_right', 'primary_diagnosis_left', 'primary_diagnosis_right'], axis=1) #Drop repetitive columns
    transformed = transformed.drop_duplicates() #drop duplicates

    result = result[result["case_id"].isin(transformed["case_id"]) == False] #remove patients in the transformed df from maaster df before merging them together.

    final = pd.concat([result, transformed], axis=0, join='outer', join_axes=None, ignore_index=True, sort=False) #combine master df with transformed df
    final = final[['case_id', 'Gene', 'KRAS', 'TP53', 'primary_diagnosis', 'tissue_or_organ_of_origin']] #reorder columns

    final['pancreatic_cancer'] = 0 #add pancreatic cancer boolean column
    list = (final[final['tissue_or_organ_of_origin'].str.contains('pancreas') == True].index.tolist()) #get list of rows that are associated with pancreatic cancer
    final.loc[list, 'pancreatic_cancer'] = 1 #set pancreatic cancer column to 1 for rows found above
    final = final.drop(['case_id', 'primary_diagnosis', 'tissue_or_organ_of_origin'], axis=1) #now that we know who has pancreatic cancer, we can drop the columns that are no longer needed
    final.to_excel("output.xlsx") #Save to excel


    #Make sure columns are strings for encoding
    final['Gene'] = final['Gene'].astype(str)
    final['KRAS'] = final['KRAS'].astype(str)
    final['TP53'] = final['TP53'].astype(str)

    #encode categorical data and save to new excel file
    le_gene = LabelEncoder()
    le_kras = LabelEncoder()
    le_tp53 = LabelEncoder()
    final['Gene'] = le_gene.fit_transform(final['Gene'])
    final['TP53'] = le_tp53.fit_transform(final['TP53'])
    final['KRAS'] = le_kras.fit_transform(final['KRAS'])
    final.to_excel("output_encoded.xlsx")


def agg_genes(data):
    list =[]
    list.append(data.Gene_left)
    list.append(data.Gene_right)
    list.sort()
    data['Gene'] = str(list[0]) + '->' + str(list[1])
    return data

if __name__ == "__main__":
    main()
