#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 21 11:47:43 2025

@ author: L-F-S

various utilities to run once before (por during) running the pipeline.


"""
import os
if not os.getcwd().endswith('modules'):
    os.chdir('modules')
import sys
import numpy as np
import pandas as pd
import pickle
import pyreadr
import time
from conf import BASE_DIR, DISEASE,MITH_IN_DRUG, MITH_OUT_DRUG,\
 TSR_OUT_DRUG, TSR_OUT_DISEASE, TSR_OUT_CSCORE, alias_2geneid

#%% Create a dictionary that maps gene symbols to gene ids using
 # the updated map provided by Alaimo gene_map.tsv
# symbol_to_id_data =pd.read_csv(BASE_DIR+'other_data/gene_map.tsv', sep='\t', dtype={'NCBI.Gene.ID':str}, usecols=[0,3])
# #symbol_to_id.dropna(inplace=True)
# symbol_to_id_transformed=symbol_to_id_data.set_index('Approved.symbol').to_dict()
# symbol_to_id_dict=symbol_to_id_transformed['NCBI.Gene.ID']

# symbol_to_id_data[symbol_to_id_data['Approved.symbol']==gene]
# #%%
# with open(BASE_DIR+'other_data/symbol_2geneid.pkl', 'wb') as f:
#     pickle.dump(symbol_to_id_dict, f)
#%% filter DEG drug signatures by only keeping non duplicate gene ids 
#  To run AFTER  creation of LINCS DEG data

def load_unfiltered_single_drug_signature(drug, mith=False):
    '''
    Loads UNFILTERED DEG drug wise drug signature

    inputs:
        DISEASE: str disease symbol (directory name)

        mith: bool
            def. False: loading DEG FC signature. if True, loading MIthRIL 
            perturbation signature.

    returns a pandas.Dataframe where rows are genes.
    '''
    
    if not mith:
        result = pyreadr.read_r(TSR_OUT_DRUG+'/LINCS/metanalysis_drug_wise/'+drug+'_metanalysis.Rds')
        return result[None].reset_index(drop=True)
    
    else:
        return pd.read_csv(TSR_OUT_DRUG+'/LINCS/metanalysis_mith3_drug_wise/'+drug+'_metanalysis.csv', sep='\t')
    
def convert_and_remove_duplicates(DEG_drug_signature):
    '''
    Used to map gene symbol to gene id for DEG FC data for LINCS
    drug wise signature data.
    Works with both a single drug, or all drugs
    Output:
        pd.DataFrame of FC signatures, with unique gene ids.
    '''
    
    def map_name_to_id(genename):
        if genename in alias_2geneid.keys():
            return str(alias_2geneid[genename])
        return np.nan

    DEG_drug_signature['gene_id']=DEG_drug_signature['gene'].apply(lambda x : map_name_to_id(x))
    
    nbefore=DEG_drug_signature.shape[0]
    DEG_drug_signature.dropna(subset='gene_id', inplace=True)
    nafter=DEG_drug_signature.shape[0]
    print(nafter-nbefore, 'genes not found in symbol: gene id mapping')

    DEG_drug_signature.drop_duplicates(subset='gene_id', keep='first', inplace=True)
    print(DEG_drug_signature.shape[0]-nafter,'duplicate gene ids removed')

    
    return DEG_drug_signature

def get_drugs_list():
    drugs_list=os.listdir(MITH_OUT_DRUG)
    
    drugs_list=[x.split('.')[0] for x in drugs_list]
    drugs_list=[x.split('_')[0] for x in drugs_list]
    drugs_list=list(np.unique(drugs_list))
    return drugs_list

def get_common_genes(disease_signature, drug_signature):
    missing_genes=set.difference(set(disease_signature.gene_id),set(drug_signature.gene_id))
    # print('missing genes from mithril', len(missing_genes))
    addional_drug_items=set.difference(set(drug_signature.gene_id),set(disease_signature.gene_id))
    # print('addional_mith_items', len(addional_drug_items))
    
    
    common_genes=set.intersection(set(disease_signature.gene_id),set(drug_signature.gene_id))
    # print('genes in common', len(common_genes))
    if not len(common_genes)+len(addional_drug_items)==drug_signature.shape[0]:
        raise ValueError('common genes+ additional mithril items != len(mithril_gene)')
    if not len(common_genes)+len(missing_genes)==disease_signature.shape[0]:
        raise ValueError('riprova: common genes+ missing mithril items != len(DEG_genes)')

    # Get subset of signatures with common genes
    drug_common_genes_signatures=drug_signature[drug_signature['gene_id'].isin(common_genes)].sort_values(by='gene_id')
    disease_common_genes_signatures=disease_signature[disease_signature.gene_id.isin(common_genes)].sort_values(by='gene_id')
    
    if not np.unique(disease_signature.iloc[disease_common_genes_signatures.index].reset_index()['gene_id']==drug_signature.iloc[drug_common_genes_signatures.index].reset_index()['gene_id']):
        raise ValueError('common gene indexes actually different between drug and disease!')
    return disease_common_genes_signatures.index, drug_common_genes_signatures.index

#%%

if __name__=='__main__':
    
    # called by 7_filter_all_drugs
    drugs_list =get_drugs_list()
    i1=int(sys.argv[1])
    i2=int(sys.argv[2])
    
    start=time.time()
    for h, drug in enumerate(drugs_list[i1:i2]):
        print(drug)
        output_filename=TSR_OUT_DRUG+'/LINCS/metanalysis_drug_wise_filtered/'+drug+'_metanalysis.csv'  
    
        if not os.path.isfile(output_filename): 
            
            DEG_drug_signature=load_unfiltered_single_drug_signature(drug, mith=False)
            filtered_df = convert_and_remove_duplicates(DEG_drug_signature)
            
            filtered_df.to_csv(output_filename, sep='\t', index=False)
