#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 13 09:50:20 2024

@ author: L-F-S

Maps  Mithril 3 batch output into
tsr metanalysis-like input matrix.
"""

import os
if not os.getcwd().endswith('modules'):
    os.chdir('modules')
import sys
import numpy as np
import pandas as pd
import time
from conf import MITH_OUT_DRUG, TSR_OUT_DRUG
from preprocessing_utils import get_drugs_list


def remove_special_characters(df):
    special_chars=';/:*?\"|'
    for char in special_chars:
        indexes=df[df['Gene Name'].str.contains('\\'+char)].index
        if len(indexes)>0:
            df.drop(index=indexes, inplace=True)
    return df


#%%
print('loading drugs list')

drugs_list=get_drugs_list()

tot_drugs=len(drugs_list)

        





# i1=int(sys.argv[1])
# i2=int(sys.argv[2])
print('Mapping mith3 output into tsr connectivity input: \n\
          Removing pathway duplicates, Merging three timepoint datasets in a single file and filtering ')
start=time.time()
for h, drug in enumerate(drugs_list[0:20]):
    
    output_filename=TSR_OUT_DRUG+'LINCS/metanalysis_mith3_drug_wise/'+drug+'_metanalysis.csv'   
    if not os.path.isfile(output_filename): 
    
    
        drug_start=time.time()
        print(drug, h+1, 'of', tot_drugs)
        #6h
        mith_perturb_signature_file_6h=drug+'_6h.perturbations.txt'
        mith_perturb_signature_6h = pd.read_csv(MITH_OUT_DRUG+mith_perturb_signature_file_6h, sep='\t', index_col=False, engine='python', usecols=['Gene Id', 'Gene Name', 'Perturbation', 'pValue', 'adj_pValue'])
        mith_perturb_signature_6h.drop_duplicates(inplace=True, keep='first') # Duplicate gene ids for pathways (reduces row number from 250k to 14k)
        mith_perturb_signature_6h.rename(columns={'Gene Id':'gene_id', 'Perturbation':'Perturbation_6h', 'pValue':'p.value_6h', 'adj_pValue':"adj.p.value_6h"}, inplace=True)
        mith_perturb_signature_6h['drug'] = drug
        mith_perturb_signature_6h['t.value_like_statistic_6h']=mith_perturb_signature_6h['p.value_6h'].apply(lambda p_value : 2001*(1-p_value))
    
        
        #24h
        mith_perturb_signature_file_24h=drug+'_24h.perturbations.txt'
        mith_perturb_signature_24h = pd.read_csv(MITH_OUT_DRUG+mith_perturb_signature_file_24h, sep='\t', index_col=False, engine='python', usecols=['Gene Id', 'Perturbation', 'pValue', 'adj_pValue'])
        mith_perturb_signature_24h.drop_duplicates(inplace=True, keep='first')
        mith_perturb_signature_24h.rename(columns={'Gene Id':'gene_id', 'Perturbation':'Perturbation_24h', 'pValue':'p.value_24h', 'adj_pValue':"adj.p.value_24h"}, inplace=True)
        mith_perturb_signature_24h['t.value_like_statistic_24h']=mith_perturb_signature_24h['p.value_24h'].apply(lambda p_value : 2001*(1-p_value))
    
    
        mith_perturb_signature_meta=mith_perturb_signature_6h.merge(mith_perturb_signature_24h, on='gene_id', how='left')
       
        #6h 24h
        mith_perturb_signature_file_6h_24h=drug+'_6h_24h.perturbations.txt'
        mith_perturb_signature_6h_24h = pd.read_csv(MITH_OUT_DRUG+mith_perturb_signature_file_6h_24h, sep='\t', index_col=False, engine='python', usecols=['Gene Id','Perturbation', 'pValue', 'adj_pValue'])
        mith_perturb_signature_6h_24h.drop_duplicates(inplace=True, keep='first')
        mith_perturb_signature_6h_24h.rename(columns={'Gene Id':'gene_id', 'Perturbation':'Perturbation_6h_24h', 'pValue':'p.value_6h_24h', 'adj_pValue':"adj.p.value_6h_24h"}, inplace=True)
        mith_perturb_signature_6h_24h['t.value_like_statistic_6h_24h']=mith_perturb_signature_6h_24h['p.value_6h_24h'].apply(lambda p_value : 2001*(1-p_value))
        
        mith_perturb_signature_meta=mith_perturb_signature_meta.merge(mith_perturb_signature_6h_24h, on='gene_id', how='left')
    
        mith_perturb_signature_meta=remove_special_characters(mith_perturb_signature_meta)
        print('drug', drug, 'processed in', np.round(time.time()-drug_start, 1))
        
        print('writing mit3 connectivity input metanalysis files for drug', drug)
        mith_perturb_signature_meta.to_csv(output_filename, sep='\t', index=False)



