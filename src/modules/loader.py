#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 29 10:24:39 2025

@ author: L-F-S
"""
import os
if not os.getcwd().endswith('modules'):
    os.chdir('modules')
import pandas as pd
import pickle
from conf import TSR_OUT_DRUG, TSR_OUT_DISEASE, TSR_OUT_CSCORE
import numpy as np
from preprocessing_utils import get_drugs_list

def load_disease_signature(DISEASE, mith=False):
    '''
    inputs:
        DISEASE: str disease symbol (directory name)

        mith: bool
            def. False: loading DEG FC signature. if True, loading MIthRIL 
            perturbation signature.

    returns a pandas.Dataframe where rows are genes.
    '''
    
    if not mith:
        filename=TSR_OUT_DISEASE+DISEASE+'/'+DISEASE+'_signature_gene_id.csv'
        return pd.read_csv(filename, sep=';',decimal=',', dtype={'gene_id':'str'})
   
    else:
        filename=TSR_OUT_DISEASE+DISEASE+'/'+DISEASE+'_mith3_signature.csv'
        return pd.read_csv(filename, sep=';')

def load_single_drug_signature(drug, mith=False, pkl=True):
    '''
    Loads DEG drug wise drug signature (with unique gene id rows)

    inputs:
        DISEASE: str disease symbol (directory name)

        mith: bool
            def. False: loading DEG FC signature. if True, loading MIthRIL 
            perturbation signature.

    returns a pandas.Dataframe where rows are genes.
    '''
    # print('pikle',pkl)
    if not mith:
        
        filename=TSR_OUT_DRUG+'LINCS/metanalysis_drug_wise_filtered/'+drug+'_metanalysis'
        
        if not pkl:
            return pd.read_csv(filename+'.csv', sep='\t',  dtype={'gene_id':'str'})
        else:
            # print('comepossibile',filename)
            with open(filename+'.pkl', 'rb') as f:
                data=pickle.load(f)
            # print('\nAOOOOOOOOOOOOOOOOOOOOOOOO',data,'wewew')
            return data
                
    else:
        
        filename=TSR_OUT_DRUG+'LINCS/metanalysis_mith3_drug_wise/'+drug+'_metanalysis'
        
        if not pkl:
            return pd.read_csv(filename+'.csv', sep='\t')
        else:
            with open(filename+'.pkl', 'rb') as f:
                data=pickle.load(f)
            return data
            

def load_drug_signatures(mith=False, pkl=True):
    '''
    inputs:

        mith: bool
            def. False: loading DEG FC signature. if True, loading MIthRIL 
            perturbation signature.
    returns a pandas.Dataframe which is a concatenation by row
    of all drug signature data. rows are genes signature for given drug 
    (genes will appear multiple time, once for every drug)
    '''
    
    print('loading all drug signatures...')
    
    drugs_list=get_drugs_list()
    
    if not mith:
        DEG_drug_list=[]
        for n, drug_file in enumerate(drugs_list): 
            drug=drug_file.split('_')[0]
            DEG_drug_list.append(load_single_drug_signature(drug, mith, pkl))
    
    else:
        DEG_drug_list=[]
        for n, drug_file in enumerate(drugs_list):
            drug=drug_file.split('_')[0]
            DEG_drug_list.append(load_single_drug_signature(drug, mith, pkl))
    
    print(n, 'drugs loaded')
    return pd.concat(DEG_drug_list)
