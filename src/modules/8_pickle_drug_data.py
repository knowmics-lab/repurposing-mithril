#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  7 09:31:28 2025

@ author: L-F-S

convert drug metanalyss data into pickled  binaryformat for faster loading 

"""
import os
if not os.getcwd().endswith('modules'):
    os.chdir('modules')
import sys
import numpy as np
import pandas as pd
import pickle
from loader import load_single_drug_signature
from preprocessing_utils import get_drugs_list
from conf import TSR_OUT_DRUG
drugs_list=get_drugs_list()
i1=int(sys.argv[1])
i2=int(sys.argv[2])
for n, drug in enumerate(drugs_list[i1:i2]):
    print(drug)
    deg_filename=TSR_OUT_DRUG+'LINCS/metanalysis_drug_wise_filtered/'+drug+'_metanalysis.csv'
    deg_pickled_filename=TSR_OUT_DRUG+'LINCS/metanalysis_drug_wise_filtered/'+drug+'_metanalysis.pkl'
    
    data=load_single_drug_signature(drug, mith=False, pkl=False)
    with open(deg_pickled_filename, 'wb') as f:
        pickle.dump(data, f)
        
    mith_filename=TSR_OUT_DRUG+'LINCS/metanalysis_mith3_drug_wise/'+drug+'_metanalysis.csv'
    mith_pickled_filename=TSR_OUT_DRUG+'LINCS/metanalysis_mith3_drug_wise/'+drug+'_metanalysis.pkl'
    data=load_single_drug_signature(drug, mith=True, pkl=False)
    with open(mith_pickled_filename, 'wb') as f:
        pickle.dump(data, f)

    
    

