"""
Created on Wed Nov 13 09:55:46 2024

@ author: L-F-S
"""
import os
import sys
import numpy as np
import pandas as pd
from conf import BASE_DIR, MITH_OUT_DISEASE, MITH_OUT, TSR_OUT_DISEASE, mith_to_tsr, TSR_OUT_DRUG, TSR_OUT   
#%% Disease mithril2
# Open MITHRIL condition OutPut
condition_name='ALS_NYGC_mock_FC'#sys.argv[1]#
mith_file_condition=MITH_OUT_DISEASE+condition_name+'.PF.Acc.Pval'
mith_data_condition=pd.read_csv(mith_file_condition, header=0, sep='\t')

# 
condition_pathway='als_NYGC_lorenzo/'# sys.argv[2]# 

#mith_file_LINCS

# Convert toutput as input for connectivity score calculatioon
mith_data_condition[['Gene Id','Gene Name','Perturbation']]

# should look like tsr/output/disease_signature/als_NYGC, which looks like:
tsr_output_name=condition_name+'mithril_signature.csv'
mith_to_tsr(mith_data_condition, TSR_OUT_DISEASE+condition_pathway, tsr_output_name)

#%%drug MITHRIL 3
condition_name='zuclopenthixol_24h'

mith_file=MITH_OUT+condition_name+'perturbations.txt'
mith_data=pd.read_csv(mith_file_condition, header=0, sep='\t')
