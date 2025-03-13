#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 21 15:55:25 2024

@ author: L-F-S
"""
import os
import subprocess
from conf import HOME_DIR, BASE_DIR, MITH_IN_DRUG, MITH_OUT_DRUG
#%%
def run_MITHRIL(mith_drug_input_filenames,i,j):
    for x in range(i,j):
        xth_drug_signature_name=mith_drug_input_filenames[x]
        drug_signature_name=xth_drug_signature_name.rstrip('.mi')
        print('running mithril for', drug_signature_name)  
        
        # launch MITHrIL (if not already launched)
#        if not os.path.isfile(MITH_OUT_DRUG+drug_signature_name+'.mo'):
#            command=['java', '-jar', 'MITHrIL2.jar', 'mithril',\
#                     '-o', MITH_OUT_DRUG+drug_signature_name+'.mo',\
#                     '-p', MITH_OUT_DRUG+drug_signature_name+'.PF.Acc.Pval',\
#                     '-i', MITH_IN_DRUG+ith_drug_signature_name,\
#                     '-m', '-verbose']
#    
#            subprocess.Popen(command)
#        
#            ran_drugs.append(ith_drug_signature_name)

def get_indexes(x, n_cores, n_drugs):
    #n_drugs=3222*3
    chunk_size=int(n_drugs/n_cores)
    if x ==0:
        print(0, chunk_size)
        return 0, chunk_size
    print( (x-1)*chunk_size,x*chunk_size)
    return (x-1)*chunk_size,x*chunk_size
    

#%%

#make sure you make a folder with ONLY drug based inputs:  MITH_IN/DRUG_SIGNATURE/
mith_drug_input_filenames=os.listdir(MITH_IN_DRUG)

os.chdir(HOME_DIR)

ran_drugs=[]
n_cores=30

for x in range(n_cores):
    if x!=0:
        i,j = get_indexes(x, n_cores, len(mith_drug_input_filenames))
        run_MITHRIL(mith_drug_input_filenames,i,j, )

#%%
f=open(MITH_OUT_DRUG+'drugs_already_ran_with_mithril_through_parallel_wrapper_launch_multiple_mithril_drug.txt','a')
f.write('\n'.join(ran_drugs))
f.close()
    