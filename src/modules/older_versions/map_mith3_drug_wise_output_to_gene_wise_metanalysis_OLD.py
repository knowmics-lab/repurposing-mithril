#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 13 09:50:20 2024

@ author: L-F-S

Maps row names of genes of Mithril 3 input matrix from gene names to gene ids
2: parallelize over genes because otherwise it takes 7 hours
"""
import os
if not os.getcwd().endswith('modules'):
    os.chdir('modules')
import sys
import numpy as np
import pandas as pd
import pickle
from conf import HOME_DIR, BASE_DIR, MITH_IN_DRUG, MITH_OUT_DRUG, TSR_OUT_DRUG
import collections
import time

#%%#%%
ex_mith=pd.read_csv(MITH_OUT_DRUG+'oligomycin-a_24h.perturbations.txt', sep='\t')
#%%
#print(ex_mith.columns)
#
##'# Pathway Id', 'Pathway Name', 'Gene Id', 'Gene Name', 'Perturbation', 'Accumulator', 'pValue'
#print(ex_mith.shape)
#print(len(ex_mith['Gene Id'].unique()))
## genes appear several times , let's check that they always have the same perturbation value
#
#for gene in ex_mith['Gene Id'].unique():
#    print(gene, len(ex_mith[ex_mith['Gene Id']==gene]), 'unique pert', \
#          len(ex_mith[ex_mith['Gene Id']==gene]['Perturbation'].unique()),\
#          'unique pval', len(ex_mith[ex_mith['Gene Id']==gene]['pValue'].unique()),\
#          'unique acc', len(ex_mith[ex_mith['Gene Id']==gene]['Accumulator'].unique()),\
#          'unique path', len(ex_mith[ex_mith['Gene Id']==gene]['# Pathway Id'].unique()))
#    break

# so this just prints the PATHWAYS in which a given gene appears!
          # important 8write in log!!|!
# every line is A GENE IN A PATHWAY, but its value of acc,pval and pert is the same in every pathway.
          # but other question
#%%

drugs_list=os.listdir(MITH_OUT_DRUG)

drugs_list=[x.split('.')[0] for x in drugs_list]
drugs_list=[x.split('_')[0] for x in drugs_list]
drugs_list=list(np.unique(drugs_list))
tot_drugs=len(drugs_list)
def filter_tsr_genes(genes_list):
    print('TODO decidi se usare tutti  i geni di mithril o tornare ai tsr only??')
    return genes_list

genes_list=ex_mith['Gene Id'].unique()
genes_list=filter_tsr_genes(genes_list)
tot_genes=len(genes_list)

def name_from_id(gene):
    unique_genes_list=ex_mith[ex_mith['Gene Id']==gene]['Gene Name'].unique()
    if len(unique_genes_list)==1:
        return unique_genes_list[0]
    raise ValueError('More than one gene name for gene id', gene,':',unique_genes_list)

#%%initialize dictionary of lists to fill for every gene:
columns=["gene","drug" ,"DE_log2_FC_6h","std.error_6h" ,"t.value_6h",
         "p.value_6h" , "adj.p.value_6h" , "DE_log2_FC_24h", "std.error_24h", "t.value_24h",\
         "p.value_24h"  , "adj.p.value_24h", "DE_log2_FC_6h_24h", "std.error_6h_24h", "t.value_6h_24h","p.value_6h_24h"]

i1=int(sys.argv[1])
i2=int(sys.argv[2])
for i, gene in enumerate(genes_list[i1:i2]):
    start=time.time()
    gene_name=name_from_id(gene)
    gene_wise_filename=TSR_OUT_DRUG+'LINCS_lorenzo/metanalysis_mith3_gene_wise/'+gene_name+'_'+gene+'_metanalysis.txt'
    print('Processing', gene, gene_name, 'gene', i, 'of', len(genes_list[i1:i2]), 'elapsed time', np.round(time.time()-start, 1))

    if not os.path.isfile(gene_wise_filename):
        
        print('writing file for gene', gene, gene_name, 'gene', i, 'of', len(genes_list[i1:i2]))
        df=[]
        for h, drug in enumerate(drugs_list):
            print('gene', gene, drug, h+1, 'of', tot_drugs, 'elapsed time', np.round(time.time()-start, 1))
            #6h
            mith_perturb_signature_file_6h=drug+'_6h.perturbations.txt'
            mith_perturb_signature_6h = pd.read_csv(MITH_OUT_DRUG+mith_perturb_signature_file_6h, sep='\t', index_col=False, engine='python')
            #24h
            mith_perturb_signature_file_24h=drug+'_24h.perturbations.txt'
            mith_perturb_signature_24h = pd.read_csv(MITH_OUT_DRUG+mith_perturb_signature_file_24h, sep='\t', index_col=False, engine='python')
            #6h 24h
            mith_perturb_signature_file_6h_24h=drug+'_6h_24h.perturbations.txt'
            mith_perturb_signature_6h_24h = pd.read_csv(MITH_OUT_DRUG+mith_perturb_signature_file_6h_24h, sep='\t', index_col=False, engine='python')
    
            #    if len(mith_perturb_signature['Gene Id'].unique()) == len(genes_list): #doublecheck they all have the same genes per pathway list       
            df_line = [gene, drug, mith_perturb_signature_6h[mith_perturb_signature_6h['Gene Id']==gene]['Perturbation'].unique()[0]\
                       ,'NA', 'NA',\
                       mith_perturb_signature_6h[mith_perturb_signature_6h['Gene Id']==gene]['pValue'].unique()[0], \
                       mith_perturb_signature_6h[mith_perturb_signature_6h['Gene Id']==gene]['adj_pValue'].unique()[0] ]                    
            df_line += [mith_perturb_signature_24h[mith_perturb_signature_24h['Gene Id']==gene]['Perturbation'].unique()[0],\
                        'NA','NA',\
                        mith_perturb_signature_24h[mith_perturb_signature_24h['Gene Id']==gene]['pValue'].unique()[0],\
                        mith_perturb_signature_6h[mith_perturb_signature_24h['Gene Id']==gene]['adj_pValue'].unique()[0]]
            df_line += [mith_perturb_signature_6h_24h[mith_perturb_signature_6h_24h['Gene Id']==gene]['Perturbation'].unique()[0],\
                        'NA','NA',\
                        mith_perturb_signature_6h_24h[mith_perturb_signature_6h_24h['Gene Id']==gene]['pValue'].unique()[0]]
            
            #  print(df_line)
            df.append(df_line)
            
        df= pd.DataFrame(df,columns=columns)
        print('DONE', gene, gene_name, 'gene', i, 'of', len(genes_list[i1:i2]), 'in ', np.round(time.time()-start, 1), 'seconds')
        df.to_csv(gene_wise_filename, sep='\t')
    
