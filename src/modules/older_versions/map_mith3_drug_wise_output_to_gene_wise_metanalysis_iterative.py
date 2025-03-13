#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 13 09:50:20 2024

@ author: L-F-S

Maps row names of genes of Mithril 3 input matrix from gene names to gene ids
"""
import os
if not os.getcwd().endswith('modules'):
    os.chdir('modules')
import sys
import numpy as np
import pandas as pd
import time
import pickle
from conf import HOME_DIR, BASE_DIR, MITH_IN_DRUG, MITH_OUT_DRUG, TSR_OUT_DRUG
import collections

#%%#%%
print('loading genes list')
ex_mith=pd.read_csv(MITH_OUT_DRUG+'oligomycin-a_24h.perturbations.txt', sep='\t')

#print(ex_mith.columns)
#
##'# Pathway Id', 'Pathway Name', 'Gene Id', 'Gene Name', 'Perturbation', 'Accumulator', 'pValue'
#print(ex_mith.shape)
#print(len(ex_mith['Gene Id'].unique()))
# genes appear several times , let's check that they always have the same perturbation value

def remove_special_characters(df):
    special_chars=';/:*?\"|'
    for char in special_chars:
        indexes=df[df['Gene Name'].str.contains('\\'+char)].index
        if len(indexes)>0:
            df.drop(index=indexes, inplace=True)
    return df

ex_mith=remove_special_characters(ex_mith)

#print('after removing special chars')
#print(ex_mith.shape)
#print(len(ex_mith['Gene Id'].unique()))

# printing  stuff to see wat's inside ex_mith
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
print('loading drugs list')

drugs_list=os.listdir(MITH_OUT_DRUG)

drugs_list=[x.split('.')[0] for x in drugs_list]
drugs_list=[x.split('_')[0] for x in drugs_list]
drugs_list=list(np.unique(drugs_list))
tot_drugs=len(drugs_list)

        
def filter_tsr_genes(genes_list):
    print('se vuoi usare solo i tsr genes, e filtrare i geni/metaboliti in piu trovati da mithril, fallo in questa funzione')
    
# 53     i1=400#sys.argv[1]                                                                                                                                                              
# 54     i2=400+int(len(genes_list)/30)#sys.argv[2]                                                                                                                                      
# 55     return genes_list[i1:i2]   
    if len(sys.argv)>1:
        i1=int(sys.argv[1])
        i2=int(sys.argv[2])
    else:
        i1=0
        i2=10
    return genes_list[i1:i2]

def name_from_id(gene):
    unique_genes_list=ex_mith[ex_mith['Gene Id']==gene]['Gene Name'].unique()
    if len(unique_genes_list)==1:
        return unique_genes_list[0]
    raise ValueError('More than one gene name for gene id', gene,':',unique_genes_list)

genes_list=np.sort(list(ex_mith['Gene Id'].unique()))




genes_list=filter_tsr_genes(genes_list)
tot_genes=len(genes_list)
print('current tot genes', tot_genes)
#%%initialize dictionary of lists to fill for every gene:
df_of_gene=collections.defaultdict(list)
columns=["gene","drug" ,"DE_log2_FC_6h","std.error_6h" ,"t.value_6h",
         "p.value_6h" , "adj.p.value_6h" , "DE_log2_FC_24h", "std.error_24h", "t.value_24h",\
         "p.value_24h"  , "adj.p.value_24h", "DE_log2_FC_6h_24h", "std.error_6h_24h", "t.value_6h_24h","p.value_6h_24h"]

#i1=0#sys.argv[1]
#i2=10#sys.argv[2]
print('CREATING gene wise dataframes for all genes all drugs')
start=time.time()
for h, drug in enumerate(drugs_list):
    drug_start=time.time()
    print(drug, h+1, 'of', tot_drugs)
    #6h
    mith_perturb_signature_file_6h=drug+'_6h.perturbations.txt'
    mith_perturb_signature_6h = pd.read_csv(MITH_OUT_DRUG+mith_perturb_signature_file_6h, sep='\t', index_col=False, engine='python')
    #24h
    mith_perturb_signature_file_24h=drug+'_24h.perturbations.txt'
    mith_perturb_signature_24h = pd.read_csv(MITH_OUT_DRUG+mith_perturb_signature_file_24h, sep='\t', index_col=False, engine='python')
    #6h 24h
    mith_perturb_signature_file_6h_24h=drug+'_6h_24h.perturbations.txt'
    mith_perturb_signature_6h_24h = pd.read_csv(MITH_OUT_DRUG+mith_perturb_signature_file_6h_24h, sep='\t', index_col=False, engine='python')

    start_gene=time.time()
    for gi, gene in enumerate(genes_list):            
        
        # CHECK if file exists
        gene_name=name_from_id(gene)
        gene_wise_filename=TSR_OUT_DRUG+'LINCS_lorenzo/metanalysis_mith3_gene_wise/'+gene_name+'_'+gene+'_metanalysis.txt'
        if not os.path.isfile(gene_wise_filename): 
            
            gene_start=time.time()

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
            
    #        print(df_line)
            df_of_gene[gene].append(df_line)
            
    print('drug', drug, 'processed in', np.round(time.time()-drug_start, 1))

        
#%%
print('writing gene wise mith3 files')

for i, gene in enumerate(genes_list):
    gene_name=name_from_id(gene)
    gene_wise_filename=TSR_OUT_DRUG+'LINCS_lorenzo/metanalysis_mith3_gene_wise/'+gene_name+'_'+gene+'_metanalysis.txt'
    print(gene)
    if not os.path.isfile(gene_wise_filename): 
        print('writing')
        gene_name=name_from_id(gene)
        df= pd.DataFrame(df_of_gene[gene],columns=columns)
        print(gene, gene_name, 'gene', i, 'of',tot_genes, 'elapsed time', np.round(time.time()-start, 1))
        try:
            df.to_csv(gene_wise_filename, sep='\t')
        except:
            print('failed to write file for gene', gene)
            continue
    
