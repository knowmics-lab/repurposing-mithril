#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 18 20:37:14 2025

@ author: L-F-S
calculate RGES connectivity score between disease signature and drug sinatures.

"""
import os
if not os.getcwd().endswith('modules'):
    os.chdir('modules')
import sys
import numpy as np
import pandas as pd
import time
from scipy import stats
from conf import DISEASE, CS_OUT
from loader import load_disease_signature, load_single_drug_signature
from preprocessing_utils import get_drugs_list

# FUNCTIONS
def get_common_genes(disease_signature, drug_signature):
    '''
    Identifies common genes between drug and disease signature. Returns indexes
    of common genes in both signatures.
    Input:
        - disease_signature: pd.DataFrame with columns ['gene_id', ...]
        - drug_signature: pd.DataFrame with columns ['gene_id', ...]
    Output:
        - Tuple of DataFrame.index of common genes in disease and drug signatures
    
    '''

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
 

def rank_genes(drug_signature, disease_signature, drug_col_name, disease_col_name):
    '''
    Compute the array V of sorted drug gene indexes, indexed by sorted disease
    indexes.
   Input:
       - drug_signature: pd.DataFrame with columns ['gene_id', drug_col_name, ...]
       - disease_signature: pd.DataFrame with columns ['gene_id', disease_col_name, ...]
       - drug_col_name: str, column name to sort drug_signature by
       - disease_col_name: str, column name to sort disease_signature by
    Output:
       NumPy array where the indices are sorted disease gene indexes 
       and the values are indexes of corresponding sorted drug genes.
   '''
    
    # Sort genes in both datasets:
    sorted_drug_genes=drug_signature.sort_values(by=drug_col_name).reset_index(drop=True)
    sorted_disease_genes=disease_signature.sort_values(by=disease_col_name)
    
    # Merge dataset and keep disease index sorting
    merged_data_on_disease_indexes = pd.merge(sorted_disease_genes, sorted_drug_genes.reset_index(), on='gene_id')
    
    # Extract drug index of sorted genes
    V=merged_data_on_disease_indexes['index'].to_numpy() 
    return V

def compute_KS(disease_disregulated_genes, V, drug_genes):
    '''
    Compute Kolmogorov-Smirnov (KS) statistic.
    From Lamb et al., 2006 supplementary
    Input:
       - disease_disregulated_genes: list or array of gene IDs
       - V: NumPy array of sorted drug gene indices
       - drug_genes: list or array of gene IDs in the drug signature
    Output:
       - a: float, maximum positive difference
       - b: float, maximum negative difference
       - s: int, number of disregulated genes
       - r: int, total number of genes in reference drug expression data
    '''
    
    # number of disregulated genes:
    s=len(disease_disregulated_genes) 

    # total number of genes in reference drug expression data:
    r=len(drug_genes)
    
    V_over_r = V/r
    a = np.max(np.arange(1, s+1)/s - V_over_r)
    b = np.max(V_over_r - (np.arange(s)/s))
    
    return a, b, s, r

def compute_KS_evil_twin(disease_disregulated_genes, V, drug_genes):
    '''
    Compute Kolmogorov-Smirnov (KS) statistic.
    From Lamb et al., 2006 supplementary
    Input:
       - disease_disregulated_genes: list or array of gene IDs
       - V: NumPy array of sorted drug gene indices
       - drug_genes: list or array of gene IDs in the drug signature
    Output:
       - a: float, maximum positive difference
       - b: float, maximum negative difference
       - s: int, number of disregulated genes
       - r: int, total number of genes in reference drug expression data
    '''
    
    # number of disregulated genes:
    s=len(disease_disregulated_genes) 

    # total number of genes in reference drug expression data:
    r=len(drug_genes)
    
    V_over_r = V/r
    a = np.max(np.abs(np.arange(1, s+1)/s - V_over_r))
    b = np.max(np.abs(V_over_r - (np.arange(s)/s)))
    ks=np.max(a,b)
    
    return ks, s, r



def lamb_normalize(cs_list):
    '''
    CS normalization following Lamb et al., 2006
    Input:
        cs_list: list, connectivity scores for drugs of interest
    Output:
        list, normalized connectivity scores
    '''
    
    p=max(cs_list)
    q=min(cs_list)
    
    def extr(x):
        if x>0:
            return p
        return -q
    
    return [x/extr(x) for x in cs_list]

def calculate_CS(a_up, a_down, b_up, b_down):
    '''
    Calculate connectivity score.
    from Lamb et al., 2006
    Input:
    - a_up: float, maximum positive difference for upregulated genes
    - a_down: float, maximum positive difference for downregulated genes
    - b_up: float, maximum negative difference for upregulated genes
    - b_down: float, maximum negative difference for downregulated genes
    Output:
        - CS: float, CS value in the interval [-1, 1]
    '''
    
    if a_up>b_up:
        ks_up=a_up
    else:
        ks_up=-b_up
        
    if a_down>b_down:
        ks_down=a_down
    else:
        ks_down=-b_down
    
    if ks_up*ks_down>0:
        return 0
    return ks_up-ks_down   

def calculate_CS_evil_twin(ks_up, ks_down):
    '''
    Calculate connectivity score.
    from Lamb et al., 2006
    Input:
    - a_up: float, maximum positive difference for upregulated genes
    - a_down: float, maximum positive difference for downregulated genes
    - b_up: float, maximum negative difference for upregulated genes
    - b_down: float, maximum negative difference for downregulated genes
    Output:
        - CS: float, CS value in the interval [-1, 1]
    '''
    return ks_up-ks_down 

def calculate_RGES(a_up, a_down, b_up, b_down):
    '''
    Calculate Reverse Gene Expression Score (RGES).
    from Bin Chen et al., 2017
    Input:
    - a_up: float, maximum positive difference for upregulated genes
    - a_down: float, maximum positive difference for downregulated genes
    - b_up: float, maximum negative difference for upregulated genes
    - b_down: float, maximum negative difference for downregulated genes
    Output:
        - RGES: float, RGES value in the interval [-2, 2]
    '''
    
    if a_up>b_up:
        ks_up=a_up
    else:
        ks_up=-b_up
        
    if a_down>b_down:
        ks_down=a_down
    else:
        ks_down=-b_down
    
    return ks_up-ks_down   

def montecarlo_connectivity(s_up, s_down, r, n_iterations=1000, score_type='bin_chen'):
    '''
    Randomly sample RGES
    Input:
        - s_up: int, number of upregulated genes
        - s_down: int, number of downregulated genes
        - r: int, total number of genes in reference drug expression data
        - n_iterations: int, number of Monte Carlo iterations (default: 1000)
        - score_type: str, optional, which calculation to use. 
            options: ['bin_chen', 'sirota', 'lamb'], default: 'bin_chen'
    Output:
        - list of float, sampled RGES values
    '''
    
    random_RGES_list=[]
    drug_genes=np.arange(r)
    
    for i in range(n_iterations):
        
        # Generate random indexes for up and downregulated genes
        random_up_and_down_indexes=np.random.choice(r, s_up + s_down, replace=False)
         
        # Create random index dictionary:
        random_V_up=random_up_and_down_indexes[:s_up]
        random_V_down=random_up_and_down_indexes[-s_down:]
        
        # Compute random KS stats:
        random_a_up, random_b_up, _, _ = compute_KS(random_up_and_down_indexes[:s_up], random_V_up, drug_genes)
        random_a_down , random_b_down, _, _ = compute_KS(random_up_and_down_indexes[-s_down:], random_V_down, drug_genes)

        # calculate RGES evil twin:
        if score_type=='evil_twin':
            ks_up,  _, _ = compute_KS(random_up_and_down_indexes[:s_up], random_V_up, drug_genes)
            ks_down, _, _ = compute_KS(random_up_and_down_indexes[-s_down:], random_V_down, drug_genes)
            random_RGES_list.append(calculate_CS_evil_twin(ks_up, ks_down, random_b_up, random_b_down))
        
        # Calculate random RGES:
        if score_type=='bin_chen':
            random_RGES_list.append(calculate_RGES(random_a_up, random_a_down, random_b_up, random_b_down))
        if (score_type=='lamb') or(score_type=='sirota') :
            random_RGES_list.append(calculate_CS(random_a_up, random_a_down, random_b_up, random_b_down))

    if score_type=='lamb':
        random_RGES_list=lamb_normalize(random_RGES_list)

    return random_RGES_list

def bin_chen_connectivity(disease_signature, drug_signature, rank_on='p_value'):
    '''calculates the Reverse Gene Expression Score (RGES), a
   connectivity score, as defined in Bin, Chen, 2017
    Input:
        - disease_signature: pd.DataFrame with columns ['gene_id', signature data, p value]
        - drug_signature: pd.DataFrame with columns ['gene_id', signature data, p value]
    Output:
        -   measured_RGES: float, calculated RGES value
        -   p_value: float, p-value from Monte Carlo simulation
    '''
    
    disease_p_val_col_name=disease_signature.columns[2]
    disease_signature_col_name=disease_signature.columns[1]
    
    if rank_on=='p_value':
        ranking_col_name=drug_signature.columns[2]
    if rank_on=='magnitude':
        ranking_col_name=drug_signature.columns[1]
    
    # # Get lists of up (down) regulated genes: older with 2vs
    disease_signature_up = disease_signature[disease_signature[disease_signature_col_name]>0]
    disease_signature_down = disease_signature[disease_signature[disease_signature_col_name]<0]
    
    # Calculate rank map V for up and down regulated genes:
    V_up = rank_genes(drug_signature, disease_signature_up  , ranking_col_name, disease_p_val_col_name)
    V_down = rank_genes(drug_signature, disease_signature_down, ranking_col_name, disease_p_val_col_name)
    
    # Compute KS statistic for up and down regulated genes:
    a_up, b_up, s_up, r = compute_KS(disease_signature_up['gene_id'], V_up, drug_signature['gene_id'])
    a_down , b_down, s_down, r = compute_KS(disease_signature_down['gene_id'], V_down, drug_signature['gene_id'])
    
    # Compute RGES:
    measured_RGES = calculate_RGES(a_up, a_down, b_up, b_down)
    
    # Calculate two tailed p-value (no assumption on the direction of RGES 
    # between disease and drug) for measured RGES, using random sampling:
    n_iterations=1000 
    random_RGES_list = montecarlo_connectivity(s_up, s_down, r, n_iterations, score_type='bin_chen')
    p_value=np.sum(np.abs(np.array(random_RGES_list))>np.abs(measured_RGES))/n_iterations
    return measured_RGES, p_value


#%% 
# if __name__=='__main__':
# ###############################################################################
# #             CALCULATE CONNECTIVITY FOR DEG/MITH DATA
# ###############################################################################
#     start=time.time()
#     print(    '''calculates connectivity score between a disease and all LINCS database drugs, for given drugs
#         and given drug perturbation times. 
#         Output: a dataframe of connectivity scores and  pvalues, for all combinations of conditions''')
    
#     # PARAMETERS:
        
#     mith=bool(int(sys.argv[3]))#False #loads DEG data
    
#     pert_times=['6h','24h','6h_24h'] #which drug experimental perturbation times (not to be confused  with mithril's perturbation value) to use

#     singature_uom = 'DE_log2_FC' if not mith else 'Perturbation'
    

#     print('mith tag:', mith, '\ndisease:', DISEASE)

#     # load disease signature:
#     disease_signature=load_disease_signature(DISEASE, mith=mith)
#     disease_signature=disease_signature[['gene_id', singature_uom, 'adj.p.value']]
    
    
#     # Discard disease genes that are 
#     # not found in drug genes
#     # DO NOT simply select common genes, as this would impact
#     # connectrivity calculations
#         #load one drug signature (all durg signatures have the same genes in the same order)
#     drug_signature=load_single_drug_signature('ibuprofen', mith=mith)
    
#     disease_common_index, drug_common_index = get_common_genes(disease_signature, drug_signature)  
#     # Note: indexing every time with pandas .iloc is slower than writing
#     # and loading binaries for filtered dataframes, for each drug
#     # however this is not relevant, as we only need to filter
#     # disease data, while we can take all drug data.
#     # Question for later: with mithril how does metabolite data impact the result?
    
#         # filter disease signature with common genes:
#     disease_signature=disease_signature.iloc[disease_common_index].reset_index(drop=True)
#     print('common disease genes', len(disease_signature))
    
#         # Optional: filter for 150 most significant:
#             # n_genes=150
#         #TODO disease_signature=filter_most_significant(disease_signature, n_genes)
#         # e volendo aggiungici il nuovo drug_common_index se voglio fare i pearsons vari..
#     # load all drugs list:
#     drugs_list=get_drugs_list()
    
#     # Initialize dataframe:
#     data = []
    
#     #get drug indexes for parallelization:
#     i1=int(sys.argv[1])
#     i2=int(sys.argv[2])
    
#     #%% Calculate connectivity score between disease and drugs:     
#     for drug in drugs_list[i1:i2]:
#         print(drug)
        
#         # load drug signature
        
#         drug_signature=load_single_drug_signature(drug, mith=mith)
        
#         for pert_time in pert_times:
#             print(pert_time)
            
#             # Calculate connectivity score between drug and disease:
            
#             # select columns of interest:
            
#             p_val_str= 'adj.p.value' if not pert_time=='6h_24h' else 'p.value'
#             columns_of_interest=['gene_id', singature_uom+'_'+pert_time, p_val_str+'_'+pert_time]
#             rank_on=sys.argv[4]
#             cscore, p_value = bin_chen_connectivity(disease_signature, drug_signature[columns_of_interest], rank_on=rank_on)
            
#             # Calculate other correlations BETWEEN COMMON GENES
#             disease_signature_signature=disease_signature[singature_uom]
#             drug_signature_signature=drug_signature[singature_uom+'_'+pert_time].iloc[drug_common_index]
#             pearson=stats.pearsonr(disease_signature_signature, drug_signature_signature)
#             spearman=stats.spearmanr(disease_signature_signature, drug_signature_signature)
#             cosine_sim=np.dot(np.array(disease_signature_signature),np.array(drug_signature_signature))/(np.linalg.norm(np.array(disease_signature_signature))*np.linalg.norm(np.array(drug_signature_signature)))
            
#             data.append([DISEASE, drug, pert_time, cscore, p_value, pearson[0], pearson[1], spearman[0], spearman[1], cosine_sim]) # other correlations, genes subset
                
    
#     connectivity_data= pd.DataFrame(data, columns=["disease","drug","perturbation_time","connectivity_score","cs_p_value",'pearson','pearson_p_value','spearman','spearman_p_value','cos_sim']) #add other correlations, genes subset
#     print(connectivity_data)
    
#     # save connectivity scores
#     connectivity_dataset_filename=CS_OUT+str(i1)+'_'+str(i2)+'_DEG_connectivity_score.tsv' if not mith else  CS_OUT+str(i1)+'_'+str(i2)+'_mith_connectivity_score.tsv'

#     connectivity_data.to_csv(connectivity_dataset_filename, sep='\t', index=False)

#     print('total elapsed time for', len(drugs_list[i1:i2]),' drugs: ', time.time()-start)

#%%     CALCULATE CONNECTIVITY FOR ONE DRUG DATA

if __name__=='__main__':
    mith=True
    pert_time='6h_24h' 
    drug='ibuprofen'
    singature_uom = 'DE_log2_FC' if not mith else 'Perturbation'
    # load disease signature:
    disease_signature=load_disease_signature(DISEASE, mith=mith)
    disease_signature=disease_signature[['gene_id', singature_uom, 'adj.p.value']]
    drug_signature=load_single_drug_signature(drug, mith=mith)
    
    disease_common_index, drug_common_index = get_common_genes(disease_signature, drug_signature)  
    
    # filter disease signature with common genes:
    disease_signature=disease_signature.iloc[disease_common_index].reset_index(drop=True)
    print('common disease genes', len(disease_signature))
    # load drug signature
    p_val_str= 'adj.p.value' if not pert_time=='6h_24h' else 'p.value'
    columns_of_interest=['gene_id', singature_uom+'_'+pert_time, p_val_str+'_'+pert_time]
    cscore, p_value = bin_chen_connectivity(disease_signature, drug_signature[columns_of_interest], rank_on='magnitude')
    print(cscore, p_value)
