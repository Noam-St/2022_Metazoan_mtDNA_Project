"""
Created on Wed Dec 11 12:29:04 2019
This code finds all possible gene clusters in a df of gene orders and returns the frequency of each gene cluster.

@author: noam1
"""
#%%
import pandas as pd
import os
from ast import literal_eval
import time
from sys import platform
from functools import wraps
path = os.getcwd()
if platform != 'win32':
    dash = '/'
else:
    dash = '\\'

def timer_deco(func):
    """ Timer decorator, adds a timer to function run in wrapper"""
    @wraps (func)
    def timer_deco_wrapper(*args, **kwargs):
        start = time.perf_counter()
        returned = func(*args, **kwargs)
        end = time.perf_counter()
        runtime = end - start
        print(f'Finished running {func.__name__} in {round(runtime,2)} seconds!')
        return returned
    return timer_deco_wrapper

#%%
def is_a_in_x_circ(A, X):
    """ Checks if a list is a sub-componenet of another list,
    considers lists as a round object. Returns True if A is inside X and False if it isnt"""
    for i in range(len(X)):
        if A == X[i:i+len(A)]: return True
        if i+len(A)>len(X):
            if A == X[i::] + X[0:len(A)-len(X[i::])]: return True
    return False
def is_a_in_x_lin(A,X):
    """ Checks if a list is a sub-component of another list,
    considers list as linear object, returns True if A is inside X and False if it isnt"""
    for i in range(len(X)):
        if A == X[i:i+len(A)]: return True
    return False 

@timer_deco
def seq_common_to_all(series,topologies):
    """ Recives a series of gene orders and topologies, checks which sub-sequences of genes repeat them selves
    in over 50% of the organisms and returns a dict of sequence:% of repeats"""
    seq_prec={}
    all_seq={}
    for x in range(0,len(series)): # Iterate over organism's gene order indices
        try: series[x] = literal_eval(series[x])
        except ValueError: pass
        for i in range(len(series[x]),0,-1): # Iterate over the current gene_order list indices in reverse
            for j in range(0,i): # Iterate over the current gene order from 0 to i index
                test_seq = series[x][j:i] # use that slice as current test_seq
                if str(test_seq) in all_seq.keys(): # Remove existing sequences from the test_seq
                    continue
                if len(test_seq)==1: # If the test_seq is a single gene, skip it
                    continue
                t=0
                f=0
                print(f'Test sequence is : {test_seq}')
                topology = topologies[x]
                for seq in series:
                    try:seq=literal_eval(seq)
                    except ValueError: pass
                    if topology == 'circular': # If the topology is circular, use is_a_in_x_circ
                        if (is_a_in_x_circ(test_seq,seq)) or (is_a_in_x_circ([f'-{i}' for i in test_seq],seq)):
                            t+=1
                        else:
                            f+=1
                    else: # If the topology is linear, use is_a_in_x_lin
                        if (is_a_in_x_lin(test_seq,seq)) or (is_a_in_x_lin([f'-{i}' for i in test_seq],seq)):
                            t+=1
                        else:
                            f+=1
                all_seq[str(test_seq)]=round((t/(f+t))*100,2) # Calculate the prevalence of the current test_seq
                if (round((t/(f+t))*100,2)>50) or (len(test_seq)<=3): # Keep only sequences that are more than 50% common
                    seq_prec[str(test_seq)]=round((t/(f+t))*100,2)
    return all_seq


# %%
if __name__ == '__main__':
    # Load the ref dataframe
    all_org_comb = pd.read_csv(os.path.join(os.path.dirname(path), 'for_dmatrix.csv'))
    all_repeats = seq_common_to_all(all_org_comb.Gene_order, all_org_comb.topology)
    all_gene_grps_df = pd.DataFrame.from_dict(all_repeats,orient='index')\
    .reset_index().rename(columns={'index':'Gene_Group',0:'Prevalence'})
    all_gene_grps_df.to_csv('gene_groups_post-fix.csv' ,index=False)
