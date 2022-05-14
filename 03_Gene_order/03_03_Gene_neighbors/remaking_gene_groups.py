# -*- coding: utf-8 -*-
"""
Created on Fri Mar 20 18:00:55 2020

@author: noam1
"""

import pandas as pd
import os,sys
from ast import literal_eval
path = os.getcwd()
from consts import REP_DICT
from utils import Which_tRNA

if sys.platform != 'win32':
    dash = '/'
else:
    dash = '\\'



#%%
for k,v in REP_DICT.items():
    REP_DICT[k] = Which_tRNA(v)

def replace_func(to_rep,by_dict):
    """ recives a str and dict, replaces str with any matching values in dict and returns dict reversed to replace string annotations
    with numerical annotations"""
    counter=0
    bools=False

    if '*' in to_rep: #checks if this gene is a duplicate marked by *
        counter=to_rep.count('*')
    to_rep=to_rep.replace('*','')
    if '-' in to_rep[0]: #checks if the gene is in the complementary strand
        bools=True
        to_rep=to_rep.replace('-','',1)
    for k,v in by_dict.items(): #iterates over the replacement dictionary and looks for matching gene strings
        if to_rep==k:
            to_rep=v
            break
    to_rep=to_rep+'*'*counter
    if bools: to_rep='-'+to_rep #If gene is in complementry, if the answer is yes 
    return to_rep


def replace_func_on_list(listy):
    """ Apply replace function on each item in list"""
    if type(listy)==str:
        listy=literal_eval(listy)
    if type(listy)==list:
        for i in range(len(listy)):
#            if '-' in listy[i]:
#                listy[i]='-' + replace_func(listy[i].replace('-',''),REP_DICT)
#            else:
            listy[i]=replace_func(listy[i],REP_DICT)
    return listy

def is_a_in_x(A, X):
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
      
  

def seq_common_to_all(series,patterns):
    """ Recives a series of gene orders, checks which sub-sequences of genes repeat them selves
    in over 50% of the organisms and returns a dict of sequence:% of repeats"""
    seq_prec={}
    all_seq={}
    for test_seq in patterns:
        t=0
        f=0
        for seq in series:
            if (is_a_in_x(test_seq,seq)) or (is_a_in_x([f'-{i}' for i in test_seq],seq)):
                t+=1
            else:
                f+=1
        all_seq[str(test_seq)]=round((t/(f+t))*100,2)
        if round((t/(f+t))*100,2)>50:
            seq_prec[str(test_seq)]=round((t/(f+t))*100,2)
    return all_seq

if __name__ == '__main__':
    full_taxa=pd.read_csv(f'{path}{dash}Gene_clusters{dash}for_dmatrix.csv',index_col=0)
    all_repeats2=pd.read_csv(f'{path}{dash}Gene_clusters{dash}gene_groups_post-fix.csv')
    all_repeats2_real=all_repeats2.copy()
    all_repeats2_real.Gene_Group=all_repeats2.Gene_Group.apply(lambda x : replace_func_on_list(x))
    full_taxa.Gene_order=[literal_eval(i) for i in full_taxa.Gene_order]
    all_repeats2.Gene_Group=[literal_eval(i) for i in all_repeats2.Gene_Group]
    all_repeats2_real.Gene_Group=[str(i) for i in all_repeats2_real.Gene_Group]
    all_repeats2_real.set_index(keys='Gene_Group',inplace=True)
    types=list(full_taxa['phylum'].unique())
    ## This code calculates gene groups percentages for all types in the types list, and based on a DF with all possible repeats,
    # and returns a dataframe with prevalence for each group  

    for i in types:
        x=seq_common_to_all(full_taxa[full_taxa['phylum']==i].Gene_order,all_repeats2.Gene_Group)
        x_df=pd.DataFrame.from_dict(x,orient='index',columns=['Prevalence']).reset_index().rename(columns={'index':'Gene_Group'})
        x_df.Gene_Group=x_df.Gene_Group.apply(lambda x: replace_func_on_list(x))
        x_df.Gene_Group=[str(i) for i in x_df.Gene_Group]
        x_df.set_index(keys='Gene_Group',inplace=True)
        x_df.rename(columns={'Prevalence':f'{i}  Prevalence\n N='+str(full_taxa.phylum[full_taxa['phylum']==i].count())},inplace=True)
        all_repeats2_real=pd.concat([all_repeats2_real,x_df],axis=1)
        
    #CODE CURRENTLY SET FOR PHYLUM$$$$$$$$$$$$$$$$$$$$$

    all_repeats2_real.to_csv('gene_groups_post-fix_by_phylum.csv')