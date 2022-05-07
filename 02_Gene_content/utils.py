from consts import TRNA_DICT, REPL_DICT
import re
import os
from Bio import Seq
import pandas as pd
PATH = os.getcwd()

def count_appearances(df, ph, level = 'phylum', col = 'Gene_set'):
    """
    Count the amount of times each gene content set appears in a certain phylum
    Return a dictionary where keys are frozensets and the values are the percentage of organisms within a given phylum with
    that gene content set.
    """
    set_dict = {}
    if ph == 'all':
        chosen_ph = df
    else:
        chosen_ph = df.loc[df[level] == ph, col]
    size = len(chosen_ph)
    for s in chosen_ph:
        if s in set_dict.keys():
            set_dict[s] += 1
        else:
            set_dict[s] = 1
    set_dict = {k : (v/size) * 100 for k,v in sorted(set_dict.items(), key = lambda x : x[1], reverse = True)}
    return set_dict

def Which_tRNA(name):
    """
    This code recieves a gene name, tries to match it with common tRNA motifs, and then with any matching values in the replacement dict, if nothing matches, returns None
    which is later removed
    """
    if name in REPL_DICT.keys(): return name 
    c=0
    name=str(name)
    if 'transfer' in name:
        name = name.replace('transfer ','t')
    codon = re.search(r'\(\w{3}\)',name) # Grab codon if it exists
    if codon: codon = codon.group(0)
    else: codon=''
    name = re.sub(pattern=r'\(\w{3}\)',repl = '',string = name)
    if name in TRNA_DICT.keys():  # if the tRNA is in tRNA-X format, replace with my format and return
        name = TRNA_DICT[name]
        return name+codon

    #Following lines try to match with each possible tRNA.
    if re.match(r't[r,R]\D{1,2}(H|His)',name):
        name='tRNA-His'
        
    if re.match(r't[r,R]\D{1,2}(K|Lys)',name):
        name='tRNA-Lys'
        
    if re.match(r't[r,R]\D{1,2}(R|Arg)',name):
        name='tRNA-Arg'
        
    if re.match(r't[r,R]\D{1,2}(D|Asp)',name):
        name='tRNA-Asp'
        
    if re.match(r't[r,R]\D{1,2}(E|Glu)',name):
        name='tRNA-Glu'
        
    if re.match(r't[r,R]\D{1,2}S',name):
        name='tRNA-Ser'
        
    if re.match(r't[r,R]\D{1,2}T',name):
        name='tRNA-Thr'
        
    if re.match(r't[r,R]\D{1,2}(N|Asn)',name):
        name='tRNA-Asn'
        
    if re.match(r't[r,R]\D{1,2}(Q|Gln)',name):
        name='tRNA-Gln'
        
    if re.match(r't[r,R]\D{1,2}V',name):
        name='tRNA-Val'
        
    if re.match(r't[r,R]\D{1,2}L',name):
        name='tRNA-Leu'
        
    if re.match(r't[r,R]\D{1,2}I',name):
        name='tRNA-Ile'
        
    if re.match(r't[r,R]\D{1,2}M',name):
        name='tRNA-Met'
        
    if re.match(r't[r,R]\D{1,2}(F|Phe)',name):
        name='tRNA-Phe'
        
    if re.match(r't[r,R]\D{1,2}(Y|Tyr)',name):
        name='tRNA-Tyr'
        
    if re.match(r't[r,R]\D{1,2}(W|Trp)',name):
        name='tRNA-Trp'
        
    if re.match(r't[r,R]\D{1,2}P',name):
        name='tRNA-Pro'
        
    if re.match(r't[r,R]\D{1,2}G',name):
        name='tRNA-Gly'
        
    if re.match(r't[r,R]\D{1,2}C',name):
        name='tRNA-Cys'
    
    if re.match(r't[r,R]\D{1,2}A',name):
        name='tRNA-Ala'
        
    for k,v in REPL_DICT.items(): # if not tRNA, this loop matches according to replacement dict
        if name in v:
            name=k
            return name
        if c==len(REPL_DICT)-1:
            name=None
            return name
        c=+1
    try:
        if name in TRNA_DICT.keys(): # if matched with any of the if's above (for tRNAs), replace with my annotation accordingly.
            name = TRNA_DICT[name]
            return name+codon
    except KeyError:
        pass
    print(name)
    if 'RNA' in name or 'rna' in name:
        return None
    elif 'orf' in name or 'ORF' in name:
        return 'ORFX'
    else:
        return None

def reverse_complement(seq):
    """
    """
    return str(Seq.Seq(seq).reverse_complement())

def row_iter_to_df(df, func, *args, **kwargs):
    """
    Iterates over rows of a dataframe and applies a function to each row.
    Returns a dataframe with the results.

    Parameters
    ----------
    df : pandas.DataFrame
        Dataframe to iterate over.
    func : function
        Function to apply to each row.
    *args : list
        Arguments to pass to func.
    **kwargs : dict 
        Keyword arguments to pass to func.
    
    Returns
    -------
    pandas.DataFrame
        Dataframe with results of func applied to each row.
    
    """
    total_df = pd.DataFrame()
    for _, row in df.iterrows():
        temp = pd.DataFrame(func(row, *args, **kwargs))
        if type(temp) != pd.DataFrame: continue
        if total_df.size == 0: total_df = temp
        else: total_df = total_df.append(temp)
    return total_df

def move_legend(ax, new_loc, **kws):
    old_legend = ax.legend_
    handles = old_legend.legendHandles
    labels = [t.get_text() for t in old_legend.get_texts()]
    title = old_legend.get_title().get_text()
    ax.legend(handles, labels, loc=new_loc, title=title, **kws)
    