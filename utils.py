import pandas as pd
from Bio import Seq
REP_DICT={
          '1':'tRNA-Phe',
          '17':'tRNA-His',
          '14':'tRNA-Lys',
          '16':'tRNA-Arg',
          '13':'tRNA-Asp',
          '19':'tRNA-Glu',
          '12':'tRNA-Ser',
          '20':'tRNA-Thr',
          '9':'tRNA-Asn',
          '5':'tRNA-Gln',
          '8':'tRNA-Ala',
          '2':'tRNA-Val',
          '18':'tRNA-Leu',
          '11':'tRNA-Tyr',
          '7':'tRNA-Trp',
          '21':'tRNA-Pro',
          '15':'tRNA-Gly',
          '10':'tRNA-Cys',
          '6':'tRNA-Met',
          '4':'tRNA-Ile',
          '22':'12S ribosomal RNA',
          '23':'16S ribosomal RNA',
          '24':'ATP6',
          '25':'ATP8',
          '26':'COX1',
          '27':'COX2',
          '28':'COX3',
          '29':'CYTB',
          '30':'ND1',
          '31':'ND2',
          '32':'ND3',
          '33':'ND4',
          '34':'ND4L',
          '35':'ND5',
          '36':'ND6',
          '37':'ATP9',
          '38':'mutS',
          '39':'heg',
          '40':'secY',
          '41':'Reph',
          '42':'ORFX',
          '43':'RNAX'}
        
def gene_symbol_to_number(to_rep,number2symbol = False):
    """ recives a str and dict, replaces str with any matching values in dict and returns dict reversed to replace string annotations
    with numerical annotations
    
    Parameters
    ----------
    to_rep : str
        string to replace
    number2symbol : bool
        if true, returns dict reversed
    
    Returns
    -------
    str
        string with replaced values
    dict
        dict with replaced values
    """
    counter=0
    bools=False
    if number2symbol:
        temp_rep = {v:k for k,v in REP_DICT.items()}
    else:
        temp_rep = REP_DICT
    if '*' in to_rep: #checks if this gene is a duplicate marked by *
        counter=to_rep.count('*')
    to_rep=to_rep.replace('*','')
    if '-' in to_rep[0]: #checks if the gene is in the complementary strand
        bools=True
        to_rep=to_rep.replace('-','',1)
    for k,v in temp_rep.items(): #iterates over the replacement dictionary and looks for matching gene strings
        if to_rep==v:
            to_rep=k
            break
    to_rep=to_rep+'*'*counter
    if bools: to_rep='-'+to_rep #If gene is in complementry, if the answer is yes 
    return to_rep

def list_of_genes_to_number(genes, number2symbol = False):
    return [gene_symbol_to_number(i, number2symbol = number2symbol) for i in genes]

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

def reverse_complement(seq):
    """
    """
    return str(Seq.Seq(seq).reverse_complement())

def parse_trnascan(file):
    """
    """
    pass

def move_legend(ax, new_loc, **kws):
    old_legend = ax.legend_
    handles = old_legend.legendHandles
    labels = [t.get_text() for t in old_legend.get_texts()]
    title = old_legend.get_title().get_text()
    ax.legend(handles, labels, loc=new_loc, title=title, **kws)
    