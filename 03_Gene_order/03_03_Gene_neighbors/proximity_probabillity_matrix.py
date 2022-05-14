#%%
import pandas as pd
import numpy as np
import re,os
from ast import literal_eval
from sys import platform
from functools import wraps
path = os.getcwd()
if platform != 'win32':
    dash = '/'
else:
    dash = '\\'


###############
#CURRENTLY SET IN NO TRNA MODE!
###############
def adjacent_genes(gorder, gene1, gene2):
  """
  Recieves a single gorder list and two gene names, returns True if the genes
  are adjascent and False if not, takes into account the circularity of mtDNA

  Parameters
  ----------
  gorder : list
    An organism's gene order
  gene1 : str
    A gene name str
  gene2 : str
    A gene name str
  
  Returns
  -------
  bool
    True if the genes are adj, False otherwise (the genes are adj if the absolute value
    of the diff between their indices equals 1 or if they are at both ends
    of the gene order because the mtDNA is circular).
  """
  gindex1 = [i for i,v in enumerate(gorder) if v == gene1] #iterate over all occurances of gene1 and grab indices
  gindex2 = [i for i,v in enumerate(gorder) if v == gene2] #iterate over all occurances of gene1 and grab indices
  for i in gindex1:
    for j in gindex2:
      result = (abs(i - j) == 1) or (set([gorder[-1], gorder[0]]) == set([gene1, gene2]))
      if result: return result
  return False



def pair_prev(df, gene1, gene2):
  """
  Runs the adjacent_genes function on every single gene order in df, looking for
  adjacency of a specific gene pair, returns the prevalence of that specific pair.

  Parameters
  ----------
  df : pd.DataFrame
    A dataframe of all my organisms
  gene1 : str
    The first gene in the pair
  gene2 : str
    The second gene in the pair
  
  Returns
  -------
  prev : int
    Prevalence of gene pair adjacency in df
  """
  if gene1 == gene2:
    return 1
  prev = np.mean(df.Gene_order.apply(adjacent_genes, args = [gene1, gene2]))
  return prev

def gorder_fixer(gorder, minimal = False):
  """
  Quick function to fix the gorder (remove codons and duplicate markers - THIS IS PROBLEMATIC)

  Parameters
  ----------
  gorder : list
    List of gorder to 'fix'
  
  Returns
  -------
  gorder : list
    Fixed gorder
  """
  gorder = [re.sub('\(\w{3}\)','',i) for i in gorder]
  if minimal:
    gorder = [i.replace('*', '') for i in gorder]
  else:
    gorder = [i.replace('*', '').replace('-','') for i in gorder]


  return gorder

def prox_mat(df, ph, genelist):
  """
  Recieve the full organism df, a specific phylum and a genelist, return an N x N
  (N = number of organisms in phylum)matrix of all genes and their proximity 
  prevalence, ordered by genelist.

  Parameters
  ----------
  df : pd.DataFrame
    A dataframe of all organisms
  ph : str
    A phlyum to work on
  genelist : list
    A list of gene names to order the returned matrix by
  
  Returns
  -------
  gmatrix : np.Array
    N x N array (N = number of organisms in phylum) of the prevalence of 
    gene x with gene y
  """
  df = df.loc[df.phylum == ph, :]
  df.Gene_order = df.Gene_order.apply(gorder_fixer)
  gmatrix = pd.DataFrame(index = genelist, columns = genelist)
  for ref_gene in genelist:
    for gene in genelist:
      cur_pair = pair_prev(df, ref_gene, gene)
      gmatrix.loc[ref_gene, gene] = cur_pair
      gmatrix.loc[gene, ref_gene] = cur_pair
  return gmatrix

def main():
  """
  Run the prob_mat on all phyla
  """
  for ph in df.phylum.unique():
    mat = prox_mat(df, ph, genelist)
    mat.to_csv(os.path.join(PATH, 'prob_matrices', 'no_trna', f'{ph}_prob.csv'))
# %%
if __name__ == '__main__':
  PATH = os.getcwd()
  df = pd.read_csv(os.path.join(os.path.dirname(os.path.dirname(PATH)), 'DB_csvs', 'final.csv'), index_col = 0)
  df.Gene_order = df.Gene_order.apply(literal_eval)
  df.Gene_order = df.Gene_order.apply(lambda x: [i for i in x if 'trn' not in i])
  genelist = [i.replace('-','').replace('*','') for i in df.loc['Homo sapiens', 'Gene_order']]
  genelist = [re.sub('\(\w{3}\)','',i) for i in genelist] 
  genelist = [i for i in genelist if 'trn' not in i] #TRNA FILTERING MODE
  genelist = list(set(genelist))
  main()