"""
FLANKING REGIONS CODE
This code filters the output of annotation.py with the desired flanking regions
TODO: Make this code a stand-alone unix pipeline step.
"""

import numpy as np
import pandas as pd
import os
from importlib import reload
import utils
PATH = os.getcwd()
reload (utils)

def num_sim(n1, n2):
  """ calculates a similarity score between 2 numbers """
  return 1 - abs(n1 - n2) / (n1 + n2)

def genomic_ranges(sample):
    """
    Based on a sample df made by the annotation.py script, return a list in the following fromat:
    [start, end, gene_name]

    Parameters
    ----------
    sample : pd.DataFrame
      A dataframe made using the annotation.py code. MUST HAVE COLUMNS:
      Gene, Position

    Returns
    -------
    granges : list
      List of all the granges in the following format:
      [start, end, gene_name]
    """
    genelist = list(sample.Gene.dropna().unique())
    try: genelist.remove(np.nan)
    except ValueError: pass
    sample.fillna
    granges = []
    for gene in genelist:
        cur_gene = sample.loc[sample.Gene == gene, :]
        cur_pos = cur_gene.Position.to_list()

        zero_ind = cur_pos[0]
        last_ind = cur_pos[-1]
        granges.append([zero_ind, last_ind, gene])
    return granges
  
def wrap(start, end, pos, flank, size):
    """
    If the current position + the flanking region or the current position - the flanking region exceeds the 
    boundaries of the mtDNA, wrap around (because the mtDNA is circular)

    Parameters
    ----------
    start : int
      Gene start position
    end : int
      Gene ending position
    pos : int
      Current position (in the loop)
    flank : int
      The size of the flanking region we want
    size : int
      The mtDNA genome size
    
    Returns 
    -------
    bool
      A boolean that is True if the current position wraps around the mtDNA and is flanking.
    """
    if start - flank < 0:
        return pos >= start - flank + size
    elif start + flank > size:
        return pos < start + flank - size
    elif end - flank < 0:
        return pos >= end - flank + size
    elif end + flank > size:
        return pos < end + flank - size
    else: return False

def downstream_wrap(end, pos, flank, size):
  if end + flank > size:
    return pos < end + flank - size   
  else: return False

def flanking(pos, genelist, flank, gsize, direction):
    """
    Receive a genelist, return True if the current position is within flanking range (including wraps)
    
    Parameters
    ----------
    pos : int
      The current position to check
    genelist : list
      The list to compare the position with, MUST be in the following format: [start, end, gene_name]
    flank : int
      The flanking region size we are interested in.
    gsize : int
      The size of the mtDNA genome.
    
    Returns
    -------
    bool
      True if the pos is within flanking range.
    """
    if direction == 'both':
      return any([(grange[0] - flank <= pos < grange[0] + flank) or (grange[1] - flank <= pos < grange[1] + flank) or (wrap(grange[0], grange[1], pos, flank, gsize)) for grange in genelist])
    elif direction == 'downstream':
      return any([(grange[1] < pos < grange[1] + flank) or (downstream_wrap(grange[1], pos, flank, gsize)) for grange in genelist])
    elif direction == 'upstream':pass #TODO(Noam): add upstream option.

def main(sample, flank, org, filtered = True, direction = 'both', only_code = False): #TODO - nogenes
  """
  The main function, receives a sample and flank window, returns a modified dataframe that only contains rows that are within the flank range of all genes

  Parameters
  ----------
  sample : pd.DataFrame
    The dataframe to work on (MUST be in the annotation.py format)
  flank : int
    The wanted flank size
  
  Returns
  -------
  pd.DataFrame
    The filtered dataframe that only contains rows that are within the flank range.
  """
  try:
    org_df = pd.read_csv('final.csv', index_col = 'organism')
  except FileNotFoundError:
    raise FileNotFoundError('final.csv must be in the same folder as this code!\n')
  phylum = org_df.loc[org, 'phylum']
  granges = genomic_ranges(sample)
  flanks = sample.Position.apply(flanking, args = [granges, flank, len(sample), direction])
  if filtered:
    if not only_code:
      return sample.loc[flanks, :]
    else:
      return sample.loc[flanks & sample['Gene'].isnull(), :]
  else:
    sample2  = sample.copy()
    sample2.loc[~flanks, sample2.columns != 'Position'] = 0
    return sample2

